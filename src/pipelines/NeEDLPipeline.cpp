//
// Created by juli on 15.07.22.
//

#include "NeEDLPipeline.hpp"

#include "../../src/jobs/InstanceLoader.hpp"
#include "../../src/jobs/SequentialJob.hpp"
#include "../../src/jobs/DbSNPAnnotator.hpp"
#include "../../src/jobs/SameAnnotationConnector.hpp"
#include "../../src/jobs/BioGridConnector.hpp"
#include "../../src/jobs/NetworkStatsPrinter.hpp"
#include "../../../src/jobs/InstanceSaver.hpp"
#include "../../../src/jobs/MMAFilter.hpp"
#include "../../../src/jobs/SeedingRandomConnected.hpp"
#include "../../../src/jobs/WriteSets.hpp"
#include "../../../src/jobs/SaveNetwork.hpp"
#include "../../../src/jobs/SeedingQuantumComputing.hpp"
#include "../../../src/jobs/LocalSearch.hpp"
#include "../jobs/ShuffleNetwork.hpp"
#include "../util/RepeaterList.hpp"
#include "../jobs/JointDegreeAnalyser.hpp"
#include "../jobs/MultiNetworkAggregator.hpp"
#include "../jobs/SnpCsvAnnotator.hpp"
#include "../jobs/ShinyAppLauncher.hpp"


NeEDLPipeline::NeEDLPipeline(std::string input_path, std::string input_format, std::string input_phenotype,
                                  size_t input_num_categories, std::string covariates_file, bool disable_save_network, bool calculate_advanced_network_stats, bool joint_degree_analysis,
                                  bool no_additional_scores, std::string output_directory,
                                  int num_threads, std::string data_directory) {

    this->data_directory = data_directory;
    this->calculate_advanced_network_stats = calculate_advanced_network_stats;
    this->joint_degree_analysis = joint_degree_analysis;
    this->disable_save_network = disable_save_network;

    if (!no_additional_scores) {
        additional_scores = epi::options::get_all_epistasis_scores(!covariates_file.empty());
    }

    if (output_directory.empty()) {
        data = std::make_shared<epi::DataModel>(true);
    } else {
        data = std::make_shared<epi::DataModel>(output_directory, true, true);
    }

    if (num_threads < 0) {
        throw epi::Error("Invalid thread number provided.");
    }
    if (num_threads == 0) {
        num_threads = omp_get_num_procs();
    }
    omp_set_num_threads(num_threads);
    Logger::logLine("OpenMP threads: " + std::to_string(omp_get_num_procs()) + " available, " + std::to_string(num_threads) + " selected");

    // create preparation pipeline (everything before the multi-network stuff)

    instance_loader =
            std::make_shared<epi::InstanceLoader>(
                    input_path,
                    input_format,
                    input_phenotype,
                    input_num_categories,
                    covariates_file
            );

    preparationPipeline = {
            instance_loader
    };

    snpAnnotationPipeline = {};
    multiNetworkAggregator = std::make_shared<epi::MultiNetworkAggregator>();
    final_search_pipeline = {};
}

void NeEDLPipeline::run() {
    if (snpAnnotationPipeline.empty()) {
        throw epi::Error("No annotation source specified.");
    }

    if (multiNetworkAggregator->num_networks() == 0) {
        throw epi::Error("No networks are specified.");
    }

    std::vector<std::shared_ptr<epi::Job>> all_combined;
    all_combined.insert(all_combined.end(), preparationPipeline.begin(), preparationPipeline.end());

    if (mmaFilter != nullptr) all_combined.push_back(mmaFilter);

    all_combined.insert(all_combined.end(), snpAnnotationPipeline.begin(), snpAnnotationPipeline.end());
    all_combined.push_back(std::make_shared<epi::SameAnnotationConnector>());
    all_combined.push_back(std::make_shared<epi::NetworkStatsPrinter>(calculate_advanced_network_stats));
    all_combined.push_back(multiNetworkAggregator);
    all_combined.insert(all_combined.end(), final_search_pipeline.begin(), final_search_pipeline.end());

    epi::SequentialJob all(all_combined);
    all.add(std::make_shared<epi::ShinyAppLauncher>(instance_loader));

    // save pipeline config
    if(data->outputDirectory != nullptr) {
        epi::Job::save_pipeline_config(std::make_shared<epi::SequentialJob>(all), data->outputDirectory->get_output_directory() + "pipeline_config.json");
    }

    // actually run the pipeline
    all.run(data);

    Logger::logLine("NeEDL-Pipeline finished successfully.");
    Logger::logResourceUsage();
}

void NeEDLPipeline::add_snp_annotation_source_dbSNP() {
    snpAnnotationPipeline.push_back(std::make_shared<epi::DbSNPAnnotator>(data_directory));
}

void NeEDLPipeline::add_snp_annotation_source(const epi::SnpCsvAnnotator& annotation_source) {
    snpAnnotationPipeline.push_back(std::make_shared<epi::SnpCsvAnnotator>(annotation_source));
}

void NeEDLPipeline::add_network(std::string name, std::shared_ptr<epi::Job> network_annotator, std::shared_ptr<epi::Job> seeding, epi::LocalSearch local_search) {
    std::vector<std::shared_ptr<epi::Job>> joblist = {
            network_annotator,
            std::make_shared<epi::NetworkStatsPrinter>(calculate_advanced_network_stats)
    };
    if (do_network_shuffling) {
        // add shuffle method
        joblist.push_back(std::make_shared<epi::ShuffleNetwork>(network_shuffling));
        joblist.push_back(std::make_shared<epi::NetworkStatsPrinter>(calculate_advanced_network_stats));
    }
    if (joint_degree_analysis) {
        joblist.push_back(std::make_shared<epi::JointDegreeAnalyser>(name + "_joint_degree"));
    }

    // seeding before save network --> leiden clusters can be saved as well
    joblist.push_back(seeding);

    if (!disable_save_network) {
        // joblist.push_back(std::make_shared<epi::SaveNetwork>("NODE_EDGE_LIST", name + "_network"));
        joblist.push_back(std::make_shared<epi::SaveNetwork>("SQLITE", name + "_network"));
    }
    joblist.push_back(std::make_shared<epi::WriteSets>(local_search.get_model_name(), additional_scores, name + "_seeds"));
    joblist.push_back(std::make_shared<epi::LocalSearch>(local_search));
    joblist.push_back(std::make_shared<epi::WriteSets>(local_search.get_model_name(), additional_scores, name + "_results"));
    joblist.push_back(std::make_shared<epi::WriteSets>(local_search.get_model_name(), additional_scores, name + "_ind_SNP_scores", true));



    multiNetworkAggregator->add(std::make_shared<epi::SequentialJob>(joblist), name);
}


void NeEDLPipeline::activate_network_shuffling(std::string method) {
    network_shuffling = method;
    do_network_shuffling = true;
}

void NeEDLPipeline::activate_MMA_filter(double cutoff, bool correct_BH) {
    mmaFilter = std::make_shared<epi::MMAFilter>(cutoff, correct_BH);
}

