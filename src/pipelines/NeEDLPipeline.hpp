//
// Created by juli on 15.07.22.
//

#ifndef GENEPISEEKER_NEEDL_DEFAULT_HPP
#define GENEPISEEKER_NEEDL_DEFAULT_HPP

#include <string>
#include "../jobs/Job.hpp"
#include "../jobs/MultiNetworkAggregator.hpp"
#include "../jobs/SnpCsvAnnotator.hpp"
#include "../jobs/NetworkCsvConnector.hpp"
#include "../jobs/LocalSearch.hpp"
#include "../jobs/BioGridConnector.hpp"
#include "../jobs/JointDegreeAnalyser.hpp"
#include "../jobs/SaveNetwork.hpp"
#include "../jobs/NetworkStatsPrinter.hpp"
#include "../jobs/ShuffleNetwork.hpp"
#include "../jobs/WriteSets.hpp"
#include "../jobs/MMAFilter.hpp"
#include "../jobs/InstanceLoader.hpp"

/**
 * Represents the default local search pipeline of NeEDL
 */
class NeEDLPipeline {
public:
    NeEDLPipeline(
                             // instance stuff
                             std::string input_path,
                             std::string input_format,
                             std::string input_phenotype,
                             size_t input_num_categories,

                             // network stats
                             bool disable_save_network,
                             bool calculate_advanced_network_stats,
                             bool joint_degree_analysis,

                             // other
                             bool no_additional_scores,
                             std::string output_directory,
                             int num_threads,
                             std::string data_directory = "./data/"
                             );

    // shuffle methods
    void activate_network_shuffling(std::string method);

    // filtering
    void activate_MMA_filter(double cutoff, bool correct_BH);

    // snp annotation
    void add_snp_annotation_source_dbSNP();
    void add_snp_annotation_source_eQTL(const std::vector<std::string>& tissue_selection, double pvalue_cutoff, bool bh_correction);
    void add_snp_annotation_source(const epi::SnpCsvAnnotator& annotation_source);

    // networks to search (eg. BIOGRID)
    template<class SeedingRoutine>
    void add_network_BIOGRID(SeedingRoutine seeding, epi::LocalSearch local_search) {
        auto network_source = std::make_shared<epi::BioGridConnector>(data_directory);
        add_network(network_source->get_name(), network_source, std::make_shared<SeedingRoutine>(seeding), local_search);
    }

    template<class SeedingRoutine>
    void add_network(const epi::NetworkCsvConnector network_source, SeedingRoutine seeding, epi::LocalSearch local_search) {
        add_network(network_source.get_name(), std::make_shared<epi::NetworkCsvConnector>(network_source), std::make_shared<SeedingRoutine>(seeding), local_search);
    }

    // seeding settings
    template<class SeedingRoutine>
    void set_final_local_search(SeedingRoutine seeding, epi::LocalSearch local_search) {
        final_search_pipeline = std::vector<std::shared_ptr<epi::Job>>({
            std::make_shared<epi::NetworkStatsPrinter>(calculate_advanced_network_stats),
        });

        if (do_network_shuffling) {
            // add shuffle method
            final_search_pipeline.push_back(std::make_shared<epi::ShuffleNetwork>(network_shuffling));
            final_search_pipeline.push_back(std::make_shared<epi::NetworkStatsPrinter>(calculate_advanced_network_stats));
        }

        if (joint_degree_analysis) {
            final_search_pipeline.push_back(std::make_shared<epi::JointDegreeAnalyser>("result_joint_degree"));
        }

        final_search_pipeline.push_back(std::make_shared<SeedingRoutine>(seeding));
        // if (disable_save_network) {
            // final_search_pipeline.push_back(std::make_shared<epi::SaveNetwork>("NODE_EDGE_LIST", "result_network"));
            final_search_pipeline.push_back(std::make_shared<epi::SaveNetwork>("SQLITE", "result_network"));
        // }
        final_search_pipeline.push_back(std::make_shared<epi::WriteSets>(local_search.get_model_name(), additional_scores, "result_seeds"));
        final_search_pipeline.push_back(std::make_shared<epi::LocalSearch>(local_search));
        final_search_pipeline.push_back(std::make_shared<epi::WriteSets>(local_search.get_model_name(), additional_scores, "result_results"));
        final_search_pipeline.push_back(std::make_shared<epi::WriteSets>(local_search.get_model_name(), additional_scores, "result_ind_SNP_scores", true));
    }

    // run pipeline
    void run();

private:
    void add_network(std::string name, std::shared_ptr<epi::Job> network_annotator, std::shared_ptr<epi::Job> seeding, epi::LocalSearch local_search);

    std::shared_ptr<epi::DataModel> data;
    std::vector<std::shared_ptr<epi::Job>> preparationPipeline;
    std::vector<std::shared_ptr<epi::Job>> snpAnnotationPipeline;
    std::vector<std::shared_ptr<epi::Job>> final_search_pipeline;
    std::shared_ptr<epi::MultiNetworkAggregator> multiNetworkAggregator;
    std::shared_ptr<epi::MMAFilter> mmaFilter = nullptr;

    std::string data_directory;
    bool disable_save_network;
    bool calculate_advanced_network_stats;
    bool joint_degree_analysis;

    std::string network_shuffling;
    bool do_network_shuffling = false;

    std::vector<std::string> additional_scores{};

    std::shared_ptr<epi::InstanceLoader> instance_loader;
};

#ifdef HEADER_ONLY
#include "NeEDLPipeline.cpp"
#endif

#endif //GENEPISEEKER_NEEDL_DEFAULT_HPP
