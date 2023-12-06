//
// Created by juli on 25.08.22.
//

#include "HelperPipelines.hpp"
#include "../jobs/SequentialJob.hpp"
#include "../jobs/InstanceLoader.hpp"
#include "../jobs/InstanceSaver.hpp"
#include "../jobs/WriteSets.hpp"
#include "../jobs/CreateRandomSets.hpp"
#include "../jobs/CreateSetKMers.hpp"
#include "../jobs/ShufflePhenotypes.hpp"
#include "../jobs/WriteSNPSetsLD.hpp"


void
HelperPipelines::convert_dataset_to_binary(std::string input_path, std::string input_format, std::string input_phenotype,
                                          size_t input_num_categories, std::string output_path) {
    epi::SequentialJob seq{{
                                   std::make_shared<epi::InstanceLoader>(
                                           input_path,
                                           input_format,
                                           input_phenotype,
                                           input_num_categories
                                   ),
                                   std::make_shared<epi::InstanceSaver>(output_path)
                           }};

    auto data = std::make_shared<epi::DataModel>(true);
    seq.run(data);
}

void HelperPipelines::calculate_scores(
        std::string input_path,
        std::string input_format,
        std::string input_phenotype,
        size_t input_num_categories,
        std::shared_ptr<epi::Job> snp_reader,
        std::vector<std::string> epistasis_models,
        std::string rank_model,
        std::string out_path,
        std::vector<std::shared_ptr<epi::Job>> snp_annotation_pipeline,
        bool shuffle_phenotypes,
        int num_threads,
        bool random_sets,
        size_t num_random_sets,
        bool create_k_mers,
        size_t k_min,
        size_t k_max,
        std::string ld_directory,
        bool ld_only
        ) {

    if (num_threads < 0) {
        throw epi::Error("Invalid thread number provided.");
    }
    if (num_threads == 0) {
        num_threads = omp_get_num_procs();
    }
    omp_set_num_threads(num_threads);
    Logger::logLine("OpenMP threads: " + std::to_string(omp_get_num_procs()) + " available, " + std::to_string(num_threads) + " selected");

    if (epistasis_models.empty()) throw epi::Error("At least one epistasis model needed to perform actions");

    epi::SequentialJob seq {{
        std::make_shared<epi::InstanceLoader>(
                input_path,
                input_format,
                input_phenotype,
                input_num_categories
        )
    }};

    if (shuffle_phenotypes) {
        seq.add(std::make_shared<epi::ShufflePhenotypes>());
    }

    seq.add(snp_annotation_pipeline);
    seq.add(snp_reader);

    if (random_sets) {
        seq.add(std::make_shared<epi::CreateRandomSets>(num_random_sets));
    }

    if(create_k_mers) {
        seq.add(std::make_shared<epi::CreateSetKMers>(k_min, k_max));
    }

    if (!ld_only) {
        auto writer = std::make_shared<epi::WriteSets>(rank_model, epistasis_models);
        writer->outfile_path(out_path);
        seq.add(writer);
    }

    if (!ld_directory.empty()) {
        auto write_ld = std::make_shared<epi::WriteSNPSetsLD>("results_LD", rank_model);
        write_ld->set_root_path(ld_directory);
        seq.add(write_ld);
    }

    auto data = std::make_shared<epi::DataModel>(true);
    seq.run(data);
}
