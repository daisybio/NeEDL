//
// Created by juli on 25.08.22.
//

#ifndef GENEPISEEKER_HELPERPIPELINES_HPP
#define GENEPISEEKER_HELPERPIPELINES_HPP

#include <string>
#include <memory>
#include "../jobs/Job.hpp"

class HelperPipelines {
public:

    static void convert_dataset_to_binary(
              std::string input_path,
              std::string input_format,
              std::string input_phenotype,
              size_t input_num_categories,
              std::string output_path
          );

    static void calculate_scores(
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
            int num_threads = 0,
            bool random_sets = false,
            size_t num_random_sets = 1000,
            bool create_k_mers = false,
            size_t k_min = 0,
            size_t k_max = 0,
            std::string LD_directory = "",
            bool ld_only = false
    );
};


#ifdef HEADER_ONLY
#include "HelperPipelines.cpp"
#endif

#endif //GENEPISEEKER_HELPERPIPELINES_HPP
