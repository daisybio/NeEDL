//
// Created by juli on 25.05.22.
//

#ifndef GENEPISEEKER_DATAMODEL_HPP
#define GENEPISEEKER_DATAMODEL_HPP

#include <memory>
#include "SNPStorage.hpp"
#include "SNPNetwork.hpp"
#include "../util/OutputDirectory.hpp"

namespace epi {

    class DataModel {
        typedef std::mt19937 random_device_type;
    public:
        explicit DataModel(bool log_to_stdout);
        DataModel(const std::string& output_directory, bool log_to_stdout, bool log_to_file);

        std::shared_ptr<SNPStorage> snpStorage = nullptr;
        std::shared_ptr<SNPNetwork> snpNetwork = nullptr;
        std::vector<SNPSet> snpSetStorage{};
        std::shared_ptr<OutputDirectory> outputDirectory = nullptr;

        const std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

        std::vector<random_device_type> random_device;


        // currently only used by preprocessing pipelines
        std::vector<std::string> disease_snps;

        struct dataset_stat {
            std::string name{};
            size_t num_variants = 0;
            size_t num_samples = 0;
            size_t num_samples_cases = 0;
            size_t num_samples_controls = 0;
        };
        std::vector<struct dataset_stat> dataset_stats;
    };

} // epi

#ifdef HEADER_ONLY
#include "DataModel.cpp"
#endif

#endif //GENEPISEEKER_DATAMODEL_HPP
