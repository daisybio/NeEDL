//
// Created by juli on 25.05.22.
//

#include "DataModel.hpp"

namespace epi {
    DataModel::DataModel(bool log_to_stdout) {
        std::random_device rd;
        size_t num_threads = omp_get_max_threads();
        for (size_t i = 0; i < num_threads; i++) {
            std::random_device::result_type seed_data[(random_device_type::state_size * sizeof(typename random_device_type::result_type)) / sizeof(rd()) +1];
            std::generate(std::begin(seed_data), std::end(seed_data), std::ref(rd));
            std::seed_seq seeds(std::begin(seed_data), std::end(seed_data));

            random_device.emplace_back(seeds);
            random_device.back().discard(700000);
        }

        Logger::setConsoleLogging(log_to_stdout);
    }

    DataModel::DataModel(const std::string& output_directory, bool log_to_stdout, bool log_to_file) : DataModel(log_to_stdout) {
        outputDirectory = std::make_shared<OutputDirectory>(output_directory, random_device[omp_get_thread_num()]);
        Logger::setFileLogging(log_to_file, outputDirectory->get_output_directory());
    }
} // epi