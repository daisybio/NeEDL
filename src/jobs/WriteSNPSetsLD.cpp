//
// Created by juli on 29.11.23.
//

#include "WriteSNPSetsLD.hpp"
#include "../util/TimeLogger.hpp"
#include <filesystem>

namespace epi {
    WriteSNPSetsLD::WriteSNPSetsLD(const std::string &directory_name, const std::string &rank_model) {
        this->directory_name = directory_name;
        if (!rank_model.empty()) {
            this->rank_model = options::epistasis_score_from_string(rank_model);
        } else {
            this->rank_model = {};
        }
    }

    void WriteSNPSetsLD::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("Creating LD matrices for each SNP set");

        std::filesystem::path path;
        if (root_path.has_value()) {
            path = std::filesystem::path(root_path.value()) / directory_name;
        } else if (data->outputDirectory != nullptr) {
            path = data->outputDirectory->get_free_filepath(directory_name, "/");
        } else {
            throw epi::Error("No output directory specified where the LD matrices should be stored.");
        }

        // create directory
        std::filesystem::create_directory(path);


        auto all_sets = data->snpSetStorage;
        // rank if requested
        if (rank_model.has_value()) {
            bool need_minimizing = data->snpStorage->need_minimizing_score(rank_model.value());
            std::sort(all_sets.begin(), all_sets.end(), [&, need_minimizing] (auto a, auto b) -> bool {
                double score_a = a.calculate_score(rank_model.value());
                double score_b = b.calculate_score(rank_model.value());
                return (need_minimizing && score_a < score_b) || (!need_minimizing && score_a > score_b);
            });
        }

        for (size_t x = 0; x < all_sets.size(); ++x) {
            auto &snp_set = all_sets[x];
            std::ofstream matrix_file (path / ("SNP_set_" + std::to_string(x) + ".csv"));

            for (size_t i = 0; i < snp_set.size(); ++i) {
                matrix_file << '\t';
                matrix_file << data->snpStorage->snp_get_name(snp_set[i]);
            }
            matrix_file << '\n';

            for (size_t i = 0; i < snp_set.size(); ++i) {
                matrix_file << data->snpStorage->snp_get_name(snp_set[i]);

                for (size_t j = 0; j < i; ++j) {
                    matrix_file << " \t";
                }
                for (size_t j = i; j < snp_set.size(); ++j) {
                    matrix_file << '\t';
                    matrix_file << Logger::to_string(data->snpStorage->calculate_LD(snp_set[i], snp_set[j]));
                }
                matrix_file << '\n';
            }
        }

        logger.stop();
    }

    void WriteSNPSetsLD::set_root_path(const std::string &root_path_) {
        this->root_path = root_path_;
    }
} // epi