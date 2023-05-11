//
// Created by juli on 25.04.23.
//

#include "PlinkMergeDatasets.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkMergeDatasets::PlinkMergeDatasets(std::vector<std::string> input_paths, std::string output_path,
                                           std::string ext_path, int num_threads) {

        this->input_paths = input_paths;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
    }

    void PlinkMergeDatasets::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("merge input files");

        std::ofstream mergelist (output_path + ".mergelist");
        for (auto & input_path : input_paths) {
            mergelist << input_path << '\n';
        }
        mergelist.close();

        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--merge-list",
                output_path + ".mergelist",
                "--noweb",
                "--allow-no-sex",
                "--make-bed",
                "--out",
                output_path
        });

        remove_plink_temp_files(output_path, { ".log", ".hh", ".mergelist", ".nosex" });

        logger.stop();
    }
} // epi