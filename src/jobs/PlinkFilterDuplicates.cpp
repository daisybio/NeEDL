//
// Created by juli on 20.04.23.
//

#include "PlinkFilterDuplicates.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkFilterDuplicates::PlinkFilterDuplicates(std::string input_path, std::string output_path, std::string ext_path,
                                                 int num_threads) {

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;

    }

    void PlinkFilterDuplicates::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("filter duplicates");

        Logger::logLine("Run duplicates test");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--list-duplicate-vars",
                "ids-only",
                "suppress-first",
                "--noweb",
                "--allow-no-sex",
                "--out",
                output_path
        });

        Logger::logLine("Remove reported duplicates from dataset");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--exclude",
                output_path + ".dupvar",
                "--noweb",
                "--make-bed",
                "--out",
                output_path
        });

        remove_plink_temp_files(output_path, { ".log", ".hh", ".dupvar", ".nosex" });

        logger.stop();
    }
} // epi