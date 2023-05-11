//
// Created by juli on 20.04.23.
//

#include "PlinkSetHHMissing.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkSetHHMissing::PlinkSetHHMissing(std::string input_path, std::string output_path, std::string ext_path,
                                         int num_threads) {

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;

    }

    void PlinkSetHHMissing::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("filtering out heterozygous haploid and nonmale Y chromosome genotype calls");

        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--set-hh-missing",
                "--noweb",
                "--make-bed",
                "--out",
                output_path
        });

        remove_plink_temp_files(output_path, { ".log", ".hh", ".nosex" });

        logger.stop();

    }
} // epi