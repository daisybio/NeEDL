//
// Created by juli on 20.04.23.
//

#include "PlinkFilterRSIds.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkFilterRSIds::PlinkFilterRSIds(std::string input_path, std::string output_path, std::string ext_path, int num_threads) {
        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
    }

    void PlinkFilterRSIds::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("filter for RS-IDs");

        Logger::logLine("find entries that do not start with \"rs\"");
        CSVParser geno_parser;
        geno_parser.parse(input_path + ".bim", '\t');

        std::ofstream exclude_file (output_path + ".notrsids");
        // find all entries that do not start with rs...
        for (size_t i = 0; i < geno_parser.num_rows(); ++i) {
            if (geno_parser.cell(i, 1).substr(0, 2) != "rs") {
                exclude_file << geno_parser.cell(i, 1) << '\n';
            }
        }
        exclude_file.close();

        Logger::logLine("Exclude those variants from the data");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--exclude",
                output_path + ".notrsids",
                "--noweb",
                "--make-bed",
                "--out",
                output_path
        });

        remove_plink_temp_files(output_path, { ".log", ".hh", ".notrsids", ".nosex" });

        logger.stop();
    }
} // epi