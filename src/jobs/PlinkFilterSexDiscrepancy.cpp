//
// Created by juli on 17.04.23.
//

#include "PlinkFilterSexDiscrepancy.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkFilterSexDiscrepancy::PlinkFilterSexDiscrepancy(std::string input_path, std::string output_path,
                                                         std::string ext_path, int num_threads) {

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;

    }

    void PlinkFilterSexDiscrepancy::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("filtering out samples with sex discrepancy");

        Logger::logLine("Impute sex based on genome.");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--check-sex",
                "--noweb",
                "--out",
                output_path
        });

        Logger::logLine("Find problematic samples");

        CSVParser parser;
        parser.parse(output_path + ".sexcheck", ' ');

        std::ofstream filtered_file (output_path + ".sexcheck_filtered");
        for (size_t i = 0; i < parser.num_rows(); ++i) {
            if (parser.cell(i, 4) == "OK") {
                filtered_file << parser.cell(i, 0) << " " << parser.cell(i, 1) << '\n';
            }
        }
        filtered_file.close();

        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--keep",
                output_path + ".sexcheck_filtered",
                "--noweb",
                "--make-bed",
                "--out",
                output_path
        });

        remove_plink_temp_files(output_path, { ".log", ".hh", ".sexcheck", ".sexcheck_filtered", ".nosex" });

        logger.stop();

    }
} // epi