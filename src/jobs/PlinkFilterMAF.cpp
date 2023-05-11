//
// Created by juli on 17.04.23.
//

#include "PlinkFilterMAF.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkFilterMAF::PlinkFilterMAF(std::string input_path, std::string output_path, std::string ext_path,
                                   int num_threads) {

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;

    }

    void PlinkFilterMAF::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("MAF filtering");

        Logger::logLine("Rule for MAF threshold: " + maf_threshold_small + " IF num_samples <= 50000 ELSE " + maf_threshold_large);

        CSVParser parser;
        parser.parse(input_path + ".fam");

        std::string maf_threshold = maf_threshold_small;
        // if(parser.num_rows() > 50000)
        //    maf_threshold = maf_threshold_large;

        Logger::logLine("Selected MAF threshold " + maf_threshold);

        Logger::logLine("Filter based on selected MAF with plink");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--maf",
                maf_threshold,
                "--noweb",
                "--make-bed",
                "--out",
                output_path
        });

        remove_plink_temp_files(output_path, { ".hh", ".log", ".nosex"});

        logger.stop();
    }
} // epi