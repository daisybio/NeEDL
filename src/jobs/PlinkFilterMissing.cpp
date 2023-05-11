//
// Created by juli on 03.04.23.
//

#include "PlinkFilterMissing.hpp"
#include "../util/helper_functions.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    PlinkFilterMissing::PlinkFilterMissing(std::string input_path, std::string output_path, std::string ext_path, std::string filter_cutoff, int num_threads) {
        /*
        if (!boost::filesystem::exists(input_path + ".bim")) {
            throw epi::Error("Cannot find bim file at \"" + input_path + ".bim\"");
        }

        if (!boost::filesystem::exists(input_path + ".bed")) {
            throw epi::Error("Cannot find bed file at \"" + input_path + ".bed\"");
        }

        if (!boost::filesystem::exists(input_path + ".fam")) {
            throw epi::Error("Cannot find fam file at \"" + input_path + ".fam\"");
        }
         */

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
        this->filter_cutoff = filter_cutoff;
    }

    void PlinkFilterMissing::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("filtering out missing information in dataset");

        Logger::logLine("filter out loci with missing data");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--geno",
                filter_cutoff,
                "--noweb",
                "--make-bed",
                "--allow-no-sex",
                "--out",
                output_path + "_temp"
        });

        Logger::logLine("filter out individuals with missing genotype data or missing phenotype");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                output_path + "_temp",
                "--mind",
                filter_cutoff,
                "--prune",
                "--noweb",
                "--make-bed",
                "--allow-no-sex",
                "--out",
                output_path
        });

        remove_plink_temp_files(output_path + "_temp", { ".log", ".bim", ".fam", ".bed", ".hh", ".nosex" });
        remove_plink_temp_files(output_path, { ".hh", ".log", ".nosex"});

        logger.stop();
    }
} // epi