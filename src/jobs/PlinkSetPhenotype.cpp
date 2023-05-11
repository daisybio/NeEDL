//
// Created by juli on 26.04.23.
//

#include "PlinkSetPhenotype.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkSetPhenotype::PlinkSetPhenotype(std::string phenotype, std::string input_path, std::string output_path,
                                         std::string ext_path, int num_threads) {

        char *pos;
        strtod(phenotype.c_str(), &pos);
        if (pos == 0) {
            throw epi::Error("Phenotype '" + phenotype + "' cannot be parsed.");
        }

        this->phenotype = phenotype;
        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
    }

    void PlinkSetPhenotype::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("set phenotype");

        Logger::logLine("Read fam file and create phenotype file");
        CSVParser ind_parser;
        ind_parser.parse(input_path + ".fam", ' ');
        if (ind_parser.num_columns() < 6) ind_parser.parse(input_path + ".fam", '\t');

        std::ofstream pheno_file(output_path + ".pheno");

        for (size_t i = 0; i < ind_parser.num_rows(); ++i) {
            pheno_file << ind_parser.cell(i, 0) << '\t' << ind_parser.cell(i, 1) << '\t' << phenotype << '\n';
        }
        pheno_file.close();

        Logger::logLine("Set phenotype with plink");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--pheno",
                output_path + ".pheno",
                "--noweb",
                "--make-bed",
                "--out",
                output_path
        });

        remove_plink_temp_files(output_path, {".log", ".hh", ".pheno", ".nosex"});

        logger.stop();

    }
} // epi