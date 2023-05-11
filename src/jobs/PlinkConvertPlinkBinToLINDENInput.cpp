//
// Created by juli on 01.05.23.
//

#include "PlinkConvertPlinkBinToLINDENInput.hpp"
#include "../util/helper_functions.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    PlinkConvertPlinkBinToLINDENInput::PlinkConvertPlinkBinToLINDENInput(std::string input_path,
                                                                         std::string output_path, std::string ext_path,
                                                                         std::string phenotype, int num_threads) {

        if (phenotype != "DICHOTOMOUS") {
            throw epi::Error("LINDEN files can only be created for DICHOTOMOUS phenotypes.");
        }

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;

    }

    void PlinkConvertPlinkBinToLINDENInput::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("convert binary plink data to macoed format");

        Logger::logLine("Recode additive phenotype and convert binary plink data to raw text format");
        run_plink_command(ext_path, {
            "--threads",
            std::to_string(num_threads),
            "--bfile",
            input_path,
            "--recode",
            "A-transpose",
            "--noweb",
            "--out",
            output_path + "_temp"
        });


        Logger::logLine("Create macoed file from traw file");
        std::ifstream plink_raw (output_path + "_temp.traw");
        if (!plink_raw.good()) {
            throw epi::Error("Cannot open file \"" + output_path + "_temp.traw\"");
        }


        std::string line;
        std::getline(plink_raw, line);
        auto header = string_split(line, '\t', false);

        if (header.size() < 7) throw epi::Error("Raw plink file does not have enough columns.");

        std::vector<bool> phenotype_list;
        {
            CSVParser ind_parser;
            ind_parser.parse(input_path + ".fam", ' ');
            if (ind_parser.num_columns() < 6) ind_parser.parse(input_path + ".fam", '\t');

            std::unordered_map<std::string, bool> phenotype_map;
            for (size_t i = 0; i < ind_parser.num_rows(); ++i) {
                phenotype_map.insert({ ind_parser.cell(i, 0) + '_' + ind_parser.cell(i, 1), ind_parser.cell(i, 5) == "2"});
            }

            for (size_t i = 6; i < header.size(); ++i) {
                auto entry = phenotype_map.find(header[i]);
                if (entry == phenotype_map.end()) throw epi::Error("Individual id '" + header[i] + "' not found in .fam file.");
                phenotype_list.push_back(entry->second);
            }
        }


        std::ofstream cases_file (output_path + ".linden.cases");
        std::ofstream controls_file (output_path + ".linden.controls");
        std::ofstream loci_file (output_path + ".linden.loci");




        size_t row_index = 1;

        size_t num_lines_per_thread = num_geno_cells_per_thread / header.size();
        if (num_lines_per_thread == 0) num_lines_per_thread = 1;
        std::vector<std::string> lines_buffer;
        std::vector<std::vector<std::string>> splits_buffer (num_lines_per_thread * num_threads, std::vector<std::string>());
        while(!plink_raw.eof()) {
            lines_buffer.clear();
            // splits_buffer must not be cleared here

            // multi_threaded splitting
            for (size_t i = 0; i < num_lines_per_thread * num_threads && !plink_raw.eof(); ++i) {
                std::getline(plink_raw, line);
                lines_buffer.push_back(line);
            }

#pragma omp parallel for default(none) shared(lines_buffer, splits_buffer) schedule(static)
            for (size_t i = 0; i < lines_buffer.size(); ++i) {
                auto splits = string_split(lines_buffer[i], '\t', false);
                splits_buffer[i] = splits;
            }

            for (size_t j = 0; j < lines_buffer.size(); ++j) {
                auto &splits = splits_buffer[j];

                if (splits.empty()) continue;
                if (splits.size() != header.size()) throw epi::Error(
                            "Raw plink file contains rows with varying number of columns: header has " +
                            std::to_string(header.size()) + " cols and line " + std::to_string(row_index) + " has " +
                            std::to_string(splits.size()) + " cols");


                for (size_t i = 6; i < splits.size(); ++i) {
                    if (phenotype_list[i - 6]) {
                        cases_file << splits[i];
                    } else {
                        controls_file << splits[i];
                    }
                }
                cases_file << '\n';
                controls_file << '\n';

                // add phenotype as last column
                loci_file << splits[1] << '\t' << splits[0] << '\t' << splits[3] << '\n';

                ++row_index;
            }
        }


        loci_file.close();
        cases_file.close();
        controls_file.close();


        plink_raw.close();
        logger.stop();

        remove_plink_temp_files(output_path + "_temp", { ".log", ".traw", ".hh", ".nosex" });
    }
} // epi