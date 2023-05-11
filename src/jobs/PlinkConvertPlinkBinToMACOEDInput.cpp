//
// Created by juli on 01.05.23.
//

#include "PlinkConvertPlinkBinToMACOEDInput.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkConvertPlinkBinToMACOEDInput::PlinkConvertPlinkBinToMACOEDInput(std::string input_path, std::string output_path,
                                                                         std::string ext_path, std::string phenotype, int num_threads) {

        if (phenotype != "DICHOTOMOUS") {
            throw epi::Error("MACOED files can only be created for DICHOTOMOUS phenotypes.");
        }

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
    }

    void PlinkConvertPlinkBinToMACOEDInput::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("convert binary plink data to macoed format");

        Logger::logLine("Recode additive phenotype and convert binary plink data to raw text format");
        run_plink_command(ext_path, {
            "--threads",
            std::to_string(num_threads),
            "--bfile",
            input_path,
            "--recode",
            "A",
            "--noweb",
            "--out",
            output_path + "_temp"
        });

        Logger::logLine("Create macoed file from raw file");
        std::ifstream plink_raw (output_path + "_temp.raw");

        // std::ofstream json_file (output_path);
        std::ofstream macoed_file;
        const size_t buffer_size = 524'288'000;
        char *macoed_buffer = new char[buffer_size];
        macoed_file.rdbuf()->pubsetbuf(macoed_buffer, buffer_size);
        macoed_file.open(output_path);

        if (!macoed_file.good()) {
            throw epi::Error("Cannot open file \"" + output_path + "\" for writing macoed dataset");
        }

        if (!plink_raw.good()) {
            throw epi::Error("Cannot open file \"" + output_path + "_temp.raw\"");
        }


        std::string line;
        std::getline(plink_raw, line);
        auto header = string_split(line, ' ', false);

        if (header.size() < 7) throw epi::Error("Raw plink file does not have enough columns.");

        // write header line with snp names
        for (size_t i = 6; i < header.size(); ++i) {
            size_t end_pos = header[i].size();
            for (int j = header[i].size(); j >= 0; --j) {
                if (header[i][j] == '_') {
                    end_pos = j;
                    break;
                }
            }
            macoed_file << header[i].substr(0, end_pos) << ',';
        }
        macoed_file << "Class\n";

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
                auto splits = string_split(lines_buffer[i], ' ', false);
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
                    macoed_file << splits[i] << ',';
                }

                // add phenotype as last column
                macoed_file << std::stoi(splits[5]) - 1 << '\n';

                ++row_index;
            }
        }


        macoed_file.close();
        delete[] macoed_buffer;


        plink_raw.close();
        logger.stop();

        remove_plink_temp_files(output_path + "_temp", { ".log", ".raw", ".hh", ".nosex" });
    }
} // epi