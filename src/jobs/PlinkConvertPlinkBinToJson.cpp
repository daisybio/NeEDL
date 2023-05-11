//
// Created by juli on 03.04.23.
//

#include "PlinkConvertPlinkBinToJson.hpp"
#include "../util/helper_functions.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    PlinkConvertPlinkBinToJson::PlinkConvertPlinkBinToJson(std::string input_path, std::string output_path, std::string ext_path, std::string phenotype, size_t num_categories, int num_threads) {
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

        if (phenotype != "CATEGORICAL" && phenotype != "DICHOTOMOUS" && phenotype != "QUANTITATIVE") {
            throw epi::Error("Unknown phenotype option. Allowed is one of CATEGORICAL, DICHOTOMOUS, QUANTITATIVE.");
        }

        if (num_categories <= 1) {
            throw epi::Error("At least 2 categories are necessary for categorical phenotype.");
        }

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
        this->num_categories = num_categories;
        this->phenotype = phenotype;
    }

    void PlinkConvertPlinkBinToJson::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("convert binary plink data to json format");

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

        Logger::logLine("Create json file from raw file");
        std::ifstream plink_raw (output_path + "_temp.traw");

        // std::ofstream json_file (output_path);
        std::ofstream json_file;
        const size_t buffer_size = 524'288'000;
        char *json_buffer = new char[buffer_size];
        json_file.rdbuf()->pubsetbuf(json_buffer, buffer_size);
        json_file.open(output_path);

        if (!json_file.good()) {
            throw epi::Error("Cannot open file \"" + output_path + "\" for writing json dataset");
        }

        if (!plink_raw.good()) {
            throw epi::Error("Cannot open file \"" + output_path + "_temp.traw\"");
        }


        // start with writing genotype data while parsing the raw file. save phenotype and snp information in memory and write them afterwards
        std::string line;
        std::getline(plink_raw, line);
        auto header = string_split(line, '\t', false);

        if (header.size() < 7) throw epi::Error("Raw plink file does not have enough columns.");

        struct snp_data {
            std::string name;
            std::string chromosome;
            std::string position;
            std::string allele_from;
            std::string allele_to;
            unsigned long long genotype_sum;
        };
        std::vector<snp_data> snp_data;
        std::unordered_map<std::string, size_t> snp_name_map;

        // start of genotype matrix
        json_file << "{\"genotype\":[";
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

                if (row_index > 1) json_file << ',';
                json_file << '[';


                unsigned long long genotype_sum = 0;
                for (size_t i = 6; i < splits.size(); ++i) {
                    if (i > 6) json_file << ',';
                    json_file << splits[i];
                    genotype_sum += std::stoull(splits[i]);
                }

                snp_name_map.insert({splits[1], snp_data.size()});
                snp_data.push_back({
                                           splits[1],
                                           splits[0],
                                           splits[3],
                                           splits[4],
                                           splits[5],
                                           genotype_sum
                                   });

                json_file << ']';

                ++row_index;
            }
        }
        // end of genotype matrix
        json_file << ']';

        // begin of snps list
        json_file << ",\"snps\":[";

        bool first_snp_data = true;
        for (auto & snp : snp_data) {
            if (first_snp_data) first_snp_data = false;
            else json_file << ',';

            json_file << "[\"" << snp.name << "\",\"" << snp.chromosome << "\",\"" << snp.position << "\",\"" << snp.allele_from << "\",\"" << snp.allele_to << "\"]";
        }

        // end of snps list
        json_file << ']';

        // read phenotype information from fam file (need map as individuals might be sorted differently than in traw file
        CSVParser ind_parser;
        ind_parser.parse(input_path + ".fam", ' ');
        if (ind_parser.num_columns() < 6) ind_parser.parse(input_path + ".fam", '\t');

        std::unordered_map<std::string, std::string> ind_phenotype_map;
        for (size_t i = 0; i < ind_parser.num_rows(); ++i) {
            ind_phenotype_map.insert({
                                             ind_parser.cell(i, 0) + "_" + ind_parser.cell(i, 1),
                                             ind_parser.cell(i, 5)
                                     });
        }

        // begin of phenotype vector
        json_file << ",\"phenotype\":[";
        bool is_dichotomous = phenotype == "DICHOTOMOUS";
        for (size_t i = 6; i < header.size(); ++i) {
            if (i > 6) json_file << ',';

            auto pheno = ind_phenotype_map.find(header[i]);
            if (pheno == ind_phenotype_map.end()) {
                throw epi::Error("Indivdual '" + header[i] + "' not found in .fam file.");
            }

            if (is_dichotomous) {
                json_file << std::stoi(pheno->second) - 1;
            } else {
                json_file << pheno->second;
            }
        }

        // end of phenotype vector
        json_file << ']';

        // begin of disease snps list
        json_file << ",\"disease_snps\":[";
        bool first_disease_snp = true;
        for (auto & snp : data->disease_snps) {
            auto snp_id = snp_name_map.find(snp);
            if (snp_id == snp_name_map.end()) {
                throw epi::Error("Disease SNP " + snp + " cannot be found in provided dataset.");
            } else {
                if (first_disease_snp) first_disease_snp = false;
                else json_file << ',';
                json_file << snp_id->second;
            }
        }

        // end of disease snps list
        json_file << ']';

        // begin of maf list
        json_file << ",\"mafs\":[";
        bool first_maf = true;
        auto default_precision = json_file.precision();
        for (auto & snp : snp_data) {
            if (first_maf) first_maf = false;
            else json_file << ',';

            double maf = .5 * double(snp.genotype_sum) / double(ind_parser.num_rows() - 1);
            json_file << std::setprecision(4) << maf;
        }
        json_file << std::setprecision(default_precision);

        // end of maf list
        json_file << ']';


        // additional info
        json_file << ",\"model_type\":\"" << phenotype << "\"";
        if (phenotype != "QUANTITATIVE") json_file << ",\"num_categories\":" << num_categories;
        json_file << ",\"num_inds\":" << ind_parser.num_rows();
        json_file << ",\"num_snps\":" << row_index - 1;


        // end of json file
        json_file << '}';

        json_file.close();
        delete[] json_buffer;


        plink_raw.close();
        logger.stop();

        remove_plink_temp_files(output_path + "_temp", { ".log", ".traw", ".hh", ".nosex" });
    }
} // epi