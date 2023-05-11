//
// Created by juli on 03.04.23.
//

#include "DiseaseSNPReader.hpp"
#include "../util/FileFinder.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    DiseaseSNPReader::DiseaseSNPReader(std::string disease_snp_file, std::string plink_filter_file) {
        FileFinder ff;
        this->path = find_file_by_ending(std::move(disease_snp_file), ff);

        this->plink_filter_file = plink_filter_file;
    }

    void DiseaseSNPReader::run(std::shared_ptr<DataModel> data) {
        Logger::logLine("Reading disease SNPs");
        std::ifstream file(path);
        std::string line;
        if(file.is_open()){
            getline(file, line);
            boost::replace_all(line, "\"", "");
            auto disease_snps = string_split(line, ',');

            // remove duplicates
            std::unordered_set<std::string> disease_snps_set (disease_snps.begin(), disease_snps.end());

            data->disease_snps = { disease_snps_set.begin(), disease_snps_set.end() };
        } else {
            throw epi::Error("Cannot open disease SNPs file.");
        }
        file.close();

        Logger::logLine("found " + std::to_string(data->disease_snps.size()) + " disease SNPs");

        if (!plink_filter_file.empty()) {
            // read the missing file and select loci without missing entries
            Logger::logLine("Filter disease SNPs");
            CSVParser loci_parser;
            loci_parser.parse(plink_filter_file + ".bim", '\t');

            std::unordered_set<std::string> selected_loci;
            for (size_t i = 1; i < loci_parser.num_rows(); ++i) {
                if (std::stoi(loci_parser.cell(i, 2)) == 0) {
                    selected_loci.insert(loci_parser.cell(i, 1));
                    selected_loci.insert(loci_parser.cell(i, 0) + ":" + loci_parser.cell(i, 1));
                }
            }

            // filter disease SNPs
            std::vector<std::string> disease_snps_filtered;
            for (auto &snp: data->disease_snps) {
                if (selected_loci.find(snp) != selected_loci.end()) {
                    disease_snps_filtered.push_back(snp);
                }
            }
            Logger::logLine("filtered disease SNPs: " + std::to_string(disease_snps_filtered.size()) + " of " +
                            std::to_string(data->disease_snps.size()) + " SNPs retained");
            data->disease_snps = disease_snps_filtered;
        }
    }
} // epi