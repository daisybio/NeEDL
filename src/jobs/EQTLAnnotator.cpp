//
// Created by juli on 28.01.24.
//

#include "EQTLAnnotator.hpp"
#include "../util/FileFinder.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    EQTLAnnotator::EQTLAnnotator(const std::vector<std::string> &tissue_types, double pvalue_cutoff, bool bh_correction, const std::string &data_directory)
    : SnpCsvAnnotator(true, '\t',  -1, -1) {
        this->tissue_types = { tissue_types.begin(), tissue_types.end() };
        this-> pvalue_cutoff = pvalue_cutoff;
        this->bh_correction = bh_correction;

        FileFinder ff;
        ff.add_ends_with(".csv");
        ff.add_contains("eqtl_mapping");
        path = find_file_by_ending(data_directory + "eqtl_mapping/", ff);

        check_columns(0, 1);

        if (!this->tissue_types.empty()) {
            Logger::logLine("Checking if given tissue types for eQTL mapping are valid.");
            // check tissue types and num columns
            CSVParser parser;
            parser.parse(path, '\t');
            if (parser.num_columns() != 4) throw epi::Error("EQTL annotation file does not contain exactly 4 columns");

            std::unordered_set<std::string> tissues;
            for (size_t i = 1; i < parser.num_rows(); ++i) {
                tissues.insert(parser.cell(i, 3));
            }
            bool check_failed = false;
            for (const auto &tissue: this->tissue_types) {
                if (!tissues.contains(tissue)) {
                    Logger::logLine("tissue " + tissue + " does not exist in eQTL dataset.");
                    check_failed = true;
                }
            }

            if (check_failed) {
                Logger::logLine("Allowed tissues are: ");
                for (const auto &tissue: tissues) {
                    Logger::logLine("     " + tissue);
                }
                throw epi::Error(
                        "Some tissue types were selected for eQTL annotation that do not exist in the dataset.");
            }
        }
    }

    std::vector<bool> EQTLAnnotator::filter_entries(const CSVParser &parser) {
        if (parser.num_columns() != 4) throw epi::Error("EQTL annotation file does not contain exactly 4 columns");

        Logger::logLine("Filter eQTL entries based on selected tissues.");

        // filter by tissue
        std::vector<size_t> indices;
        std::vector<double> pvalues;
        for (size_t i = 1; i < parser.num_rows(); ++i) {
            if (tissue_types.empty() || tissue_types.contains(parser.cell(i, 3))) {
                indices.push_back(i);
                pvalues.push_back(std::stod(parser.cell(i, 2)));
            }
        }

        // do bh filtering
        if (bh_correction) {
            Logger::logLine("Apply Benjamini-Hochberg correction to p-values of eQTL entries");
            benjamini_hochberg_correction(pvalues);
        }

        // create filter mask
        Logger::logLine("Filter eQTL entries with p-value cutoff " + Logger::to_string(pvalue_cutoff));
        std::vector<bool> result (parser.num_rows(), false);
        size_t num_selected = 0;
        for (size_t i = 0; i < indices.size(); ++i) {
            if (pvalues[i] < pvalue_cutoff) {
                result[indices[i]] = true;
                ++num_selected;
            }
        }
        Logger::logLine("selected " + Logger::to_string(num_selected) + " of " + Logger::to_string(indices.size()) + " eQTL entries");

        return result;
    }

} // epi