//
// Created by juli on 28.01.24.
//

#include "EQTLAnnotator.hpp"
#include "../util/FileFinder.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    EQTLAnnotator::EQTLAnnotator(const std::vector<std::string> &tissue_types, const std::string &data_directory) {
        this->tissue_types = { tissue_types.begin(), tissue_types.end() };

        FileFinder ff;
        ff.add_ends_with(".csv");
        ff.add_contains("eqtl_mapping_p0.01");
        path = find_file_by_ending(data_directory + "dbSNP/eqtl_mapping/", ff);

        check_columns(0, 1);
    }

    std::vector<bool> EQTLAnnotator::filter_entries(const CSVParser &parser) {
        if (parser.num_columns() < 3) {
            throw epi::Error("Tissue filter column missing in eqtl mapping file.");
        }
        std::vector<bool> result;
        result.reserve(parser.num_rows());
        for (size_t i = 0; i < parser.num_rows(); ++i) {
            result[i] = tissue_types.contains(parser.cell(i, 2));
        }

        return result;
    }

} // epi