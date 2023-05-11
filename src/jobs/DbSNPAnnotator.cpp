//
// Created by juli on 26.05.22.
//

#include "DbSNPAnnotator.hpp"

#include <utility>
#include "../util/helper_functions.hpp"
#include "../util/TimeLogger.hpp"
#include "../data_model/SNP_t.hpp"

namespace epi {
    DbSNPAnnotator::DbSNPAnnotator(const std::string& data_directory)
    : SnpCsvAnnotator(false, '\t',  ';', -1) {
        FileFinder ff;
        ff.add_ends_with(".csv");
        ff.add_contains("snps_restruc_full");
        path = find_file_by_ending(data_directory + "dbSNP/inc_pseudogenes/", ff);

        check_columns(1, 0);
    }

} // epi