//
// Created by juli on 26.05.22.
//

#include "BioGridConnector.hpp"
#include "../util/FileFinder.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"
#include <unordered_set>

namespace epi {

    BioGridConnector::BioGridConnector(const std::string& data_directory)
     : NetworkCsvConnector("BIOGRID", true, '\t', -1, -1){

        FileFinder ff;
        ff.add_ends_with(".txt");
        ff.add_contains("BIOGRID");
        path = find_file_by_ending(data_directory + "BIOGRID", ff);

        check_columns("Official Symbol Interactor A", "Official Symbol Interactor B");
    }
} // epi