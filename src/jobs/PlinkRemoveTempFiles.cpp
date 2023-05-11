//
// Created by juli on 17.04.23.
//

#include "PlinkRemoveTempFiles.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkRemoveTempFiles::PlinkRemoveTempFiles(std::string path, std::vector<std::string> endings) {
        this->path = path;
        this->endings = endings;
    }

    void PlinkRemoveTempFiles::run(std::shared_ptr<DataModel> data) {
        remove_plink_temp_files(path, endings);
    }
} // epi