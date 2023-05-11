//
// Created by juli on 24.04.23.
//

#include "InternalBimLoader.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    InternalBimLoader::InternalBimLoader(std::string bim_path) {
        this->bim_path = bim_path;
    }

    void InternalBimLoader::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("loading GWAS data");
        Logger::logLine("input file: " + bim_path);

        auto storage = std::make_shared<SNPStorage_WithoutGeno>(bim_path);
        auto snpStorage = std::static_pointer_cast<SNPStorage>(storage);
        data->snpStorage = snpStorage;
        SNPStorage::currentSnpStorage = snpStorage;

        logger.stop();
    }
} // epi