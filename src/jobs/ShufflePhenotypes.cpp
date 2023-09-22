//
// Created by juli on 22.09.23.
//

#include "ShufflePhenotypes.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    void ShufflePhenotypes::run(std::shared_ptr<DataModel> data) {
        if (data->snpStorage == nullptr) {
            throw epi::Error("Cannot shuffle phenotypes because no instance was loaded yet.");
            return;
        }

        TimeLogger logger("shuffling phenotypes");

        data->snpStorage->shuffle_phenotypes();

        logger.stop();
    }
} // epi