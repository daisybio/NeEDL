//
// Created by juli on 27.05.22.
//

#include "InstanceSaver.hpp"
#include "../util/TimeLogger.hpp"

#include <utility>

namespace epi {
    InstanceSaver::InstanceSaver(std::string path) {
        this->path = std::move(path);
    }

    void InstanceSaver::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("saving GWAS data");
        if (data->snpStorage == nullptr) {
            throw epi::Error("Cannot save the instance as no instance was loaded yet.");
            return;
        }

        data->snpStorage->save(path);
        logger.stop();
    }

    rapidjson::Value InstanceSaver::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("InstanceSaver"), doc.GetAllocator());
        obj.AddMember("path", rapidjson::Value().SetString(path.c_str(), path.size(), doc.GetAllocator()), doc.GetAllocator());
        return obj;
    }
} // epi