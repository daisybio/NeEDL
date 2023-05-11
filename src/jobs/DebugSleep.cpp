//
// Created by juli on 05.10.22.
//

#include "DebugSleep.hpp"
#include "../util/TimeLogger.hpp"
#include <thread>

namespace epi {
    DebugSleep::DebugSleep(unsigned long sleep_ms) {
        this->sleep_ms = sleep_ms;
    }

    void DebugSleep::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("sleeping " + std::to_string(sleep_ms) + " ms");
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
        logger.stop();
    }

    rapidjson::Value DebugSleep::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("DebugSleep"), doc.GetAllocator());
        obj.AddMember("sleep_ms", rapidjson::Value().SetUint64(sleep_ms), doc.GetAllocator());
        return obj;
    }
} // epi