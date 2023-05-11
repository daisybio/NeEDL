//
// Created by juli on 26.05.22.
//

#include "TimeLogger.hpp"

namespace epi {
    void TimeLogger::stop() {
        std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
        float secs = (std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count() /
                      1000000000.0f);

        std::string logline = "Finished " + desc + " in ";
        if (secs > 120.f) {
            logline += std::to_string(secs / 60.f) + " min";
        } else {
            logline += std::to_string(secs) + " sec";
        }
        Logger::logLine(logline);
    }

    TimeLogger::TimeLogger(const std::string& desc_) {
        this->desc = desc_;
        Logger::logLine("Begin " + desc_);
        startTime = std::chrono::high_resolution_clock::now();
    }
} // epi