//
// Created by juli on 26.05.22.
//

#ifndef GENEPISEEKER_TIMELOGGER_HPP
#define GENEPISEEKER_TIMELOGGER_HPP

#include "Logger.hpp"
#include <string>

namespace epi {

    class TimeLogger {
    public:
        TimeLogger(const std::string& desc);
        void stop();
    private:
        std::string desc;
        std::chrono::high_resolution_clock::time_point startTime;
    };

} // epi


#ifdef HEADER_ONLY
#include "TimeLogger.cpp"
#endif

#endif //GENEPISEEKER_TIMELOGGER_HPP
