//
// Created by juli on 05.10.22.
//

#ifndef GENEPISEEKER_DEBUGSLEEP_HPP
#define GENEPISEEKER_DEBUGSLEEP_HPP

#include "Job.hpp"

namespace epi {

    class DebugSleep : public Job {
    public:
        DebugSleep(unsigned long sleep_ms);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        unsigned long sleep_ms;
    };

} // epi

#ifdef HEADER_ONLY
#include "DebugSleep.cpp"
#endif

#endif //GENEPISEEKER_DEBUGSLEEP_HPP
