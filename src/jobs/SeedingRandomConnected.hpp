//
// Created by juli on 21.06.22.
//

#ifndef GENEPISEEKER_SEEDINGRANDOMCONNECTED_HPP
#define GENEPISEEKER_SEEDINGRANDOMCONNECTED_HPP

#include "Job.hpp"

namespace epi {

    class SeedingRandomConnected : public Job {
    public:
        SeedingRandomConnected(uint_fast32_t num_seeds);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        uint_fast32_t num_seeds;
    };

} // epi

#ifdef HEADER_ONLY
#include "SeedingRandomConnected.cpp"
#endif


#endif //GENEPISEEKER_SEEDINGRANDOMCONNECTED_HPP
