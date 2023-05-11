//
// Created by juli on 21.10.22.
//

#ifndef GENEPISEEKER_CREATERANDOMSETS_HPP
#define GENEPISEEKER_CREATERANDOMSETS_HPP

#include "Job.hpp"

namespace epi {

    class CreateRandomSets : public Job {
    public:
        CreateRandomSets(size_t num_sets);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;

    private:
        size_t num_sets;
    };

} // epi

#ifdef HEADER_ONLY
#include "CreateRandomSets.cpp"
#endif

#endif //GENEPISEEKER_CREATERANDOMSETS_HPP
