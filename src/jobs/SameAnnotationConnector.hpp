//
// Created by juli on 26.05.22.
//

#ifndef GENEPISEEKER_SAMEANNOTATIONCONNECTOR_HPP
#define GENEPISEEKER_SAMEANNOTATIONCONNECTOR_HPP

#include "Job.hpp"

namespace epi {

    class SameAnnotationConnector : public Job {
    public:
        SameAnnotationConnector() = default;
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    };

} // epi

#ifdef HEADER_ONLY
#include "SameAnnotationConnector.cpp"
#endif


#endif //GENEPISEEKER_SAMEANNOTATIONCONNECTOR_HPP
