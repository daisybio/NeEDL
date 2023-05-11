//
// Created by juli on 04.08.22.
//

#ifndef GENEPISEEKER_JOINTDEGREEANALYSER_HPP
#define GENEPISEEKER_JOINTDEGREEANALYSER_HPP

#include "Job.hpp"

namespace epi {

    class JointDegreeAnalyser : public Job {
    public:
        explicit JointDegreeAnalyser(std::string name);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        std::string name;
    };

} // epi

#ifdef HEADER_ONLY
#include "JointDegreeAnalyser.cpp"
#endif


#endif //GENEPISEEKER_JOINTDEGREEANALYSER_HPP
