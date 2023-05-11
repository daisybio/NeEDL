//
// Created by juli on 27.05.22.
//

#ifndef GENEPISEEKER_INSTANCESAVER_HPP
#define GENEPISEEKER_INSTANCESAVER_HPP

#include "Job.hpp"

namespace epi {

    class InstanceSaver : public Job {
    public:
        InstanceSaver(std::string path);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        std::string path;
    };

} // epi

#ifdef HEADER_ONLY
#include "InstanceSaver.cpp"
#endif


#endif //GENEPISEEKER_INSTANCESAVER_HPP
