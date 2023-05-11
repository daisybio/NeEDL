//
// Created by juli on 25.05.22.
//

#ifndef GENEPISEEKER_JOB_HPP
#define GENEPISEEKER_JOB_HPP


#include "../../ext/rapidjson/document.h"
#include "../data_model/DataModel.hpp"

/**
 * @brief This is an interface that can is implemented for every job in a possible epistasis pipeline
 */
namespace epi {
    class Job {
    public:
        /**
         * @brief this method will be called to run this job
         */
        virtual void run(std::shared_ptr<DataModel> data) = 0;
        virtual rapidjson::Value getConfig(rapidjson::Document & doc);
        virtual ~Job() = default;

        static void save_pipeline_config(std::shared_ptr<Job> job, std::string path);
    };
}

#ifdef HEADER_ONLY
#include "Job.cpp"
#endif

#endif //GENEPISEEKER_JOB_HPP
