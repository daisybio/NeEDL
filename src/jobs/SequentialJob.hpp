//
// Created by juli on 25.05.22.
//

#ifndef GENEPISEEKER_SEQUENTIALJOB_HPP
#define GENEPISEEKER_SEQUENTIALJOB_HPP

#include <vector>
#include <memory>

#include "Job.hpp"

/**
 * @brief this class runs a list of given jobs in sequence
 */
 namespace epi {
     class SequentialJob : public Job {
     public:
         SequentialJob() = default;
         SequentialJob(std::vector<std::shared_ptr<Job>> job_list);

         void add(const std::shared_ptr<Job>& job);
         void add(std::vector<std::shared_ptr<Job>> job_list_);

         void run(std::shared_ptr<DataModel> data) override;
         rapidjson::Value getConfig(rapidjson::Document &doc) override;

     protected:
         std::vector<std::shared_ptr<Job>> job_list;
     };
 }

#ifdef HEADER_ONLY
#include "SequentialJob.cpp"
#endif

#endif //GENEPISEEKER_SEQUENTIALJOB_HPP
