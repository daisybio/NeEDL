//
// Created by juli on 09.08.22.
//

#ifndef GENEPISEEKER_MULTINETWORKAGGREGATOR_HPP
#define GENEPISEEKER_MULTINETWORKAGGREGATOR_HPP

#include "Job.hpp"
#include "SequentialJob.hpp"


namespace epi {

    /*!
     * Accepts pipelines of multiple runs of local search (including network construction and seeding) and combines the result sets and their structure within each current network into a new network.
     * Every network run of the multi-network approach need to be passed as one (Sequential)Job
     */
    class MultiNetworkAggregator : public Job {
    public:
        MultiNetworkAggregator(std::vector<std::shared_ptr<Job>> job_list, std::vector<std::string> names);
        MultiNetworkAggregator() = default;

        void add(const std::shared_ptr<Job> &job, std::string name);
        void add(std::vector<std::shared_ptr<Job>> job_list_, std::vector<std::string> names_);

        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;

        size_t num_networks();

    protected:
        std::vector<std::shared_ptr<Job>> job_list;
        std::vector<std::string> names;

        // intermediate results
        std::vector<std::vector<SNPSet>> result_sets;
        std::vector<std::unordered_map<SNP_t, std::vector<SNP_t>, SNP_t::SNPHash>> adjacency_data;
    };

} // epi

#ifdef HEADER_ONLY
#include "MultiNetworkAggregator.cpp"
#endif


#endif //GENEPISEEKER_MULTINETWORKAGGREGATOR_HPP
