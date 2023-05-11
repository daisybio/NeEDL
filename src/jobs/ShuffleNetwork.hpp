//
// Created by juli on 25.07.22.
//

#ifndef GENEPISEEKER_SHUFFLENETWORK_HPP
#define GENEPISEEKER_SHUFFLENETWORK_HPP

#include "Job.hpp"

namespace epi {

    class ShuffleNetwork : public Job {
    public:
        ShuffleNetwork(const std::string &algorithm);

        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        void shuffle_topology_preserving(const std::shared_ptr<DataModel> &data, bool preserve_snp_degree);

        void shuffle_expected_degree(const std::shared_ptr<DataModel> &data, bool preserve_individual_snp_degree);

        void shuffle_expected_degree_ind2(const std::shared_ptr<DataModel> &data, size_t num_iterations);

        void shuffle_expected_degree_ind3(const std::shared_ptr<DataModel> &data);

        void analyze_randomization_step(const std::shared_ptr<DataModel> &data, epi::SNPNetwork &initial_network);

        enum {
            TOPOLOGY_PRESERVING_WITH_SNP_DEGREE = 0,
            TOPOLOGY_PRESERVING_WITHOUT_SNP_DEGREE = 1,
            EXPECTED_DEGREE_KEEP_DEGREE_DISTRIBUTION = 2,
            EXPECTED_DEGREE_KEEP_INDIVIDUAL_DEGREE = 3
        } shuffle_algorithm;
        std::array<std::string, 4> shuffle_algorithm_names = {
            "TOPOLOGY_PRESERVING_WITH_SNP_DEGREE",
            "TOPOLOGY_PRESERVING_WITHOUT_SNP_DEGREE",
            "EXPECTED_DEGREE_KEEP_DEGREE_DISTRIBUTION",
            "EXPECTED_DEGREE_KEEP_INDIVIDUAL_DEGREE"
        };

        const size_t expected_degree_max_retries = 1;
    };
} // epi

#ifdef HEADER_ONLY
#include "ShuffleNetwork.cpp"
#endif


#endif //GENEPISEEKER_SHUFFLENETWORK_HPP
