//
// Created by juli on 24.08.22.
//

#ifndef GENEPISEEKER_SEEDINGCOMMUNITYWISE_HPP
#define GENEPISEEKER_SEEDINGCOMMUNITYWISE_HPP

#include "Job.hpp"

namespace epi {

    class SeedingCommunityWise : public Job {
    public:
        explicit SeedingCommunityWise(
                const std::string & scoring_model = "PENETRANCE",
                double quantile = .25,
                unsigned long max_cluster_size = 1000,
                unsigned long num_sets_per_cluster = 5,
                unsigned long num_snps_per_set = 2
                        );

        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;

    protected:

        std::vector<std::vector<SNP_t>> leiden_with_size_constraint(const std::shared_ptr<SNPNetwork>& snpNetwork, unsigned long maximal_allowed_cluster_size, bool optimize_min, unsigned long &min_cluster_size, unsigned long &max_cluster_size);
        void refine_leiden_clustering(const std::shared_ptr<DataModel>& data, std::vector<std::vector<SNP_t>> & clusters, size_t maximal_allowed_cluster_size);
        void set_clusters_as_attributes(const std::shared_ptr<DataModel>& data, std::string attribute_name, std::vector<std::vector<SNP_t>> clusters);
        std::vector<std::vector<SNPSet>> generate_random_sets(const std::shared_ptr<DataModel> & data, const std::vector<std::vector<SNP_t>> & clusters);
        std::vector<SNPSet> select_start_seeds (std::vector<std::vector<SNPSet>> qc_sets);

        options::EpistasisScore epi_model;
        double selection_quantile;
        unsigned long max_cluster_size;
        unsigned long num_sets_per_cluster;
        unsigned long num_snps_per_set;


        // leiden algorithm
        const float leiden_forward_search_speed = .5f;
        const size_t leiden_num_binary_search_steps = 4;
        const float leiden_beta = 0.01f;
        const size_t leiden_max_leiden_steps = 4;
    };

} // epi

#ifdef HEADER_ONLY
#include "SeedingCommunityWise.cpp"
#endif


#endif //GENEPISEEKER_SEEDINGCOMMUNITYWISE_HPP
