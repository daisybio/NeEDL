//
// Created by juli on 22.06.22.
//

#ifndef GENEPISEEKER_SEEDINGQUANTUMCOMPUTING_HPP
#define GENEPISEEKER_SEEDINGQUANTUMCOMPUTING_HPP

#include "SeedingCommunityWise.hpp"

namespace epi {

    class SeedingQuantumComputing : public SeedingCommunityWise {
    public:

        explicit SeedingQuantumComputing(const std::string& qc_mode,
                                const std::string & scoring_model = "PENETRANCE",
                                double quantile = .25,
                                 unsigned long max_cluster_size = 1000,
                                 unsigned long num_sets_per_cluster = 5,
                                 unsigned long num_snps_per_set = 2,
                                unsigned long min_qc_cluster_size = 100,
                                int qc_n_clique = 2,
                                int qc_k = 3,
                                double qc_nu = 0.2,
                                double qc_lambda0 = 5.,
                                double qc_lambda1 = 1.,
                                double qc_lambda2 = 1.
                                        );

        void set_simulated_annealing_params(int num_samples = 10,
                                            int num_sweeps = 1000,
                                            uint64_t seed = 1234);

        void set_quantum_annealing_params(std::string token = "",
                                          int num_reads = 10,
                                          int solver_idx = 4,
                                          double fw_annealing_ramp_time = 1.,
                                          double fw_annealing_pause_time = 1.,
                                          double rev_annealing_ramp_time = 1.,
                                          double rev_annealing_pause_time = 1.,
                                          double rev_annealing_s_target = 1.);

        void set_qaoa_params(std::string vendor = "",
                             std::string azure_subscription_id = "",
                             std::string azure_resource_group = "",
                             std::string azure_name = "",
                             std::string azure_location = "",
                             std::string azure_backend = "",
                             std::string optimizer = "",
                             int maxiter = 10,
                             int reps = 5,
                             int n_shots = 5,
                             int is_recursive_qaoa = 0);

        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        std::vector<std::vector<SNPSet>> generate_qc_sets(const std::shared_ptr<DataModel> & data, const std::vector<std::vector<SNP_t>> & clusters);

        // ATTENTION: only set defaults in above functions!

        // quantum computing
        unsigned long min_qc_cluster_size;
        int qc_n_clique;
        int qc_k;
        double qc_nu, qc_lambda0, qc_lambda1, qc_lambda2;

        enum {
            SIMULATED_ANNEALING = 0,
            QUANTUM_ANNEALING = 1,
            QAOA = 2
        } QCMode = SIMULATED_ANNEALING;
        std::array<std::string, 3> QCMode_names = { "SIMULATED_ANNEALING", "QUANTUM_ANNEALING", "QAOA" };

        // qc - simulated annealing
        uint64_t qc_sa_seed;
        int qc_sa_num_samples;
        int qc_sa_num_sweeps;

        // qc - quantum annealing
        std::string qc_qa_token;
        int qc_qa_num_reads;
        int qc_qa_solver_idx;
        double qc_qa_fw_annealing_ramp_time;
        double qc_qa_fw_annealing_pause_time;
        double qc_qa_rev_annealing_ramp_time;
        double qc_qa_rev_annealing_pause_time;
        double qc_qa_rev_annealing_s_target;

        // qc - qaoa
        std::string qc_qaoa_vendor;
        std::string qc_qaoa_azure_subscription_id;
        std::string qc_qaoa_azure_resource_group;
        std::string qc_qaoa_azure_name;
        std::string qc_qaoa_azure_location;
        std::string qc_qaoa_azure_backend;
        std::string qc_qaoa_optimizer;
        int qc_qaoa_maxiter;
        int qc_qaoa_reps;
        int qc_qaoa_n_shots;
        int qc_qaoa_is_recursive_qaoa;



    };

} // epi

#ifdef HEADER_ONLY
#include "SeedingQuantumComputing.cpp"
#endif


#endif //GENEPISEEKER_SEEDINGQUANTUMCOMPUTING_HPP
