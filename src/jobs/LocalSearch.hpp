//
// Created by juli on 06.07.22.
//

#ifndef GENEPISEEKER_LOCALSEARCH_HPP
#define GENEPISEEKER_LOCALSEARCH_HPP

#include "Job.hpp"

namespace epi {

    class LocalSearch : public Job {
    public:
        explicit LocalSearch(
               std::string epistasis_model = "PENETRANCE",
               bool collapse_identical_results = true,
               size_t max_rounds = 300,
               double search_time_limit_minutes = 0.,
               double per_seed_time_limit_minutes = 0.,
               std::string annealing_type = "SIMULATED_ANNEALING",
               double cooling_factor = 1.,
               double annealing_start_prob = .8,
               double annealing_end_prob = .01,
               bool output_score_development = true,
               std::string score_development_file_name = "score_over_time",
               size_t min_set_size = 2,
               size_t max_set_size = 10,
               bool calculate_monte_carlo_pvalues = false,
               size_t monte_carlo_permutations = 10000
            );

        void activate_LD_check(
                std::string ld_file,
                std::string ld_mode,
                double cutoff
                );

        void activate_LD_check(
                std::string ld_file,
                std::string ld_mode,
                size_t mc_min_set,
                size_t mc_max_set,
                size_t mc_sample_size
                );

        void run(std::shared_ptr<DataModel> data) override;

        std::string get_model_name();

        rapidjson::Value getConfig(rapidjson::Document &doc) override;

    private:
        SNPSet process_start_seed(const SNPSet& start_seed, const std::shared_ptr<DataModel> &data);
        bool get_annealing_decision(
                const std::shared_ptr<DataModel> &data,
                int rounds,
                double score_now,
                double score_before,
                double delta_sum,
                double temperature,
                int iterations_without_improvement
            );


        double annealing_start_prob;
        double annealing_end_prob;
        double cooling_factor;

        size_t max_rounds;
        size_t min_set;
        size_t max_set;
        double search_time_limit;
        double time_per_start_seed;

        std::chrono::high_resolution_clock::time_point search_start_time;

        bool calculate_monte_carlo;
        size_t monte_carlo_num_permutations;

        bool collapse_identical_results;

        bool output_scores_over_time;
        std::string scores_over_time_name;

        options::EpistasisScore model;

        enum {
            RANDOM_ANNEALING,
            HYPERBOLIC_TAN_ANNEALING,
            SIMULATED_ANNEALING
        } annealing_type;
        std::string annealing_type_str;


        // LD check
        enum {
            DISABLED = 0,
            CUTOFF_MODE = 1,
            MC_MODE = 2
        } ld_check = DISABLED;
        std::array<std::string, 3> ld_check_str = { "DISABLED", "CUTOFF_MODE", "MC_MODE" };

        std::shared_ptr<LDTester> ld_tester = nullptr;
        std::string ld_file, ld_mode;
        double ld_cutoff;
        size_t ld_mc_min_set, ld_mc_max_set, ld_mc_sample_size;
    };

} // epi

#ifdef HEADER_ONLY
#include "LocalSearch.cpp"
#endif


#endif //GENEPISEEKER_LOCALSEARCH_HPP
