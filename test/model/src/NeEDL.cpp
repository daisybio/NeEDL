//
// Created by juli on 25.05.22.
//

#ifdef CMAKE_RELEASE
#pragma message ("CMAKE_RELEASE set")
    #define HEADER_ONLY
#endif

#include <CLI11.hpp>
#include <omp.h>

#include "../../src/pipelines/NeEDLPipeline.hpp"
#include "../../src/util/helper_functions.hpp"
#include "../../../src/jobs/SeedingRandomConnected.hpp"
#include "../../../src/jobs/SeedingCommunityWise.hpp"
#include "../../../src/jobs/SeedingQuantumComputing.hpp"






int main(int argc, char** argv) {
    Logger::logCmdArgs(argc, argv);
    // Function to check that argument is an integer greater than 1.
    std::function<std::string(std::string &)> check_integer_greater_1 = [](std::string &arg) -> std::string {
        int arg_as_int;
        bool failed{false};
        try {
            arg_as_int = std::stoi(arg);
        }
        catch (...) {
            failed = true;
        }
        if (not failed) {
            failed = (arg_as_int < 2);
        }
        if (failed) {
            return "arguments must be integers greater than 1";
        }
        return "";
    };


    CLI::App app("NeEDL - Scalable network-based epistasis detection via local search");


    // environmental stuff
    std::string output_directory = "";
    app.add_option("--output-directory", output_directory,
                   "The output directory. If no directory is set, no output files will be created.");

    int num_threads = 0;
    app.add_option("--num-threads", num_threads,
                   "The number of threads. 0 means all available cores. DEFAULT: use all available cores (" +
                   std::to_string(omp_get_num_procs()) + " on this machine)")->check(CLI::NonNegativeNumber);

    std::string data_directory = "./data/";
    app.add_option("--data-directory", data_directory,
                   "Path where additional data like dbSNP and BIOGRID is located. DEFAULT: ./data/");


    // instance related
    std::string input_path;
    app.add_option("--input-path", input_path,
                   "Either a directory with exactly one input file ending with json/csv or a full path to the input file.")->required();

    std::string input_format;
    app.add_option("--input-format", input_format, "The format of the input file.")->required();

    std::string input_phenotype;
    app.add_option("--phenotype", input_phenotype, "The type of the phenotypes.")->required();

    size_t input_num_categories = 2;
    app.add_option("--num-categories", input_num_categories,
                   "The number of categories for categorical phenotypes. DEFAULT: 2.")->check(
            CLI::Validator(check_integer_greater_1, "GREATER1"));

    std::string covariates_file;
    app.add_option("--covariates-file", covariates_file, "CSV file containing covariate data.");


    /*
    app.add_option("--input-SNP-format", options.input_SNP_format,
                    "The format of the SNPs, either rsIDs (RSID) or position based (POSITION_BASED)")->required();
    */


    // network statistics
    bool do_joint_degree_analysis = false;
    app.add_flag("--do-joint-degree-analysis", do_joint_degree_analysis, "If set, the individual and joint degree distribution is saved to a file.");

    bool calculate_advanced_network_statistics = false;
    app.add_flag("--calculate-advanced-network-statistics", calculate_advanced_network_statistics, "If set, metrics like the network diameter are calculated that are more computationally intensive than other metrics. This option should not be used on large networks.");

    bool disable_save_network = false;
    app.add_flag("--disable-save-network", disable_save_network, "if set, no network structure will be saved into the output directory. This can save time and disk space.");

    bool no_additional_scores = false;
    app.add_flag("--no-additional-scores", no_additional_scores, "If set, only the score for the model used in the local search is saved to the output files. This flag can be used to speed up the writing of SNP set and individual SNP scores if the other scores are not needed.");

    // network shuffling
    std::string network_shuffle_method = "";
    app.add_option("--network-shuffle-method", network_shuffle_method, "This option allows to shuffle the network before local search to evaluate the advantage gained through the particular network. Options are TOPOLOGY_PRESERVING_WITH_SNP_DEGREE, TOPOLOGY_PRESERVING_WITHOUT_SNP_DEGREE, EXPECTED_DEGREE_KEEP_DEGREE_DISTRIBUTION, EXPECTED_DEGREE_KEEP_INDIVIDUAL_DEGREE.");

    // filters
    double mma_cutoff = -1;
    app.add_option("--mma-filter", mma_cutoff, "Set the cutoff of the MMA filter. If not specified MMA filter will not be applied.");
    bool mma_use_bh_correction = false;
    app.add_flag("--mma-use-BH-correction", mma_use_bh_correction,
                 "Use Benjamini-Hochberg correction for the maximum marginal association filter p-values. The correction is done according to the R method p.adjust(method = 'BH').");

    // snp annotation methods
    bool use_dbSNP = false;
    app.add_flag("--snp-annotate-dbSNP", use_dbSNP, "This option annotates the SNPs with the help of the included dbSNP database. This results in a gene annotation.");

    std::vector<std::string> snp_annotation_methods;
    app.add_option("--snp-annotate", snp_annotation_methods, "Declare annotation sources with this option. Option can be used multiple times to add multiple sources. The file must be a csv file with one column for the snp names (as in the input file) and one column for the annotations. Both columns can be multi-columns that contain multiple entries per row separated with a separating character. If a column is not a multi-column, give -1 as separator. If both column descriptors can be casted to a number this is used as a null-based index of the column Otherwise, the two values are searched for in the header line. The format of this option is <path>|<has-header ? 'yes' : 'no'>|<snp-col>|<annotation-col>|[<csv-separator>]|[<snp-separator>]|[<annotation-separator>].");


    // multi-network networks
    bool use_BIOGRID = false;
    app.add_flag("--network-BIOGRID", use_BIOGRID, "Use the supplied BIOGRID as a network for local search.");

    std::vector<std::string> networks;
    app.add_option("--network", networks, "With this option a network can be added. Multiple networks are possible. The network file needs to be a csv file where one row represents an edge between two SNPs based on their annotation. Both columns can be multi-columns that contain more than one entry separated by a separating character. If a column is not a multi-column, give -1 as separator. If both column descriptors can be casted to a number this is used as a null-based index of the column Otherwise, the two values are searched for in the header line. The name is used for logging purposes and should not contain any special characters. Format: <name>|<path>|<has-header ? 'yes' : 'no'>|<col1>|<col1>|[<csv-separator>]|[<col1-separator>]|[<col2-separator>]");


    // multi-network local search
    std::vector<std::string> ms_model;
    app.add_option("--ms-model", ms_model,
                   "The model used as the score function in the final local search. Options are: BAYESIAN, VARIANCE, PENETRANCE, REGRESSION")->required();

    std::vector<unsigned long> ms_max_rounds;
    app.add_option("--ms-max-rounds", ms_max_rounds,
                   "The size of maximal internal rounds of local search. DEFAULT: 300.");

    std::vector<unsigned long> ms_min_set_size;
    app.add_option("--ms-min-set", ms_min_set_size, "The size of minimal disease SNPs in set for the final local search. DEFAULT: 2.");

    std::vector<unsigned long> ms_max_set_size;
    app.add_option("--ms-max-set", ms_max_set_size, "The size of maximal disease SNPs in set for the final local search. DEFAULT: 10.");

    std::vector<std::string> ms_per_seed_time_limit_str;
    app.add_option("--ms-per-seed-time-limit", ms_per_seed_time_limit_str,
                   "If this option is set, the local search is stopped for a particular seed if processing that seed takes longer than the given max time. Use m/h/d to give a number in minutes, hours or days (eg. '5d' = 5 days). Default: no time limit");

    std::vector<std::string> ms_search_time_limit_str;
    app.add_option("--ms-search-time-limit", ms_search_time_limit_str,
                   "If this option is set, the local search is stopped if it takes longer than the given max time. Use m/h/d to give a number in minutes, hours or days (eg. '5d' = 5 days). Default: no time limit");

    std::vector<std::string> ms_annealing_type;
    app.add_option("--ms-annealing-type", ms_annealing_type,
                   "Used annealing type for the final local search: RANDOM_ANNEALING, HYPERBOLIC_TAN_ANNEALING, SIMULATED_ANNEALING. DEFAULT: SIMULATED_ANNEALING.");

    std::vector<double> ms_cooling_factor;
    app.add_option("--ms-cooling-factor", ms_cooling_factor, "Cooling factor for SIMULATED_ANNEALING. DEFAULT: 1.0");

    std::vector<double> ms_annealing_start_prob;
    app.add_option("--ms-annealing-start-prob", ms_annealing_start_prob,
                   "Start probability for SIMULATED_ANNEALING. DEFAULT: 0.8");
    std::vector<double> ms_annealing_end_prob;
    app.add_option("--ms-annealing-end-prob", ms_annealing_end_prob,
                   "End probability for SIMULATED_ANNEALING. DEFAULT: 0.01");

    bool ms_ld_check = false;
    app.add_flag("--ms-ld-check", ms_ld_check, "Activates linkage disequilibrium checking during local search");

    std::string ms_ld_matrix = "";
    app.add_option("--ms-ld-matrix", ms_ld_matrix, "Path to an LD matrix in CSV format that has the same order of SNPs as the input file");

    std::string ms_ld_mode = "MEAN";
    app.add_option("--ms-ld-mode", ms_ld_mode,"The LD mode specifies how individual LD values are aggregated. OPTIONS are  MEAN, MAX. DEFAULT: MEAN");

    double ms_ld_cutoff = -std::numeric_limits<double>::infinity();
    app.add_option("--ms-ld-cutoff", ms_ld_cutoff, "if this LD cutoff is not set, monte carlo sampling is performed");

    size_t ms_ld_mc_min_set = 2;
    app.add_option("--ms-ld-mc-min-set", ms_ld_mc_min_set, "If no LD cutoff is set, monte carlo sampling is performed which samples random SNP sets with given min size. DEFAULT: 2");

    size_t ms_ld_mc_max_set = 10;
    app.add_option("--ms-ld-mc-max-set", ms_ld_mc_max_set, "If no LD cutoff is set, monte carlo sampling is performed which samples random SNP sets with given max size. DEFAULT: 10");

    size_t ms_ld_mc_num_samples = 1000;
    app.add_option("--ms-ld-mc-num-samples", ms_ld_mc_max_set, "If no LD cutoff is set, monte carlo sampling is performed. This options specifies how many samples are generated during the monte carlo process DEFAULT: 1000");


    // multi-network local search seeding
    std::vector<std::string> ms_seeding_routine;
    app.add_option("--ms-seeding-routine", ms_seeding_routine,
                   "Local search start seed type: RANDOM_CONNECTED, COMMUNITY_WISE, QUANTUM_COMPUTING. DEFAULT: RANDOM_CONNECTED.");

    std::vector<long> ms_rc_start_seeds;
    app.add_option("--ms-rc-start-seeds", ms_rc_start_seeds,
                   "Number of start seeds when using seeding routine RANDOM_CONNECTED for the final local search. DEFAULT: 300.")->check(
            CLI::Validator(check_integer_greater_1, "GREATER1"));

    std::vector<std::string> ms_qc_mode;
    app.add_option("--ms-qc-mode", ms_qc_mode,
                   "Mode for quantum computing seeding of the final local search. Options are: SIMULATED_ANNEALING, QUANTUM_ANNEALING, QAOA. DEFAULT: SIMULATED_ANNEALING");

    std::vector<double> ms_cw_quantile;
    app.add_option("--ms-cw-quantile", ms_cw_quantile, "Selects how many of the best start seeds are selected for COMMUNITY_WISE or QUANTUM_COMPUTING for the final local search. DEFAULT: 0.25");

    std::vector<unsigned long> ms_cw_max_cluster_size;
    app.add_option("--ms-cw-max-cluster-size", ms_cw_max_cluster_size, "Maximum cluster size that is allowed for the seeding COMMUNITY_WISE or QUANTUM_COMPUTING for the final local search. DEFAULT: 1000");

    std::vector<unsigned long> ms_cw_num_sets_per_cluster;
    app.add_option("--ms-cw-num-sets-per-cluster", ms_cw_num_sets_per_cluster, "How many start seeds candidates are cosnidered for every cluster for the seeding COMMUNITY_WISE or QUANTUM_COMPUTING for the final local search. DEFAULT: 5");

    std::vector<unsigned long> ms_cw_num_snps_per_set;
    app.add_option("--ms-cw-num-snps-per-set", ms_cw_num_snps_per_set, "How many SNPs each start seed should contain if it was not generated with quantum computing for the seeding COMMUNITY_WISE or QUANTUM_COMPUTING for the final local search. DEFAULT: 2");

    std::vector<unsigned long> ms_qc_min_cluster_size;
    app.add_option("--ms-qc-min-cluster-size", ms_qc_min_cluster_size,"Quantum computing is only used to find seeds for large enough clusters. Seeds for smaller clusters are picked randomly. This parameter defines the threshold at what size to pick quantum computing. DEFAULT 100");

    std::vector<int> ms_qc_n_clique;
    app.add_option("--ms-qc-n-clique", ms_qc_n_clique,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 2");

    std::vector<int> ms_qc_k;
    app.add_option("--ms-qc-k", ms_qc_k,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 3");

    std::vector<double> ms_qc_nu;
    app.add_option("--ms-qc-nu", ms_qc_nu,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 3.0");

    std::vector<double> ms_qc_lambda0;
    app.add_option("--ms-qc-lambda0", ms_qc_lambda0,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 5.0");

    std::vector<double> ms_qc_lambda1;
    app.add_option("--ms-qc-lambda1", ms_qc_lambda1,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 1.0");

    std::vector<double> ms_qc_lambda2;
    app.add_option("--ms-qc-lambda2", ms_qc_lambda2,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 1.0");

    std::vector<int> ms_qc_sa_num_samples;
    app.add_option("--ms-qc-sa-num-samples", ms_qc_sa_num_samples,"Parameter for seeding QUANTUM_COMPUTING, mode SIMULATED_ANNEALING, DEFAULT 10");

    std::vector<int> ms_qc_sa_num_sweeps;
    app.add_option("--ms-qc-sa-num-sweeps", ms_qc_sa_num_sweeps,"Parameter for seeding QUANTUM_COMPUTING, mode SIMULATED_ANNEALING, DEFAULT 1000");

    std::vector<uint64_t> ms_qc_sa_seed;
    app.add_option("--ms-qc-sa-seed", ms_qc_sa_seed,"Parameter for seeding QUANTUM_COMPUTING, mode SIMULATED_ANNEALING, DEFAULT 1234");

    std::vector<std::string> ms_qc_qa_token;
    app.add_option("--ms-qc-qa-token", ms_qc_qa_token,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT ''");

    std::vector<int> ms_qc_qa_num_reads;
    app.add_option("--ms-qc-qa-num-reads", ms_qc_qa_num_reads,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 10");

    std::vector<int> ms_qc_qa_solver_idx;
    app.add_option("--ms-qc-qa-solver-idx", ms_qc_qa_solver_idx,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 4");

    std::vector<double> ms_qc_qa_fw_annealing_ramp_time;
    app.add_option("--ms-qc-qa-fw-annealing-ramp-time", ms_qc_qa_fw_annealing_ramp_time,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 1.");

    std::vector<double> ms_qc_qa_fw_annealing_pause_time;
    app.add_option("--ms-qc-qa-fw-annealing-pause-time", ms_qc_qa_fw_annealing_pause_time,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 1.");

    std::vector<double> ms_qc_qa_rev_annealing_ramp_time;
    app.add_option("--ms-qc-qa-rev-annealing-ramp-time", ms_qc_qa_rev_annealing_ramp_time,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 1.");

    std::vector<double> ms_qc_qa_rev_annealing_pause_time;
    app.add_option("--ms-qc-qa-rev-annealing-pause-time", ms_qc_qa_rev_annealing_pause_time,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 1.");

    std::vector<double> ms_qc_qa_rev_annealing_s_target;
    app.add_option("--ms-qc-qa-rev-annealing-s-target", ms_qc_qa_rev_annealing_s_target,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 1.");

    std::vector<std::string> ms_qc_qaoa_vendor;
    app.add_option("--ms-qc-qaoa-vendor", ms_qc_qaoa_vendor,"Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::vector<std::string> ms_qc_qaoa_azure_subscription_id;
    app.add_option("--ms-qc-qaoa-azure-subscription-id", ms_qc_qaoa_azure_subscription_id,"Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::vector<std::string> ms_qc_qaoa_azure_resource_group;
    app.add_option("--ms-qc-qaoa-azure-resource-group", ms_qc_qaoa_azure_resource_group, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::vector<std::string> ms_qc_qaoa_azure_name;
    app.add_option("--ms-qc-qaoa-azure-name", ms_qc_qaoa_azure_name, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::vector<std::string> ms_qc_qaoa_azure_location;
    app.add_option("--ms-qc-qaoa-azure-location", ms_qc_qaoa_azure_location, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::vector<std::string> ms_qc_qaoa_azure_backend;
    app.add_option("--ms-qc-qaoa-azure-backend", ms_qc_qaoa_azure_backend, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::vector<std::string> ms_qc_qaoa_optimizer;
    app.add_option("--ms-qc-qaoa-optimizer", ms_qc_qaoa_optimizer, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::vector<int> ms_qc_qaoa_maxiter;
    app.add_option("--ms-qc-qaoa-maxiter", ms_qc_qaoa_maxiter, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT 10.");

    std::vector<int> ms_qc_qaoa_reps;
    app.add_option("--ms-qc-qaoa-reps", ms_qc_qaoa_reps, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT 5.");

    std::vector<int> ms_qc_qaoa_n_shots;
    app.add_option("--ms-qc-qaoa-n-shots", ms_qc_qaoa_n_shots, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT 5.");

    std::vector<int> ms_qc_qaoa_is_recursive_qaoa;
    app.add_option("--ms-qc-qaoa-is-recursive-qaoa", ms_qc_qaoa_is_recursive_qaoa, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT 0.");

    // final local search
    bool disable_final_search = false;
    app.add_flag("--disable-final-search", disable_final_search, "Disable the final search of the network. This flag can be used if only multiple networks are used but the results should not be aggregated.");

    std::string fs_model;
    app.add_option("--fs-model", fs_model,
                   "The model used as the score function in the final local search. Options are: BAYESIAN, VARIANCE, PENETRANCE, REGRESSION");

    unsigned long fs_max_rounds = 300;
    app.add_option("--fs-max-rounds", fs_max_rounds,
                   "The size of maximal internal rounds of local search. DEFAULT: 300.");

    unsigned long fs_min_set_size = 2;
    app.add_option("--fs-min-set", fs_min_set_size, "The size of minimal disease SNPs in set for the final local search. DEFAULT: 2.");

    unsigned long fs_max_set_size = 10;
    app.add_option("--fs-max-set", fs_max_set_size, "The size of maximal disease SNPs in set for the final local search. DEFAULT: 10.");

    std::string fs_per_seed_time_limit_str;
    app.add_option("--fs-per-seed-time-limit", fs_per_seed_time_limit_str,
                   "If this option is set, the local search is stopped for a particular seed if processing that seed takes longer than the given max time. Use m/h/d to give a number in minutes, hours or days (eg. '5d' = 5 days). Default: no time limit");

    std::string fs_search_time_limit_str;
    app.add_option("--fs-search-time-limit", fs_search_time_limit_str,
                   "If this option is set, the local search is stopped if it takes longer than the given max time. Use m/h/d to give a number in minutes, hours or days (eg. '5d' = 5 days). Default: no time limit");

    std::string fs_annealing_type = "SIMULATED_ANNEALING";
    app.add_option("--fs-annealing-type", fs_annealing_type,
                   "Used annealing type for the final local search: RANDOM_ANNEALING, HYPERBOLIC_TAN_ANNEALING, SIMULATED_ANNEALING. DEFAULT: SIMULATED_ANNEALING.");

    double fs_cooling_factor = 1.0;
    app.add_option("--fs-cooling-factor", fs_cooling_factor, "Cooling factor for SIMULATED_ANNEALING. DEFAULT: 1.0");

    double fs_annealing_start_prob = 0.8;
    app.add_option("--fs-annealing-start-prob", fs_annealing_start_prob,
                   "Start probability for SIMULATED_ANNEALING. DEFAULT: 0.8");
    double fs_annealing_end_prob = 0.01;
    app.add_option("--fs-annealing-end-prob", fs_annealing_end_prob,
                   "End probability for SIMULATED_ANNEALING. DEFAULT: 0.01");

    bool fs_ld_check = false;
    app.add_flag("--fs-ld-check", fs_ld_check, "Activates linkage disequilibrium checking during local search");

    std::string fs_ld_matrix = "";
    app.add_option("--fs-ld-matrix", fs_ld_matrix, "Path to an LD matrix in CSV format that has the same order of SNPs as the input file");

    std::string fs_ld_mode = "MEAN";
    app.add_option("--fs-ld-mode", fs_ld_mode,"The LD mode specifies how individual LD values are aggregated. OPTIONS are  MEAN, MAX. DEFAULT: MEAN");

    double fs_ld_cutoff = -std::numeric_limits<double>::infinity();
    app.add_option("--fs-ld-cutoff", fs_ld_cutoff, "if this LD cutoff is not set, monte carlo sampling is performed");

    size_t fs_ld_mc_min_set = 2;
    app.add_option("--fs-ld-mc-min-set", fs_ld_mc_min_set, "If no LD cutoff is set, monte carlo sampling is performed which samples random SNP sets with given min size. DEFAULT: 2");

    size_t fs_ld_mc_max_set = 10;
    app.add_option("--fs-ld-mc-max-set", fs_ld_mc_max_set, "If no LD cutoff is set, monte carlo sampling is performed which samples random SNP sets with given max size. DEFAULT: 10");

    size_t fs_ld_mc_num_samples = 1000;
    app.add_option("--fs-ld-mc-num-samples", fs_ld_mc_max_set, "If no LD cutoff is set, monte carlo sampling is performed. This options specifies how many samples are generated during the monte carlo process DEFAULT: 1000");


    // final local search seeding
    std::string fs_seeding_routine = "RANDOM_CONNECTED";
    app.add_option("--fs-seeding-routine", fs_seeding_routine,
                   "Local search start seed type: RANDOM_CONNECTED, COMMUNITY_WISE, QUANTUM_COMPUTING. DEFAULT: RANDOM_CONNECTED.");

    long fs_rc_start_seeds = 300;
    app.add_option("--fs-rc-start-seeds", fs_rc_start_seeds,
                   "Number of start seeds when using seeding routine RANDOM_CONNECTED for the final local search. DEFAULT: 300.")->check(
            CLI::Validator(check_integer_greater_1, "GREATER1"));

    std::string fs_qc_mode = "SIMULATED_ANNEALING";
    app.add_option("--fs-qc-mode", fs_qc_mode,
                   "Mode for quantum computing seeding of the final local search. Options are: SIMULATED_ANNEALING, QUANTUM_ANNEALING, QAOA. DEFAULT: SIMULATED_ANNEALING");

    double fs_cw_quantile = .25;
    app.add_option("--fs-cw-quantile", fs_cw_quantile, "Selects how many of the best start seeds are selected for COMMUNITY_WISE or QUANTUM_COMPUTING for the final local search. DEFAULT: 0.25");

    unsigned long fs_cw_max_cluster_size = 1000;
    app.add_option("--fs-cw-max-cluster-size", fs_cw_max_cluster_size, "Maximum cluster size that is allowed for the seeding COMMUNITY_WISE or QUANTUM_COMPUTING for the final local search. DEFAULT: 1000");

    unsigned long fs_cw_num_sets_per_cluster= 5;
    app.add_option("--fs-cw-num-sets-per-cluster", fs_cw_num_sets_per_cluster, "How many start seeds candidates are cosnidered for every cluster for the seeding COMMUNITY_WISE or QUANTUM_COMPUTING for the final local search. DEFAULT: 5");

    unsigned long fs_cw_num_snps_per_set = 2;
    app.add_option("--fs-cw-num-snps-per-set", fs_cw_num_snps_per_set, "How many SNPs each start seed should contain if it was not generated with quantum computing for the seeding COMMUNITY_WISE or QUANTUM_COMPUTING for the final local search. DEFAULT: 2");

    unsigned long fs_qc_min_cluster_size = 100;
    app.add_option("--fs-qc-min-cluster-size", fs_qc_min_cluster_size,"Quantum computing is only used to find seeds for large enough clusters. Seeds for smaller clusters are picked randomly. This parameter defines the threshold at what size to pick quantum computing. DEFAULT 100");

    int fs_qc_n_clique = 2;
    app.add_option("--fs-qc-n-clique", fs_qc_n_clique,"Parameter for seeding QUANTUM_COMPUTING. DEFAULT 2");

    int fs_qc_k = 3;
    app.add_option("--fs-qc-k", fs_qc_k,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 3");

    double fs_qc_nu = 3.;
    app.add_option("--fs-qc-nu", fs_qc_nu,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 3.0");

    double fs_qc_lambda0 = 5.;
    app.add_option("--fs-qc-lambda0", fs_qc_lambda0,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 5.0");

    double fs_qc_lambda1 = 1.;
    app.add_option("--fs-qc-lambda1", fs_qc_lambda1,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 1.0");

    double fs_qc_lambda2 = 1.;
    app.add_option("--fs-qc-lambda2", fs_qc_lambda2,"Parameter for seeding QUANTUM_COMPUTING, DEFAULT 1.0");

    int fs_qc_sa_num_samples = 10;
    app.add_option("--fs-qc-sa-num-samples", fs_qc_sa_num_samples,"Parameter for seeding QUANTUM_COMPUTING, mode SIMULATED_ANNEALING, DEFAULT 10");

    int fs_qc_sa_num_sweeps = 1000;
    app.add_option("--fs-qc-sa-num-sweeps", fs_qc_sa_num_sweeps,"Parameter for seeding QUANTUM_COMPUTING, mode SIMULATED_ANNEALING, DEFAULT 1000");

    uint64_t fs_qc_sa_seed = 1234;
    app.add_option("--fs-qc-sa-seed", fs_qc_sa_seed,"Parameter for seeding QUANTUM_COMPUTING, mode SIMULATED_ANNEALING, DEFAULT 1234");

    std::string fs_qc_qa_token = "";
    app.add_option("--fs-qc-qa-token", fs_qc_qa_token,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT ''");

    int fs_qc_qa_num_reads = 10;
    app.add_option("--fs-qc-qa-num-reads", fs_qc_qa_num_reads,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 10");

    int fs_qc_qa_solver_idx = 4;
    app.add_option("--fs-qc-qa-solver-idx", fs_qc_qa_solver_idx,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 4");

    double fs_qc_qa_fw_annealing_ramp_time = 1.;
    app.add_option("--fs-qc-qa-fw-annealing-ramp-time", fs_qc_qa_fw_annealing_ramp_time,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 1.");

    double fs_qc_qa_fw_annealing_pause_time = 1.;
    app.add_option("--fs-qc-qa-fw-annealing-pause-time", fs_qc_qa_fw_annealing_pause_time,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 1.");

    double fs_qc_qa_rev_annealing_ramp_time = 1.;
    app.add_option("--fs-qc-qa-rev-annealing-ramp-time", fs_qc_qa_rev_annealing_ramp_time,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 1.");

    double fs_qc_qa_rev_annealing_pause_time = 1.;
    app.add_option("--fs-qc-qa-rev-annealing-pause-time", fs_qc_qa_rev_annealing_pause_time,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 1.");

    double fs_qc_qa_rev_annealing_s_target = 1.;
    app.add_option("--fs-qc-qa-rev-annealing-s-target", fs_qc_qa_rev_annealing_s_target,"Parameter for seeding QUANTUM_COMPUTING, mode QUANTUM_ANNEALING, DEFAULT 1.");

    std::string fs_qc_qaoa_vendor = "";
    app.add_option("--fs-qc-qaoa-vendor", fs_qc_qaoa_vendor,"Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::string fs_qc_qaoa_azure_subscription_id = "";
    app.add_option("--fs-qc-qaoa-azure-subscription-id", fs_qc_qaoa_azure_subscription_id,"Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::string fs_qc_qaoa_azure_resource_group = "";
    app.add_option("--fs-qc-qaoa-azure_resource_group", fs_qc_qaoa_azure_resource_group,"Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::string fs_qc_qaoa_azure_name = "";
    app.add_option("--fs-qc-qaoa-azure-name", fs_qc_qaoa_azure_name,"Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::string fs_qc_qaoa_azure_location = "";
    app.add_option("--fs-qc-qaoa-azure-location", fs_qc_qaoa_azure_location, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::string fs_qc_qaoa_azure_backend = "";
    app.add_option("--fs-qc-qaoa-azure-backend", fs_qc_qaoa_azure_backend, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    std::string fs_qc_qaoa_optimizer = "";
    app.add_option("--fs-qc-qaoa-optimizer", fs_qc_qaoa_optimizer, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT ''.");

    int fs_qc_qaoa_maxiter = 10;
    app.add_option("--fs-qc-qaoa-maxiter", fs_qc_qaoa_maxiter, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT 10.");

    int fs_qc_qaoa_reps = 5;
    app.add_option("--fs-qc-qaoa-reps", fs_qc_qaoa_reps, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT 5.");

    int fs_qc_qaoa_n_shots = 5;
    app.add_option("--fs-qc-qaoa-n-shots", fs_qc_qaoa_n_shots, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT 5.");

    int fs_qc_qaoa_is_recursive_qaoa = 0;
    app.add_option("--fs-qc-qaoa-is-recursive-qaoa", fs_qc_qaoa_is_recursive_qaoa, "Parameter for seeding QUANTUM_COMPUTING, mode QAOA, DEFAULT 0.");


    // Optional options.
    /*
    options.input_LD_matrix = "-1";
    app.add_option("--input-LD-matrix", options.input_LD_matrix, "if this path is set LD function is activated");
    options.cutoff_LD = 0.0;
    app.add_option("--cutoff-LD", options.cutoff_LD, "if this LD cutoff is not set, monte carlo sampling is performed");
    options.mode_LD = "MEAN";
    app.add_option("--mode-LD", options.mode_LD,
                   "sets the mode of LD, can have the values  MEAN, MAX, or SUM. DEFAULT: MEAN");

    options.input_SNP_configuration = "SNP_REGION";
    app.add_option("--input-SNP-configuration", options.input_SNP_configuration,
                   "Configures, if SNP format is position based, which position mapping should be performed. KNOWN_SNPS: only known SNPs of dbSNP, SNP_REGIONS: known SNPs of dbSNPs and regions of known SNPs, CLOSEST_GENE: known SNPs of dbSNP and closest genes, GENE_SNP_REGIONS: known SNPs of dbSNP and closest gene and regions of known SNPs. DEFAULT: SNP_REGIONS");
    options.input_chromosome_information = "-1";
    app.add_option("--input-chromosome-information", options.input_chromosome_information,
                   "Additional chromosome information for CSV input");
    options.input_mafs_information_csv = "-1";
    app.add_option("--input-maf-information", options.input_mafs_information_csv,
                   "Additional MAF information for CSV input");
    options.maf_filter = -1.0;
    app.add_option("--maf-filter", options.maf_filter, "Activate MAF filter");
    options.maximum_marginal_association_filter = -1.0;

     options.max_seconds_per_start_point = -1;
    app.add_option("--max-seconds-per-start-point", options.max_seconds_per_start_point,
                   "The size of maximal time used per start points in local search. -1 deactivates this feature. DEFAULT: -1.")->check(
            CLI::Validator(check_integer_greater_1, "GREATER1"));
    options.annealing_type = "SIMULATED_ANNEALING";
    app.add_option("--annealing-type", options.annealing_type,
                   "Used annealing type: RANDOM_ANNEALING, HYPERBOLIC_TAN_ANNEALING, SIMULATED_ANNEALING. DEFAULT: SIMULATED_ANNEALING.");
    options.cooling_factor = 1.0;
    app.add_option("--cooling-factor", options.cooling_factor, "Cooling factor for SIMULATED_ANNEALING. DEFAULT: 1.0");
    options.annealing_start_prob = 0.8;
    app.add_option("--annealing-start-prob", options.annealing_start_prob,
                   "Start probability for SIMULATED_ANNEALING. DEFAULT: 0.8");
    options.annealing_end_prob = 0.01;
    app.add_option("--annealing-end-prob", options.annealing_end_prob,
                   "End probability for SIMULATED_ANNEALING. DEFAULT: 0.01");
    app.add_option("--model-cache-size", model_cache_size,
                   "Number of model scores that are cached for each round of the local search. Every thread uses its own cache! DEFAULT: 10000");
    num_pareto_iterations = 10;
    app.add_option("--pareto-num-iterations", num_pareto_iterations,
                   "Number of iterations of pareto for the majority vote. DEFAULT 10");

    calculate_monte_carlo = false;
    app.add_flag("--calculate-monte-carlo", calculate_monte_carlo,
                 "activate or deactivate the calculation of monte carlo p-values for the results. DEFAULT: deactivated");
    monte_carlo_num_permutations = 10000;
    app.add_option("--monte-carlo-num-permutations", monte_carlo_num_permutations,
                   "number of permutations used to calculate the monte carlo p-values. DEFAULT: 10000");

    app.add_flag("--silent-mode", silentMode, "Silent mode: Deactivates all console logging");
    app.add_flag("--file-logging", fileLogging, "Write logging to a logfile.");

    app.add_option("--log-directory", log_directory,
                   "Path to the directories where logfiles should be created if the log files should not be saved in the same directory as the output files.");
    app.add_option("--best-results-per-file", bestResultsPerFile,
                   "Number of best results that should be printed. DEFAULT: 1");


    std::string global_time_limit_str;
    app.add_option("--global-time-limit", global_time_limit_str,
                   "Can be used to limit time used for local search. If local search was not completed within the limit, it will be stopped and the preliminary result will be returned. Use m/h/d to give a number in minutes, hours or days (eg. '5d' = 5 days). Default: no time limit");

    options.getGraph_flag = "";
    bool get_graph_json = false;
    bool get_graph_csv = false;
    app.add_flag("--get-graph-json", get_graph_json,
                 "if set no local search will be performed. Only adjacency list of SNP-SNP interaction graph will be written to a json file.");
    app.add_flag("--get-graph-csv", get_graph_csv,
                 "if set no local search will be performed. Only adjacency matrix of SNP-SNP interaction graph will be written to a csv file.");
    app.add_flag("--score-over-time", score_over_time,
                 "If set the best score will be stored for every minute during the local search.");
    */

    // Parse the options.
    try {
        app.parse((argc), (argv));
    } catch(const CLI::ParseError &e) {
        return app.exit(e);
    }

    try {
        NeEDLPipeline pipeline(
                input_path,
                input_format,
                input_phenotype,
                input_num_categories,
                covariates_file,
                disable_save_network,
                calculate_advanced_network_statistics,
                do_joint_degree_analysis,
                no_additional_scores,
                output_directory,
                num_threads,
                data_directory
        );

        if (mma_cutoff != -1) pipeline.activate_MMA_filter(mma_cutoff, mma_use_bh_correction);

        if (!network_shuffle_method.empty()) pipeline.activate_network_shuffling(network_shuffle_method);

        if (use_dbSNP) pipeline.add_snp_annotation_source_dbSNP();
        for (auto &method: snp_annotation_methods) {
            pipeline.add_snp_annotation_source(epi::SnpCsvAnnotator::parse_from_source_string(method));
        }

        // multi-network stuff
        size_t num_networks = networks.size();
        if (use_BIOGRID) num_networks += 1;

        // expand all parameters for thw multi-network configuration
        epi::expand_multi_param("--ms-model", ms_model, std::string(""), num_networks, true);
        epi::expand_multi_param("--ms-max-rounds", ms_max_rounds, 300UL, num_networks, false);
        std::vector<double> ms_per_seed_search_limit_minutes;
    for (auto & val : ms_per_seed_time_limit_str) ms_per_seed_search_limit_minutes.push_back(epi::parseTimespanString(val));
        epi::expand_multi_param("--ms-per-seed-time_limit", ms_per_seed_search_limit_minutes, 0., num_networks, false);
        std::vector<double> ms_search_time_limit_minutes;
        for (auto &val: ms_search_time_limit_str) ms_search_time_limit_minutes.push_back(epi::parseTimespanString(val));
        epi::expand_multi_param("--ms-search-time-limit", ms_search_time_limit_minutes, 0., num_networks, false);
    epi::expand_multi_param("--ms-annealing-type", ms_annealing_type, std::string("SIMULATED_ANNEALING"), num_networks, false);
        epi::expand_multi_param("--ms-cooling-factor", ms_cooling_factor, 1.0, num_networks, false);
        epi::expand_multi_param("--ms-annealing-start-prob", ms_annealing_start_prob, 0.8, num_networks, false);
        epi::expand_multi_param("--ms-annealing-end-prob", ms_annealing_end_prob, 0.01, num_networks, false);
        epi::expand_multi_param("--ms-min-set-size", ms_min_set_size, 2UL, num_networks, false);
        epi::expand_multi_param("--ms-max-set-size", ms_max_set_size, 10UL, num_networks, false);
    epi::expand_multi_param("--ms-seeding-routine", ms_seeding_routine, std::string("RANDOM_CONNECTED"), num_networks, false);
        epi::expand_multi_param("--ms-rc-start-seeds", ms_rc_start_seeds, 300L, num_networks, false);
        epi::expand_multi_param("--ms-qc-mode", ms_qc_mode, std::string("SIMULATED_ANNEALING"), num_networks, false);
        epi::expand_multi_param("--ms-cw-quantile", ms_cw_quantile, .25, num_networks, false);
        epi::expand_multi_param("--ms-cw-max-cluster-size", ms_cw_max_cluster_size, 1000UL, num_networks, false);
        epi::expand_multi_param("--ms-cw-num-sets-per-cluster", ms_cw_num_sets_per_cluster, 5UL, num_networks, false);
        epi::expand_multi_param("--ms-cw-num-snps-per-set", ms_cw_num_snps_per_set, 2UL, num_networks, false);
        epi::expand_multi_param("--ms-qc-min-cluster-size", ms_qc_min_cluster_size, 100UL, num_networks, false);
        epi::expand_multi_param("--ms-qc-n-clique", ms_qc_n_clique, 2, num_networks, false);
        epi::expand_multi_param("--ms-qc-k", ms_qc_k, 3, num_networks, false);
        epi::expand_multi_param("--ms-qc-nu", ms_qc_nu, 3., num_networks, false);
        epi::expand_multi_param("--ms-qc-lambda0", ms_qc_lambda0, 5., num_networks, false);
        epi::expand_multi_param("--ms-qc-lambda1", ms_qc_lambda1, 1., num_networks, false);
        epi::expand_multi_param("--ms-qc-lambda2", ms_qc_lambda2, 1., num_networks, false);
        epi::expand_multi_param("--ms-qc-sa-num-samples", ms_qc_sa_num_samples, 10, num_networks, false);
        epi::expand_multi_param("--ms-qc-sa-num-sweeps", ms_qc_sa_num_sweeps, 1000, num_networks, false);
        epi::expand_multi_param("--ms-qc-sa-seed", ms_qc_sa_seed, uint64_t(1234), num_networks, false);
        epi::expand_multi_param("--ms-qc-qa-token", ms_qc_qa_token, std::string(""), num_networks, false);
        epi::expand_multi_param("--ms-qc-qa-num-reads", ms_qc_qa_num_reads, 10, num_networks, false);
        epi::expand_multi_param("--ms-qc-qa-solver-idx", ms_qc_qa_num_reads, 4, num_networks, false);
    epi::expand_multi_param("--ms-qc-qa-fw-annealing-ramp-time", ms_qc_qa_fw_annealing_ramp_time, 1., num_networks, false);
    epi::expand_multi_param("--ms-qc-qa-fw-annealing-pause-time", ms_qc_qa_fw_annealing_pause_time, 1., num_networks, false);
    epi::expand_multi_param("--ms-qc-qa-rev-annealing-ramp-time", ms_qc_qa_rev_annealing_ramp_time, 1., num_networks, false);
    epi::expand_multi_param("--ms-qc-qa-rev-annealing-pause-time", ms_qc_qa_rev_annealing_pause_time, 1., num_networks, false);
    epi::expand_multi_param("--ms-qc-qa-rev-annealing-s-target", ms_qc_qa_rev_annealing_s_target, 1., num_networks, false);
        epi::expand_multi_param("--ms-qc-qaoa-vendor", ms_qc_qaoa_vendor, std::string(""), num_networks, false);
    epi::expand_multi_param("--ms-qc-qaoa-azure-subscription-id", ms_qc_qaoa_azure_subscription_id, std::string(""), num_networks, false);
    epi::expand_multi_param("--ms-qc-qaoa-azure-resource-group", ms_qc_qaoa_azure_resource_group, std::string(""), num_networks, false);
        epi::expand_multi_param("--ms-qc-qaoa-azure-name", ms_qc_qaoa_azure_name, std::string(""), num_networks, false);
    epi::expand_multi_param("--ms-qc-qaoa-azure-backend", ms_qc_qaoa_azure_backend, std::string(""), num_networks, false);
        epi::expand_multi_param("--ms-qc-qaoa-optimizer", ms_qc_qaoa_optimizer, std::string(""), num_networks, false);
        epi::expand_multi_param("--ms-qc-qaoa-maxiter", ms_qc_qaoa_maxiter, 10, num_networks, false);
        epi::expand_multi_param("--ms-qc-qaoa-reps", ms_qc_qaoa_reps, 5, num_networks, false);
        epi::expand_multi_param("--ms-qc-qaoa-n-shots", ms_qc_qaoa_n_shots, 5, num_networks, false);
        epi::expand_multi_param("--ms-qc-qaoa-is-recursive-qaoa", ms_qc_qaoa_is_recursive_qaoa, 0, num_networks, false);


        auto escapedCharMap = epi::getEscapedCharMap();
        for (size_t i = 0; i < num_networks; i++) {
            bool is_biogrid = true;
            std::string name = "BIOGRID";
            epi::NetworkCsvConnector *connector = nullptr;
            if (!use_BIOGRID || i >= 1) {
                is_biogrid = false;

                auto splits = epi::string_split(networks[use_BIOGRID ? (i - 1) : i], '|');
                if (splits.size() < 5) throw epi::Error("Not enough arguments to specify a network.");

                name = splits[0];
                std::string path = splits[1];
                bool has_header = epi::toUpperCase(splits[2]) == "YES";
                std::string column1 = splits[3];
                std::string column2 = splits[4];

                bool col1_cast = false, col2_cast = false;
                size_t col1_num = epi::try_parse_number(column1, col1_cast);
                size_t col2_num = epi::try_parse_number(column2, col2_cast);

                char csv_sep = ';';
                if (splits.size() > 5 && !splits[5].empty()) {
                    if (escapedCharMap.find(splits[5]) != escapedCharMap.end()) {
                        splits[5][0] = escapedCharMap[splits[5]];
                    } else throw epi::Error("Separator can only be a single character");
                    if (!splits[5].empty()) csv_sep = splits[5][0];
                }

                char col1_sep = ';';
                if (splits.size() > 6 && !splits[6].empty()) {
                    if (escapedCharMap.find(splits[6]) != escapedCharMap.end()) {
                        splits[6][0] = escapedCharMap[splits[6]];
                    } else throw epi::Error("Separator can only be a single character");
                    if (!splits[6].empty()) col1_sep = splits[6][0];
                }

                char col2_sep = ';';
                if (splits.size() > 7 && !splits[7].empty()) {
                    if (escapedCharMap.find(splits[7]) != escapedCharMap.end()) {
                        splits[7][0] = escapedCharMap[splits[7]];
                    } else throw epi::Error("Separator can only be a single character");
                    if (!splits[7].empty()) col2_sep = splits[7][0];
                }

                if (col1_cast && col2_cast) {
                    connector = new epi::NetworkCsvConnector(name, path, col1_num, col2_num, has_header, csv_sep,
                                                             col1_sep, col2_sep);
                } else if (has_header) {
                    connector = new epi::NetworkCsvConnector(name, path, column1, column2, csv_sep, col1_sep, col2_sep);
                } else {
                throw epi::Error("Error in network source: A csv without header was given but the columns are referenced by column name.");
                }
            }

       epi::LocalSearch network_local_search(
                    ms_model[i],
                    true, // collapse identical results
                    ms_max_rounds[i],
                    ms_search_time_limit_minutes[i],
                    ms_per_seed_search_limit_minutes[i],
                    ms_annealing_type[i],
                    ms_cooling_factor[i],
                    ms_annealing_start_prob[i],
                    ms_annealing_end_prob[i],
                    true, // score development
                    name + "_search_score_over_time",
                    ms_min_set_size[i],
                    ms_max_set_size[i]
            );

            if (ms_ld_check) {
                if (ms_ld_cutoff == -std::numeric_limits<double>::infinity()) {
                    network_local_search.activate_LD_check(ms_ld_matrix, ms_ld_mode, ms_ld_mc_min_set, ms_ld_mc_max_set,
                                                           ms_ld_mc_num_samples);
                } else {
                    network_local_search.activate_LD_check(ms_ld_matrix, ms_ld_mode, ms_ld_cutoff);
                }
            }

            // parse seeding
            if (ms_seeding_routine[i] == "RANDOM_CONNECTED") {
                epi::SeedingRandomConnected seeding(ms_rc_start_seeds[i]);
                if (is_biogrid) pipeline.add_network_BIOGRID(seeding, network_local_search);
                else pipeline.add_network(*connector, seeding, network_local_search);
            } else if (ms_seeding_routine[i] == "COMMUNITY_WISE") {
                epi::SeedingCommunityWise seeding(
                        ms_model[i],
                        ms_cw_quantile[i],
                        ms_cw_max_cluster_size[i],
                        ms_cw_num_sets_per_cluster[i],
                        ms_cw_num_snps_per_set[i]
                );
                if (is_biogrid) pipeline.add_network_BIOGRID(seeding, network_local_search);
                else pipeline.add_network(*connector, seeding, network_local_search);
            } else if (ms_seeding_routine[i] == "QUANTUM_COMPUTING") {
                epi::SeedingQuantumComputing seeding(
                        ms_qc_mode[i],
                        ms_model[i],
                        ms_cw_quantile[i],
                        ms_cw_max_cluster_size[i],
                        ms_cw_num_sets_per_cluster[i],
                        ms_cw_num_snps_per_set[i],
                        ms_qc_min_cluster_size[i],
                        ms_qc_n_clique[i],
                        ms_qc_k[i],
                        ms_qc_nu[i],
                        ms_qc_lambda0[i],
                        ms_qc_lambda1[i],
                        ms_qc_lambda2[i]
                );
                if (ms_qc_mode[i] == "SIMULATED_ANNEALING") {
                    seeding.set_simulated_annealing_params(
                            ms_qc_sa_num_samples[i],
                            ms_qc_sa_num_sweeps[i],
                            ms_qc_sa_seed[i]
                    );
                } else if (ms_qc_mode[i] == "QUANTUM_ANNEALING") {
                    seeding.set_quantum_annealing_params(
                            ms_qc_qa_token[i],
                            ms_qc_qa_num_reads[i],
                            ms_qc_qa_solver_idx[i],
                            ms_qc_qa_fw_annealing_ramp_time[i],
                            ms_qc_qa_fw_annealing_pause_time[i],
                            ms_qc_qa_rev_annealing_ramp_time[i],
                            ms_qc_qa_rev_annealing_pause_time[i],
                            ms_qc_qa_rev_annealing_s_target[i]
                    );
                } else if (ms_qc_mode[i] == "QAOA") {
                    seeding.set_qaoa_params(
                            ms_qc_qaoa_vendor[i],
                            ms_qc_qaoa_azure_subscription_id[i],
                            ms_qc_qaoa_azure_resource_group[i],
                            ms_qc_qaoa_azure_name[i],
                            ms_qc_qaoa_azure_location[i],
                            ms_qc_qaoa_azure_backend[i],
                            ms_qc_qaoa_optimizer[i],
                            ms_qc_qaoa_maxiter[i],
                            ms_qc_qaoa_reps[i],
                            ms_qc_qaoa_n_shots[i],
                            ms_qc_qaoa_is_recursive_qaoa[i]
                    );
                } else {
                    throw epi::Error("Unknown value für --ms-qc-mode");
                }

                if (is_biogrid) pipeline.add_network_BIOGRID(seeding, network_local_search);
                else pipeline.add_network(*connector, seeding, network_local_search);
            } else {
                throw epi::Error("Unknown value for --ms-seeding-routine");
            }

            if (!is_biogrid) delete connector;
        }


        // final local search
        if (!disable_final_search && num_networks > 1) {
            if (fs_model.empty()) {
                throw epi::Error("Parameter --fs-model missing.");
            }
            double fs_per_seed_time_limit_minutes = 0.;
            if (!fs_per_seed_time_limit_str.empty())
                fs_per_seed_time_limit_minutes = epi::parseTimespanString(fs_per_seed_time_limit_str);

            double fs_search_time_limit_minutes = 0.;
            if (!fs_search_time_limit_str.empty())
                fs_search_time_limit_minutes = epi::parseTimespanString(fs_search_time_limit_str);

            epi::LocalSearch final_local_search(
                    fs_model,
                    true, // collapse identical results
                    fs_max_rounds,
                    fs_search_time_limit_minutes,
                    fs_per_seed_time_limit_minutes,
                    fs_annealing_type,
                    fs_cooling_factor,
                    fs_annealing_start_prob,
                    fs_annealing_end_prob,
                    true, // score development
                    "result_search_score_over_time",
                    fs_min_set_size,
                    fs_max_set_size
            );

            if (fs_ld_check) {
                if (fs_ld_cutoff == -std::numeric_limits<double>::infinity()) {
                    final_local_search.activate_LD_check(fs_ld_matrix, fs_ld_mode, fs_ld_mc_min_set, fs_ld_mc_max_set,
                                                         fs_ld_mc_num_samples);
                } else {
                    final_local_search.activate_LD_check(fs_ld_matrix, fs_ld_mode, fs_ld_cutoff);
                }
            }
            if (fs_seeding_routine == "RANDOM_CONNECTED") {
                pipeline.set_final_local_search(epi::SeedingRandomConnected(fs_rc_start_seeds), final_local_search);
            } else if (fs_seeding_routine == "COMMUNITY_WISE") {
                pipeline.set_final_local_search(epi::SeedingCommunityWise(
                        fs_model,
                        fs_cw_quantile,
                        fs_cw_max_cluster_size,
                        fs_cw_num_sets_per_cluster,
                        fs_cw_num_snps_per_set
                ), final_local_search);
            } else if (fs_seeding_routine == "QUANTUM_COMPUTING") {
                epi::SeedingQuantumComputing seeding_qc(
                        fs_qc_mode,
                        fs_model,
                        fs_cw_quantile,
                        fs_cw_max_cluster_size,
                        fs_cw_num_sets_per_cluster,
                        fs_cw_num_snps_per_set,
                        fs_qc_min_cluster_size,
                        fs_qc_n_clique,
                        fs_qc_k,
                        fs_qc_nu,
                        fs_qc_lambda0,
                        fs_qc_lambda1,
                        fs_qc_lambda2
                );

                if (fs_qc_mode == "SIMULATED_ANNEALING") {
                    seeding_qc.set_simulated_annealing_params(
                            fs_qc_sa_num_samples,
                            fs_qc_sa_num_sweeps,
                            fs_qc_sa_seed
                    );
                } else if (fs_qc_mode == "QUANTUM_ANNEALING") {
                    seeding_qc.set_quantum_annealing_params(
                            fs_qc_qa_token,
                            fs_qc_qa_num_reads,
                            fs_qc_qa_solver_idx,
                            fs_qc_qa_fw_annealing_ramp_time,
                            fs_qc_qa_fw_annealing_pause_time,
                            fs_qc_qa_rev_annealing_ramp_time,
                            fs_qc_qa_rev_annealing_pause_time,
                            fs_qc_qa_rev_annealing_s_target
                    );
                } else if (fs_qc_mode == "QAOA") {
                    seeding_qc.set_qaoa_params(
                            fs_qc_qaoa_vendor,
                            fs_qc_qaoa_azure_subscription_id,
                            fs_qc_qaoa_azure_resource_group,
                            fs_qc_qaoa_azure_name,
                            fs_qc_qaoa_azure_location,
                            fs_qc_qaoa_azure_backend,
                            fs_qc_qaoa_optimizer,
                            fs_qc_qaoa_maxiter,
                            fs_qc_qaoa_reps,
                            fs_qc_qaoa_n_shots,
                            fs_qc_qaoa_is_recursive_qaoa
                    );
                } else {
                    throw epi::Error("Unknown value for --fs-qc-mode");
                }

                pipeline.set_final_local_search(seeding_qc, final_local_search);
            } else {
                throw epi::Error("Unknown value for --fs-seeding-routine");
            }
        }

// run the created pipeline
        pipeline.run();
    } catch (...) {
        // capture the exception
        auto eptr = std::current_exception();
        Logger::logError(eptr);
        std::rethrow_exception(eptr);
    }
}