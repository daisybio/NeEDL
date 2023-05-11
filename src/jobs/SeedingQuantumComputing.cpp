//
// Created by juli on 22.06.22.
//

#include "SeedingQuantumComputing.hpp"

#include <utility>
#include "../util/TimeLogger.hpp"
#include "../../quepistasis/header/snps_optimization.h"

namespace epi {
    SeedingQuantumComputing::SeedingQuantumComputing(const std::string& qc_mode, const std::string &scoring_model,
                                                     double quantile,
                                                     unsigned long max_cluster_size,
                                                     unsigned long num_sets_per_cluster,
                                                     unsigned long num_snps_per_set,
                                                     unsigned long min_qc_cluster_size,
                                                     int qc_n_clique, int qc_k, double qc_nu, double qc_lambda0,
                                                     double qc_lambda1, double qc_lambda2)
    : SeedingCommunityWise(scoring_model, quantile, max_cluster_size, num_sets_per_cluster, num_snps_per_set)
    {
        if (qc_mode == "SIMULATED_ANNEALING") {
            QCMode = SIMULATED_ANNEALING;
            set_simulated_annealing_params();
        } else if (qc_mode == "QUANTUM_ANNEALING") {
            QCMode = QUANTUM_ANNEALING;
            set_quantum_annealing_params();
        } else if (qc_mode == "QAOA") {
            QCMode = QAOA;
            set_qaoa_params();
        } else {
            throw epi::Error("Unknown quantum computing mode " + qc_mode);
        }

        if (min_qc_cluster_size < 2) {
            throw epi::Error("min_qc_cluster_size must be >= 2");
        }

        this->min_qc_cluster_size = min_qc_cluster_size;
        this->qc_n_clique = qc_n_clique;
        this->qc_k = qc_k;
        this->qc_nu = qc_nu;
        this->qc_lambda0 = qc_lambda0;
        this->qc_lambda1 = qc_lambda1;
        this->qc_lambda2 = qc_lambda2;
    }

    void SeedingQuantumComputing::set_simulated_annealing_params(int num_samples, int num_sweeps, uint64_t seed) {
        if (QCMode != SIMULATED_ANNEALING) {
            throw epi::Error("Setting attributes for SIMULATED_ANNEALING only allowed if selected mode is SIMULATED_ANNEALING");
        }
        qc_sa_num_samples = num_samples;
        qc_sa_num_sweeps = num_sweeps;
        qc_sa_seed = seed;
    }

    void SeedingQuantumComputing::set_quantum_annealing_params(std::string token, int num_reads, int solver_idx,
                                                               double fw_annealing_ramp_time,
                                                               double fw_annealing_pause_time,
                                                               double rev_annealing_ramp_time,
                                                               double rev_annealing_pause_time,
                                                               double rev_annealing_s_target) {
        if (QCMode != QUANTUM_ANNEALING) {
            throw epi::Error("Setting attributes for QUANTUM_ANNEALING only allowed if selected mode is QUANTUM_ANNEALING");
        }

        qc_qa_token = std::move(token);
        qc_qa_num_reads = num_reads;
        qc_qa_solver_idx = solver_idx;
        qc_qa_fw_annealing_ramp_time = fw_annealing_ramp_time;
        qc_qa_fw_annealing_pause_time = fw_annealing_pause_time;
        qc_qa_rev_annealing_ramp_time = rev_annealing_ramp_time;
        qc_qa_rev_annealing_pause_time = rev_annealing_pause_time;
        qc_qa_rev_annealing_s_target = rev_annealing_s_target;
    }


    void SeedingQuantumComputing::set_qaoa_params(std::string vendor, std::string azure_subscription_id,
                                                  std::string azure_resource_group, std::string azure_name,
                                                  std::string azure_location, std::string azure_backend,
                                                  std::string optimizer, int maxiter, int reps, int n_shots,
                                                  int is_recursive_qaoa) {

        if (QCMode != QAOA) {
            throw epi::Error("Setting attributes for QAOA only allowed if selected mode is QAOA");
        }

        qc_qaoa_vendor = std::move(vendor);
        qc_qaoa_azure_subscription_id = std::move(azure_subscription_id);
        qc_qaoa_azure_resource_group = std::move(azure_resource_group);
        qc_qaoa_azure_name = std::move(azure_name);
        qc_qaoa_azure_location = std::move(azure_location);
        qc_qaoa_azure_backend = std::move(azure_backend);
        qc_qaoa_optimizer = std::move(optimizer);
        qc_qaoa_maxiter = maxiter;
        qc_qaoa_reps = reps;
        qc_qaoa_n_shots = n_shots;
        qc_qaoa_is_recursive_qaoa = is_recursive_qaoa;
    }

    void SeedingQuantumComputing::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("seeding with method QUANTUM_COMPUTING");
        // create a clustering with respect to maximal cluster size for QC

        TimeLogger clusteringLogger("creating optimal clustering");
        unsigned long min_size, max_size;
        auto clusters = leiden_with_size_constraint(data->snpNetwork, max_cluster_size, true, min_size, max_size);
        set_clusters_as_attributes(data, "leiden_size_constraint", clusters);
        Logger::logLine("clustering: " + std::to_string(clusters.size()) + " clusters with size in [" + std::to_string(min_size) + ", " + std::to_string(max_size) + "]");

        // refine clustering
        refine_leiden_clustering(data, clusters, max_cluster_size);
        set_clusters_as_attributes(data, "after_refinement", clusters);
        Logger::logLine("after cluster refinement: " + std::to_string(clusters.size()) + " clusters");
        clusteringLogger.stop();

        TimeLogger qcLogger("separating clusters into small and large ones");
        std::vector<std::vector<SNP_t>> small_clusters, large_clusters;
        for (const auto & cluster : clusters) {
            if (cluster.size() < min_qc_cluster_size) small_clusters.push_back(cluster);
            else large_clusters.push_back(cluster);
        }

        // seed generation
        Logger::logLine("generate random seeds for small communities (< " + std::to_string(min_qc_cluster_size) + ")");
        auto candidate_snp_sets_small = generate_random_sets(data, small_clusters);

        Logger::logLine("generate start seeds with quantum computing for large communities (>= " + std::to_string(min_qc_cluster_size) + ")");
        auto candidate_snp_sets_large = generate_qc_sets(data, large_clusters);

        Logger::logLine("join found results");
        for (const auto & cluster : candidate_snp_sets_large) candidate_snp_sets_small.push_back(cluster);

        // seed selection
        Logger::logLine("selecting start seeds based on selected quantile");
        data->snpSetStorage = select_start_seeds(candidate_snp_sets_small);
        Logger::logLine("Selected " + std::to_string(data->snpSetStorage.size()) + " start seeds for local search.");

        logger.stop();
    }


    std::vector<std::vector<SNPSet>>
    SeedingQuantumComputing::generate_qc_sets(const std::shared_ptr<DataModel> & data, const std::vector<std::vector<SNP_t>> &clusters) {
        std::vector<std::vector<SNPSet>> result_sets(clusters.size());

        size_t cluster_processed = 1;
        size_t num_sets_exceeded_max_size = 0;
#pragma omp parallel for default(none) shared(clusters, cluster_processed, data, result_sets) reduction(+:num_sets_exceeded_max_size)
        for (size_t cluster_i = 0; cluster_i < clusters.size(); cluster_i++) {
            auto &cluster = clusters[cluster_i];
            std::unordered_set<SNPSet, SNPSetHash> cluster_result_sets;

            Logger::logProgress("Processing cluster " + std::to_string(cluster_processed) + " of " +
                                std::to_string(clusters.size()));

#pragma omp critical
            {
                cluster_processed++;
            }

            int N_SNPS = cluster.size();
            auto *stat_corr = new double[N_SNPS * N_SNPS]{};
            auto *bio_corr = new double[N_SNPS * N_SNPS]{};

            for (size_t i = 0; i < cluster.size(); i++) {
                // set diagonal to 0
                stat_corr[i + i * cluster.size()] = 0;

                // get all connected snps for cluster[i]
                auto connected_snps = data->snpNetwork->get_adjacent_snps(cluster[i]);
                std::unordered_set<SNP_t, SNP_t::SNPHash> connected_map(connected_snps.begin(),
                                                                        connected_snps.end());

                // create upper triangle matrix
                for (size_t j = i + 1; j < cluster.size(); j++) {
                    // pairwise score --> statistical correlation
                    double score = data->snpStorage->calculate_score({{cluster[i], cluster[j]}}, epi_model);
                    stat_corr[j + i * cluster.size()] = score;

                    // edge in SSI-network --> biological correlation
                    double connect_val = connected_map.find(cluster[j]) == connected_map.end() ? 0. : 1.;
                    bio_corr[j + i * cluster.size()] = connect_val;
                }
            }

            matrix stat_corr_mat(stat_corr, N_SNPS);
            matrix bio_corr_mat(bio_corr, N_SNPS);

            snps_qubo_matrix qubo(qc_n_clique, N_SNPS);

            qubo.fill_n_max_weighted_k_clique_qubo(stat_corr_mat, bio_corr_mat, qc_k, qc_nu, qc_lambda0, qc_lambda1,
                                                   qc_lambda2);

            std::vector<std::vector<int>> snps_set_list;

            std::string temp_file = data->outputDirectory->get_free_filepath("qc_data_" + std::to_string(cluster_i), ".txt");

            switch (QCMode) {
                case SIMULATED_ANNEALING:
                    snps_set_list = qubo.solve_cpu_simulated_annealing(qc_sa_num_samples,
                                                                   qc_sa_num_sweeps,
                                                                   qc_sa_seed,
                                                                   nullptr, temp_file.c_str());
                    break;
                case QUANTUM_ANNEALING:
                    snps_set_list = qubo.solve_dwave_quantum_annealing(qc_qa_token.c_str(),
                                                                 qc_qa_num_reads,
                                                                 qc_qa_solver_idx,
                                                                 qc_qa_fw_annealing_ramp_time,
                                                                 qc_qa_fw_annealing_pause_time,
                                                                 qc_qa_rev_annealing_ramp_time,
                                                                 qc_qa_rev_annealing_pause_time,
                                                                 qc_qa_rev_annealing_s_target,
                                                                 temp_file.c_str());
                    break;
                case QAOA:
                    snps_set_list = qubo.solve_azure_QAOA(qc_qaoa_vendor.c_str(),
                                                    qc_qaoa_azure_subscription_id.c_str(),
                                                    qc_qaoa_azure_resource_group.c_str(),
                                                    qc_qaoa_azure_name.c_str(),
                                                    qc_qaoa_azure_location.c_str(),
                                                    qc_qaoa_azure_backend.c_str(),
                                                    qc_qaoa_optimizer.c_str(),
                                                    qc_qaoa_maxiter,
                                                    qc_qaoa_reps,
                                                    qc_qaoa_n_shots,
                                                    qc_qaoa_is_recursive_qaoa,
                                                    temp_file.c_str());
                    break;
            }

            for (auto &set: snps_set_list) {
                std::vector<SNP_t> set_converted;
                for (auto &snp: set) set_converted.push_back(cluster[snp]);

                // TODO: fix qc step to only return sets with correct size
                if (set_converted.size() <= MAXIMUM_SNP_SET_SIZE) {
                    auto converted = SNPSet(set_converted);
                    converted.set_attribute("SEED_ORIGIN", "QUANTUM_COMPUTING");
                    cluster_result_sets.insert(converted);
                } else ++num_sets_exceeded_max_size;
            }

            result_sets[cluster_i] = std::vector<SNPSet>(cluster_result_sets.begin(), cluster_result_sets.end());
        }

        Logger::logLine(std::to_string(num_sets_exceeded_max_size) + " SNP-sets returned by the QC method exceeded the maximum SNP set size and where thus omitted.");

        return result_sets;
    }

    rapidjson::Value SeedingQuantumComputing::getConfig(rapidjson::Document &doc) {
        auto obj= SeedingCommunityWise::getConfig(doc);
        obj.RemoveMember("JOB");
        obj.AddMember("JOB", rapidjson::Value().SetString("SeedingQuantumComputing"), doc.GetAllocator());
        obj.AddMember("n_clique", rapidjson::Value().SetInt(qc_n_clique), doc.GetAllocator());
        obj.AddMember("k", rapidjson::Value().SetInt(qc_k), doc.GetAllocator());
        obj.AddMember("nu", rapidjson::Value().SetDouble(qc_nu), doc.GetAllocator());
        obj.AddMember("lambda0", rapidjson::Value().SetDouble(qc_lambda0), doc.GetAllocator());
        obj.AddMember("lambda1", rapidjson::Value().SetDouble(qc_lambda1), doc.GetAllocator());
        obj.AddMember("lambda2", rapidjson::Value().SetDouble(qc_lambda2), doc.GetAllocator());
        obj.AddMember("qc_mode", rapidjson::Value().SetString(QCMode_names[QCMode].c_str(), QCMode_names[QCMode].size(), doc.GetAllocator()), doc.GetAllocator());

        if (QCMode == SIMULATED_ANNEALING) {
            obj.AddMember("num_samples", rapidjson::Value().SetInt(qc_sa_num_samples), doc.GetAllocator());
            obj.AddMember("num_sweeps", rapidjson::Value().SetInt(qc_sa_num_sweeps), doc.GetAllocator());
            obj.AddMember("seed", rapidjson::Value().SetUint64(qc_sa_seed), doc.GetAllocator());
        } else if (QCMode == QUANTUM_ANNEALING) {
            obj.AddMember("token", rapidjson::Value().SetString(qc_qa_token.c_str(), qc_qa_token.size(), doc.GetAllocator()), doc.GetAllocator());
            obj.AddMember("num_reads", rapidjson::Value().SetInt(qc_qa_num_reads), doc.GetAllocator());
            obj.AddMember("solver_idx", rapidjson::Value().SetInt(qc_qa_solver_idx), doc.GetAllocator());
            obj.AddMember("fw_annealing_ramp_time", rapidjson::Value().SetDouble(qc_qa_fw_annealing_ramp_time), doc.GetAllocator());
            obj.AddMember("fw_annealing_pause_time", rapidjson::Value().SetDouble(qc_qa_fw_annealing_pause_time), doc.GetAllocator());
            obj.AddMember("rev_annealing_ramp_time", rapidjson::Value().SetDouble(qc_qa_rev_annealing_ramp_time), doc.GetAllocator());
            obj.AddMember("rev_annealing_pause_time", rapidjson::Value().SetDouble(qc_qa_rev_annealing_pause_time), doc.GetAllocator());
            obj.AddMember("rev_annealing_s_target", rapidjson::Value().SetDouble(qc_qa_rev_annealing_s_target), doc.GetAllocator());
        } else if (QCMode == QAOA) {
            obj.AddMember("vendor", rapidjson::Value().SetString(qc_qaoa_vendor.c_str(), qc_qaoa_vendor.size(), doc.GetAllocator()), doc.GetAllocator());
            obj.AddMember("azure_subscription_id", rapidjson::Value().SetString(qc_qaoa_azure_subscription_id.c_str(), qc_qaoa_azure_subscription_id.size(), doc.GetAllocator()), doc.GetAllocator());
            obj.AddMember("azure_resource_group", rapidjson::Value().SetString(qc_qaoa_azure_resource_group.c_str(), qc_qaoa_azure_resource_group.size(), doc.GetAllocator()), doc.GetAllocator());
            obj.AddMember("azure_name", rapidjson::Value().SetString(qc_qaoa_azure_name.c_str(), qc_qaoa_azure_name.size(), doc.GetAllocator()), doc.GetAllocator());
            obj.AddMember("azure_location", rapidjson::Value().SetString(qc_qaoa_azure_location.c_str(), qc_qaoa_azure_location.size(), doc.GetAllocator()), doc.GetAllocator());
            obj.AddMember("azure_backend", rapidjson::Value().SetString(qc_qaoa_azure_backend.c_str(), qc_qaoa_azure_backend.size(), doc.GetAllocator()), doc.GetAllocator());
            obj.AddMember("optimizer", rapidjson::Value().SetString(qc_qaoa_optimizer.c_str(), qc_qaoa_optimizer.size(), doc.GetAllocator()), doc.GetAllocator());
            obj.AddMember("maxiter", rapidjson::Value().SetInt(qc_qaoa_maxiter), doc.GetAllocator());
            obj.AddMember("reps", rapidjson::Value().SetInt(qc_qaoa_reps), doc.GetAllocator());
            obj.AddMember("n_shots", rapidjson::Value().SetInt(qc_qaoa_n_shots), doc.GetAllocator());
            obj.AddMember("is_recursive_qaoa", rapidjson::Value().SetInt(qc_qaoa_is_recursive_qaoa), doc.GetAllocator());
        }
        return obj;
    }


} // epi