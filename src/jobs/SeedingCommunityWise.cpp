//
// Created by juli on 24.08.22.
//

#include "SeedingCommunityWise.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    SeedingCommunityWise::SeedingCommunityWise(
            const std::string &scoring_model,
            double quantile,
            unsigned long max_cluster_size,
            unsigned long num_sets_per_cluster,
            unsigned long num_snps_per_set
    ) {
        epi_model = options::epistasis_score_from_string(scoring_model);
        selection_quantile = quantile;
        this->max_cluster_size = max_cluster_size;
        this->num_sets_per_cluster = num_sets_per_cluster;
        this->num_snps_per_set = num_snps_per_set;
    }


    std::vector<std::vector<SNP_t>>
    SeedingCommunityWise::leiden_with_size_constraint(const std::shared_ptr<SNPNetwork>& snpNetwork, unsigned long maximal_allowed_cluster_size, bool optimize_min, unsigned long &best_min_size, unsigned long &best_max_size) {
        // this stuff here is binary search without sub-clustering
        bool reached_max = false;
        float resolution_min = 0.f; // std::numeric_limits<float>::min();
        float resolution_max = 0.f; // std::numeric_limits<float>::min();
        size_t binary_search_steps = 0;

        std::vector<std::vector<SNP_t>> best_clustering;
        float best_resolution = 0.f;
        best_min_size = 0;
        best_max_size = 0;

        size_t cluster_min_size, cluster_max_size, previous_max_size = 0, previous_num_clusters; // , cluster_largest_clique = 0;
        do {
            cluster_min_size = -1;
            cluster_max_size = 0;
            float resolution;
            if (reached_max) {
                // do binary search to find optimal clustering
                resolution = (resolution_min + resolution_max) / 2;
            } else {
                // do the forward pass
                resolution = resolution_max;
            }

            auto clusters = snpNetwork->cluster_leiden(resolution, leiden_beta, leiden_max_leiden_steps);

            for (auto &cluster: clusters) {
                cluster_min_size = std::min(cluster_min_size, cluster.size());
                cluster_max_size = std::max(cluster_max_size, cluster.size());
            }

            if (reached_max) {
                // binary search
                if (optimize_min) {
                    if (cluster_max_size <= maximal_allowed_cluster_size) {
                        resolution_max = resolution;
                        best_clustering = clusters;
                        best_resolution = resolution;
                        best_min_size = cluster_min_size;
                        best_max_size = cluster_max_size;
                    } else {
                        resolution_min = resolution;
                    }
                } else {
                    if (cluster_min_size <= maximal_allowed_cluster_size) {
                        resolution_max = resolution;
                        best_clustering = clusters;
                        best_resolution = resolution;
                        best_min_size = cluster_min_size;
                        best_max_size = cluster_max_size;
                    } else {
                        resolution_min = resolution;
                    }
                }
                binary_search_steps ++;
            } else {
                // forward pass
                if ((optimize_min && cluster_max_size <= maximal_allowed_cluster_size) ||
                    (!optimize_min && cluster_min_size <= maximal_allowed_cluster_size)) {
                    reached_max = true;

                    best_clustering = clusters;
                    best_resolution = resolution;
                    best_min_size = cluster_min_size;
                    best_max_size = cluster_max_size;
                } else {
                    best_clustering = clusters;
                    resolution_min = resolution_max;
                    resolution_max += leiden_forward_search_speed;
                }
            }

            if (!reached_max && previous_max_size == cluster_max_size && previous_num_clusters == clusters.size()) break;
            previous_max_size = cluster_max_size;
            previous_num_clusters = clusters.size();
        } while (!reached_max || (binary_search_steps < leiden_num_binary_search_steps && resolution_min != resolution_max));

        return best_clustering;
    }

    void SeedingCommunityWise::refine_leiden_clustering(const std::shared_ptr<DataModel>& data, std::vector<std::vector<SNP_t>> & clusters, size_t maximal_allowed_cluster_size) {
        // create node cluster map so connections can be found easily
        std::unordered_map<SNP_t, size_t, SNP_t::SNPHash> node_cluster_map;
        for (size_t i = 0; i < clusters.size(); i++) {
            for (auto& node : clusters[i]) node_cluster_map.insert({node, i });
        }

        for (size_t i = 0; i < clusters.size(); i++) {
            auto &cluster = clusters[i];
            if (cluster.size() < maximal_allowed_cluster_size && !cluster.empty()) {
                // try to find any other cluster that is connected to this one
                std::unordered_set<size_t> other_clusters;
                for (auto& node : cluster) {
                    auto adjacent_nodes = data->snpNetwork->get_adjacent_snps(node);
                    for (auto & adjacent_node : adjacent_nodes) {
                        auto map_item = node_cluster_map.find(adjacent_node);
                        if (map_item != node_cluster_map.end()) {
                            size_t cluster_index = map_item->second;
                            if (cluster_index != i) other_clusters.insert(cluster_index);
                        }
                    }
                }

                if (!other_clusters.empty()) {
                    // there is a connection to at least one other cluster with connection
                    std::vector<std::pair<size_t,size_t>> other_clusters_with_size;
                    for (auto& c : other_clusters) other_clusters_with_size.emplace_back( c, clusters[c].size() );
                    // sort by cluster size
                    std::sort(other_clusters_with_size.begin(), other_clusters_with_size.end(), [] (std::pair<size_t,size_t> a, std::pair<size_t,size_t> b) {
                        return a.second < b.second;
                    });

                    // add other clusters to this one until the maximal size is reached or no other clusters are left
                    for (auto &other : other_clusters_with_size) {
                        // no other cluster can be added
                        if (cluster.size() + other.second > maximal_allowed_cluster_size) break;

                        // add this cluster
                        auto &other_cluster = clusters[other.first];
                        for (auto& other_node : other_cluster) node_cluster_map[other_node] = i;
                        cluster.insert(cluster.end(), other_cluster.begin(), other_cluster.end());
                        other_cluster.clear();
                    }
                }
            }
        }

        // remove empty clusters
        std::vector<std::vector<SNP_t>> refined;
        for (auto &cluster : clusters) {
            if (!cluster.empty()) refined.push_back(cluster);
        }
        clusters = refined;
    }

    void SeedingCommunityWise::set_clusters_as_attributes(const std::shared_ptr<DataModel> &data,
                                                             std::string attribute_name,
                                                             std::vector<std::vector<SNP_t>> clusters) {

        unsigned long i = 0;
        for (auto & cluster : clusters) {
            for (auto & snp : cluster) {
                data->snpStorage->snp_set_variable_attribute(snp, attribute_name, std::to_string(i));
            }
            i++;
        }
    }

    std::vector<std::vector<SNPSet>> SeedingCommunityWise::generate_random_sets(const std::shared_ptr<DataModel> &data,
                                                                                const std::vector<std::vector<SNP_t>> &clusters) {
        std::vector<std::vector<SNPSet>> result_sets(clusters.size());

        size_t cluster_processed = 1;
#pragma omp parallel for default(none) shared(clusters, cluster_processed, num_sets_per_cluster, num_snps_per_set, data, result_sets)
        for (size_t cluster_i = 0; cluster_i < clusters.size(); cluster_i++) {
            auto &cluster = clusters[cluster_i];
            std::unordered_set<SNPSet, SNPSetHash> cluster_result_sets;

            Logger::logProgress("Processing cluster " + std::to_string(cluster_processed) + " of " +
                                std::to_string(clusters.size()));

#pragma omp critical
            {
                cluster_processed++;
            }

            // add a start seed for that cluster
            if (cluster.size() <= num_snps_per_set) {
                // insert full cluster as start seed
                auto set = SNPSet(cluster);
                set.set_attribute("SEED_ORIGIN", "COMMUNITY_WISE");
                cluster_result_sets.insert(set);
            } else {
                for (size_t j = 0; j < num_sets_per_cluster; j++) {
                    // pick random snp selection that is connected
                    std::uniform_int_distribution<size_t> dist(0, cluster.size() - 1);
                    size_t start_snp_pos = dist(data->random_device[omp_get_thread_num()]);
                    const SNP_t &start_snp = cluster[start_snp_pos];

                    std::vector<SNP_t> selected_snps({start_snp});
                    std::set<SNP_t> add_options;
                    auto adj_snps_start = data->snpNetwork->get_adjacent_snps(start_snp);
                    // checks for every snp if it is in the cluster and is not yet selected
                    for (auto &snp: adj_snps_start) {
                        if (std::find(cluster.begin(), cluster.end(), snp) != cluster.end() &&
                            std::find(selected_snps.begin(), selected_snps.end(), snp) == selected_snps.end()) {

                            add_options.insert(snp);
                        }
                    }

                    for (size_t i = 1; i < num_snps_per_set; i++) {
                        if (add_options.empty()) break;

                        // pick random option from the add_options set and remove it from the set
                        dist = std::uniform_int_distribution<size_t>(0, add_options.size() - 1);
                        size_t pos = dist(data->random_device[omp_get_thread_num()]);
                        auto snp_it = std::next(add_options.begin(), pos);
                        const SNP_t &next_snp = *snp_it;
                        add_options.erase(next_snp);
                        selected_snps.push_back(next_snp);

                        auto adj_snps = data->snpNetwork->get_adjacent_snps(start_snp);
                        // checks for every snp if it is in the cluster and is not yet selected
                        for (auto &snp: adj_snps) {
                            if (std::find(cluster.begin(), cluster.end(), snp) != cluster.end() &&
                                std::find(selected_snps.begin(), selected_snps.end(), snp) == selected_snps.end()) {

                                add_options.insert(snp);
                            }
                        }
                    }

                    auto set = SNPSet(selected_snps);
                    set.set_attribute("SEED_ORIGIN", "COMMUNITY_WISE");
                    cluster_result_sets.insert(set);
                }
            }

            result_sets[cluster_i] = std::vector<SNPSet>(cluster_result_sets.begin(), cluster_result_sets.end());
        }

        return result_sets;
    }

    std::vector<SNPSet> SeedingCommunityWise::select_start_seeds(std::vector<std::vector<SNPSet>> qc_sets) {
        // calculate all set scores first and order them by score (best score at the back)
        auto score_direction = SNPStorage::currentSnpStorage->need_minimizing_score(epi_model);

#pragma omp parallel for default(none) shared(qc_sets, score_direction)
        for (size_t i = 0; i < qc_sets.size(); i++) {
            auto & set_list = qc_sets[i];
            std::sort(set_list.begin(), set_list.end(), [&score_direction, this] (auto & left, auto & right) -> bool {
                return score_direction == (left.calculate_score(epi_model) > right.calculate_score(epi_model));
            });
        }

        std::unordered_set<SNPSet, SNPSetHash> selected_sets;

        // select the best score from every set
        for (auto & set_list : qc_sets) {
            if (!set_list.empty()) selected_sets.insert(set_list.back());
        }

        // join all scores together and select the quantile of best scores
        std::vector<SNPSet> all_candidates;
        for (auto & set_list : qc_sets) {
            all_candidates.insert(all_candidates.end(), set_list.begin(), set_list.end());
        }

        // sort again --> this time best scores at the front
        std::sort(all_candidates.begin(), all_candidates.end(), [&score_direction, this] (auto & left, auto & right) -> bool {
            return score_direction == (left.calculate_score(epi_model) < right.calculate_score(epi_model));
        });

        double num_remaining_exact = selection_quantile * double(all_candidates.size());
        auto num_remaining = long(ceil(num_remaining_exact));
        selected_sets.insert(all_candidates.begin(), all_candidates.begin() + num_remaining);

        return { selected_sets.begin(), selected_sets.end() };
    }

    void SeedingCommunityWise::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("seeding with method COMMUNITY_WISE");
        // create a clustering with respect to maximal cluster size for QC

        TimeLogger clusteringLogger("creating optimal clustering");
        unsigned long min_size, max_size;
        auto clusters = leiden_with_size_constraint(data->snpNetwork, max_cluster_size, true, min_size, max_size);
        set_clusters_as_attributes(data, "leiden_cluster", clusters);
        Logger::logLine("clustering: " + std::to_string(clusters.size()) + " clusters with size in [" + std::to_string(min_size) + ", " + std::to_string(max_size) + "]");

        // refine clustering
        refine_leiden_clustering(data, clusters, max_cluster_size);
        set_clusters_as_attributes(data, "leiden_cluster_after_refinement", clusters);
        Logger::logLine("after cluster refinement: " + std::to_string(clusters.size()) + " clusters left");
        clusteringLogger.stop();

        // seed generation
        Logger::logLine("generate random seeds for every community");
        auto candidate_snp_sets = generate_random_sets(data, clusters);

        // seed selection
        Logger::logLine("selecting start seeds based on selected quantile");
        data->snpSetStorage = select_start_seeds(candidate_snp_sets);
        Logger::logLine("Selected " + std::to_string(data->snpSetStorage.size()) + " start seeds for local search.");

        logger.stop();
    }

    rapidjson::Value SeedingCommunityWise::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("SeedingCommunityWise"), doc.GetAllocator());
        auto epiModelStr = options::epistasis_score_to_string(epi_model);
        obj.AddMember("model", rapidjson::Value().SetString(epiModelStr.c_str(), epiModelStr.size(), doc.GetAllocator()), doc.GetAllocator());

        obj.AddMember("quantile", rapidjson::Value().SetDouble(selection_quantile), doc.GetAllocator());
        obj.AddMember("max_cluster_size", rapidjson::Value().SetUint64(max_cluster_size), doc.GetAllocator());
        obj.AddMember("num_sets_per_cluster", rapidjson::Value().SetUint64(num_sets_per_cluster), doc.GetAllocator());
        obj.AddMember("num_snps_per_set", rapidjson::Value().SetUint64(num_snps_per_set), doc.GetAllocator());

        obj.AddMember("leiden_forward_search_speed", rapidjson::Value().SetFloat(leiden_forward_search_speed), doc.GetAllocator());
        obj.AddMember("leiden_num_binary_search_steps", rapidjson::Value().SetUint64(leiden_num_binary_search_steps), doc.GetAllocator());
        obj.AddMember("leiden_beta", rapidjson::Value().SetFloat(leiden_beta), doc.GetAllocator());
        obj.AddMember("leiden_max_leiden_steps", rapidjson::Value().SetUint64(leiden_max_leiden_steps), doc.GetAllocator());
        return obj;
    }

} // epi