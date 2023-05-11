//
// Created by juli on 25.07.22.
//

#include "ShuffleNetwork.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/RepeaterList.hpp"

namespace epi {
    void ShuffleNetwork::run(std::shared_ptr<DataModel> data) {
        TimeLogger *logger;
        if (data->snpNetwork == nullptr) throw epi::Error("No network available to apply shuffle method.");
        auto initial_network = *data->snpNetwork;

        switch(shuffle_algorithm) {
            case TOPOLOGY_PRESERVING_WITHOUT_SNP_DEGREE:
                logger = new TimeLogger("shuffle network with algorithm TOPOLOGY PRESERVING not preserving SNP degree");
                shuffle_topology_preserving(data, false);
                break;
            case TOPOLOGY_PRESERVING_WITH_SNP_DEGREE:
                logger = new TimeLogger("shuffle network with algorithm TOPOLOGY PRESERVING preserving SNP degree");
                shuffle_topology_preserving(data, true);
                break;
            case EXPECTED_DEGREE_KEEP_DEGREE_DISTRIBUTION:
                logger = new TimeLogger("shuffle network with algorithm EXPECTED DEGREE preserving the degree distribution");
                // shuffle_expected_degree(data, false);
                shuffle_topology_preserving(data, false);
                shuffle_expected_degree_ind3(data);
                break;
            case EXPECTED_DEGREE_KEEP_INDIVIDUAL_DEGREE:
                logger = new TimeLogger("shuffle network with algorithm EXPECTED DEGREE preserving each SNP's degree");
                // shuffle_expected_degree(data, true);
                // shuffle_expected_degree_ind2(data, 4);
                shuffle_expected_degree_ind3(data);
                break;
        }

        analyze_randomization_step(data, initial_network);

        logger->stop();
        delete logger;
    }

    void ShuffleNetwork::shuffle_topology_preserving(const std::shared_ptr<DataModel>& data, bool preserve_snp_degree) {
        // put snps into bins based on their degree
        std::unordered_map<size_t, std::vector<SNP_t>> degree_bins;
        auto snps_iter = data->snpNetwork->all();
        for (auto snp: snps_iter) {
            size_t degree = 0;
            if (preserve_snp_degree) {
                degree = data->snpNetwork->get_degree(snp);
            }

            auto bin_it = degree_bins.find(degree);
            if (bin_it == degree_bins.end()) {
                // insert new bin
                degree_bins.insert({degree, {snp}});
            } else {
                bin_it->second.push_back(snp);
            }
        }

        // shuffle graph nodes
        size_t num_changed = 0;
        size_t num_same = 0;

        std::vector<std::pair<SNP_t, SNP_t>> replacement_pairs;
        for (auto bin_it : degree_bins) {
            std::vector<SNP_t> new_indices (bin_it.second.begin(), bin_it.second.end() );
            std::shuffle(new_indices.begin(), new_indices.end(), data->random_device[omp_get_thread_num()]);
            for (size_t i = 0; i < bin_it.second.size(); i++) {
                if (bin_it.second[i] == new_indices[i])
                    num_same++;
                else {
                    num_changed++;
                    // change both snps
                    replacement_pairs.emplace_back(bin_it.second[i], new_indices[i]);
                }
            }
        }
        data->snpNetwork->replace_nodes(replacement_pairs);

        Logger::logLine("Shuffled SNPs: " + Logger::to_string(num_changed) + " changed, " + Logger::to_string(num_same) +
                        " unchanged");
    }

    void ShuffleNetwork::shuffle_expected_degree(const std::shared_ptr<DataModel>& data, bool preserve_individual_snp_degree) {
        // create list of snps and their degrees
        std::vector<std::pair<SNP_t, size_t>> node_degree_list;
        size_t num_edges_before = 0;
        // this block is used to delete helper vectors
        {
            // std::unordered_map<size_t, std::vector<SNP_t>> degrees;
            std::vector<SNP_t> snp_list;
            std::vector<size_t> degree_list;
            auto snps_iter = data->snpNetwork->all();
            for (auto snp: snps_iter) {
                snp_list.push_back(snp);
                degree_list.push_back(data->snpNetwork->get_degree(snp));
            }

            // shuffle degrees_list if requested
            if (!preserve_individual_snp_degree) {
                std::shuffle(degree_list.begin(), degree_list.end(), data->random_device[omp_get_thread_num()]);
            }

            // copy back into pairwise list and skip nodes without connection
            for (size_t i = 0; i < snp_list.size(); i++) {
                if (degree_list[i] > 0) node_degree_list.emplace_back(snp_list[i], degree_list[i]);
                num_edges_before += degree_list[i];
            }
        }

        // shuffle nodes in list as well to get better distr later
        std::shuffle(node_degree_list.begin(), node_degree_list.end(), data->random_device[omp_get_thread_num()]);

        // delete old network and re-insert all previous nodes
        data->snpNetwork->clear();
        for (auto & n : node_degree_list) {
            data->snpNetwork->add_node(n.first);
        }

        // create new edges
        size_t total_removed = 0;
        size_t last_reported = 0;
        size_t edges_done = 0;

        // convert to optimized RepeaterList
        RepeaterList<SNP_t> node_list (node_degree_list);

        struct backtracking_data {
            SNPEdge edge;           // edge (to remove from the network)
            size_t id1;             // group id 1 to restore count in node_list
            size_t id2;             // group id 2 to restore count in node_list
            size_t snp1_pos;        // to restore main search loop
            size_t snp2_pos;        // to restore main search loop
            size_t search_position;  // to restore main search loop
        };
        // in this list every decision is stored to be able to backtrack if we get a deadlock
        std::vector<struct backtracking_data> backtracking_list;

        size_t snp1_pos, snp2_pos, search_start;
        bool do_backtracking = false;

        while(node_list.size() > 1) {
            if (last_reported - 10 > edges_done) {
                last_reported = edges_done;
                Logger::logProgress("Processing edges: " + Logger::to_string(edges_done / 2) + " of " + Logger::to_string(num_edges_before / 2) + " done");
            }

            SNPEdge current_edge;

            // get two nodes that are not equal and do not already have a connection
            std::uniform_int_distribution<size_t> distr_edges(0, node_list.size() - 1);

            if (!do_backtracking) {
                // select random first SNP
                snp1_pos = distr_edges(data->random_device[omp_get_thread_num()]);

                // select random second SNP
                snp2_pos = node_list.get_group_start(distr_edges(data->random_device[omp_get_thread_num()]));
                search_start = 0;
            }

            auto snp1 = node_list[snp1_pos];
            bool found_pair = false;

            // walk through data with SNP2 until new combination is found or SNP1 can be discarded from search
            size_t snp2_real, i;
            for (i = search_start; i < node_list.size(); i++) {
                snp2_real = (snp2_pos + i) % node_list.size();
                // move to the end to skip all items of the same group
                size_t advanced = node_list.get_group_end(snp2_real);
                i += advanced - snp2_real;

                // std::cout << i << ", " << advanced << ", " << snp2_real << ", " << node_list.size() << std::endl;

                auto snp2 = node_list[snp2_real];

                // std::cout << SNPSet({snp1, snp2}).get_snp_string() << std::endl;

                // if both are equal or edge already exists continue search
                if (snp1 == snp2) continue;

                current_edge = { snp1, snp2 };
                if (data->snpNetwork->edge_exists(current_edge)) continue;

                // new combination found --> insert
                found_pair = true;
                data->snpNetwork->add_edge(current_edge);

                break;
            }


            // delete entry for snp1 even though not all necessary connections are made --> this can theoretically happen
            if (found_pair) {
                auto id1 = node_list.erase(std::max(snp1_pos, snp2_real));
                auto id2 = node_list.erase(std::min(snp1_pos, snp2_real));
                backtracking_list.push_back({ current_edge, id1, id2, snp1_pos, snp2_pos, i });
                if (do_backtracking) {
                    do_backtracking = false;
                }
                std::cout << '[' << backtracking_list.size() << "] selected: " << current_edge.get_snp_string() << std::endl;
                edges_done += 2;
            } else {
                if (!backtracking_list.empty()) {
                    // do backtracking
                    do_backtracking = true;
                    for (size_t j = 0; j < 1000 && !backtracking_list.empty(); j ++) {
                        auto last_decision = backtracking_list.back();
                        backtracking_list.pop_back();
                        // restore network
                        data->snpNetwork->remove_edge(last_decision.edge);
                        // restore node_list data structure
                        node_list.restore_item_of_group(last_decision.id1);
                        node_list.restore_item_of_group(last_decision.id2);
                        // restore loop params
                        snp1_pos = last_decision.snp1_pos;
                        snp2_pos = last_decision.snp2_pos;
                        search_start = last_decision.search_position + 1;
                        std::cout << '[' << backtracking_list.size() << "] do backtracking for SNP "
                                  << data->snpStorage->snp_get_name(node_list[snp1_pos]) << std::endl;
                    }
                } else {
                    std::cout << '[' << backtracking_list.size() << "] no backtracking possible for " << data->snpStorage->snp_get_name(snp1) << std::endl;
                    total_removed += node_list.erase_group(snp1_pos);
                    edges_done += 2;
                }
            }


        }
        Logger::logLine("Processing edges: " + Logger::to_string(edges_done / 2) + " of " + Logger::to_string(num_edges_before / 2) + " done");


        Logger::logLine(Logger::to_string(total_removed / 2) + " of " + Logger::to_string(num_edges_before / 2) + " edges were discarded.");
        // Logger::logLine(Logger::to_string(to_remove.size()) + " of " + Logger::to_string(num_nodes) + " SNPs were discarded.");

    }

    void ShuffleNetwork::shuffle_expected_degree_ind2(const std::shared_ptr<DataModel>& data, size_t num_iterations) {
        if (data->snpNetwork->num_edges() == 0) return;

        for (size_t iter = 1; iter <= num_iterations; iter++) {
            Logger::logLine("Shuffle iteration " + std::to_string(iter) + " of " + std::to_string(num_iterations));
            // create list of snps and their degrees
            std::vector<std::pair<SNP_t, size_t>> node_degree_list;
            std::unordered_map<SNP_t, std::vector<SNP_t>, SNP_t::SNPHash> adjacency_data;

            auto snps_iter = data->snpNetwork->all();
            for (auto snp: snps_iter) {
                node_degree_list.emplace_back(snp, data->snpNetwork->get_degree(snp));
                adjacency_data.insert({snp, data->snpNetwork->get_adjacent_snps(snp)});
            }

            // shuffle nodes in list as well to get better distr later
            std::shuffle(node_degree_list.begin(), node_degree_list.end(), data->random_device[omp_get_thread_num()]);

            // convert to optimized RepeaterList
            RepeaterList<SNP_t> node_list(node_degree_list);

            size_t num_rounds_failed = 0;
            // size_t num_retries = 0;
            size_t num_done = 0;

            std::unordered_set<SNPEdge, SNPEdgeHash> processed_snps;
            // size_t num_recurring = 0;

            size_t num_rounds = data->snpNetwork->num_edges(); // / 4;
            for (size_t i1 = 1; i1 <= num_rounds; i1++) {
                bool performed_change = false;
                // size_t retry_i;
                // for (retry_i = 0; retry_i < expected_degree_max_retries; retry_i++) {
                if (node_list.size() == 0) break;
                // select a random edge
                std::uniform_int_distribution<size_t> distr_snp1(0, node_list.size() - 1);
                auto snp1_pos = distr_snp1(data->random_device[omp_get_thread_num()]);
                auto snp1 = node_list[snp1_pos];

                // auto adjacent1 = data->snpNetwork->get_adjacent_snps(snp1);
                auto &adjacent1 = adjacency_data[snp1];

                std::uniform_int_distribution<size_t> distr_snp2(0, adjacent1.size() - 1);
                auto snp2_initial = distr_snp2(data->random_device[omp_get_thread_num()]);
                size_t snp2_pos;
                SNP_t snp2;
                for (size_t i2 = 0; i2 < adjacent1.size() && !performed_change; i2++) {
                    snp2_pos = (snp2_initial + i2) % adjacent1.size();
                    snp2 = adjacent1[snp2_pos];

                    if (processed_snps.find({snp1, snp2}) != processed_snps.end()) {
                        adjacent1.erase(adjacent1.begin() + snp2_pos);
                        i2 --;
                        continue;
                    }

                    // select a node that is neither connected to snp1 nor to snp2
                    auto snp3_initial = distr_snp1(data->random_device[omp_get_thread_num()]);
                    size_t snp3_pos;
                    SNP_t snp3;
                    for (size_t i3 = 0; i3 < node_list.size() && !performed_change; i3++) {
                        snp3_pos = (snp3_initial + i3) % node_list.size();
                        // move to the end to skip all items of the same group
                        size_t advanced = node_list.get_group_end(snp3_pos);
                        i3 += advanced - snp3_pos;

                        snp3 = node_list[snp3_pos];

                        if (!data->snpNetwork->edge_exists({snp1, snp3})
                            && !data->snpNetwork->edge_exists({snp2, snp3})) {
                            // this snp3 is valid

                            // find a snp4 that is connected to max. one of snp1 and snp2
                            // auto adjacent3 = data->snpNetwork->get_adjacent_snps(snp3);
                            auto &adjacent3 = adjacency_data[snp3];
                            std::uniform_int_distribution<size_t> distr_snp4(0, adjacent3.size() - 1);
                            auto snp4_initial = distr_snp4(data->random_device[omp_get_thread_num()]);
                            size_t snp4_pos;
                            SNP_t snp4;
                            for (size_t i4 = 0; i4 < adjacent3.size() && !performed_change; i4++) {
                                snp4_pos = (snp4_initial + i4) % adjacent3.size();
                                snp4 = adjacent3[snp4_pos];

                                if (processed_snps.find({snp3, snp4}) != processed_snps.end()) {
                                    adjacent3.erase(adjacent3.begin() + snp4_pos);
                                    i4 --;
                                    continue;
                                }

                                bool not_adjacent_to_1 = !data->snpNetwork->edge_exists({snp1, snp4});
                                bool not_adjacent_to_2 = !data->snpNetwork->edge_exists({snp2, snp4});

                                if (not_adjacent_to_1 || not_adjacent_to_2) {
                                    auto for_snp1 = snp3;
                                    auto for_snp2 = snp4;
                                    if (not_adjacent_to_1) {
                                        for_snp1 = snp4;
                                        for_snp2 = snp3;
                                    }

                                    SNPEdge edge_1(snp1, for_snp1);
                                    SNPEdge edge_2(snp2, for_snp2);

                                    bool edge_exists = processed_snps.find(edge_1) != processed_snps.end()
                                                       || processed_snps.find(edge_2) != processed_snps.end();

                                    if (!edge_exists) {
                                        processed_snps.insert(edge_1);
                                        processed_snps.insert(edge_2);
                                        processed_snps.insert({snp3, snp4});
                                        processed_snps.insert({snp1, snp2});

                                        adjacent1.erase(adjacent1.begin() + snp2_pos);
                                        adjacent3.erase(adjacent3.begin() + snp4_pos);

                                        node_list.erase(std::max(snp1_pos, snp3_pos));
                                        node_list.erase(std::min(snp1_pos, snp3_pos));

                                        // this snp4 is valid --> finally change edges
                                        data->snpNetwork->remove_edge({snp1, snp2});
                                        data->snpNetwork->remove_edge({snp3, snp4});

                                        data->snpNetwork->add_edge(edge_1);
                                        data->snpNetwork->add_edge(edge_2);

                                        performed_change = true;
                                    }
                                }
                            }
                        }

                    }

                    if (!performed_change) {
                        adjacent1.erase(adjacent1.begin() + snp2_pos);
                        i2 --;
                    }
                }

                num_done++;
                Logger::logProgress(
                        "exchanging edges: " + std::to_string(num_done) + " of " +
                        std::to_string(num_rounds) +
                        " done");

                // std::cout << "node_list: " << node_list.size() << ", processed_snps: " << processed_snps.size() << std::endl;

                if (!performed_change) {
                    num_rounds_failed++;
                    node_list.erase_group(snp1_pos);
                }
                // }

                // num_retries += retry_i;

            }

            Logger::logLine("exchanging edges: " + std::to_string(num_rounds) + " of " +
                            std::to_string(num_rounds) + " done");
            // Logger::logLine("number of retries: " + std::to_string(num_retries) + ", completely failed iterations: " +
            //                std::to_string(num_iterations_failed));
            Logger::logLine("# failed exchanges: " + std::to_string(num_rounds_failed));
            // Logger::logLine("number of recurring: " + std::to_string(num_recurring));
        }
    }

    void ShuffleNetwork::analyze_randomization_step(const std::shared_ptr<DataModel> &data, epi::SNPNetwork &initial_network) {
        Logger::logLine("Analyzing resulting network");
        // check how many edges changed
        auto snps_new_iter = data->snpNetwork->all();
        size_t num_degree_changed = 0;
        size_t edges_found = 0;
        size_t edges_unchanged = 0;
        double total_degree_change = 0, total_degree_change_abs = 0;
        for (auto snp: snps_new_iter) {
            auto adjacent_old = initial_network.get_adjacent_snps(snp);
            auto adjacent_new = data->snpNetwork->get_adjacent_snps(snp);

            if (adjacent_old.size() != adjacent_new.size()) {
                num_degree_changed++;
                auto diff = long(adjacent_new.size()) - long(adjacent_old.size());
                total_degree_change += diff;
                total_degree_change_abs += std::abs(diff);
            }

            std::sort(adjacent_old.begin(), adjacent_old.end());
            std::sort(adjacent_new.begin(), adjacent_new.end());
            std::vector<SNP_t> intersect;
            std::set_intersection(adjacent_old.begin(), adjacent_old.end(), adjacent_new.begin(), adjacent_new.end(),
                                  std::back_inserter(intersect));

            edges_found += adjacent_new.size();
            edges_unchanged += intersect.size();
        }

        total_degree_change /= data->snpNetwork->num_nodes();
        total_degree_change_abs /= data->snpNetwork->num_nodes();

        Logger::logLine("degree change: num. nodes: " + std::to_string(num_degree_changed) + ", avg. change: " +
                        Logger::to_string(total_degree_change) + ", avg. abs. change: " +
                        Logger::to_string(total_degree_change_abs));
        Logger::logLine(
                "unchanged edges: " + std::to_string(edges_unchanged / 2) + " of " + std::to_string(edges_found / 2));
    }

    void ShuffleNetwork::shuffle_expected_degree_ind3(const std::shared_ptr<DataModel> &data) {
        size_t num_iterations = data->snpNetwork->num_edges();

        // create list of snps and their degrees
        std::vector<std::pair<SNP_t, size_t>> node_degree_list;
        auto snps_iter = data->snpNetwork->all();
        for (auto snp: snps_iter) {
            node_degree_list.emplace_back(snp, data->snpNetwork->get_degree(snp));
        }
        // shuffle nodes in list as well to get better distr later
        std::shuffle(node_degree_list.begin(), node_degree_list.end(), data->random_device[omp_get_thread_num()]);
        // convert to optimized RepeaterList
        RepeaterList<SNP_t> node_list(node_degree_list);

        data->snpNetwork->clear_edges();

        size_t num_failed = 0, num_successful = 0;
        std::uniform_int_distribution<size_t> distr(0, node_list.size() - 1);
        while (num_successful < num_iterations && num_failed < num_iterations) {
            auto snp1_pos = distr(data->random_device[omp_get_thread_num()]);
            auto snp1 = node_list[snp1_pos];

            auto group_start = node_list.get_group_start(snp1_pos);
            auto group_end = node_list.get_group_end(snp1_pos);

            std::uniform_int_distribution<size_t> distr(0, node_list.size() - 1 - (group_end - group_start + 1));
            auto snp2_pos = distr(data->random_device[omp_get_thread_num()]);
            if (snp2_pos >= group_start) snp2_pos += group_end - group_start + 1;
            auto snp2 = node_list[snp2_pos];

            if (snp1 == snp2) {
                throw epi::Error("This is not supposed to happen!!");
            }

            SNPEdge edge (snp1, snp2);

            if (!data->snpNetwork->edge_exists(edge)) {
                data->snpNetwork->add_edge(edge, "NET_SHUFFLE");
                ++num_successful;
            } else ++num_failed;

            Logger::logProgress("Shuffle iteration " + std::to_string(num_successful) + " of " + std::to_string(num_iterations) + " (" + std::to_string(num_failed) + " failed)");
        }
        Logger::logLine("Shuffle iteration " + std::to_string(num_successful) + " of " + std::to_string(num_iterations) + " (" + std::to_string(num_failed) + " failed)");
    }

    ShuffleNetwork::ShuffleNetwork(const std::string& algorithm) {
        if (algorithm == "TOPOLOGY_PRESERVING_WITH_SNP_DEGREE") {
            shuffle_algorithm = TOPOLOGY_PRESERVING_WITH_SNP_DEGREE;
        } else if (algorithm == "TOPOLOGY_PRESERVING_WITHOUT_SNP_DEGREE") {
            shuffle_algorithm = TOPOLOGY_PRESERVING_WITHOUT_SNP_DEGREE;
        } else if (algorithm == "EXPECTED_DEGREE_KEEP_DEGREE_DISTRIBUTION") {
            shuffle_algorithm = EXPECTED_DEGREE_KEEP_DEGREE_DISTRIBUTION;
        } else if (algorithm == "EXPECTED_DEGREE_KEEP_INDIVIDUAL_DEGREE") {
            shuffle_algorithm = EXPECTED_DEGREE_KEEP_INDIVIDUAL_DEGREE;
        } else {
            throw epi::Error("Unknown shuffle algorithm " + algorithm);
        }
    }

    rapidjson::Value ShuffleNetwork::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("ShuffleNetwork"), doc.GetAllocator());
        obj.AddMember("algorithm", rapidjson::Value().SetString(shuffle_algorithm_names[shuffle_algorithm].c_str(), shuffle_algorithm_names[shuffle_algorithm].size(), doc.GetAllocator()), doc.GetAllocator());
        return obj;
    }
} // epi