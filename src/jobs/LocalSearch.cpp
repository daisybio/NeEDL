//
// Created by juli on 06.07.22.
//

#include "LocalSearch.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/LDTester.hpp"

namespace epi {
    void LocalSearch::run(std::shared_ptr<DataModel> data) {
        if (ld_check == CUTOFF_MODE) {
            ld_tester = std::make_shared<LDTester>(ld_file, ld_mode, data, ld_cutoff);
        } else if (ld_check == MC_MODE) {
            ld_tester = std::make_shared<LDTester>(ld_file, ld_mode, data, ld_mc_min_set, ld_mc_max_set, ld_mc_sample_size);
        }

        std::shared_ptr<TimeLogger> logger;
        search_start_time = std::chrono::high_resolution_clock::now();
        if (search_time_limit > 0.) {
            logger = std::make_shared<TimeLogger>("local search with time limit " + Logger::to_string(search_time_limit) + " min");
        } else {
            logger = std::make_shared<TimeLogger>("local search without time limit");
        }

        // statistics
        //long model_values_cached = 0, model_values_seen = 0;
        size_t finished_seeds = 0;
        size_t total_seeds = data->snpSetStorage.size();
        size_t num_skipped = 0;

        // for score development over time
        bool need_minimizing = data->snpStorage->need_minimizing_score(model);
        double current_best_score = std::numeric_limits<double>::max();
        if (!need_minimizing) current_best_score = -current_best_score;
        for (auto &seed : data->snpSetStorage) {
            if (need_minimizing) {
                current_best_score = std::min(current_best_score, seed.calculate_score(model));
            } else {
                current_best_score = std::max(current_best_score, seed.calculate_score(model));
            }
        }
        std::vector<std::pair<double, long>> best_score_timepoints;
        best_score_timepoints.emplace_back(current_best_score, 0);

        // resulting seeds
        std::vector<SNPSet> result_sets (total_seeds);
        std::chrono::high_resolution_clock::time_point searchStart = std::chrono::high_resolution_clock::now();

#pragma omp parallel for default(none) shared(data, need_minimizing, searchStart, finished_seeds, total_seeds, num_skipped, current_best_score, best_score_timepoints, result_sets) schedule(dynamic)
        for (std::size_t i = 0; i < data->snpSetStorage.size(); i++) {
            bool allow_search = true;

            if (search_time_limit > 0.) {
                std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
                double current_duration =
                        std::chrono::duration_cast<std::chrono::nanoseconds>(now - search_start_time).count() /
                        60000000000.0f;
                allow_search = current_duration < search_time_limit;
            }

            SNPSet result_set;
            bool result_retrieved = false;

            if (allow_search) {
                result_set = process_start_seed(data->snpSetStorage[i], data);
                result_retrieved = true;

                if (calculate_monte_carlo) {
                    result_set.set_attribute("MONTE_CARLO_SCORE", Logger::to_string(
                            data->snpStorage->calculate_monte_carlo_score(result_set, model,
                                                                          monte_carlo_num_permutations)));
                }
            }


#pragma omp critical
            {
                if (result_retrieved && result_set.size() > 0) {
                    // model_values_cached += current_cache->get_values_cached();
                    // model_values_seen += current_cache->get_values_seen();
                    result_sets[i] = result_set;

                    double set_score = result_set.calculate_score(model);
                    if ((need_minimizing && set_score < current_best_score) ||
                        (!need_minimizing && set_score > current_best_score)) {
                        current_best_score = set_score;
                        std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
                        long time_point = std::chrono::duration_cast<std::chrono::microseconds>(
                                now - searchStart).count();
                        best_score_timepoints.emplace_back(set_score, time_point);
                    }
                } else {
                    num_skipped++;
                }

                finished_seeds++;
                std::string line = "Local search running: " + Logger::to_string(finished_seeds) + " of " +
                                   Logger::to_string(total_seeds) + " start points done. Current best score is " +
                                   Logger::to_string(current_best_score);
                if (finished_seeds < total_seeds) Logger::logProgress(line);
                else Logger::logLine(line);
            }

        }


        if (search_time_limit > 0.) {
            Logger::logLine("Runs skipped due to exceeded search time limit: " + Logger::to_string(num_skipped));
        }

        // ModelCache<PhenoType>::print_usage(model_values_seen, model_values_cached);

        if (output_scores_over_time) {
            // write best score every minute to file
            std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
            long end = std::chrono::duration_cast<std::chrono::milliseconds>(now - searchStart).count();

            if (data->outputDirectory == nullptr) {
                throw epi::Error("Output directory need to be specificed to save scores over time data.");
            }

            std::ofstream time_file (data->outputDirectory->get_output_directory() + scores_over_time_name + ".csv");
            time_file << "time (ms)\tscore\n";

            for (auto & time_point : best_score_timepoints) {
                auto time_str = std::to_string(time_point.second);
                while (time_str.size() < 4) time_str = '0' + time_str;
                time_str = time_str.substr(0, time_str.size() - 3) + '.' + time_str.substr(time_str.size() - 3);
                time_file << time_str << '\t' << Logger::to_string(time_point.first) << '\n';
            }

            // size_t scores_list_pos = 0;
            // for (long t = 0; t < end; t += 60000) {
            //     long nextTimestamp = best_score_timepoints[scores_list_pos].second;
            //     double score = best_score_timepoints[scores_list_pos].first;
            //     while (t > nextTimestamp) {
            //         score = best_score_timepoints[scores_list_pos].first;
            //         if (scores_list_pos >= best_score_timepoints.size() - 1) break;
            //         scores_list_pos++;
            //         nextTimestamp = best_score_timepoints[scores_list_pos].second;
            //     }
            //     if (scores_list_pos > 0) scores_list_pos --;

            //     time_file << long(t / 60000) << '\t' << score << '\n';
            // }
        }

        // remove empty results
        std::vector<SNPSet> all_sets (result_sets.begin(), result_sets.end());
        result_sets.clear();
        for (auto & set : all_sets) {
            if (set.size() > 0) result_sets.push_back(set);
        }

        if (collapse_identical_results) {
            // make distinct set
            std::unordered_map<SNPSet, std::vector<SNPSet>, SNPSetHash> results_map;

            for (auto & res : result_sets) {
                auto val = results_map.find(res);
                if (val == results_map.end()) {
                    // value not found --> insert it
                    results_map.insert({ res, std::vector<SNPSet>({ res }) });
                } else {
                    // add information to already existing element
                    val->second.push_back(res);
                }
            }

            result_sets.clear();
            for (auto & group : results_map) {
                auto merged_set = group.first.clone();
                merged_set.clear_attributes();

                merged_set.set_attribute("NUM_MERGED", Logger::to_string(group.second.size()));

                // collect all attribute keys
                std::set<std::string> attrib_keys;
                for (auto & res : group.second) {
                    auto keys = res.get_attribute_keys();
                    attrib_keys.insert(keys.begin(), keys.end());
                }

                // collect all attribute values
                for (auto & key : attrib_keys) {
                    bool is_int = true;
                    std::vector<long long> int_vals;
                    bool is_float = true;
                    std::vector<double> float_vals;

                    std::vector<std::string> str_vals;

                    for (auto & set : group.second) {
                        auto val = set.get_attribute(key);

                        // check if it is int
                        if (is_int) {
                            char *pos;
                            long long v = strtoll(val.c_str(), &pos, 10);
                            if (pos == 0) int_vals.push_back(v);
                            else is_int = false;
                        }

                        if (is_float) {
                            char *pos;
                            double v = strtod(val.c_str(), &pos);
                            if (pos == 0) float_vals.push_back(v);
                            else is_float = false;
                        }

                        str_vals.push_back(val);
                    }

                    if (is_int) {
                        long long sum = 0, min = int_vals[0], max = int_vals[0];
                        std::string concatenated;
                        bool first = true;
                        for (auto & val : int_vals) {
                            sum += val;
                            if (first) first = false;
                            else concatenated += ';';
                            concatenated += Logger::to_string(val);

                            min = std::min(val, min);
                            max = std::max(val, max);
                        }

                        std::set<long long> vals_distinct (int_vals.begin(), int_vals.end());
                        first = true;
                        std::string concatenated_distinct;
                        for (auto & val : vals_distinct) {
                            if (first) first = false;
                            else concatenated_distinct += ';';
                            concatenated_distinct += Logger::to_string(val);
                        }

                        merged_set.set_attribute(key + "_AVG", Logger::to_string(double(sum) / double(int_vals.size())));
                        merged_set.set_attribute(key + "_MIN", Logger::to_string(min));
                        merged_set.set_attribute(key + "_MAX", Logger::to_string(max));
                        merged_set.set_attribute(key + "_DISTINCT", concatenated_distinct);
                        merged_set.set_attribute(key + "_ALL", concatenated);
                    } else if (is_float) {
                        double sum = 0, min = float_vals[0], max = float_vals[0];
                        std::string concatenated;
                        bool first = true;
                        for (auto & val : float_vals) {
                            sum += val;
                            if (first) first = false;
                            else concatenated += ';';
                            concatenated += Logger::to_string(val);

                            min = std::min(val, min);
                            max = std::max(val, max);
                        }

                        std::set<double> vals_distinct (float_vals.begin(), float_vals.end());
                        first = true;
                        std::string concatenated_distinct;
                        for (auto & val : vals_distinct) {
                            if (first) first = false;
                            else concatenated_distinct += ';';
                            concatenated_distinct += Logger::to_string(val);
                        }

                        merged_set.set_attribute(key + "_AVG", Logger::to_string(double(sum) / double(int_vals.size())));
                        merged_set.set_attribute(key + "_MIN", Logger::to_string(min));
                        merged_set.set_attribute(key + "_MAX", Logger::to_string(max));
                        merged_set.set_attribute(key + "_DISTINCT", concatenated_distinct);
                        merged_set.set_attribute(key + "_ALL", concatenated);
                    } else {
                        std::string concatenated;
                        bool first = true;
                        for (auto & val : str_vals) {
                            if (first) first = false;
                            else concatenated += ';';
                            concatenated += val;
                        }

                        std::set<std::string> vals_distinct (str_vals.begin(), str_vals.end());
                        first = true;
                        std::string concatenated_distinct;
                        for (auto & val : vals_distinct) {
                            if (first) first = false;
                            else concatenated_distinct += ';';
                            concatenated_distinct += val;
                        }


                        merged_set.set_attribute(key + "_DISTINCT", concatenated_distinct);
                        merged_set.set_attribute(key + "_ALL", concatenated);
                    }
                }

                result_sets.push_back(merged_set);
            }
        }

        data->snpSetStorage = result_sets;

        logger->stop();
    }

    SNPSet LocalSearch::process_start_seed(const SNPSet& start_seed, const std::shared_ptr<DataModel> &data) {
        std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

        bool need_minimizing = data->snpStorage->need_minimizing_score(model);

        enum {
            SIMULATED_ANNEALING,
            CONVERGENCE,
            SEED_TIMEOUT,
            SEARCH_TIMEOUT
        } stopping_reason = CONVERGENCE;

        SNPSet result = start_seed.clone();

        double delta_sum = 0;
        double temperature = 1.0 / std::log(annealing_start_prob);


        // add adjacent snps until min_set criterion is fulfilled or no adjacent nodes are found anymore
        while (result.size() < min_set) {
            std::vector<SNP_t> adjacent_snps = data->snpNetwork->get_adjacent_snps(result.vec());
            if (adjacent_snps.empty()) break;

            // add a random snp to the set
            std::uniform_int_distribution<size_t> dist(0, adjacent_snps.size() - 1);
            result += adjacent_snps[dist(data->random_device[omp_get_thread_num()])];
        }

        if (result.size() < min_set) {
            // return empty set with information
            SNPSet empty;
            empty.set_attribute("STOPPING_REASON", "MIN_SET CRITERION VIOLATED");
            empty.set_attribute("NUM_ROUNDS", std::to_string(0));
            Logger::logLine("min_set criterion was violated with start seed " + start_seed.get_snp_string());
            return empty;
        }

        // create subnetwork
        SNPNetwork local_subnetwork (true);
        local_subnetwork.add_nodes(result.begin(), result.end());
        // add start seeds to result set
        // result += {{ start_seed.begin(), start_seed.end() }};

        // fully connect all nodes in subnetwork
        std::vector<SNPEdge> start_edges;
        for (size_t i = 0; i < result.size(); i++) {
            for (size_t j = i + 1; j < result.size(); j++) {
                start_edges.emplace_back(result[i], result[j]);
            }
        }
        local_subnetwork.add_edges(start_edges.begin(), start_edges.end());

        int iterations_without_improvement = 0;
        SNPSet best_result = result.clone();

        size_t current_round = 1;
        for (; current_round <= max_rounds; current_round++) {
            // check time limits
            std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();

            if (search_time_limit > 0.) {
                // search time limit
                double current_duration =
                        std::chrono::duration_cast<std::chrono::nanoseconds>(now - search_start_time).count() /
                        60000000000.0;
                if (current_duration >= search_time_limit) {
                    stopping_reason = SEARCH_TIMEOUT;
                    break;
                }
            }
            if (time_per_start_seed > 0.) {
                // per seed timeout
                double current_duration =
                        std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_time).count() / 60000000000.0;
                if (current_duration >= time_per_start_seed) {
                    stopping_reason = SEED_TIMEOUT;
                    break;
                }
            }

            // current score
            double current_score = result.calculate_score(model);
            double previous_score = current_score;
            SNPSet current_set = result.clone();
            SNPSet second_set;
            double second_score = 0.; // this will be set as soon as first set is selected
            enum {
                DO_NOTHING,
                ADD,
                DELETE,
                SUBSTITUTE
            } first_change_option = DO_NOTHING, second_change_option = DO_NOTHING;
            SNP_t first_change_add_snp, first_change_delete_snp, second_change_add_snp, second_change_delete_snp;

            // get articulation points
            std::vector<SNP_t> set_items (result.begin(), result.end());
            std::sort(set_items.begin(), set_items.end());
            auto articulation_points = local_subnetwork.get_articulation_points();
            std::sort(articulation_points.begin(), articulation_points.end());

            // get only snps that are not articulation points and in current result set
            std::vector<SNP_t> allowed_to_delete;
            std::set_difference(set_items.begin(), set_items.end(),
                                articulation_points.begin(), articulation_points.end(),
                                std::back_inserter(allowed_to_delete));


            // test all possible add options
            if (result.size() < max_set) {
                // get add candidates
                std::vector<SNP_t> adjacent_snps = data->snpNetwork->get_adjacent_snps(result.vec());

                for (auto & snp : adjacent_snps) {
                    // check for LD
                    if (ld_check != DISABLED && ld_tester->test(result, snp)) continue;

                    // create new set with additional snp
                    SNPSet new_set = result + snp;
                    double new_score = new_set.calculate_score(model);
                    if ((need_minimizing && new_score < current_score) || (!need_minimizing && new_score > current_score)) {
                        current_set = new_set;
                        current_score = new_score;
                        first_change_option = ADD;
                        first_change_add_snp = snp;
                    } else if (second_change_option == DO_NOTHING || (need_minimizing && new_score < second_score) || (!need_minimizing && new_score > second_score)) {
                        second_set = new_set;
                        second_score = new_score;
                        second_change_option = ADD;
                        second_change_add_snp = snp;
                    }
                }
            }

            //test all delete-options
            if (result.size() > min_set) {
                for (auto &snp: allowed_to_delete) {
                    // create new set by removing snp
                    SNPSet new_set = result - snp;
                    double new_score = new_set.calculate_score(model);
                    if ((need_minimizing && new_score < current_score) ||
                        (!need_minimizing && new_score > current_score)) {
                        current_set = new_set;
                        current_score = new_score;
                        first_change_option = DELETE;
                        first_change_delete_snp = snp;
                    } else if (second_change_option == DO_NOTHING || (need_minimizing && new_score < second_score) ||
                               (!need_minimizing && new_score > second_score)) {
                        second_set = new_set;
                        second_score = new_score;
                        second_change_option = DELETE;
                        second_change_delete_snp = snp;
                    }
                }
            }

            //test all substitution-options
            if (!allowed_to_delete.empty()) {
                size_t num_combinations_checked = 0;
                for (auto & delete_snp : allowed_to_delete) {
                    // delete snp from set
                    SNPSet new_set_del = result - delete_snp;

                    // get add candidates
                    std::vector<SNP_t> adjacent_snps = data->snpNetwork->get_adjacent_snps(new_set_del.vec());

                    for (auto &add_snp: adjacent_snps) {
                        num_combinations_checked ++;

                        // check for LD
                        if (ld_check != DISABLED && ld_tester->test(new_set_del, add_snp)) continue;

                        SNPSet new_set = new_set_del + add_snp;
                        double new_score = new_set.calculate_score(model);

                        if ((need_minimizing && new_score < current_score) || (!need_minimizing && new_score > current_score)) {
                            current_set = new_set;
                            current_score = new_score;
                            first_change_option = SUBSTITUTE;
                            first_change_add_snp = add_snp;
                            first_change_delete_snp = delete_snp;

                        } else if (second_change_option == DO_NOTHING || (need_minimizing && new_score < second_score) ||
                                   (!need_minimizing && new_score > second_score)) {
                            second_set = new_set;
                            second_score = new_score;
                            second_change_option = SUBSTITUTE;
                            second_change_add_snp = add_snp;
                            second_change_delete_snp = delete_snp;
                        }
                    }
                }
            }

            // check whether we have at least one next turn to do
            if (first_change_option == DO_NOTHING) {
                if (second_change_option == DO_NOTHING) {
                    // out of options
                    break;
                } else {
                    // update delta
                    double delta = fabs(second_score - previous_score);
                    delta_sum += delta;

                    if (get_annealing_decision(data, current_round, second_score, previous_score, delta_sum, temperature, iterations_without_improvement)) {
                        // set second option to first option
                        first_change_option = second_change_option;
                        first_change_add_snp = second_change_add_snp;
                        first_change_delete_snp = second_change_delete_snp;
                        current_set = second_set;
                        current_score = second_score;
                    } else {
                        stopping_reason = SIMULATED_ANNEALING;
                        break;
                    }
                }

                iterations_without_improvement ++;
            } else {
                // update delta
                double delta = fabs(current_score - previous_score);
                delta_sum += delta;
            }

            // update subnetwork
            if (first_change_option == DELETE || first_change_option == SUBSTITUTE) {
                // delete node from network
                local_subnetwork.remove_node(first_change_delete_snp);
            }
            if (first_change_option == ADD || first_change_option == SUBSTITUTE) {
                // add new node to the network
                local_subnetwork.add_node(first_change_add_snp);

                // add all edges from the main network
                auto connected_snps = data->snpNetwork->get_adjacent_snps(first_change_add_snp);
                std::vector<SNPEdge> edges;
                for (auto &adj_snp : connected_snps) {
                    if (local_subnetwork.contains_node(adj_snp)) edges.emplace_back(first_change_add_snp, adj_snp);
                }

                local_subnetwork.add_edges(edges.begin(), edges.end());
            }

            result = current_set;
            if ((need_minimizing && current_score < best_result.calculate_score(model)) || (!need_minimizing && current_score > best_result.calculate_score(model))) {
                best_result = result;
            }

            temperature *= cooling_factor;
        }

        best_result.set_attribute("NUM_ROUNDS", std::to_string(current_round));
        switch (stopping_reason) {
            case CONVERGENCE:
                best_result.set_attribute("STOPPING_REASON", "CONVERGENCE");
                break;
            case SEARCH_TIMEOUT:
                best_result.set_attribute("STOPPING_REASON", "SEARCH_TIMEOUT");
                break;
            case SEED_TIMEOUT:
                best_result.set_attribute("STOPPING_REASON", "SEED_TIMEOUT");
                break;
            case SIMULATED_ANNEALING:
                best_result.set_attribute("STOPPING_REASON", "SIMULATED_ANNEALING");
                break;
        }

        return best_result;
    }


    bool LocalSearch::get_annealing_decision(const std::shared_ptr<DataModel> &data, int rounds, double score_now, double score_before, double delta_sum, double temperature, int iterations_without_improvement) {
        bool do_it{false};

        if (annealing_type == RANDOM_ANNEALING) {
            std::uniform_real_distribution<double> distr(0, 1);
            //random annealing
            int do_it_val = distr(data->random_device[omp_get_thread_num()]);
            if (do_it_val > .5 and rounds < max_rounds - 1) {
                do_it = true;
            }

        } else if (annealing_type == HYPERBOLIC_TAN_ANNEALING) {
            std::uniform_real_distribution<double> distr(0, 1);
            //Simulated annealing
            //gets bigger if comes to end
            double score_anneal = (score_now - score_before) / (max_rounds - rounds);
            //hyperbolic tan -> strictly increasing f'(0) = 1, f(0) = 0 and f(eternity)=1
            double score_normalized = 1 - (2 / (exp(2 * score_anneal) + 1));

            double r = distr(data->random_device[omp_get_thread_num()]);

            if (r >= score_normalized) {
                do_it = true;
            }
        } else if (annealing_type == SIMULATED_ANNEALING) {
            std::uniform_real_distribution<double> distr(0, 1);

            double delta{std::fabs(score_now - score_before)};
            double delta_avg{delta_sum / rounds};
            double random_number = distr(data->random_device[omp_get_thread_num()]);
            double condition = std::exp((-delta) / (delta_avg * temperature)) - 1;
            double div_iterations_no_improvement =
                    static_cast<double>(iterations_without_improvement) / static_cast<double>(rounds);

            if ((random_number > condition and delta != 0) or
                (random_number < div_iterations_no_improvement and delta != 0)) {
                do_it = true;
            }
        }

        return do_it;
    }

    LocalSearch::LocalSearch(std::string epistasis_model,
                             bool collapse_identical_results,
                             size_t max_rounds,
                             double search_time_limit_minutes,
                             double per_seed_time_limit_minutes,
                             std::string annealing_type,
                             double cooling_factor,
                             double annealing_start_prob,
                             double annealing_end_prob,
                             bool output_score_development,
                             std::string score_development_file_name,
                             size_t min_set_size,
                             size_t max_set_size,
                             bool calculate_monte_carlo_pvalues,
                             size_t monte_carlo_permutations) {


        this->model = options::epistasis_score_from_string(epistasis_model);

        this->collapse_identical_results = collapse_identical_results;
        this->max_rounds = max_rounds;
        this->search_time_limit = search_time_limit_minutes;
        this->time_per_start_seed = per_seed_time_limit_minutes;

        if (annealing_type == "RANDOM_ANNEALING") {
            this->annealing_type = RANDOM_ANNEALING;
        } else if (annealing_type == "HYPERBOLIC_TAN_ANNEALING") {
            this->annealing_type = HYPERBOLIC_TAN_ANNEALING;
        } else if (annealing_type == "SIMULATED_ANNEALING") {
            this->annealing_type = SIMULATED_ANNEALING;
        } else {
            throw epi::Error("Unknown annealing_type " + annealing_type);
        }
        annealing_type_str = annealing_type;


        this->cooling_factor = cooling_factor;
        this->annealing_start_prob = annealing_start_prob;
        this->annealing_end_prob = annealing_end_prob;
        this->output_scores_over_time = output_score_development;
        this->scores_over_time_name = score_development_file_name;
        this->min_set = min_set_size;
        this->max_set = max_set_size;
        this->calculate_monte_carlo = calculate_monte_carlo_pvalues;
        this->monte_carlo_num_permutations = monte_carlo_permutations;

        // set the cooling factor if we have a max round limit
        double temperature = 1.0 / std::log(annealing_start_prob);
        if (max_rounds > 1) {
            this->cooling_factor = std::pow((1.0 / std::log(annealing_end_prob)) / temperature,
                                      1.0 / static_cast<double>(max_rounds - 1));
        }
    }

    std::string LocalSearch::get_model_name() {
        return options::epistasis_score_to_string(model);
    }

    void LocalSearch::activate_LD_check(std::string ld_file_, std::string ld_mode_, double cutoff) {
        ld_file = ld_file_;
        ld_mode = ld_mode_;
        ld_cutoff = cutoff;
        ld_check = CUTOFF_MODE;
    }

    void
    LocalSearch::activate_LD_check(std::string ld_file_, std::string ld_mode_, size_t mc_min_set, size_t mc_max_set, size_t mc_sample_size) {
        ld_file = ld_file_;
        ld_mode = ld_mode_;
        ld_mc_min_set = mc_min_set;
        ld_mc_max_set = mc_max_set;
        ld_mc_sample_size = mc_sample_size;
        ld_check = MC_MODE;
    }

    rapidjson::Value LocalSearch::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("LocalSearch"), doc.GetAllocator());
        auto epiModelStr = options::epistasis_score_to_string(model);
        obj.AddMember("model", rapidjson::Value().SetString(epiModelStr.c_str(), epiModelStr.size(), doc.GetAllocator()), doc.GetAllocator());
        obj.AddMember("collapse_identical_results", rapidjson::Value().SetBool(collapse_identical_results), doc.GetAllocator());
        obj.AddMember("max_rounds", rapidjson::Value().SetUint64(max_rounds), doc.GetAllocator());
        obj.AddMember("search_time_limit_minutes", rapidjson::Value().SetDouble(search_time_limit), doc.GetAllocator());
        obj.AddMember("per_seed_search_time_limit_minutes", rapidjson::Value().SetDouble(time_per_start_seed), doc.GetAllocator());
        obj.AddMember("annealing_type", rapidjson::Value().SetString(annealing_type_str.c_str(), annealing_type_str.size(), doc.GetAllocator()), doc.GetAllocator());
        obj.AddMember("annealing_cooling_factor", rapidjson::Value().SetDouble(cooling_factor), doc.GetAllocator());
        obj.AddMember("annealing_start_prob", rapidjson::Value().SetDouble(annealing_start_prob), doc.GetAllocator());
        obj.AddMember("annealing_end_prob", rapidjson::Value().SetDouble(annealing_end_prob), doc.GetAllocator());
        obj.AddMember("create_score_development_file", rapidjson::Value().SetBool(output_scores_over_time), doc.GetAllocator());
        if (output_scores_over_time) {
            obj.AddMember("score_development_file_name",
                          rapidjson::Value().SetString(scores_over_time_name.c_str(), scores_over_time_name.size(), doc.GetAllocator()),
                          doc.GetAllocator());
        }
        obj.AddMember("min_set_size", rapidjson::Value().SetUint64(min_set), doc.GetAllocator());
        obj.AddMember("max_set_size", rapidjson::Value().SetUint64(max_set), doc.GetAllocator());
        obj.AddMember("calculate_monte_carlo_pvalues", rapidjson::Value().SetBool(calculate_monte_carlo), doc.GetAllocator());
        if (calculate_monte_carlo) {
            obj.AddMember("monte_carlo_num_permutations", rapidjson::Value().SetUint64(monte_carlo_num_permutations), doc.GetAllocator());
        }

        obj.AddMember("ld_check", rapidjson::Value().SetString(ld_check_str[ld_check].c_str(), ld_check_str[ld_check].size(), doc.GetAllocator()), doc.GetAllocator());
        if (ld_check != DISABLED) {
            obj.AddMember("ld_file", rapidjson::Value().SetString(ld_file.c_str(), ld_file.size(), doc.GetAllocator()), doc.GetAllocator());
            obj.AddMember("ld_mode", rapidjson::Value().SetString(ld_mode.c_str(), ld_mode.size(), doc.GetAllocator()), doc.GetAllocator());
            if (ld_check == CUTOFF_MODE) {
                obj.AddMember("ld_cutoff", rapidjson::Value().SetDouble(ld_cutoff), doc.GetAllocator());
            } else if (ld_check == MC_MODE) {
                obj.AddMember("ld_mc_min_set", rapidjson::Value().SetUint64(ld_mc_min_set), doc.GetAllocator());
                obj.AddMember("ld_mc_max_set", rapidjson::Value().SetUint64(ld_mc_max_set), doc.GetAllocator());
                obj.AddMember("ld_mc_sample_size", rapidjson::Value().SetUint64(ld_mc_sample_size), doc.GetAllocator());
            }
        }

        return obj;
    }
} // epi