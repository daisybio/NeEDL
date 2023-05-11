//
// Created by juli on 09.08.22.
//

#include "MultiNetworkAggregator.hpp"
#include <algorithm>

namespace epi {
    MultiNetworkAggregator::MultiNetworkAggregator(std::vector<std::shared_ptr<Job>> job_list, std::vector<std::string> names) {
        this->job_list = std::move(job_list);
        this->names = std::move(names);
        if (this->job_list.size() != this->names.size()) throw epi::Error("Unequal number of multi-network jobs and names");
    }

    void MultiNetworkAggregator::add(const std::shared_ptr<Job> &job, std::string name) {
        this->job_list.push_back(job);
        this->names.push_back(name);
    }

    void MultiNetworkAggregator::add(std::vector<std::shared_ptr<Job>> job_list_, std::vector<std::string> names_) {
        this->job_list.insert(this->job_list.end(), job_list_.begin(), job_list_.end());
        this->names.insert((this->names.end()), names_.begin(), names_.end());
        if (this->job_list.size() != this->names.size()) throw epi::Error("Unequal number of multi-network jobs and names");
    }


    void MultiNetworkAggregator::run(std::shared_ptr<DataModel> data) {
        SNPNetwork initial_network = *data->snpNetwork;
        std::vector<SNPSet> initial_result_set(data->snpSetStorage.begin(), data->snpSetStorage.end());

        size_t job_i = 0;
        for (auto &job: this->job_list) {
            Logger::logLine(
                    "Multi-Network pipeline " + names[job_i] + " " + std::to_string(job_i + 1) + " of " +
                    std::to_string(this->job_list.size()));

            job_i++;

            // reset initial state
            *data->snpNetwork = initial_network;
            data->snpSetStorage = initial_result_set;

            job->run(data);

            if (data->snpNetwork == nullptr) {
                throw epi::Error("Pipeline need to create an SNP-network in order to do multi-network aggregation.");
            }

            // extract and store result data
            result_sets.emplace_back(data->snpSetStorage.begin(), data->snpSetStorage.end());
            std::unordered_map<SNP_t, std::vector<SNP_t>, SNP_t::SNPHash> adjacency;
            for (auto &set: data->snpSetStorage) {
                std::vector<SNP_t> ordered_set(set.begin(), set.end());
                std::sort(ordered_set.begin(), ordered_set.end());

                for (auto &snp: set) {
                    if (adjacency.find(snp) == adjacency.end()) {
                        auto snp_adjacency = data->snpNetwork->get_adjacent_snps(snp);
                        std::sort(snp_adjacency.begin(), snp_adjacency.end());

                        std::vector<SNP_t> set_adjacency;
                        std::set_intersection(ordered_set.begin(), ordered_set.end(),
                                              snp_adjacency.begin(), snp_adjacency.end(),
                                              std::back_inserter(set_adjacency));

                        adjacency[snp] = set_adjacency;
                    }
                }
            }
            adjacency_data.push_back(adjacency);
        }

        Logger::logLine("Construct new network from the preliminary results");
        // clear data from previous runs
        data->snpSetStorage.clear();
        data->snpNetwork->clear();

        for (size_t i = 0; i < result_sets.size(); i++) {
            for (auto &set: result_sets[i]) {
                // insert vertices
                data->snpNetwork->add_nodes(set.begin(), set.end());

                // insert edges
                for (auto &snp: set) {
                    data->snpStorage->snp_set_or_add_variable_attribute(snp, "ms_source", names[i], ';');

                    auto &adjacent_snps = adjacency_data[i][snp];
                    std::for_each(adjacent_snps.begin(), adjacent_snps.end(), [&data, &snp, this, i](const SNP_t &s) {
                        data->snpNetwork->add_edge({s, snp}, names[i]);
                    });
                }
            }
        }
    }

    size_t MultiNetworkAggregator::num_networks() {
        return job_list.size();
    }

    rapidjson::Value MultiNetworkAggregator::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("MultiNetworkAggregator"), doc.GetAllocator());

        rapidjson::Value job_li(rapidjson::kArrayType);
        for (size_t i = 0; i < job_list.size(); ++i) {
            rapidjson::Value job_obj(rapidjson::kObjectType);
            job_obj.AddMember("name", rapidjson::Value().SetString(names[i].c_str(), names[i].size(), doc.GetAllocator()), doc.GetAllocator());
            job_obj.AddMember("config", job_list[i]->getConfig(doc), doc.GetAllocator());

            job_li.PushBack(job_obj, doc.GetAllocator());
        }
        obj.AddMember("jobs", job_li, doc.GetAllocator());

        return obj;
    }
} // epi
