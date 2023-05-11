//
// Created by juli on 21.10.22.
//

#include "CreateRandomSets.hpp"

namespace epi {

    void CreateRandomSets::run(std::shared_ptr<DataModel> data) {
        // count the distribution of snp set lengths and generate random snps with same distribution

        std::vector<unsigned int> count_list;

        for (auto &snp_set: data->snpSetStorage) {
            size_t num_snps = snp_set.size();

            if (num_snps > count_list.size()) {
                // need to extend vector
                for (size_t i = count_list.size(); i <= num_snps; i++) count_list.push_back(0);
            }

            count_list[num_snps - 1]++;
        }

        // print counts
        Logger::logLine("SNP-set lengths found:");
        for (unsigned int i = 0; i < count_list.size(); i++) {
            Logger::logLine("   [" + std::to_string(i + 1) + "]: " +
                            Logger::to_string(float(count_list[i]) / float(data->snpSetStorage.size())) + " (" +
                            std::to_string(count_list[i]) + ")");
        }

        size_t total_snp_sets = data->snpSetStorage.size();

        std::uniform_int_distribution<size_t> distr_size(1, total_snp_sets);
        std::uniform_int_distribution<size_t> distr_snps(0, data->snpStorage->num_snps() - 1);

        // generate random sets
        std::vector<unsigned int> generated_count(count_list.size());
        data->snpSetStorage.clear();
        auto thread_num = omp_get_thread_num();
        for (size_t i = 0; i < this->num_sets; i++) {
            // get number of SNPs
            size_t num_snps = 0;
            size_t pos = distr_size(data->random_device[thread_num]);
            size_t curr = 0;
            for (size_t j = 0; j < count_list.size(); j++) {
                curr += count_list[j];
                if (pos <= curr) {
                    num_snps = j + 1;
                    break;
                }
            }

            generated_count[num_snps - 1]++;

            // generate random snp set of that size
            SNPSet snp_set;
            for (size_t j = 0; j < num_snps; j++) {
                snp_set += data->snpStorage->by_instance_id(distr_snps(data->random_device[thread_num]));
            }

            data->snpSetStorage.push_back(snp_set);
        }

        // print counts
        Logger::logLine("SNP-set lengths generated:");
        for (unsigned int i = 0; i < count_list.size(); i++) {
            Logger::logLine("   [" + std::to_string(i + 1) + "]: " +
                            Logger::to_string(float(generated_count[i]) / float(data->snpSetStorage.size())) + " (" +
                            std::to_string(generated_count[i]) + ")");
        }
    }

    CreateRandomSets::CreateRandomSets(size_t num_sets) {
        this->num_sets = num_sets;
    }

    rapidjson::Value CreateRandomSets::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("CreateRandomSets"), doc.GetAllocator());
        obj.AddMember("num_sets", rapidjson::Value().SetUint64(num_sets), doc.GetAllocator());
        return obj;
    }
} // epi