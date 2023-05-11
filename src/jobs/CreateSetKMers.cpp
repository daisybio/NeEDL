//
// Created by juli on 21.10.22.
//

#include "CreateSetKMers.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {

    long CreateSetKMers::binom_coeff (long n, long k) {
        if (k == 0 || k == n) return 1;
        return binom_coeff(n - 1, k - 1) + binom_coeff(n - 1, k);
    }

    int CreateSetKMers::permute(int *vec, int k, int max) {
        if (++vec[k] > max) {
            vec[k] --;

            if (k == 0) return 0;

            int min = permute(vec, k - 1, vec[k] - 1);
            vec[k] = min;
        }
        return vec[k] + 1;
    }

    CreateSetKMers::CreateSetKMers(size_t k_min, size_t k_max) {
        if (k_min == 0 || k_max == 0 || k_min > k_max) {
            throw epi::Error("k_min and k_max must not be null and k_min must be <= k_max.");
        }

        this->k_min = k_min;
        this->k_max = k_max;
    }

    CreateSetKMers::CreateSetKMers(size_t k) {
        if (k == 0) {
            throw epi::Error("k must not be null.");
        }

        this->k_min = this->k_max = k;
    }

    void CreateSetKMers::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("Begin creating SNP-set k-mers");

        std::vector<SNPSet> unique_sets (data->snpSetStorage.begin(), data->snpSetStorage.end());
        data->snpSetStorage.clear();
        for (auto & set : unique_sets) {
            data->snpSetStorage.push_back(set);

            if (set.size() >= k_min) {
                for(int k = k_min; k <= std::min(k_max, set.size()); k++) {
                    // calculate number of possible combinations
                    long num_combinations = binom_coeff(set.size(), k);

                    auto *current = new int[k];
                    for(int x = 0; x < k; x++) current[x] = x;
                    for (long combination = 0; combination < num_combinations; combination ++) {
                        // create subset
                        SNPSet subset;
                        for(int x = 0; x < k; x++) {
                            subset += set.at(current[x]);
                        }
                        data->snpSetStorage.push_back(subset);

                        // shift new combination
                        permute(current, k - 1, set.size() - 1);
                    }

                    delete[] current;
                }
            }
        }

        logger.stop();
    }

    rapidjson::Value CreateSetKMers::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("CreateSetKMers"), doc.GetAllocator());
        obj.AddMember("k_min", rapidjson::Value().SetUint64(k_min), doc.GetAllocator());
        obj.AddMember("k_max", rapidjson::Value().SetUint64(k_max), doc.GetAllocator());
        return obj;
    }
} // epi