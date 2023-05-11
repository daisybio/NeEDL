//
// Created by juli on 21.06.22.
//

#include "SeedingRandomConnected.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    SeedingRandomConnected::SeedingRandomConnected(uint_fast32_t num_seeds) {
        this->num_seeds = num_seeds;
    }

    void SeedingRandomConnected::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("seeding with method RANDOM_CONNECTED");
        if (data->snpNetwork == nullptr) {
            throw epi::Error("A SNP-SNP-interaction network is needed in order to use RANDOM_CONNECTED_SEED.");
        }

        data->snpSetStorage.clear();

        std::unordered_set<SNP_t, SNP_t::SNPHash> used_snps;
        auto all_network_snps = data->snpNetwork->get_network_snps();

        while (!all_network_snps.empty() && data->snpSetStorage.size() < num_seeds) {
            std::uniform_int_distribution<size_t> dist1(0, all_network_snps.size() - 1);
            size_t snp1_pos = dist1(data->random_device[omp_get_thread_num()]);
            SNP_t snp1 = all_network_snps[snp1_pos];

            // it is either already used or we will use it now
            all_network_snps.erase(all_network_snps.begin() + snp1_pos);

            if (used_snps.find(snp1) == used_snps.end()) {
                // this snp was not used yet --> let's see if it has any adjacent SNP that we can use
                auto adjacent_snps = data->snpNetwork->get_adjacent_snps(snp1);
                while (!adjacent_snps.empty()) {
                    std::uniform_int_distribution<size_t> dist2(0, adjacent_snps.size() - 1);
                    size_t snp2_pos = dist2(data->random_device[omp_get_thread_num()]);
                    SNP_t snp2 = adjacent_snps[snp2_pos];

                    // it is either used or we will use it now
                    adjacent_snps.erase(adjacent_snps.begin() + snp2_pos);

                    if (used_snps.find(snp2) == used_snps.end()) {
                        // we have a pair where both snps are new
                        SNPSet set ({ snp1, snp2 });
                        set.set_attribute("SEED_ORIGIN", "RANDOM_CONNECTED");
                        data->snpSetStorage.push_back(set);
                        used_snps.insert(snp1);
                        used_snps.insert(snp2);
                        break;
                    }
                }
            }
        }

        logger.stop();
    }

    rapidjson::Value SeedingRandomConnected::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("SeedingRandomConnected"), doc.GetAllocator());
        obj.AddMember("num_seeds", rapidjson::Value().SetUint64(num_seeds), doc.GetAllocator());
        return obj;
    }
} // epi