//
// Created by juli on 17.06.22.
//

#include "MMAFilter.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    MMAFilter::MMAFilter(double filter_cutoff, bool use_BH_correction) {
        this->filter_cutoff = filter_cutoff;
        this->use_BH_correction = use_BH_correction;
    }

    void MMAFilter::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("MMA filter");
        if (data->snpStorage == nullptr) {
            throw epi::Error("MMA filter needs a loaded instance to work on!");
        }

        std::vector<double> mma_scores(data->snpStorage->num_snps());
        auto iterable = data->snpStorage->all();

#pragma omp parallel for default(none) shared(iterable, data, mma_scores)
        for(size_t i = 0; i < mma_scores.size(); i++) {
            auto snp = *(iterable.begin() + i);
            mma_scores[i] = data->snpStorage->calculate_score({snp}, options::EpistasisScore::VARIANCE);
        }

        if (use_BH_correction) {
            // correct p values
            Logger::logLine("Correct MMA p-values with Benjamini-Hochberg");

            // clone vector with mma values
            std::vector<std::pair<double, size_t>> mma_values_copy;
            for(size_t i = 0; i < mma_scores.size(); i++) {
                mma_values_copy.emplace_back( mma_scores[i], i);
            }

            // sort pairs in ascending order by p value
            std::sort(mma_values_copy.begin(), mma_values_copy.end(), [] (std::pair<double, size_t> a, std::pair<double, size_t> b) { return a.first < b.first; });

            // correct p-value
            for (size_t i = 0; i < mma_values_copy.size(); i++) {
                mma_values_copy[i].first *= double(mma_values_copy.size()) / double(i + 1);
            }

            // check that lower p-value cannot lead to higher corrected one comparing the others
            double min = mma_values_copy.back().first;
            for (long i = mma_values_copy.size() - 1; i > 0; i--) {
                if (mma_values_copy[i].first > min) {
                    mma_values_copy[i].first = min;
                } else {
                    min = mma_values_copy[i].first;
                }
            }

            // sort back adjusted p-values
            for (auto x : mma_values_copy) {
                // std::cout << maximum_marginal_association_for_snps_[x.second] << " -> " << x.first << std::endl;
                mma_scores[x.second] = x.first;
            }
        }

        // store mma_scores and mark as removed if threshold is not reached
        auto score_it = mma_scores.begin();
        size_t num_removed = 0;
        for (const auto &snp : data->snpStorage->all()) {
            data->snpStorage->set_snp_mma(snp, *score_it);
            if (*score_it <= filter_cutoff) {
                data->snpStorage->set_removed(snp, true);
                ++num_removed;
            }
            score_it++;
        }

        Logger::logLine("MMA filter removed " + std::to_string(num_removed) + " SNPs from the dataset.");

        logger.stop();
    }

    rapidjson::Value MMAFilter::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("MMAFilter"), doc.GetAllocator());
        obj.AddMember("filter_cutoff", rapidjson::Value().SetDouble(filter_cutoff), doc.GetAllocator());
        obj.AddMember("use_BH_correction", rapidjson::Value().SetBool(use_BH_correction), doc.GetAllocator());
        return obj;
    }
} // epi