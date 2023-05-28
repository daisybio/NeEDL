//
// Created by juli on 21.10.22.
//

#include "SNP_t.hpp"
#include "SNPStorage.hpp"

/**
 * @brief methods in this source file were separated to fix the order of include when compiling with option HEADER_ONLY.
 */
namespace epi {
    std::string SNPSet::get_snp_string() const {
        std::string out;

        std::vector<std::string> sorted_snp_names;
        for (auto & snp : *snps) {
            sorted_snp_names.push_back(SNPStorage::currentSnpStorage->snp_get_name(snp));
        }

        bool is_first = true;
        for(auto & name : sorted_snp_names) {
            if (!is_first) out += ';';
            else is_first = false;

            out += name;
        }

        return out;
    }

    double SNPSet::calculate_score(options::EpistasisScore score) {
        // check if score was already calculated
        auto model_index = uint_fast32_t(score);
        if ((scores_calculated & (1 << model_index)) == 0) {
            // need to calculate score
            if (scores == nullptr) scores = std::make_shared<std::vector<double>>(NUM_EPISTASIS_SCORES);
            scores->at(model_index) = SNPStorage::currentSnpStorage->calculate_score(*this, score);
            scores_calculated |= 1 << model_index;
        }

        return scores->at(model_index);
    }

    std::string SNPEdge::get_snp_string() const {
        if (bits.snp1 == 0 && bits.snp2 == 0) return "{ INVALID EDGE }";
        return "{ " + SNPStorage::currentSnpStorage->snp_get_name(SNP_t(bits.snp1)) + ", " + SNPStorage::currentSnpStorage->snp_get_name(SNP_t(bits.snp2)) + " }";
    }

    std::vector<std::string> SNPSet::get_annotations() {
        std::set<std::string> all_annotations;
        for (auto & snp : *snps) {
            auto annotations = SNPStorage::currentSnpStorage->snp_get_annotations(snp);
            all_annotations.insert(annotations.begin(), annotations.end());
        }

        std::vector<std::string> all_annotations_vec { all_annotations.begin(), all_annotations.end() };

        return all_annotations_vec;
    }

}