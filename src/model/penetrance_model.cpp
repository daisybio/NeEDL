//
// Created by juli on 19.05.22.
//

#include "penetrance_model.hpp"

namespace epi {
    template<>
    PenetranceModel<QuantitativePhenoType>::
    PenetranceModel(Instance<QuantitativePhenoType> * instance):
            EpistasisModel<QuantitativePhenoType>(instance),
            min_cell_size_{10},
            score_("NLL"),
            global_distribution_(),
            snp_set_(),
            maximum_likelihood_distributions_() {}

    template<>
    PenetranceModel<CategoricalPhenoType>::
    PenetranceModel(Instance<CategoricalPhenoType> * instance):
            EpistasisModel<CategoricalPhenoType>(instance),
            min_cell_size_{10 * this->instance_->num_categories()},
            score_("NLL"),
            global_distribution_(),
            snp_set_(),
            maximum_likelihood_distributions_() {}

    template<>
    QuantitativePhenoType
    PenetranceModel<QuantitativePhenoType>::
    predict_(Ind ind) const {
        return maximum_likelihood_distributions_.at(this->instance_->genotype_at_snp_set(snp_set_, ind)).at(0);
    }

    template<>
    CategoricalPhenoType
    PenetranceModel<CategoricalPhenoType>::
    predict_(Ind ind) const {
        CategoricalPhenoType predicted_phenotype{0};
        CategoricalPhenoType current_phenotype{0};
        double max_probability{0};
        for (double probability : maximum_likelihood_distributions_.at(this->instance_->genotype_at_snp_set(snp_set_, ind))) {
            if (probability > max_probability) {
                predicted_phenotype = current_phenotype;
                max_probability = probability;
            }
            current_phenotype++;
        }
        return predicted_phenotype;
    }

    template<>
    void
    PenetranceModel<QuantitativePhenoType>::
    set_default_options_() {
        min_cell_size_ = 10;
        score_ = "NLL";
    }

    template<>
    void
    PenetranceModel<CategoricalPhenoType>::
    set_default_options_() {
        min_cell_size_ = 10 * this->instance_->num_categories();
        score_ = "NLL";
    }

    template<>
    std::vector<double>
    PenetranceModel<QuantitativePhenoType>::
    maximum_likelihood_distribution_(const std::vector<QuantitativePhenoType> & phenotypes) const {
        double mu{std::accumulate(phenotypes.begin(), phenotypes.end(), 0.0)};
        mu /= static_cast<double>(phenotypes.size());
        double sigma{0.0};
        for (QuantitativePhenoType phenotype : phenotypes) {
            sigma += (phenotype - mu) * (phenotype - mu);
        }
        sigma /= static_cast<double>(phenotypes.size());
        sigma = std::sqrt(sigma);
        if (sigma == 0) {
            sigma = 0.0001;
        }
        return std::vector<double>{mu, sigma};
    }

    template<>
    std::vector<double>
    PenetranceModel<CategoricalPhenoType>::
    maximum_likelihood_distribution_(const std::vector<CategoricalPhenoType> & phenotypes) const {
        std::vector<double> distribution(this->instance_->num_categories(), 0.0);
        for (CategoricalPhenoType phenotype : phenotypes) {
            distribution.at(phenotype) += 1.0;
        }
        double sample_size{static_cast<double>(phenotypes.size())};
        for (double & probability : distribution) {
            probability /= sample_size;
        }
        std::size_t size_support{0};
        for (double & probability : distribution) {
            if (probability > 0) {
                size_support++;
            }
        }
        if (size_support < distribution.size()) {
            for (double & probability : distribution) {
                if (probability > 0) {
                    probability -= 0.00001 *  static_cast<double>(distribution.size() - size_support)/ static_cast<double>(size_support);
                }
                else {
                    probability = 0.00001;
                }
            }
        }
        return distribution;
    }

    template<>
    double
    PenetranceModel<QuantitativePhenoType>::
    likelihood_(QuantitativePhenoType phenotype, const std::vector<double> & distribution) const {
        return misc::normal_pdf(phenotype, distribution.at(0), distribution.at(1));
    }

    template<>
    double
    PenetranceModel<CategoricalPhenoType>::
    likelihood_(CategoricalPhenoType phenotype, const std::vector<double> & distribution) const {
        return misc::ensure_valid_continuous_probability(distribution.at(phenotype));
    }

    template<>
    std::string
    PenetranceModel<QuantitativePhenoType>::
    meta_info_() const {
        return "phenotype = quantitative\n";
    }

    template<>
    std::string
    PenetranceModel<CategoricalPhenoType>::
    meta_info_() const {
        return std::string("phenotype = categorical\ncategories = ") + std::to_string(this->instance_->num_categories()) + "\n";
    }

    template<>
    void
    PenetranceModel<QuantitativePhenoType>::
    check_loaded_model_() const {
        if (maximum_likelihood_distributions_.back().size() != 2) {
            throw Error(std::string("Wrong number of parameters per genotype in loaded model. Expected: 2. Actual: ") + std::to_string(maximum_likelihood_distributions_.back().size()) + ".");
        }
    }

    template<>
    void
    PenetranceModel<CategoricalPhenoType>::
    check_loaded_model_() const {
        if (maximum_likelihood_distributions_.back().size() != this->instance_->num_categories()) {
            throw Error(std::string("Wrong number of parameters per genotype in loaded model. Expected: ") + std::to_string(this->instance_->num_categories())+ ". Actual: " + std::to_string(maximum_likelihood_distributions_.back().size()) + ".");
        }
    }

}