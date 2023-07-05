//
// Created by juli on 19.05.22.
//

#include "bayesian_model.hpp"

namespace epi {


    template<>
    void
    BayesianModel<QuantitativePhenoType>::
    initialize_() {
        mu_ = 0.0;
        for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
            mu_ += this->instance_->phenotype(ind);
        }
        mu_ /= static_cast<double>(this->instance_->num_inds());
        sigma_ = 0.0;
        for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
            sigma_ += (this->instance_->phenotype(ind) - mu_) * (this->instance_->phenotype(ind) - mu_);
        }
        sigma_ /= static_cast<double>(this->instance_->num_inds());
    }

    template<>
    void
    BayesianModel<CategoricalPhenoType>::
    initialize_() {
        num_bins_= instance_->num_categories();
    }

    template<>
    void
    BayesianModel<QuantitativePhenoType>::
    initialize_binned_phenotypes_(std::vector<CategoricalPhenoType> & binned_phenotypes, const Eigen::MatrixXd &residuals) const {
        CategoricalPhenoType bin;
        for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
            // bin = static_cast<CategoricalPhenoType>(std::floor(misc::normal_cdf(this->instance_->phenotype(ind), mu_, sigma_) * static_cast<double>(num_bins_)));
            // bin = static_cast<CategoricalPhenoType>(std::floor(misc::normal_cdf(residuals(ind, 0), mu_, sigma_) * static_cast<double>(num_bins_)));
            bin = static_cast<CategoricalPhenoType>(std::floor(misc::normal_cdf((this->get_cov_status() ? residuals(ind, 0) : this->instance_->phenotype(ind)), mu_, sigma_) * static_cast<double>(num_bins_)));
            if (bin == num_bins_) {
                bin--;
            }
            binned_phenotypes.emplace_back(bin);
        }
    }

    template<>
    void
    BayesianModel<CategoricalPhenoType>::
    initialize_binned_phenotypes_(std::vector<CategoricalPhenoType> & binned_phenotypes, const Eigen::MatrixXd &residuals) const {
        for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
            binned_phenotypes.emplace_back((this->get_cov_status() ? static_cast<CategoricalPhenoType>(std::round(residuals(ind, 0))) : this->instance_->phenotype(ind)));
            // binned_phenotypes.emplace_back(static_cast<CategoricalPhenoType>(std::round(residuals(ind, 0))));
        }
    }

    template<>
    bool
    BayesianModel<CategoricalPhenoType>::
    parse_option_(const std::string & option, const std::string & arg) {
        if (option == "num-bins") {
            return true;
        }
        return false;
    }

    template<>
    bool
    BayesianModel<QuantitativePhenoType>::
    parse_option_(const std::string & option, const std::string & arg) {
        if (option == "num-bins") {
            try {
                num_bins_ = std::stoul(arg);
            }
            catch (...) {
                throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to int greater equal 2>] [...]\"");
            }
            if (num_bins_ < 2) {
                throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to int greater equal 2>] [...]\"");
            }
            return true;
        }
        return false;
    }

    template<>
    void
    BayesianModel<QuantitativePhenoType>::
    set_default_options_() {
        num_bins_ = 100;
    }

    template<>
    void
    BayesianModel<CategoricalPhenoType>::
    set_default_options_() {
        num_bins_ = this->instance_->num_categories();
    }
}