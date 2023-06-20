//
// Created by juli on 19.05.22.
//

#include "variance_model.hpp"

namespace epi {
    template<>
    void
    VarianceModel<QuantitativePhenoType>::
    initialize_() {
        global_mean_ = 0.0;
        for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
            global_mean_ += this->instance_->phenotype(ind);
        }
        global_mean_ /= static_cast<double>(this->instance_->num_inds());
    }

    template<>
    void
    VarianceModel<CategoricalPhenoType>::
    initialize_() {
        num_inds_in_categories_ = std::vector<std::size_t>(this->instance_->num_categories(), 0);
        for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
            num_inds_in_categories_.at(this->instance_->phenotype(ind))++;
        }
        num_non_empty_categories_ = 0;
        for (std::size_t num_inds_in_category: num_inds_in_categories_) {
            if (num_inds_in_category > 0) {
                num_non_empty_categories_++;
            }
        }
    }

    template<>
    void
    VarianceModel<QuantitativePhenoType>::
    compute_test_statistic_(const std::vector<std::vector<QuantitativePhenoType>> &penetrance_table,
                            double &test_statistic, std::size_t &num_groups) const {

        // Compute number of groups and the group means.
        num_groups = 0;
        std::vector<double> group_means;
        for (const auto &cell: penetrance_table) {
            if (cell.empty()) {
                group_means.emplace_back(undefined_double());
                continue;
            }
            num_groups++;
            group_means.emplace_back(0.0);
            for (const auto &phenotype: cell) {
                group_means.back() += phenotype;
            }
            group_means.back() /= static_cast<double>(cell.size());
        }

        // Compute sums of squares within and between groups.
        double ss_between_groups{0.0};
        double ss_within_groups{0.0};
        std::vector<std::vector<QuantitativePhenoType>>::const_iterator cell;
        std::vector<double>::const_iterator group_mean;
        for (cell = penetrance_table.begin(), group_mean = group_means.begin();
             cell != penetrance_table.end() and group_mean != group_means.end(); cell++, group_mean++) {
            if (cell->empty()) {
                continue;
            }
            ss_between_groups +=
                    static_cast<double>(cell->size()) * (*group_mean - global_mean_) * (*group_mean - global_mean_);
            for (const double &phenotype: *cell) {
                ss_within_groups += (phenotype - *group_mean) * (phenotype - *group_mean);
            }
        }

        // Compute F test statistic.
        if (ss_between_groups <= 0.0 or num_groups == 1) {
            test_statistic = 0.0;
        } else if (ss_within_groups <= 0.0 or this->instance_->num_inds() == num_groups) {
            test_statistic = std::numeric_limits<double>::infinity();
        } else {
            test_statistic = (ss_between_groups / static_cast<double>(num_groups - 1)) /
                             (ss_within_groups / static_cast<double>(this->instance_->num_inds() - num_groups));
        }
    }

    template<>
    void
    VarianceModel<CategoricalPhenoType>::
    compute_test_statistic_(const std::vector<std::vector<CategoricalPhenoType>> &penetrance_table,
                            double &test_statistic, std::size_t &num_groups) const {

        // Compute chi-squared test statistic and the number of groups.
        num_groups = 0;
        test_statistic = 0.0;
        std::vector<std::size_t> num_cell_inds_in_categories(this->instance_->num_categories(), 0);
        double expected_cell_inds_in_category{0.0};
        double deviation_from_expectation{0.0};
        for (const auto &cell: penetrance_table) {
            if (cell.empty()) {
                continue;
            }
            num_groups++;
            for (auto &count: num_cell_inds_in_categories) {
                count = 0;
            }
            for (const auto &phenotype: cell) {
                num_cell_inds_in_categories.at(phenotype)++;
            }
            for (CategoricalPhenoType phenotype{0};
                 static_cast<std::size_t>(phenotype) < this->instance_->num_categories(); phenotype++) {
                if (num_inds_in_categories_.at(phenotype) == 0) {
                    continue;
                }
                expected_cell_inds_in_category =
                        static_cast<double>(num_inds_in_categories_.at(phenotype) * cell.size()) /
                        static_cast<double>(this->instance_->num_inds());
                deviation_from_expectation =
                        static_cast<double>(num_cell_inds_in_categories.at(phenotype)) - expected_cell_inds_in_category;
                test_statistic +=
                        (deviation_from_expectation * deviation_from_expectation) / expected_cell_inds_in_category;
            }
        }

        // Correct possible numerical errors accumulated during summation.
        if (num_groups == 1) {
            test_statistic = 0.0;
        }
        if (num_non_empty_categories_ == 1) {
            test_statistic = 0.0;
        }
    }

    template<>
    double
    VarianceModel<QuantitativePhenoType>::
    p_value_(double test_statistic, std::size_t num_groups) const {
        if (test_statistic <= 0.0) {
            return 1.0 - std::numeric_limits<double>::min();
        }
        if (test_statistic == std::numeric_limits<double>::infinity()) {
            return std::numeric_limits<double>::min();
        }
        double return_value{
                1.0 - misc::f_cdf(test_statistic, num_groups - 1, this->instance_->num_inds() - num_groups)};
        return misc::ensure_valid_continuous_probability(return_value);
    }

    template<>
    double
    VarianceModel<CategoricalPhenoType>::
    p_value_(double test_statistic, std::size_t num_groups) const {
        if (test_statistic <= 0.0) {
            return 1.0 - std::numeric_limits<double>::min();
        }
        double return_value{
                1.0 - misc::chi_squared_cdf(test_statistic, (num_groups - 1) * (num_non_empty_categories_ - 1))};
        return misc::ensure_valid_continuous_probability(return_value);
    }
}