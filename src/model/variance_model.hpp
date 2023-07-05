/*******************************************************************************
 *                                                                             *
 *   Copyright (C) 2020 by David B. Blumenthal                                 *
 *                                                                             *
 *   This file is part of GenEpiSeeker.                                        *
 *                                                                             *
 *   GenEpiSeeker is free software: you can redistribute it and/or modify it   *
 *   under the terms of the GNU General Public License as published by         *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   GenEpiSeeker is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with GenEpiSeeker. If not, see <http://www.gnu.org/licenses/>.      *
 *                                                                             *
 ******************************************************************************/


/*!
 * @file  variance_model.hpp
 * @brief Declaration of epi::VarianceModel.
 */

#ifndef SRC_MODEL_VARIANCE_MODEL_HPP_
#define SRC_MODEL_VARIANCE_MODEL_HPP_

#include "epistasis_model.hpp"

namespace epi {

/*!
 * @brief Epistasis models based on the variance of the phenotypes within and between the cells of the SNP set's penetrance table.
 * @details Computes the p-value of the ANOVA test (for quantitative phenotypes) or of the chi-squared test (for categorical phenotypes).
 *
 * Implements and extends epistasis models suggested in:
 * - Yupeng Wang, Xinyu Liu, Kelly Robbins, and Romdhane Rekaya:
 *   &ldquo;AntEpiSeeker: detecting epistatic interactions for case-control studies using a two-stage ant colony optimization algorithm&rdquo;,
 *   https://doi.org/10.1186/1756-0500-3-117
 *
 * Does not support any options.
 */
    template<class PhenoType>
    class VarianceModel : public EpistasisModel<PhenoType> {

    public:

        VarianceModel(Instance<PhenoType> *instance);

        void cov_activate();

        void cov_deactivate();

        virtual ~VarianceModel();

    private:

        double global_mean_;

        bool incl_cov_;

        std::vector<std::size_t> num_inds_in_categories_;

        std::size_t num_non_empty_categories_;

        virtual void initialize_() final;

        virtual options::ModelSense model_sense_() const final;

        virtual bool is_predictive_() const final;

        virtual double evaluate_(const std::vector<SNP> &snp_set, bool prepare_prediction) final;

        void
        compute_test_statistic_(const std::vector<std::vector<PhenoType>> &penetrance_table, double &test_statistic,
                                std::size_t &num_groups) const;

        double p_value_(double test_statistic, std::size_t num_groups) const;

        bool get_cov_status() const;

    };

}

#include "variance_model.ipp"

#ifdef HEADER_ONLY

#include "variance_model.cpp"

#endif


#endif /* SRC_MODEL_VARIANCE_MODEL_HPP_ */
