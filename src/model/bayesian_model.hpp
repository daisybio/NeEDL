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
 * @file  bayesian_model.hpp
 * @brief Declaration of epi::BayesianModel.
 */

#ifndef SRC_MODEL_BAYESIAN_MODEL_HPP_
#define SRC_MODEL_BAYESIAN_MODEL_HPP_

#include "epistasis_model.hpp"

namespace epi {

/*!
 * @brief Epistasis models based on Bayesian networks.
 * @details Computes the p-value of the ANOVA test (for quantitative phenotypes) or of the chi-squared test (for categorical phenotypes).
 *
 * Implements and extends epistasis models suggested in:
 * - Peng-Jie Jin and Hong-Bin Shen:
 *   &ldquo;MACOED: a multi-objective ant colony optimization algorithm for SNP epistasis detection in genome-wide association studies&rdquo;,
 *   https://doi.org/10.1093/bioinformatics/btu702
 * - Xia Jiang, Richard E. Neapolitan, M. Michael Barmada, and Shyam Visweswaran:
 *   &ldquo;Learning genetic epistasis using Bayesian network scoring criteria&rdquo;,
 *   https://doi.org/10.1186/1471-2105-12-89
 *
 * Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--num-bins @<convertible to int greater equal 2@></tt> | the number of bins used for discretizing quantitative phenotypes | @p 100 | has no effect unless @p PhenoType equals epi::QuantitativePhenoType |
 */
    template<class PhenoType>
    class BayesianModel : public EpistasisModel<PhenoType> {

    public:

        BayesianModel(Instance<PhenoType> * instance);

        void cov_activate();

        void cov_deactivate();

        virtual ~BayesianModel();

    private:

        double mu_;

        double sigma_;

        bool incl_cov_;

        std::size_t num_bins_;

        virtual void initialize_() final;

        virtual options::ModelSense model_sense_() const final;

        virtual bool is_predictive_() const final;

        virtual double evaluate_(const std::vector<SNP> & snp_set, bool prepare_prediction) final;

        virtual bool parse_option_(const std::string & option, const std::string & arg) final;

        virtual void set_default_options_() final;

        virtual std::string valid_options_() const final;

        void initialize_binned_phenotypes_(std::vector<CategoricalPhenoType> & binned_phenotypes,
                                           const Eigen::MatrixXd &residuals) const;

        double log_factorial_(std::size_t n);

        bool get_cov_status() const;

        bool cov_score_() const;
    };

}

#include "bayesian_model.ipp"

#ifdef HEADER_ONLY
#include "bayesian_model.cpp"
#endif

#endif /* SRC_MODEL_BAYESIAN_MODEL_HPP_ */
