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
 * @file  penetrance_model.hpp
 * @brief Declaration of epi::PenetranceModel.
 */

#ifndef SRC_MODEL_PENETRANCE_MODEL_HPP_
#define SRC_MODEL_PENETRANCE_MODEL_HPP_

#include "epistasis_model.hpp"

namespace epi {

/*!
 * @brief Epistasis models based on local maximum likelihood models for the cells of the SNP set's penetrance table.
 * @details For all cells of the penetrance table with sufficiently many entries, local maximum likelihood distributions
 * are computed (Normal distributions for quantitative phenotypes, categorical distributions for categorical phenotypes).
 * For cells with too few entries, the global maximum likelihood distribution is employed. Subsequently, a score based
 * on the likelihood of the joint model is computed.
 *
 * Implements the epistasis model suggested in:
 * - David B. Blumenthal et al.:
 *   &ldquo;Modeling epistatic interaction: what do we want to optimize?&rdquo;,
 *   work in progress
 *
 * Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--score LLH\|NLL\|AIC\|BIC</tt> | the returned score | @p NLL | @p LHH ~ likelihood, @p NLL ~ negative log-likelihood, @p AIC ~ Aikake information criterion, @p BIC ~ Bayesian information criterion |
 * | <tt>\--min-cell-size  @<convertible to int greater 0@></tt> | the minimal size a cell of the penetrance table must have to compute local distributions | 10 for quantitative phenotypes, 10 * number of categories for categorical phenotypes | for cells with fewer entries, the global distribution is used |
 */
    template<class PhenoType>
    class PenetranceModel : public EpistasisModel<PhenoType> {

    public:

        PenetranceModel(Instance<PhenoType> * instance);

        const std::string & get_score();

        void cov_activate();

        void cov_deactivate();

        virtual ~PenetranceModel();

    private:

        std::size_t min_cell_size_;

        std::string score_;

        bool incl_cov_;

        Eigen::MatrixXd beta_;

        std::vector<double> global_distribution_;

        std::vector<SNP> snp_set_;

        std::vector<std::vector<double>> maximum_likelihood_distributions_;

        virtual options::ModelSense model_sense_() const final;

        virtual bool is_predictive_() const final;

        virtual double evaluate_(const std::vector<SNP> & snp_set, bool prepare_prediction) final;

        virtual PhenoType predict_(Ind ind) const final;

        virtual void initialize_() final;

        virtual bool parse_option_(const std::string & option, const std::string & arg) final;

        virtual void set_default_options_() final;

        virtual std::string valid_options_() const final;

        virtual void save_(const std::string & filename) const final;

        virtual void load_(const std::string & filename) final;

        std::vector<double> maximum_likelihood_distribution_(const std::vector<PhenoType> & phenotypes) const;

        double likelihood_(PhenoType phenotype, const std::vector<double> & distribution) const;

        std::string meta_info_() const;

        void check_loaded_model_() const;

        bool get_cov_status() const;

        bool cov_score_() const;
    };

}

#include "penetrance_model.ipp"

#ifdef HEADER_ONLY
#include "penetrance_model.cpp"
#endif


#endif /* SRC_MODEL_PENETRANCE_MODEL_HPP_ */
