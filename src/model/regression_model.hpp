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
 * @file regression_model.hpp
 * @brief Declaration of epi::RegressionModel.
 */

#ifndef SRC_MODEL_REGRESSION_MODEL_HPP_
#define SRC_MODEL_REGRESSION_MODEL_HPP_


#include "epistasis_model.hpp"

namespace epi {

/*!
 * @brief Epistasis models based on linear regression with quadratic features.
 * @details For all individuals, feature vectors based on their genotypes at the given SNP set are computed.
 * The feature vector contains feature for all SNPs contained in the SNP set as well as for all of their
 * quadratic interactions. Subsequently, a regression model (linear regression for quantitative phenotypes,
 * logistic regression for dichotomous categorical phenotypes, and multinomial regression for non-dichotomous
 * categorical phenotyoes) is computed using block gradient descent. Finally, scores are computed based on the
 * likelihood of the computed model or based on its gain w.r.t. an additive model without features for quadratic interaction.
 *
 * Implements and extends epistasis models suggested in:
 * - Peng-Jie Jin and Hong-Bin Shen:
 *   &ldquo;MACOED: a multi-objective ant colony optimization algorithm for SNP epistasis detection in genome-wide association studies&rdquo;,
 *   https://doi.org/10.1093/bioinformatics/btu702
 * - Xiang Wan, Can Yang, Qiang Yang, Hong Xue, Xiaodan Fan, Nelson L. S. Tang, and Weichuan Yu:
 *   &ldquo;BOOST: A fast approach to detecting gene-gene interactions in genome-wide case-control studies&rdquo;,
 *   https://doi.org/10.1016/j.ajhg.2010.07.021
 *
 * Supports the following options:
 * | <tt>\--@<option@> @<arg@></tt> | modified parameter | default  | more information |
 * | ------------------------------ | ------------------ | -------- | ---------------- |
 * | <tt>\--score LLH\|NLL\|AIC\|BIC\|LLH-GAIN\|NLL-GAIN\|AIC-GAIN\|BIC-GAIN</tt> | the returned score | @p NLL | @p LHH ~ likelihood, @p NLL ~ negative log-likelihood, @p AIC ~ Aikake information criterion, @p BIC ~ Bayesian information criterion, @p *-GAIN ~ improvement w.r.t. additive model |
 * | <tt>\--max-itrs  @<convertible to int greater equal 0@></tt> | the maximal number of iterations of block gradient descent | 500 | if 0, no maximal number of iterations is enforced |
 * | <tt>\--learning-rate  @<convertible to double greater 0@></tt> | the learning rate employed by the block gradient descent | 0.1 | n.a. |
 * | <tt>\--epsilon @<convertible to double greater 0@></tt> | the convergence threshold employed by the block gradient descent | 0.01 | n.a. |
 */
template<class PhenoType>
class RegressionModel : public EpistasisModel<PhenoType> {

public:

	RegressionModel(Instance<PhenoType> * instance);

    const std::string & get_score();

	virtual ~RegressionModel();

private:

	std::string score_;

	double learning_rate_;

	double epsilon_;

	double max_itrs_;

	std::vector<SNP> snp_set_;

	Eigen::MatrixXd interaction_model_;

	virtual options::ModelSense model_sense_() const final;

	virtual bool is_predictive_() const final;

	virtual double evaluate_(const std::vector<SNP> & snp_set, bool prepare_prediction) final;

	virtual PhenoType predict_(Ind ind) const final;

	virtual bool parse_option_(const std::string & option, const std::string & arg) final;

	virtual void set_default_options_() final;

	virtual std::string valid_options_() const final;

	virtual void save_(const std::string & filename) const final;

	virtual void load_(const std::string & filename) final;

	void construct_feature_matrix_(const std::vector<SNP> & snp_set, bool construct_interaction_features, Eigen::MatrixXd & feature_matrix) const;

	void construct_features_(const std::vector<SNP> & snp_set, Ind ind, bool construct_interaction_features, Eigen::MatrixXd & feature_matrix) const;

	void fit_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis, double & sigma) const;

	void fit_linear_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis) const;

	void compute_linear_model_sigma_(const Eigen::MatrixXd & hypothesis, double & sigma) const;

	void fit_logistic_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis) const;

	void fit_multinomial_model_(const Eigen::MatrixXd & feature_matrix, Eigen::MatrixXd & model, Eigen::MatrixXd & hypothesis) const;

	double likelihood_(Ind ind, const Eigen::MatrixXd & hypothesis, double sigma) const;

	bool gain_score_() const;

	std::size_t degrees_of_freedom_per_feature_() const;

	std::string meta_info_() const;

	void check_loaded_model_() const;

};

}

#include "regression_model.ipp"

#ifdef HEADER_ONLY
#include "regression_model.cpp"
#endif



#endif /* SRC_MODEL_REGRESSION_MODEL_HPP_ */
