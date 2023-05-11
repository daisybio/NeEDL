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
 * @file  epistasis_model.hpp
 * @brief Declaration of epi::EpistasisModel.
 */

#ifndef SRC_MODEL_EPISTASIS_MODEL_HPP_
#define SRC_MODEL_EPISTASIS_MODEL_HPP_

#include "../util/types.hpp"
#include "../util/misc.hpp"
#include "instance.hpp"

namespace epi {

/*!
 * @brief Abstract class modeling objective functions for epistasis detection.
 * @tparam PhenoType The type of the phenotypes. Use epi::CategoricalPhenoType for categorical phenotypes and epi::QuantitativePhenoType for quantitative phenotypes.
 */
template<class PhenoType>
class EpistasisModel {

public:

	/*!
	 * @brief Constructs epistasis model for given instance.
	 * @param[in] instance Pointer to the instance on which the objective should be evaluated.
	 */
	EpistasisModel(Instance<PhenoType> * instance);

	/*!
	 * @brief Destructor.
	 * @note Must be implemented by all derived classes.
	 */
	virtual ~EpistasisModel() = 0;

	/*!
	 * @brief Set the instance on which the objective should be evaluated.
	 * @param[in] instance Pointer to the instance on which the objective should be evaluated.
	 * @note Invalidates calls to initialize().
	 */
	void set_instance(Instance<PhenoType> * instance);

	/*!
	 * @brief Initializes the objective, i.e., pre-computes values shared by each evaluation.
	 */
	void initialize();

	/*!
	 * @brief Sets the options of the objective.
	 * @param[in] options String of the form <tt>[\--@<option@> @<arg@>] [...]</tt>, where @p option contains neither spaces nor single quotes,
	 * and @p arg contains neither spaces nor single quotes or is of the form <tt>'[\--@<sub-option@> @<sub-arg@>] [...]'</tt>,
	 * where both @p sub-option and @p sub-arg contain neither spaces nor single quotes.
	 * @note Invalidates calls to initialize().
	 */
	void set_options(const std::string & options);

	/*!
	 * @brief Returns the model sense.
	 * @return Returns epi::options::ModelSense::MINIMIZE if the objective should be minimized and epi::options::ModelSense::MAXIMIZE if it should be maximized.
	 */
	options::ModelSense model_sense() const;

	/*!
	 * @brief Check whether the model is predictive.
	 * @return Boolean @p true if trained models can be used to predict phenotypes.
	 */
	bool is_predictive() const;

	/*!
	 * @brief Returns @p true if @p objective_value_1 is better than @p objective_value_2.
	 * @param[in] objective_value_1 First objective value.
	 * @param[in] objective_value_2 Second objective value.
	 * @return Boolean @p true if @p objective_value_1 is better than @p objective_value_2.
	 */
	bool is_better(double objective_value_1, double objective_value_2) const;

	/*!
	 * @brief Evaluates the objective for a given SNP set.
	 * @param[in] snp_set SNP set for which the objective should be evaluated.
	 * @return The objective of @p snp_set on the instance selected by the last call to set_instance().
	 */
	double evaluate(const std::vector<SNP> & snp_set);

	/*!
	 * @brief Evaluates the objective for a given SNP set and tracks the time of doing so.
	 * @param[in] snp_set SNP set for which the objective should be evaluated.
	 * @return The objective of @p snp_set on the instance selected by the last call to set_instance().
	 */
	double evaluate_track_time(const std::vector<SNP> & snp_set);

	/*!
	 * @brief The runtime of the last tracked evaluation of the objective.
	 * @return The runtime of the last call to evaluate_track_time().
	 */
	double evaluation_time() const;

	/*!
	 * @brief Returns the Monte Carlo p-value for the objective of a given SNP set.
	 * @details To compute the Monte Carlo p-value of the SNP set @p snp_set, the phenotypes of the individuals contained in
	 * the instance are permuted @p num_permutations times. For each permutation as well as for the original phenotypes, the objective
	 * is then evaluated at @p snp_set. Let @f$N_+@f$ be the number of permutations for which the obtained objective was better than
	 * the objective obtained for the original phenotypes. Then, the Monte Carlo p-value of @p snp_set is defined as
	 * @f$(N_+ + 1)/(\mathtt{num\_permutations} + 1)@f$.
	 * @param[in] snp_set The SNP set for which the Monte Carlo p-value should be computed.
	 * @param[in] num_permutations The number of permutations of the phenotypes used to compute the Monte Carlo p-value.
	 * @param[in] seed The seed used to generate the permutations. If set to undefined(), the seed it randomly generated.
	 * @return The Monte Carlo p-value of SNP set @p snp_set.
	 */
	double monte_carlo_p_value(const std::vector<SNP> & snp_set, std::size_t num_permutations, std::size_t seed = undefined_uint());

	/*!
	 * @brief Saves the model for the SNP set passed to the last call to evaluate(), evaluate_track_time(), or monte_carlo_p_value().
	 * @param[in] filename Path to the file where the model should be saved.
	 */
	void save(const std::string & filename) const;

	/*!
	 * @brief Loads a model from disc. The model must match the dimensions of the epistasis instance.
	 * @param[in] filename Path to the file where the model is stored.
	 */
	void load(const std::string & filename);

	/*!
	 * @brief Predicts the phenotype of a given individual.
	 * @param[in] ind The individual for whose phenotype should be predicted.
	 * @return The predicted phenotype of individual @p ind.
	 */
	PhenoType predict(Ind ind) const;

	/*!
	 * @brief Predicts the phenotype of a given individual and tracks the time of doing so.
	 * @param[in] ind The individual for whose phenotype should be predicted.
	 * @return The predicted phenotype of individual @p ind.
	 */
	PhenoType predict_track_time(Ind ind);

	/*!
	 * @brief The runtime of the last tracked prediction of the model.
	 * @return The runtime of the last call to predict_track_time().
	 */
	double prediction_time() const;


protected:

	/*!
	 * @brief A pointer to the instance on which the objective should be evaluated.
	 */
	Instance<PhenoType> * instance_;

	/*!
	 * @brief A flag that indicates whether or not the objective is initialized.
	 */
	bool initialized_;

private:

	std::map<std::string, std::string> options_;

	Seconds evaluation_time_;

	bool evaluation_time_available_;

	Seconds prediction_time_;

	bool prediction_time_available_;

	/*!
	 * @brief Returns the model sense.
	 * @return Returns epi::options::ModelSense::MINIMIZE if the objective should be minimized and epi::options::ModelSense::MAXIMIZE if it should be maximized.
	 * @note Must be overwritten by all derived classes.
	 */
	virtual options::ModelSense model_sense_() const;

	/*!
	 * @brief Returns @p true if the model is predictive.
	 * @return Boolean @p true if the model is predictive and @p false otherwise.
	 */
	virtual bool is_predictive_() const;

	/*!
	 * @brief Returns the objective for the SNP set @p snp_set.
	 * @details Returns the objective for the SNP set @p snp_set on the instance selected by the last call to set_instance(),
	 * given the options specified by the last call to set_options().
	 * @param[in] snp_set The SNP set on which the objective should be evaluated.
	 * @param[in] prepare_prediction If @p true and the epistasis model is predictive, this method must initialize all fields required by calls to predict_().
	 * @return The objective for the SNP set @p snp_set
	 * @note Must be overwritten by all derived classes.
	 */
	virtual double evaluate_(const std::vector<SNP> & snp_set, bool prepare_prediction);

	/*!
	 * @brief Predicts the phenotype of a given individual.
	 * @param[in] ind The individual for whose phenotype should be predicted.
	 * @return The predicted phenotype of individual @p ind.
	 * @note Must be overwritten by all predictive derived classes.
	 */
	virtual PhenoType predict_(Ind ind) const;

	/*!
	 * @brief Initializes the objective on the instance selected by the last call to set_instance().
	 * @note Must be overwritten by all derived classes that require initialization.
	 */
	virtual void initialize_();

	/*!
	 * @brief Parses one option.
	 * @param[in] option The name of the option.
	 * @param[in] arg The argument of the option.
	 * @return Boolean @p true if @p option is a valid option name for the method and @p false otherwise.
	 * @note Must be overwritten by all derived classes that have options.
	 */
	virtual bool parse_option_(const std::string & option, const std::string & arg);

	/*!
	 * @brief Sets all options to default values.
	 * @note Must be overwritten by all derived classes that have options.
	 */
	virtual void set_default_options_();

	/*!
	 * @brief Returns string of all valid options.
	 * @return String of the form <tt>[\--@<option@> @<arg@>] [...]</tt>.
	 * @note Must be overwritten by all derived classes that have options.
	 */
	virtual std::string valid_options_() const;

	/*!
	 * @brief Saves the model for the SNP set passed to the last call to evaluate(), evaluate_track_time(), or monte_carlo_p_value().
	 * @param[in] filename Path to the file where the model should be saved.
	 * @note Should be overwritten by all predictive derived classes.
	 */
	virtual void save_(const std::string & filename) const;

	/*!
	 * @brief Loads a model from disc. The model must match the dimensions of the epistasis instance.
	 * @param[in] filename Path to the file where the model is stored.
	 * @note Should be overwritten by all predictive derived classes.
	 */
	virtual void load_(const std::string & filename);

};

}

#include "epistasis_model.ipp"


#endif /* SRC_MODEL_EPISTASIS_MODEL_HPP_ */
