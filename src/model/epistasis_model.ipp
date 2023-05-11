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
 * @file  epistasis_model.ipp
 * @brief Definition of epi::EpistasisModel.
 */


#ifndef SRC_MODEL_EPISTASIS_MODEL_IPP_
#define SRC_MODEL_EPISTASIS_MODEL_IPP_

#include "epistasis_model.hpp"

namespace epi {

template<class PhenoType>
EpistasisModel<PhenoType>::
EpistasisModel(Instance<PhenoType> * instance):
instance_{instance},
initialized_{false},
options_(),
evaluation_time_(),
evaluation_time_available_{false},
prediction_time_(),
prediction_time_available_{false} {}

template<class PhenoType>
EpistasisModel<PhenoType>::
~EpistasisModel() {}

template<class PhenoType>
void
EpistasisModel<PhenoType>::
set_instance(Instance<PhenoType> * instance) {
	instance_ = instance;
	set_default_options_();
	initialized_ = false;
}

template<class PhenoType>
void
EpistasisModel<PhenoType>::
initialize() {
	initialize_();
	initialized_ = true;
}

template<class PhenoType>
void
EpistasisModel<PhenoType>::
set_options(const std::string & options) {
	misc::options_string_to_options_map(options, options_);
	set_default_options_();
	for (auto option_arg : options_) {
		if (not parse_option_(option_arg.first, option_arg.second)) {
			throw Error("Invalid option \"" + option_arg.first + "\". Usage: options = \"" + valid_options_() + "\".");
		}
	}
	initialized_ = false;
}

template<class PhenoType>
options::ModelSense
EpistasisModel<PhenoType>::
model_sense() const {
	return model_sense_();
}

template<class PhenoType>
bool
EpistasisModel<PhenoType>::
is_predictive() const {
	return is_predictive_();
}

template<class PhenoType>
bool
EpistasisModel<PhenoType>::
is_better(double objective_value_1, double objective_value_2) const {
	if (model_sense_() == epi::options::ModelSense::MAXIMIZE) {
		return objective_value_1 > objective_value_2;
	}
	return objective_value_1 < objective_value_2;
}

template<class PhenoType>
double
EpistasisModel<PhenoType>::
evaluate(const std::vector<SNP> & snp_set) {
	if (not initialized_) {
		initialize_();
		initialized_ = true;
	}
	return evaluate_(snp_set, true);
}

template<class PhenoType>
double
EpistasisModel<PhenoType>::
evaluate_track_time(const std::vector<SNP> & snp_set) {
	auto start = std::chrono::high_resolution_clock::now();
	if (not initialized_) {
		initialize_();
		initialized_ = true;
	}
	double objective_value{evaluate_(snp_set, true)};
	auto end = std::chrono::high_resolution_clock::now();
	evaluation_time_ = end - start;
	evaluation_time_available_ = true;
	return objective_value;
}

template<class PhenoType>
double
EpistasisModel<PhenoType>::
evaluation_time() const {
	if (not evaluation_time_available_) {
		throw Error("No evaluation time available. Call evaluate_track_time() before calling evaluation_time().");
	}
	return evaluation_time_.count();
}

template<class PhenoType>
double
EpistasisModel<PhenoType>::
prediction_time() const {
	if (not prediction_time_available_) {
		throw Error("No evaluation time available. Call evaluate_track_time() before calling evaluation_time().");
	}
	return prediction_time_.count();
}

template<class PhenoType>
double
EpistasisModel<PhenoType>::
monte_carlo_p_value(const std::vector<SNP> & snp_set, std::size_t num_permutations, std::size_t seed) {
	if (not initialized_) {
		initialize_();
		initialized_ = true;
	}
	if (seed == undefined_uint()) {
		std::random_device rng;
		seed = static_cast<std::size_t>(rng());
	}
	instance_->set_seed(seed);
	double objective_value{evaluate_(snp_set, true)};
	std::size_t num_better{0};
	for (std::size_t counter{0}; counter < num_permutations; counter++) {
		instance_->shuffle_phenotypes();
		if (not is_better(objective_value, evaluate_(snp_set, false))) {
			num_better++;
		}
	}
	instance_->restore_phenotypes();
	return (static_cast<double>(num_better + 1) / static_cast<double>(num_permutations + 1));
}

template<class PhenoType>
void
EpistasisModel<PhenoType>::
save(const std::string & filename) const {
	save_(filename);
}

template<class PhenoType>
void
EpistasisModel<PhenoType>::
load(const std::string & filename) {
	load_(filename);
}

template<class PhenoType>
PhenoType
EpistasisModel<PhenoType>::
predict(Ind ind) const {
	if (not is_predictive_()) {
		throw Error("The selected epistasis model does not allow prediction of phenotypes.");
	}
	return predict_(ind);
}

template<class PhenoType>
PhenoType
EpistasisModel<PhenoType>::
predict_track_time(Ind ind) {
	auto start = std::chrono::high_resolution_clock::now();
	if (not is_predictive_()) {
		throw Error("The selected epistasis model does not allow prediction of phenotypes.");
	}
	PhenoType prediction{predict_(ind)};
	auto end = std::chrono::high_resolution_clock::now();
	prediction_time_ = end - start;
	prediction_time_available_ = true;
	return prediction;
}

// Default implementations of virtual member functions.

template<class PhenoType>
options::ModelSense
EpistasisModel<PhenoType>::
model_sense_() const {
	return options::ModelSense::MINIMIZE;
}

template<class PhenoType>
bool
EpistasisModel<PhenoType>::
is_predictive_() const {
	return false;
}

template<class PhenoType>
double
EpistasisModel<PhenoType>::
evaluate_(const std::vector<SNP> & snp_set, bool prepare_prediction) {
	return 0.0;
}

template<class PhenoType>
PhenoType
EpistasisModel<PhenoType>::
predict_(Ind ind) const {
	PhenoType pheno_type{0};
	return pheno_type;
}

template<class PhenoType>
void
EpistasisModel<PhenoType>::
initialize_() {}

template<class PhenoType>
bool
EpistasisModel<PhenoType>::
parse_option_(const std::string & option, const std::string & arg) {
	return false;
}

template<class PhenoType>
void
EpistasisModel<PhenoType>::
set_default_options_() {}

template<class PhenoType>
std::string
EpistasisModel<PhenoType>::
valid_options_() const {
	return "";
}

template<class PhenoType>
void
EpistasisModel<PhenoType>::
save_(const std::string & filename) const {}

template<class PhenoType>
void
EpistasisModel<PhenoType>::
load_(const std::string & filename) {}


}


#endif /* SRC_MODEL_EPISTASIS_MODEL_IPP_ */
