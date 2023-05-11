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
 * @file  penetrance_model.ipp
 * @brief Definition of epi::PenetranceModel.
 */

#ifndef SRC_MODEL_PENETRANCE_MODEL_IPP_
#define SRC_MODEL_PENETRANCE_MODEL_IPP_

#include "penetrance_model.hpp"

namespace epi {



template<class PhenoType>
PenetranceModel<PhenoType>::
~PenetranceModel() {}

template<class PhenoType>
options::ModelSense
PenetranceModel<PhenoType>::
model_sense_() const {
	if (score_ == "LLH") {
		return options::ModelSense::MAXIMIZE;
	}
	return options::ModelSense::MINIMIZE;
}

template<class PhenoType>
bool
PenetranceModel<PhenoType>::
is_predictive_() const {
	return true;
}

template<class PhenoType>
double
PenetranceModel<PhenoType>::
evaluate_(const std::vector<SNP> & snp_set, bool prepare_prediction) {

	// Get the genotype IDs of all individuals at the given SNP set and construct the penetrance table.
	std::size_t size_penetrance_table{misc::size_penetrance_table(snp_set.size())};
	std::vector<std::vector<PhenoType>> penetrance_table(size_penetrance_table, std::vector<PhenoType>());
	std::vector<std::size_t> genotypes;
	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		genotypes.emplace_back(this->instance_->genotype_at_snp_set(snp_set, ind));
		penetrance_table.at(genotypes.back()).emplace_back(this->instance_->phenotype(ind));
	}

	// Fit the distributions for all cells of the penetrance table with sufficiently many samples.
	std::vector<std::vector<double>> maximum_likelihood_distributions;
	for (std::size_t genotype{0}; genotype < size_penetrance_table; genotype++) {
		if (penetrance_table.at(genotype).size() >= min_cell_size_) {
			maximum_likelihood_distributions.emplace_back(maximum_likelihood_distribution_(penetrance_table.at(genotype)));
		}
		else {
			maximum_likelihood_distributions.emplace_back(global_distribution_);
		}
	}

	// Store the distributions if required.
	if (prepare_prediction) {
		snp_set_ = snp_set;
		maximum_likelihood_distributions_ = maximum_likelihood_distributions;
	}

	// If the likelihood (LLH) is the selected score, compute and return it.
	if (score_ == "LLH") {
		double score{1.0};
		for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
			score *= likelihood_(this->instance_->phenotype(ind), maximum_likelihood_distributions.at(genotypes.at(ind)));
		}
		return score;
	}

	// Compute the negative log-likelihood.
	double score{0.0};
	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		score -= std::log(likelihood_(this->instance_->phenotype(ind), maximum_likelihood_distributions.at(genotypes.at(ind))));
	}

	// If the negative log-likelihood (NLL) is the selected score, return it.
	if (score_ == "NLL") {
		return score;
	}

	// Compute the degrees of freedom of the model, i.e., the number of estimated parameters.
	double degrees_of_freedom{static_cast<double>(size_penetrance_table * global_distribution_.size())};

	// If the Aikake information criterion (AIC) is the selected score, return it.
	if (score_ == "AIC") {
		return 2.0 * score + 2.0 * degrees_of_freedom;
	}

	// If the Bayesian information criterion (BIC) is the selected score, return it.
	return 2.0 * score + std::log(static_cast<double>(this->instance_->num_inds())) * degrees_of_freedom;
}






template<class PhenoType>
void
PenetranceModel<PhenoType>::
initialize_() {
	std::vector<PhenoType> phenotypes;
	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		phenotypes.emplace_back(this->instance_->phenotype(ind));
	}
	global_distribution_ = maximum_likelihood_distribution_(phenotypes);
}

template<class PhenoType>
bool
PenetranceModel<PhenoType>::
parse_option_(const std::string & option, const std::string & arg) {
	if (option == "score") {
		score_ = arg;
		if (score_ != "LLH" and score_ != "NLL" and score_ != "AIC" and score_ != "BIC") {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " LLH|NLL|AIC|BIC] [...]\"");
		}
		return true;
	}
	if (option == "min-cell-size") {
		try {
			min_cell_size_ = std::stoul(arg);
		}
		catch (...) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to int greater 0>] [...]\"");
		}
		if (min_cell_size_ <= 0) {
			throw Error(std::string("Invalid argument \"") + arg  + "\" for option \"--" + option + "\". Usage: options = \"[--" + option + " <convertible to int greater 0>] [...]\"");
		}
		return true;
	}
	return false;
}



template<class PhenoType>
std::string
PenetranceModel<PhenoType>::
valid_options_() const {
	return "[--score LLH|NLL|AIC|BIC] [--min-cell-size <convertible to int greater 0>]";
}

template<class PhenoType>
void
PenetranceModel<PhenoType>::
save_(const std::string & filename) const {
	std::ofstream file(filename.c_str());
	file << "[meta]\n";
	file << "model-type = PenetranceModel\n";
	file << meta_info_();
	file << "snp-set = ";
	std::size_t pos_next_item{1};
	for (SNP snp : snp_set_) {
		file << snp;
		if (pos_next_item++ < snp_set_.size()) {
			file << ",";
		}
	}
	file << "\n";
	file << "[parameters]\n";

	for (std::size_t genotype{0}; genotype < maximum_likelihood_distributions_.size(); genotype++) {
		file << genotype << " = ";
		pos_next_item = 1;
		for (double param : maximum_likelihood_distributions_.at(genotype)) {
			file << std::setprecision(10) << param;
			if (pos_next_item++ < maximum_likelihood_distributions_.at(genotype).size()) {
				file << ",";
			}
		}
		file << "\n";
	}
	file.close();
}

template<class PhenoType>
void
PenetranceModel<PhenoType>::
load_(const std::string & filename) {
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(filename, pt);
	if (pt.get<std::string>("meta.model-type") != "PenetranceModel") {
		throw Error(std::string("Unexpected model type ") + pt.get<std::string>("meta.model-type") + ".Expected: PenetranceModel.");
	}
	std::vector<std::string> tokens;
	misc::tokenize(pt.get<std::string>("meta.snp-set"), ',', '\'', true, tokens);
	snp_set_.clear();
	for (const auto & snp_as_string : tokens) {
		snp_set_.emplace_back(std::stoul(snp_as_string));
	}
	maximum_likelihood_distributions_.clear();
	for (std::size_t genotype{0}; genotype < misc::size_penetrance_table(snp_set_.size()); genotype++) {
		maximum_likelihood_distributions_.emplace_back(std::vector<double>());
		misc::tokenize(pt.get<std::string>(std::string("parameters.") + std::to_string(genotype)), ',', '\'', true, tokens);
		for (const auto & param_as_string : tokens) {
			maximum_likelihood_distributions_.back().emplace_back(std::stod(param_as_string));
		}
	}
	check_loaded_model_();
}


    template<class PhenoType>
    const std::string &PenetranceModel<PhenoType>::get_score() {
        return score_;
    }

}


#endif /* SRC_MODEL_PENETRANCE_MODEL_IPP_ */
