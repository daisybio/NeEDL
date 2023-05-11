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
 * @file  variance_model.ipp
 * @brief Definition of epi::VarianceModel.
 */

#ifndef SRC_MODEL_VARIANCE_MODEL_IPP_
#define SRC_MODEL_VARIANCE_MODEL_IPP_

#include "variance_model.hpp"


namespace epi {

template<class PhenoType>
VarianceModel<PhenoType>::
VarianceModel(Instance<PhenoType> * instance):
EpistasisModel<PhenoType>(instance),
global_mean_{0.0},
num_inds_in_categories_(),
num_non_empty_categories_{0} {}

template<class PhenoType>
VarianceModel<PhenoType>::
~VarianceModel() {}



template<class PhenoType>
options::ModelSense
VarianceModel<PhenoType>::
model_sense_() const {
	return options::ModelSense::MINIMIZE;
}

template<class PhenoType>
bool
VarianceModel<PhenoType>::
is_predictive_() const {
	return false;
}

template<class PhenoType>
double
VarianceModel<PhenoType>::
evaluate_(const std::vector<SNP> & snp_set, bool prepare_prediction) {

	// Construct the penetrance table.
	std::vector<std::vector<PhenoType>> penetrance_table(misc::size_penetrance_table(snp_set.size()), std::vector<PhenoType>());

	//std::vector<std::vector<double>> penetrance_table_string(misc::size_penetrance_table(snp_set.size()),std::vector<double>());

	for (Ind ind{0}; ind < this->instance_->num_inds(); ind++) {
		penetrance_table.at(this->instance_->genotype_at_snp_set(snp_set, ind)).emplace_back(this->instance_->phenotype(ind));

        //penetrance_table_string.at(this->instance_->genotype_at_snp_set(snp_set, ind)).emplace_back(this->instance_->phenotype(ind));
	}

    /**
    * ACHTUNG!! NUR KURZ WEGEN TESTS!! MUSS AUSKOMMENTIERT WERDEN
    */
    /*
    std::fstream appendFileToWorkWith;
    appendFileToWorkWith.open("/home/exbio/SNPs/final_outputfile_test/7137_variance_rs429358.csv",  std::fstream::in | std::fstream::out | std::fstream::trunc);

    for(int i = 0; i < penetrance_table_string.size(); i++)
    {
        std::vector<double> line = penetrance_table_string.at(i);
        for(int j = 0; j < line.size(); j++)
        {
            appendFileToWorkWith << line.at(j);

            if(j < line.size()-1)
            {
                appendFileToWorkWith<<"\t";
            }
        }
        appendFileToWorkWith << "\n";
    }

    appendFileToWorkWith.close();
    */
    /**
     * HIER ENDE!!
     */


	// Compute test statistic and number of groups, i.e., non-empty cells in penetrance table.
	double test_statistic{0};
	std::size_t num_groups{0};
	compute_test_statistic_(penetrance_table, test_statistic, num_groups);

	// Return the p-value.
	return p_value_(test_statistic, num_groups);
}




}


#endif /* SRC_MODEL_VARIANCE_MODEL_IPP_ */
