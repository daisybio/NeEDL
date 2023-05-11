//
// Created by juli on 19.05.22.
//

#include "instance.hpp"

namespace epi {

    template<>
    std::size_t
    Instance<CategoricalPhenoType>::
    num_categories() const {
        return num_categories_;
    }

    template<>
    bool
    Instance<QuantitativePhenoType>::
    quantitative_phenotypes() const {
        return true;
    }

    template<>
    bool
    Instance<CategoricalPhenoType>::
    quantitative_phenotypes() const {
        return false;
    }

    template<>
    CategoricalPhenoType
    Instance<CategoricalPhenoType>::
    parse_phenotype_(const std::string & phenotype_str, Ind ind) const {
        CategoricalPhenoType phenotype;
        try {
            phenotype = static_cast<CategoricalPhenoType>(std::stoul(phenotype_str));
        }
        catch (...) {
            throw Error("The input file contains invalid phenotype for individual " + std::to_string(ind) + ". Expected: convertible to int between 0 and " + std::to_string(num_categories_ - 1) + ".");
        }
        if (phenotype >= num_categories_) {
            throw Error("The input file contains invalid phenotype for individual " + std::to_string(ind) + ". Expected: convertible to int between 0 and " + std::to_string(num_categories_ - 1) + ".");
        }
        return phenotype;
    }

    template<>
    QuantitativePhenoType
    Instance<QuantitativePhenoType>::
    parse_phenotype_(const std::string & phenotype_str, Ind ind) const {
        QuantitativePhenoType phenotype;
        try {
            phenotype = static_cast<QuantitativePhenoType>(std::stod(phenotype_str));
        }
        catch (...) {
            throw Error("The input file contains invalid phenotype for individual " + std::to_string(ind) + ". Expected: convertible to double.");
        }
        return phenotype;
    }

}