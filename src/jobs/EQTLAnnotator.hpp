//
// Created by juli on 28.01.24.
//

#ifndef NEEDL_EQTLANNOTATOR_HPP
#define NEEDL_EQTLANNOTATOR_HPP

#include "SnpCsvAnnotator.hpp"

namespace epi {

    class EQTLAnnotator : public SnpCsvAnnotator {
    public:
        explicit EQTLAnnotator(const std::vector<std::string>& tissue_types, double pvalue_cutoff = 0.05, bool bh_correction = true, const std::string& data_directory = "data/");

    private:
        std::vector<bool> filter_entries(const epi::CSVParser &parser) override;

        std::unordered_set<std::string> tissue_types;
        double pvalue_cutoff;
        bool bh_correction;
    };

} // epi


#ifdef HEADER_ONLY
#include "EQTLAnnotator.cpp"
#endif

#endif //NEEDL_EQTLANNOTATOR_HPP
