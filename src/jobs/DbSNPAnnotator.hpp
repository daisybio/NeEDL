//
// Created by juli on 26.05.22.
//

#ifndef GENEPISEEKER_DBSNPANNOTATOR_HPP
#define GENEPISEEKER_DBSNPANNOTATOR_HPP

#include "SnpCsvAnnotator.hpp"

namespace epi {

    class DbSNPAnnotator : public SnpCsvAnnotator {
    public:
        explicit DbSNPAnnotator(const std::string& data_directory = "data/");
    };

} // epi

#ifdef HEADER_ONLY
#include "DbSNPAnnotator.cpp"
#endif

#endif //GENEPISEEKER_DBSNPANNOTATOR_HPP
