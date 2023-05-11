//
// Created by juli on 26.05.22.
//

#ifndef GENEPISEEKER_BIOGRIDCONNECTOR_HPP
#define GENEPISEEKER_BIOGRIDCONNECTOR_HPP

#include "NetworkCsvConnector.hpp"
#include <string>

namespace epi {

    class BioGridConnector : public NetworkCsvConnector {
    public:
        explicit BioGridConnector(const std::string& data_directory = "data/");
    };

} // epi

#ifdef HEADER_ONLY
#include "BioGridConnector.cpp"
#endif

#endif //GENEPISEEKER_BIOGRIDCONNECTOR_HPP
