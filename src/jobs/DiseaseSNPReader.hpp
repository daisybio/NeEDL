//
// Created by juli on 03.04.23.
//

#ifndef GENEPISEEKER_DISEASESNPREADER_HPP
#define GENEPISEEKER_DISEASESNPREADER_HPP

#include "Job.hpp"

namespace epi {

    class DiseaseSNPReader : public Job {
    public:
        DiseaseSNPReader(std::string disease_snp_file, std::string plink_filter_file = "");
        void run(std::shared_ptr<DataModel> data) override;
    protected:
        std::string path;
        std::string plink_filter_file;
    };

} // epi

#ifdef HEADER_ONLY
#include "DiseaseSNPReader.cpp"
#endif

#endif //GENEPISEEKER_DISEASESNPREADER_HPP
