//
// Created by juli on 17.04.23.
//

#ifndef GENEPISEEKER_PLINKFILTERMAF_HPP
#define GENEPISEEKER_PLINKFILTERMAF_HPP

#include "Job.hpp"

namespace epi {

    class PlinkFilterMAF : public Job {
    public:
        PlinkFilterMAF(std::string input_path, std::string output_path, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;

        std::string maf_threshold_small = "0.05";
        std::string maf_threshold_large = "0.01";
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkFilterMAF.cpp"
#endif

#endif //GENEPISEEKER_PLINKFILTERMAF_HPP
