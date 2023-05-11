//
// Created by juli on 17.04.23.
//

#ifndef GENEPISEEKER_PLINKFILTERHWE_HPP
#define GENEPISEEKER_PLINKFILTERHWE_HPP

#include "Job.hpp"

namespace epi {

    class PlinkFilterHWE : public Job {
    public:
        PlinkFilterHWE(std::string input_path, std::string output_path, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;

        double hwe_threshold_case = 1e-6;
        double hwe_threshold_control = 1e-10;
        double hwe_threshold_non_dichotomous = 1e-6;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkFilterHWE.cpp"
#endif

#endif //GENEPISEEKER_PLINKFILTERHWE_HPP
