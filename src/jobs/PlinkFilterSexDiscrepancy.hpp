//
// Created by juli on 17.04.23.
//

#ifndef GENEPISEEKER_PLINKFILTERSEXDISCREPANCY_HPP
#define GENEPISEEKER_PLINKFILTERSEXDISCREPANCY_HPP

#include "Job.hpp"

namespace epi {

    class PlinkFilterSexDiscrepancy : public Job {
    public:
        PlinkFilterSexDiscrepancy(std::string input_path, std::string output_path, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkFilterSexDiscrepancy.cpp"
#endif

#endif //GENEPISEEKER_PLINKFILTERSEXDISCREPANCY_HPP
