//
// Created by juli on 03.04.23.
//

#ifndef GENEPISEEKER_PLINKFILTERMISSING_HPP
#define GENEPISEEKER_PLINKFILTERMISSING_HPP

#include "Job.hpp"

namespace epi {

    class PlinkFilterMissing : public Job {
    public:
        PlinkFilterMissing(std::string input_path, std::string output_path, std::string ext_path, std::string filter_cutoff, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    protected:
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;
        std::string filter_cutoff;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkFilterMissing.cpp"
#endif

#endif //GENEPISEEKER_PLINKFILTERMISSING_HPP
