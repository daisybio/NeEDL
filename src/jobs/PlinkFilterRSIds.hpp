//
// Created by juli on 20.04.23.
//

#ifndef GENEPISEEKER_PLINKFILTERRSIDS_HPP
#define GENEPISEEKER_PLINKFILTERRSIDS_HPP

#include "Job.hpp"

namespace epi {

    class PlinkFilterRSIds : public Job {
    public:
        PlinkFilterRSIds(std::string input_path, std::string output_path, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkFilterRSIds.cpp"
#endif

#endif //GENEPISEEKER_PLINKFILTERRSIDS_HPP
