//
// Created by juli on 20.04.23.
//

#ifndef GENEPISEEKER_PLINKSETHHMISSING_HPP
#define GENEPISEEKER_PLINKSETHHMISSING_HPP

#include "Job.hpp"

namespace epi {

    class PlinkSetHHMissing : public Job {
    public:
        PlinkSetHHMissing(std::string input_path, std::string output_path, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkSetHHMissing.cpp"
#endif

#endif //GENEPISEEKER_PLINKSETHHMISSING_HPP
