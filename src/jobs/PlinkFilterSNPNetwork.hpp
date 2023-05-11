//
// Created by juli on 24.04.23.
//

#ifndef GENEPISEEKER_PLINKFILTERSNPNETWORK_HPP
#define GENEPISEEKER_PLINKFILTERSNPNETWORK_HPP

#include "Job.hpp"

namespace epi {

    class PlinkFilterSNPNetwork : public Job {
    public:
        PlinkFilterSNPNetwork(std::string input_path, std::string output_path, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;


    private:
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkFilterSNPNetwork.cpp"
#endif

#endif //GENEPISEEKER_PLINKFILTERSNPNETWORK_HPP
