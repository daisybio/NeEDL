//
// Created by juli on 25.04.23.
//

#ifndef GENEPISEEKER_PLINKMERGEDATASETS_HPP
#define GENEPISEEKER_PLINKMERGEDATASETS_HPP

#include "Job.hpp"

namespace epi {

    class PlinkMergeDatasets : public Job {
    public:
        PlinkMergeDatasets(std::vector<std::string> input_paths, std::string output_path, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    private:

        std::vector<std::string> input_paths;
        std::string output_path;
        std::string ext_path;
        int num_threads;

    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkMergeDatasets.cpp"
#endif

#endif //GENEPISEEKER_PLINKMERGEDATASETS_HPP
