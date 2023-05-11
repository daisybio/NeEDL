//
// Created by juli on 01.05.23.
//

#ifndef GENEPISEEKER_PLINKCONVERTPLINKBINTOLINDENINPUT_HPP
#define GENEPISEEKER_PLINKCONVERTPLINKBINTOLINDENINPUT_HPP

#include "Job.hpp"

namespace epi {

    class PlinkConvertPlinkBinToLINDENInput : public Job {
    public:
        PlinkConvertPlinkBinToLINDENInput(std::string input_path, std::string output_path, std::string ext_path, std::string phenotype, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;

        const size_t num_geno_cells_per_thread = 5000000;

    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkConvertPlinkBinToLINDENInput.cpp"
#endif

#endif //GENEPISEEKER_PLINKCONVERTPLINKBINTOLINDENINPUT_HPP
