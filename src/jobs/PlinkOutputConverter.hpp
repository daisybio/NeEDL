//
// Created by juli on 25.04.23.
//

#ifndef GENEPISEEKER_PLINKOUTPUTCONVERTER_HPP
#define GENEPISEEKER_PLINKOUTPUTCONVERTER_HPP

#include "Job.hpp"

namespace epi {

    class PlinkOutputConverter : public Job {
    public:
        PlinkOutputConverter(std::string output_format, std::string input_path, std::string output_path, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        enum {
            BIM_BED_FAM, // only copies files
            PED_MAP,
            TPED_TFAM,
            VCF
        } output_format;

        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;

    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkOutputConverter.cpp"
#endif

#endif //GENEPISEEKER_PLINKOUTPUTCONVERTER_HPP
