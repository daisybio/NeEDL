//
// Created by juli on 25.04.23.
//

#ifndef GENEPISEEKER_PLINKINPUTCONVERTER_HPP
#define GENEPISEEKER_PLINKINPUTCONVERTER_HPP

#include "Job.hpp"

namespace epi {

    /*!
     * converts everything to .bim/.bed/.fam
     */
    class PlinkInputConverter : public Job {
    public:
        PlinkInputConverter(std::string input_format, std::string input_path, std::string output_path, std::string ext_path, int num_threads, bool impute_sex = false);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        enum {
            BIM_BED_FAM, // only copies files
            PED_MAP,
            TPED_TFAM,
            VCF,
            BCF2
        } input_format;

        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;
        bool impute_sex; // needed for conversion of json
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkInputConverter.cpp"
#endif

#endif //GENEPISEEKER_PLINKINPUTCONVERTER_HPP
