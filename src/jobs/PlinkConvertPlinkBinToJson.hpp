//
// Created by juli on 03.04.23.
//

#ifndef GENEPISEEKER_PLINKCONVERTPLINKBINTOJSON_HPP
#define GENEPISEEKER_PLINKCONVERTPLINKBINTOJSON_HPP

#include "Job.hpp"

namespace epi {

    class PlinkConvertPlinkBinToJson : public Job {
    public:
        PlinkConvertPlinkBinToJson(std::string input_path, std::string output_path, std::string ext_path, std::string phenotype, size_t num_categories, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    protected:
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;
        size_t num_categories;
        std::string phenotype;

        const size_t num_geno_cells_per_thread = 5000000;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkConvertPlinkBinToJson.cpp"
#endif

#endif //GENEPISEEKER_PLINKCONVERTPLINKBINTOJSON_HPP
