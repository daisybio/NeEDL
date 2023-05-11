//
// Created by juli on 01.05.23.
//

#ifndef GENEPISEEKER_PLINKCONVERTJSONTOPLINK_HPP
#define GENEPISEEKER_PLINKCONVERTJSONTOPLINK_HPP

#include "Job.hpp"

namespace epi {

    class PlinkConvertJsonToPlink : public Job {
    public:
        PlinkConvertJsonToPlink(std::string input_path, std::string output_path, std::string phenotype, size_t num_categories, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string input_path;
        std::string output_path;
        std::string phenotype;
        size_t num_categories;
        std::string ext_path;
        int num_threads;
    };

} // epi


#ifdef HEADER_ONLY
#include "PlinkConvertJsonToPlink.cpp"
#endif

#endif //GENEPISEEKER_PLINKCONVERTJSONTOPLINK_HPP
