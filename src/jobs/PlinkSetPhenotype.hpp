//
// Created by juli on 26.04.23.
//

#ifndef GENEPISEEKER_PLINKSETPHENOTYPE_HPP
#define GENEPISEEKER_PLINKSETPHENOTYPE_HPP

#include "Job.hpp"

namespace epi {

    class PlinkSetPhenotype : public Job {
    public:
        PlinkSetPhenotype(std::string phenotype, std::string input_path, std::string output_path, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string phenotype;
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkSetPhenotype.cpp"
#endif

#endif //GENEPISEEKER_PLINKSETPHENOTYPE_HPP
