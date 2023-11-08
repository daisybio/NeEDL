//
// Created by juli on 08.11.23.
//

#ifndef GENEPISEEKER_PLINKSHUFFLEPHENOTYPE_HPP
#define GENEPISEEKER_PLINKSHUFFLEPHENOTYPE_HPP

#include "Job.hpp"

namespace epi {

    class PlinkShufflePhenotype : public Job {
    public:
        PlinkShufflePhenotype(std::string input_path, std::string output_path);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string input_path;
        std::string output_path;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkShufflePhenotype.cpp"
#endif

#endif //GENEPISEEKER_PLINKSHUFFLEPHENOTYPE_HPP
