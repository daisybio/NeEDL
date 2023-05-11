//
// Created by juli on 24.04.23.
//

#ifndef GENEPISEEKER_INTERNALBIMLOADER_HPP
#define GENEPISEEKER_INTERNALBIMLOADER_HPP

#include "Job.hpp"

namespace epi {

    class InternalBimLoader : public Job {
    public:
        InternalBimLoader(std::string bim_path);
        void run(std::shared_ptr<DataModel> data) override;
    private:
        std::string bim_path;
    };

} // epi

#ifdef HEADER_ONLY
#include "InternalBimLoader.cpp"
#endif

#endif //GENEPISEEKER_INTERNALBIMLOADER_HPP
