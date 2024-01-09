//
// Created by juli on 09.01.24.
//

#ifndef NEEDL_PLINKFILTERSNPANNOTATION_HPP
#define NEEDL_PLINKFILTERSNPANNOTATION_HPP

#include "Job.hpp"

namespace epi {

    class PlinkFilterSNPAnnotation : public Job {
    public:
        PlinkFilterSNPAnnotation(std::vector<std::string> annotations, bool exclude_annotations, std::string input_path, std::string output_path, std::string ext_path, int num_threads);
        void run(std::shared_ptr<DataModel> data) override;


    private:
        std::vector<std::string> annotations;
        bool exclude_annotations;
        std::string input_path;
        std::string output_path;
        std::string ext_path;
        int num_threads;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkFilterSNPAnnotation.cpp"
#endif

#endif //NEEDL_PLINKFILTERSNPANNOTATION_HPP
