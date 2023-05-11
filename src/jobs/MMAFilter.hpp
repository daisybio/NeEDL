//
// Created by juli on 17.06.22.
//

#ifndef GENEPISEEKER_MMAFILTER_HPP
#define GENEPISEEKER_MMAFILTER_HPP

#include "Job.hpp"

namespace epi {

    class MMAFilter : public Job {
    public:
        MMAFilter(double filter_cutoff, bool use_BH_correction);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;

    private:
        double filter_cutoff;
        bool use_BH_correction;
    };

} // epi

#ifdef HEADER_ONLY
#include "MMAFilter.cpp"
#endif


#endif //GENEPISEEKER_MMAFILTER_HPP
