//
// Created by juli on 17.06.22.
//

#ifndef GENEPISEEKER_MAFFILTER_HPP
#define GENEPISEEKER_MAFFILTER_HPP

#include "Job.hpp"

namespace epi {

    class MAFFilter : public Job {
    public:
        MAFFilter(double filter_cutoff);
        MAFFilter(double filter_cutoff, std::string additional_maf_values);

        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;

    private:
        double filter_cutoff;
        bool has_additional_maf_info;
        std::string additional_maf_file;
    };

} // epi

#ifdef HEADER_ONLY
#include "MAFFilter.cpp"
#endif


#endif //GENEPISEEKER_MAFFILTER_HPP
