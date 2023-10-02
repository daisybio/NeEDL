//
// Created by juli on 21.10.22.
//

#ifndef GENEPISEEKER_READNEEDLSETS_HPP
#define GENEPISEEKER_READNEEDLSETS_HPP

#include "Job.hpp"

namespace epi {

    class ReadNeEDLSets : public Job {
    public:
        ReadNeEDLSets(const std::string& path, bool ignore_unknown_snps = false);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        std::string path;
        size_t column_index;
        bool ignore_unknown_snps;
    };

} // epi

#ifdef HEADER_ONLY
#include "ReadNeEDLSets.cpp"
#endif


#endif //GENEPISEEKER_READNEEDLSETS_HPP
