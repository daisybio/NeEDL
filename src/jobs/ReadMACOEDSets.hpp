//
// Created by juli on 21.10.22.
//

#ifndef GENEPISEEKER_READMACOEDSETS_HPP
#define GENEPISEEKER_READMACOEDSETS_HPP

#include "Job.hpp"

namespace epi {

    class ReadMACOEDSets : public Job {
    public:
        ReadMACOEDSets(const std::string& path, bool ignore_unknown_snps = false);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        std::string path;
        bool ignore_unknown_snps;
    };

} // epi

#ifdef HEADER_ONLY
#include "ReadMACOEDSets.cpp"
#endif


#endif //GENEPISEEKER_READMACOEDSETS_HPP
