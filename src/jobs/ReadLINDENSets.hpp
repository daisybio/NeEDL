//
// Created by juli on 21.10.22.
//

#ifndef GENEPISEEKER_READLINDENSETS_HPP
#define GENEPISEEKER_READLINDENSETS_HPP

#include "Job.hpp"

namespace epi {

    class ReadLINDENSets : public Job {
    public:
        ReadLINDENSets(const std::string& path, bool ignore_unknown_snps = false);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        std::string path;
        size_t column1_index, column2_index;
        bool ignore_unknown_snps;
    };

} // epi

#ifdef HEADER_ONLY
#include "ReadLINDENSets.cpp"
#endif


#endif //GENEPISEEKER_READLINDENSETS_HPP
