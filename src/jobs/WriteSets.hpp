//
// Created by juli on 22.06.22.
//

#ifndef GENEPISEEKER_WRITESETS_HPP
#define GENEPISEEKER_WRITESETS_HPP

#include "Job.hpp"

namespace epi {

    class WriteSets : public Job {
    public:
        WriteSets(std::string rank_model = "PENETRANCE",  std::vector<std::string> scores = {}, std::string name = "output", bool individual_snps = false);
        void run(std::shared_ptr<DataModel> data) override;
        void outfile_path(std::string path);
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        std::string name = "output";
        std::string ext_outfile_path = "";
        std::vector<options::EpistasisScore> epi_models;
        options::EpistasisScore rank_model;
        bool has_rank_model;
        bool individual_snps;
    };

} // epi

#ifdef HEADER_ONLY
#include "WriteSets.cpp"
#endif


#endif //GENEPISEEKER_WRITESETS_HPP
