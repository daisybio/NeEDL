//
// Created by juli on 26.05.22.
//

#ifndef GENEPISEEKER_INSTANCELOADER_HPP
#define GENEPISEEKER_INSTANCELOADER_HPP

#include "Job.hpp"

namespace epi {

    class InstanceLoader : public Job {
    public:
        InstanceLoader(std::string path, std::string input_format, std::string phenotype, size_t num_categories = 2);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
        const std::string getInputFilePath() const;
    protected:
        std::string path;
        std::string inputFormatStr;
        options::InputFormat inputFormat;
        std::string phenotypeStr;
        epi::options::PhenoType phenotype;
        size_t num_categories;
    };

} // epi

#ifdef HEADER_ONLY
#include "InstanceLoader.cpp"
#endif

#endif //GENEPISEEKER_INSTANCELOADER_HPP
