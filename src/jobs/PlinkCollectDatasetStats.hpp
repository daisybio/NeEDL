//
// Created by juli on 20.04.23.
//

#ifndef GENEPISEEKER_PLINKCOLLECTDATASETSTATS_HPP
#define GENEPISEEKER_PLINKCOLLECTDATASETSTATS_HPP

#include "Job.hpp"

namespace epi {

    class PlinkCollectDatasetStats : public Job {
    public:
        PlinkCollectDatasetStats(std::string path, std::string phenotype, std::string name);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string path;
        std::string name;
        bool is_dichotomous = false;
    };

} // epi

#ifdef HEADER_ONLY
#include "PlinkCollectDatasetStats.cpp"
#endif


#endif //GENEPISEEKER_PLINKCOLLECTDATASETSTATS_HPP
