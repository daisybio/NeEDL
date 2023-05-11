//
// Created by juli on 20.04.23.
//

#ifndef GENEPISEEKER_PLINKPRINTDATASETSTATS_HPP
#define GENEPISEEKER_PLINKPRINTDATASETSTATS_HPP

#include "Job.hpp"

namespace epi {

    class PlinkPrintDatasetStats : public Job {
    public:
        PlinkPrintDatasetStats(std::string phenotype);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        bool is_dichtomous = false;
    };

} // epi


#ifdef HEADER_ONLY
#include "PlinkPrintDatasetStats.cpp"
#endif

#endif //GENEPISEEKER_PLINKPRINTDATASETSTATS_HPP
