//
// Created by juli on 27.05.22.
//

#ifndef GENEPISEEKER_NETWORKSTATSPRINTER_HPP
#define GENEPISEEKER_NETWORKSTATSPRINTER_HPP

#include "Job.hpp"

namespace epi {

    class NetworkStatsPrinter : public Job {
    public:
        NetworkStatsPrinter(bool printTimeExpensiveSats = true);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        bool printTimeExpensiveStats = true;
    };

} // epi

#ifdef HEADER_ONLY
#include "NetworkStatsPrinter.cpp"
#endif


#endif //GENEPISEEKER_NETWORKSTATSPRINTER_HPP
