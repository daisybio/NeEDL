//
// Created by juli on 29.11.23.
//

#ifndef NEEDL_WRITESNPSETSLD_HPP
#define NEEDL_WRITESNPSETSLD_HPP

#include "Job.hpp"

namespace epi {

    class WriteSNPSetsLD : public Job {
    public:
        explicit WriteSNPSetsLD(const std::string &directory_name, const std::string &rank_model = "PENETRANCE");
        void set_root_path(const std::string &root_path_);

        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string directory_name;
        std::optional<std::string> root_path {};
        std::optional<options::EpistasisScore> rank_model;
    };

} // epi

#ifdef HEADER_ONLY
#include "WriteSNPSetsLD.cpp"
#endif

#endif //NEEDL_WRITESNPSETSLD_HPP
