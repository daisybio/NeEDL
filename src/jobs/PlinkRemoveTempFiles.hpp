//
// Created by juli on 17.04.23.
//

#ifndef GENEPISEEKER_PLINKREMOVETEMPFILES_HPP
#define GENEPISEEKER_PLINKREMOVETEMPFILES_HPP

#include "Job.hpp"

namespace epi {

    class PlinkRemoveTempFiles : public Job {
    public:
        PlinkRemoveTempFiles(std::string path, std::vector<std::string> endings);
        void run(std::shared_ptr<DataModel> data) override;

    private:
        std::string path;
        std::vector<std::string> endings;

    };

} // epi


#ifdef HEADER_ONLY
#include "PlinkRemoveTempFiles.cpp"
#endif

#endif //GENEPISEEKER_PLINKREMOVETEMPFILES_HPP
