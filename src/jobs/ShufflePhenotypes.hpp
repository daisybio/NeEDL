//
// Created by juli on 22.09.23.
//

#ifndef NEEDL_SHUFFLEPHENOTYPES_HPP
#define NEEDL_SHUFFLEPHENOTYPES_HPP

#include "Job.hpp"

namespace epi {

    class ShufflePhenotypes : public Job {
    public:
        void run(std::shared_ptr<DataModel> data) override;
    };

} // epi

#ifdef HEADER_ONLY
#include "ShufflePhenotypes.cpp"
#endif

#endif //NEEDL_SHUFFLEPHENOTYPES_HPP
