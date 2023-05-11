//
// Created by juli on 21.10.22.
//

#ifndef GENEPISEEKER_CREATESETKMERS_HPP
#define GENEPISEEKER_CREATESETKMERS_HPP

#include "Job.hpp"

namespace epi {

    class CreateSetKMers : public Job {
    public:
        CreateSetKMers(size_t k_min, size_t k_max);
        CreateSetKMers(size_t k);

        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;

    private:
        long binom_coeff (long n, long k);
        int permute(int *vec, int k, int max);

        size_t k_min, k_max;
    };

} // epi

#ifdef HEADER_ONLY
#include "CreateSetKMers.cpp"
#endif

#endif //GENEPISEEKER_CREATESETKMERS_HPP
