//
// Created by juli on 28.09.22.
//

#ifndef GENEPISEEKER_LDTESTER_HPP
#define GENEPISEEKER_LDTESTER_HPP

#include <string>
#include "types.hpp"
#include "../data_model/DataModel.hpp"

namespace epi {

    class LDTester {
    public:
        LDTester (const std::string& ld_file, const std::string& mode, const std::shared_ptr<DataModel> &data, double cutoff);
        LDTester (const std::string& ld_file, const std::string& mode, const std::shared_ptr<DataModel> &data, size_t min_set, size_t max_set, size_t sample_size_mc);

        bool test(const SNPSet &snp_set, const SNP_t &test_snp);
    private:
        LDTester (const std::string& ld_file, const std::string& mode, const std::shared_ptr<DataModel> &data);

        std::string get_ld_mode_name() const;

        enum LD_MODE {
            MEAN,
            MAX
        } ld_mode;

        Eigen::MatrixXd ld_matrix{};
        double ld_cutoff{};
    };

} // epi


#ifdef HEADER_ONLY
#include "LDTester.cpp"
#endif

#endif //GENEPISEEKER_LDTESTER_HPP
