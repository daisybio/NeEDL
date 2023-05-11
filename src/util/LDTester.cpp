//
// Created by juli on 28.09.22.
//

#include "LDTester.hpp"
#include "csv_parser.hpp"
#include "TimeLogger.hpp"

namespace epi {
    LDTester::LDTester(const std::string& ld_file, const std::string& mode, const std::shared_ptr<DataModel> &data) {
        // parse mode
        if (mode == "MEAN") {
            this->ld_mode = MEAN;
        } else if (mode == "MAX") {
            this->ld_mode = MAX;
        } else {
            throw epi::Error("Unknown LD mode " + mode);
        }

        // read matrix
        CSVParser ld_matrix_parser;
        ld_matrix_parser.parse(ld_file, ' ');

        if (ld_matrix_parser.num_columns() != data->snpStorage->num_snps()) {
            throw epi::Error("Number of SNPs does not match the number of columns in the LD file");
        }

        ld_matrix.resize(ld_matrix_parser.num_columns(), ld_matrix_parser.num_columns());
        size_t header_length = 0;
        for (; header_length < ld_matrix_parser.num_rows(); header_length++) if (ld_matrix_parser.cell(header_length,0) == "1") break;

        if (ld_matrix_parser.num_rows() - header_length != data->snpStorage->num_snps()) {
            throw epi::Error("Number of SNPs does not match the number of rows in the LD file (detected " + std::to_string(header_length) + " header rows)");
        }

        for (size_t i = header_length; i < ld_matrix_parser.num_rows(); i++) {
            for (size_t j = 0; j < ld_matrix_parser.num_rows(); j++) {
                ld_matrix(i, j) = std::stod(ld_matrix_parser.cell(i, j));
            }
        }
    }

    std::string LDTester::get_ld_mode_name() const {
        switch (ld_mode) {
            case MEAN:
                return "MEAN";
            case MAX:
                return "MAX";
        }
        return {};
    }

    LDTester::LDTester(const std::string& ld_file, const std::string& mode, const std::shared_ptr<DataModel> &data, double cutoff) : LDTester(ld_file, mode, data) {
        this->ld_cutoff = cutoff;
    }

    LDTester::LDTester(const std::string& ld_file, const std::string& mode, const std::shared_ptr<DataModel> &data, size_t min_set, size_t max_set, size_t sample_size_mc) : LDTester(ld_file, mode, data) {
        TimeLogger tl ("Determining LD cutoff through monte carlo sampling");

        std::uniform_int_distribution<size_t> set_size_distr(min_set, max_set);
        std::uniform_int_distribution<size_t> mat_row_distr(0, ld_matrix.row(0).size() - 1);

        std::vector<double> ld_sampling_values;

        auto thread_nr = omp_get_thread_num();

        for (size_t i = 0; i < sample_size_mc; i++) {
            size_t set_size = set_size_distr(data->random_device[thread_nr]);

            // draw random ints between 0 and matrix_indices.size()-1 without replacement for set_size + 1 (new snp)
            std::unordered_set<size_t> selected_snps_set;
            while (selected_snps_set.size() < set_size + 1) {
                size_t random_index = mat_row_distr(data->random_device[thread_nr]);
                if (selected_snps_set.find(random_index) == selected_snps_set.end())
                    selected_snps_set.insert(random_index);
            }

            std::vector<size_t> selected_snps;
            for (auto &snp : selected_snps_set) selected_snps.push_back(snp);

            // first randomly sampled index is set to be the row index (new snp)
            size_t row_index = selected_snps.back();
            selected_snps.pop_back();

            // all other randomly sampled indices represent the snp subset
            double ld_result = 0;
            if (ld_mode == MEAN) {
                for (auto & val : selected_snps) ld_result += ld_matrix(row_index, val);
                ld_result /= double(selected_snps.size());
            } else if (ld_mode == MAX) {
                for (auto & val : selected_snps) ld_result += std::max(ld_result, ld_matrix(row_index, val));
            }

            ld_sampling_values.push_back(ld_result);
        }

        // determine LD cutoff
        std::sort(ld_sampling_values.begin(), ld_sampling_values.end());
        this->ld_cutoff = ld_sampling_values[static_cast<size_t>(std::floor(double(sample_size_mc) * 0.95) - 1)];

        tl.stop();
        Logger::logLine("determined LD cutoff: " + Logger::to_string(this->ld_cutoff) + " (mode: " + get_ld_mode_name() + ')');
    }

    bool LDTester::test(const SNPSet &snp_set, const SNP_t &test_snp) {
       double ld_result = 0;
       if (ld_mode == MEAN) {
            for (auto & snp : snp_set) ld_result += ld_matrix(test_snp.value, snp.value);
            ld_result /= double(snp_set.size());
        } else if (ld_mode == MAX) {
            for (auto & snp : snp_set) ld_result += std::max(ld_result, ld_matrix(test_snp.value, snp.value));
        }

        return ld_result >= ld_cutoff;
    }
} // epi