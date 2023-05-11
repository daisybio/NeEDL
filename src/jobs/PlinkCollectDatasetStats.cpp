//
// Created by juli on 20.04.23.
//

#include "PlinkCollectDatasetStats.hpp"

namespace epi {
    PlinkCollectDatasetStats::PlinkCollectDatasetStats(std::string path, std::string phenotype, std::string name) {
        this->path = path;
        this->is_dichotomous = phenotype == "DICHOTOMOUS";
        this->name = name;
    }

    void PlinkCollectDatasetStats::run(std::shared_ptr<DataModel> data) {
        size_t num_inds = 0, num_variants = 0;
        size_t num_cases = 0, num_controls = 0;

        // own block so ind_parser data gets deleted to save memory
        {
            CSVParser ind_parser;
            ind_parser.parse(path + ".fam", ' ');

            if (ind_parser.num_columns() < 6) ind_parser.parse(path + ".fam", '\t');

            num_inds = ind_parser.num_rows();
            if (is_dichotomous) {
                for (size_t i = 0; i < ind_parser.num_rows(); ++i) {
                    int pheno = std::stoi(ind_parser.cell(i, 5));
                    if (pheno == 2) {
                        ++num_cases;
                    } else if (pheno == 1) {
                        ++num_controls;
                    }
                }
            }
        }

        // own block so variants_parser data gets deleted directly afterwards
        {
            CSVParser variants_parser;
            variants_parser.parse(path + ".bim");

            num_variants = variants_parser.num_rows();
        }

        data->dataset_stats.push_back({
                                              name,
                                              num_variants,
                                              num_inds,
                                              num_cases,
                                              num_controls
                                      });

    }
} // epi