//
// Created by juli on 20.04.23.
//

#include "PlinkPrintDatasetStats.hpp"
#include "../util/TerminalTable.hpp"

namespace epi {
    void PlinkPrintDatasetStats::run(std::shared_ptr<DataModel> data) {
        TerminalTable table;

        if (is_dichtomous) {
            table.add_row({"stage", "# variants", "# samples", "# cases", "# controls"});
        } else {
            table.add_row({"stage", "# variants", "# samples"});
        }

        for (auto & stat : data->dataset_stats) {
            if (is_dichtomous) {
                table.add_row({
                    stat.name,
                    std::to_string(stat.num_variants),
                    std::to_string(stat.num_samples),
                    std::to_string(stat.num_samples_cases),
                    std::to_string(stat.num_samples_controls)
                });
            } else {
                table.add_row({
                    stat.name,
                    std::to_string(stat.num_variants),
                    std::to_string(stat.num_samples)
                });
            }
        }

        table.print();
    }

    PlinkPrintDatasetStats::PlinkPrintDatasetStats(std::string phenotype) {
        this->is_dichtomous = phenotype == "DICHOTOMOUS";
    }
} // epi