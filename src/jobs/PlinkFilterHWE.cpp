//
// Created by juli on 17.04.23.
//

#include "PlinkFilterHWE.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkFilterHWE::PlinkFilterHWE(std::string input_path, std::string output_path, std::string ext_path,
                                   int num_threads) {

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;

    }

    void PlinkFilterHWE::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("HWE filtering");

        Logger::logLine("Run HWE test");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--hardy",
                "midp",
                "--nonfounders",
                "--noweb",
                "--out",
                output_path
        });

        Logger::logLine("Parse plink HWE file");
        // parse hwe file
        CSVParser parser;
        parser.parse(output_path + ".hwe", ' ');

        struct loci_data {
            std::string chr;
            std::string name;
            bool passes_case_test = false;
            bool passes_control_test = false;
        };
        std::unordered_map<std::string, struct loci_data> loci_map;

        for (size_t i = 1; i < parser.num_rows(); ++i) {
            std::string key = parser.cell(i, 0) + " : " + parser.cell(i, 1);
            auto loci_ptr = loci_map.find(key);
            if (loci_ptr == loci_map.end()) {
                loci_map.insert({ key, { parser.cell(i, 0), parser.cell(i, 1) }});
            }
            loci_ptr = loci_map.find(key);

            const std::string& test_type = parser.cell(i, 2);
            double p_value = 0.0;
            try {
                p_value = std::stod(parser.cell(i, 8));
            } catch (const std::out_of_range& ex) {
                p_value = 0.0;
            } catch (const std::invalid_argument& ex) {
                p_value = 0.0;
            }

            if (test_type == "AFF") {
                // case test
                if (p_value >= hwe_threshold_case) loci_ptr->second.passes_case_test = true;
            } else if (test_type == "UNAFF") {
                // control test
                if (p_value >= hwe_threshold_control) loci_ptr->second.passes_control_test = true;
            } else if (test_type == "ALL(QT)" || test_type == "ALL(NP)") {
                // non-dichotomous case
                if (p_value >= hwe_threshold_non_dichotomous) loci_ptr->second.passes_control_test = loci_ptr->second.passes_case_test = true;
            }
        }

        // filter selected loci
        Logger::logLine("Select loci based on hwe test");
        Logger::logLine("HWE thresholds are");
        Logger::logLine("   case samples: " + std::to_string(hwe_threshold_case));
        Logger::logLine("   control samples: " + std::to_string(hwe_threshold_control));
        Logger::logLine("   non-dichotomous samples: " + std::to_string(hwe_threshold_non_dichotomous));

        std::vector<std::pair<std::string, std::string>> selected_loci;
        for (auto & loci : loci_map) {
            if (loci.second.passes_case_test && loci.second.passes_control_test) {
                selected_loci.emplace_back(loci.second.chr, loci.second.name);
            }
        }

        // write succeeding loci to file
        Logger::logLine("Selected " + std::to_string(selected_loci.size()) + " loci after applying the HWE test.");
        std::ofstream selected_loci_file (output_path + ".selected_loci");
        for (auto & loci : selected_loci) {
            selected_loci_file << loci.second << '\n';
        }
        selected_loci_file.close();


        // filter out failing loci
        Logger::logLine("Using plink to exclude failing loci");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--attrib",
                output_path + ".selected_loci",
                "--noweb",
                "--make-bed",
                "--out",
                output_path
        });




        remove_plink_temp_files(output_path, { ".hh", ".log", ".hwe", ".selected_loci", ".nosex" });

        logger.stop();

    }
} // epi