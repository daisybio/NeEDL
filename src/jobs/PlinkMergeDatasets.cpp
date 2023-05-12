//
// Created by juli on 25.04.23.
//

#include "PlinkMergeDatasets.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkMergeDatasets::PlinkMergeDatasets(std::vector<std::string> input_paths, std::string output_path,
                                           std::string ext_path, int num_threads) {

        this->input_paths = input_paths;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
    }

    void PlinkMergeDatasets::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("merge input files");

        std::ofstream mergelist (output_path + ".mergelist");
        for (auto & input_path : input_paths) {
            mergelist << input_path << '\n';
        }
        mergelist.close();

        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--merge-list",
                output_path + ".mergelist",
                "--noweb",
                "--allow-no-sex",
                "--make-bed",
                "--out",
                output_path
        });


        // merge failed because multiple alleles exist for some positions --> exclude them
        if (boost::filesystem::exists(output_path + "-merge.missnp")) {
            std::ofstream mergelist2(output_path + ".mergelist2");
            std::vector<std::string> intermediate_paths;

            for (size_t i = 0; i < input_paths.size(); ++i) {
                // filter out the failing SNPs from every input dataset
                std::string intermediate_path = output_path + "_exclude_failed_" + std::to_string(i);
                run_plink_command(ext_path, {
                        "--threads",
                        std::to_string(num_threads),
                        "--bfile",
                        input_paths[i],
                        "--noweb",
                        "--exclude",
                        output_path + "-merge.missnp",
                        "--make-bed",
                        "--out",
                        intermediate_path
                });

                mergelist2 << intermediate_path << '\n';

                intermediate_paths.push_back(intermediate_path);
            }
            mergelist2.close();

            // repeat merging
            run_plink_command(ext_path, {
                    "--threads",
                    std::to_string(num_threads),
                    "--merge-list",
                    output_path + ".mergelist2",
                    "--noweb",
                    "--allow-no-sex",
                    "--make-bed",
                    "--out",
                    output_path
            });

            // delete intermediate files
            for (auto & intermediate_path : intermediate_paths) {
                remove_plink_temp_files(intermediate_path, { ".log", ".hh", ".nosex", ".bim", ".bed", ".fam" });

            }
        }


        remove_plink_temp_files(output_path, { ".log", ".hh", ".mergelist", ".mergelist2", ".nosex", "-merge.missnp" });

        logger.stop();
    }
} // epi