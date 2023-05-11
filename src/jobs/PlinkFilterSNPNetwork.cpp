//
// Created by juli on 24.04.23.
//

#include "PlinkFilterSNPNetwork.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkFilterSNPNetwork::PlinkFilterSNPNetwork(std::string input_path, std::string output_path, std::string ext_path,
                                                 int num_threads) {

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
    }

    void PlinkFilterSNPNetwork::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("filter for network SNPs");
        if (data->snpNetwork == nullptr) throw epi::Error("SNPNetwork needed for filtering.");

        Logger::logLine("find entries that are not present in the network");
        CSVParser geno_parser;
        geno_parser.parse(input_path + ".bim", '\t');

        std::ofstream exclude_file (output_path + ".notnetwork");
        // find all entries that do not start with rs...
        for (size_t i = 0; i < geno_parser.num_rows(); ++i) {
            auto name = geno_parser.cell(i, 1);
            if (!data->snpNetwork->contains_node(data->snpStorage->by_name(name))) {
                exclude_file << name << '\n';
            }
        }
        exclude_file.close();

        Logger::logLine("Exclude those variants from the data");
        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--exclude",
                output_path + ".notnetwork",
                "--noweb",
                "--make-bed",
                "--out",
                output_path
        });

        remove_plink_temp_files(output_path, { ".log", ".hh", ".notnetwork", ".nosex" });

        logger.stop();

    }
} // epi