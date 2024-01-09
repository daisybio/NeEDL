//
// Created by juli on 09.01.24.
//

#include "PlinkFilterSNPAnnotation.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkFilterSNPAnnotation::PlinkFilterSNPAnnotation(std::vector<std::string> annotations, bool exclude_annotations,
                                                       std::string input_path, std::string output_path,
                                                       std::string ext_path, int num_threads) {

        this->annotations = annotations;
        this->exclude_annotations = exclude_annotations;

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
    }

    void PlinkFilterSNPAnnotation::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("filter for annotations");

        Logger::logLine("find entries that have at least one of the selected annotations");

        std::unordered_set<std::string> with_annotation_names;
        { // delete with_annotation after this block
            std::unordered_set<SNP_t, SNP_t::SNPHash> with_annotation;
            for (const auto &anno: annotations) {
                auto snps = data->snpStorage->by_annotation(anno);
                with_annotation.insert(snps.begin(), snps.end());
            }
            with_annotation_names.reserve(with_annotation.size());
            for (const auto &snp: with_annotation) {
                with_annotation_names.insert(data->snpStorage->snp_get_name(snp));
            }
        }

        CSVParser geno_parser;
        geno_parser.parse(input_path + ".bim", '\t');

        std::ofstream exclude_file (output_path + ".withanno");
        // find all entries that do not start with rs...
        for (size_t i = 0; i < geno_parser.num_rows(); ++i) {
            const auto& name = geno_parser.cell(i, 1);
            if (with_annotation_names.contains(name)) {
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
                exclude_annotations ? "--exclude" : "--extract",
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