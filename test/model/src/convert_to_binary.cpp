//
// Created by juli on 07.06.21.
//
/**
 * calculates scores of SNP-sets in a list of result files
 */

#include <CLI11.hpp>

#define HEADER_ONLY
#include "../../../src/jobs/ReadNeEDLSets.hpp"
#include "../../../src/pipelines/HelperPipelines.hpp"
#include "../../../src/util/helper_functions.hpp"

using namespace epi;

int main(int argc, char **argv) {
    // parse args
    CLI::App app("Convert to binary converts a dataset file (eg. in JSON_EPGIN format) into NeEDL binary format.");
    std::string gwas_input_format;
    app.add_option("--gwas-input-format", gwas_input_format, "The format of the input GWAS file.")->required();
    std::string gwas_input_phenotype;
    app.add_option("--gwas-input-phenotype", gwas_input_phenotype, "The phenotype of the GWAS data.")->required();
    std::string gwas_input_file;
    app.add_option("--gwas-input-file", gwas_input_file, "The file containing GWAS data.")->required();
    unsigned int num_categories = 2;
    app.add_option("--num-categories", num_categories, "Number of categories if GWAS file is CATEGORICAL. Default: 2");

    std::string output_file;
    app.add_option("--output-file", output_file,
                   "The binary output file")->required();


    // Parse the options.
    CLI11_PARSE(app, argc, argv);

    HelperPipelines::convert_dataset_to_binary(
            gwas_input_file,
            gwas_input_format,
            gwas_input_phenotype,
            num_categories,
            output_file
            );


    return 0;
}

