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
#include "../../../src/jobs/DbSNPAnnotator.hpp"
#include "../../../src/jobs/ReadMACOEDSets.hpp"
#include "../../../src/jobs/ReadLINDENSets.hpp"

using namespace epi;

int main(int argc, char **argv) {
    // parse args
    CLI::App app("Calculate scores - A tool to evaluate SNP-sets in result files");
    std::string gwas_input_format;
    app.add_option("--gwas-input-format", gwas_input_format, "The format of the input GWAS file.")->required();
    std::string gwas_input_phenotype;
    app.add_option("--gwas-input-phenotype", gwas_input_phenotype, "The phenotype of the GWAS data.")->required();
    std::string gwas_input_file;
    app.add_option("--gwas-input-file", gwas_input_file, "The file containing GWAS data.")->required();
    unsigned int num_categories = 2;
    app.add_option("--num-categories", num_categories, "Number of categories if GWAS file is CATEGORICAL. Default: 2");

    std::string snp_sets_output_file;
    app.add_option("--snp-sets-output-file", snp_sets_output_file,
                   "The output file containing the calculated scores.");
    std::string snp_sets_input_type;
    app.add_option("--snp-sets-input-type", snp_sets_input_type,
                   "What format the SNP sets are in. Options are: NEEDL, MACOED, LINDEN")->required();
    std::string snp_sets_input_file;
    app.add_option("--snp-sets-input-file", snp_sets_input_file,
                   "Input file containing SNP sets to analyze.")->required();

    int num_threads = 1;
    app.add_option("--num-threads", num_threads,
                   "Number of threads to be used for score calculation. Default: 1. 0 means all available");

    bool create_random_sets = false;
    app.add_flag("--create-random-sets", create_random_sets, "If activated, random SNP sets are created with the same size distribution as the input SNP sets.");
    int num_random_sets = 1000;
    app.add_option("--num-random-sets", num_random_sets,
                   "Number of random sets to be generated if --create-random-sets is chosen.");

    std::vector<std::string> models;
    app.add_option("--model", models,
                   "Add a model to the list of models that are applied on the SNP sets. If none is specified, all available models will be used.");

    std::string rank_model;
    app.add_option("--rank-model", rank_model,
                   "Add a model to the list of models that is used to rank the SNP sets. If none is specified, no ranking is performed.");


    std::string data_directory = "./data/";
    app.add_option("--data-directory", data_directory,
                   "Path where additional data like dbSNP and BIOGRID is located. DEFAULT: ./data/");

    // snp annotation methods
    bool use_dbSNP = false;
    app.add_flag("--snp-annotate-dbSNP", use_dbSNP, "This option annotates the SNPs with the help of the included dbSNP database. This results in a gene annotation.");

    std::vector<std::string> snp_annotation_methods;
    app.add_option("--snp-annotate", snp_annotation_methods, "Declare annotation sources with this option. Option can be used multiple times to add multiple sources. The file must be a csv file with one column for the snp names (as in the input file) and one column for the annotations. Both columns can be multi-columns that contain multiple entries per row separated with a separating character. If a column is not a multi-column, give -1 as separator. If both column descriptors can be casted to a number this is used as a null-based index of the column Otherwise, the two values are searched for in the header line. The format of this option is <path>|<has-header ? 'yes' : 'no'>|<snp-col>|<annotation-col>|[<csv-separator>]|[<snp-separator>]|[<annotation-separator>].");

    std::string k_mers;
    app.add_option("--k-mers", k_mers, "Give a number (eg. 4) or a range (eg. 2-4) for k to generate scores of all k-mers.");

    bool ignore_unknown_snps = false;
    app.add_flag("--ignore-unknown-snps", ignore_unknown_snps, "If this flag is set, all SNP sets that contain a SNP which is not present in the selected dataset are ignored.");

    bool shuffle_phenotypes = false;
    app.add_flag("--shuffle-phenotypes", shuffle_phenotypes, "If this flag is set, the individuals' phenotypes are shuffled prior to score calculation.");

    std::string ld_directory;
    app.add_option("--ld-output-directory", ld_directory, "If this path is set, a LD matrix file (csv-format) is created in that directory for every SNP set. This should be a path to an existing directory");

    bool ld_only = false;
    app.add_flag("--ld-only", ld_only, "If this flag is set, only LD calculation but no score calculation is performed.");


    // Parse the options.
    CLI11_PARSE(app, argc, argv);

    if (snp_sets_output_file.empty() && !ld_only) {
        throw epi::Error("Either specify an output file via --snp-sets-output-file or select --ld-only.");
    }
    if (!snp_sets_output_file.empty() && ld_only) {
        Logger::logLine("Warning: --ld-only specified but an output file was specified via --snp-sets-output-file. Output file will be ignored.");
    }
    if (ld_only && ld_directory.empty()) {
        throw epi::Error("--ld-only set but no LD output directory given with --ld-output-directory.");
    }


    snp_sets_input_type = toUpperCase(snp_sets_input_type);
    std::shared_ptr<epi::Job> reader;
    if (snp_sets_input_type == "NEEDL") {
        reader = std::make_shared<epi::ReadNeEDLSets>(snp_sets_input_file, ignore_unknown_snps);
    } else if (snp_sets_input_type == "MACOED") {
        reader = std::make_shared<epi::ReadMACOEDSets>(snp_sets_input_file, ignore_unknown_snps);
    } else if (snp_sets_input_type == "LINDEN") {
        reader = std::make_shared<epi::ReadLINDENSets>(snp_sets_input_file, ignore_unknown_snps);
    } else {
        throw epi::Error("Unknown SNP-sets input type: " + snp_sets_input_type);
    }

    // parse k-mer string
    bool create_k_mers = false;
    size_t k_min = 0, k_max = 0;
    if (!k_mers.empty()) {
        create_k_mers = true;

        int p = k_mers.find('-');
        if (p == std::string::npos) {
            k_min = k_max = std::stoi(k_mers);
        } else {
            k_min = std::stoi(k_mers.substr(0, p));
            k_max = std::stoi(k_mers.substr(p + 1));
        }
    }

    if (models.empty()) {
        // set all available models
        models = epi::options::get_all_epistasis_scores();
    }

    std::vector<std::shared_ptr<epi::Job>> snpAnnotationPipeline;
    if (use_dbSNP) snpAnnotationPipeline.push_back(std::make_shared<epi::DbSNPAnnotator>(data_directory));
    for (auto &method: snp_annotation_methods) {
        snpAnnotationPipeline.push_back(std::make_shared<epi::SnpCsvAnnotator>(epi::SnpCsvAnnotator::parse_from_source_string(method)));
    }


    HelperPipelines::calculate_scores(
            gwas_input_file,
            gwas_input_format,
            gwas_input_phenotype,
            num_categories,
            reader,
            models,
            rank_model,
            snp_sets_output_file,
            snpAnnotationPipeline,
            shuffle_phenotypes,
            num_threads,
            create_random_sets,
            num_random_sets,
            create_k_mers,
            k_min,
            k_max,
            ld_directory,
            ld_only
            );


    return 0;
}

