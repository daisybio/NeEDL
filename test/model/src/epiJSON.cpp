//
// Created by juli on 07.06.21.
//
/**
 * calculates scores of SNP-sets in a list of result files
 */

#include <CLI11.hpp>

#define HEADER_ONLY

#include "../../../src/pipelines/HelperPipelines.hpp"
#include "../../../src/util/helper_functions.hpp"
#include "../../../src/jobs/DiseaseSNPReader.hpp"
#include "../../../src/jobs/PlinkConvertPlinkBinToJson.hpp"
#include "../../../src/jobs/PlinkFilterMissing.hpp"
#include "../../../src/jobs/PlinkFilterSexDiscrepancy.hpp"
#include "../../../src/jobs/PlinkFilterMAF.hpp"
#include "../../../src/jobs/PlinkRemoveTempFiles.hpp"
#include "../../../src/jobs/PlinkFilterHWE.hpp"
#include "../../../src/jobs/PlinkSetHHMissing.hpp"
#include "../../../src/jobs/PlinkCollectDatasetStats.hpp"
#include "../../../src/jobs/PlinkPrintDatasetStats.hpp"
#include "../../../src/jobs/PlinkFilterDuplicates.hpp"
#include "../../../src/jobs/PlinkFilterRSIds.hpp"
#include "../../../src/jobs/PlinkFilterSNPNetwork.hpp"
#include "../../../src/jobs/InternalBimLoader.hpp"
#include "../../../src/jobs/DbSNPAnnotator.hpp"
#include "../../../src/jobs/BioGridConnector.hpp"
#include "../../../src/jobs/NetworkStatsPrinter.hpp"
#include "../../../src/jobs/SameAnnotationConnector.hpp"
#include "../../../src/jobs/PlinkInputConverter.hpp"
#include "../../../src/jobs/PlinkMergeDatasets.hpp"
#include "../../../src/jobs/PlinkOutputConverter.hpp"
#include "../../../src/jobs/PlinkSetPhenotype.hpp"
#include "../../../src/jobs/PlinkConvertPlinkBinToMACOEDInput.hpp"
#include "../../../src/jobs/PlinkConvertJsonToPlink.hpp"
#include "../../../src/jobs/PlinkConvertPlinkBinToLINDENInput.hpp"
#include "../../../src/jobs/PlinkShufflePhenotype.hpp"

using namespace epi;

int main(int argc, char **argv) {
    // parse args
    CLI::App app("epiJSON generates input files needed for NeEDL from common plink-readable formats.");
    std::vector<std::string> input_files;
    app.add_option("--input-file", input_files,
                   "Path and filename of the dataset (without file extension if consisting of multiple files, e.g., .bim/.bed/.fam).")->required();

    std::vector<std::string> input_formats;
    app.add_option("--input-format", input_formats, "Type of input. Can be one of JSON, BIM_BED_FAM, PED_MAP, TPED_TFAM, VCF or BCF2. Can be specified either once, if all input files have the same format or need to be specified for every input file individually. If not specified at all, epiJSON tries to detect the input format automatically.");

    std::vector<std::string> override_phenotypes;
    app.add_option("--override-phenotype", override_phenotypes, "This option can optionally be set to override the phenotype contained in the input file. This is helpful if cases and controls are in separate files. If used, it needs to be specified once for every input file. If you want to use this feature only for some input files, set to 'NO' to use the phenotype in the input file. Please use the same notation as used in .fam files (1 = control, 2 = case, -9/0 = missing, or numeric data for categorical/quantitative data");

    bool shuffle_phenotype = false;
    app.add_flag("--shuffle-phenotype", shuffle_phenotype, "If this flag is set, the phenotype of all samples is shuffled.");

    std::string output_directory;
    app.add_option("--output-directory", output_directory,
                   "Output directory were all intermediate and final results should be written to. This should be empty prior to starting epiJSON.")->required();

    // output formats
    bool make_all_formats = false;
    app.add_flag("--make-all-formats", make_all_formats, "Set this flag to create all supported output formats.");

    bool make_json = false;
    app.add_flag("--make-json", make_json, "Set this flag to create a JSON output file that can be used with NeEDL.");

    bool make_bim_bed_fam = false;
    app.add_flag("--make-bim-bed-fam", make_bim_bed_fam, "Set this flag to create a .bim/.bed/.fam output.");

    bool make_ped_map = false;
    app.add_flag("--make-ped-map", make_ped_map, "Set this flag to create a .ped/.map output.");

    bool make_tped_tfam = false;
    app.add_flag("--make-tped-tfam", make_tped_tfam, "Set this flag to create a .tped/.tfam output.");

    bool make_vcf = false;
    app.add_flag("--make-vcf", make_vcf, "Set this flag to create a VCFv4.2 file which is block-gzipped.");

    bool make_macoed = false;
    app.add_flag("--make-macoed", make_macoed, "Set this flag to create an output file that can be read by MACOED. Can only be used for dichotomous phenotypes.");

    bool make_linden = false;
    app.add_flag("--make-linden", make_linden, "Set this flag to create an output file that can be read by LINDEN. Can only be used for dichotomous phenotypes.");



    // std::string ident_type = "RSID";
    // app.add_option("--ident-type", ident_type, "Type of SNP identifiers. Options are RSID, POSITION. DEFAULT: RSID");

    std::string phenotype;
    app.add_option("--phenotype", phenotype, "The type of the phenotypes.")->required();

    size_t num_categories = 2;
    app.add_option("--num-categories", num_categories,
                   "The number of categories for categorical phenotypes. DEFAULT: 2.");

    std::string disease_snps_file;
    app.add_option("--disease-snps", disease_snps_file,
                   "Optional: a file containing comma separated disease SNPs which should be included into the final JSON. Omit if no disease SNPs should be added");

    std::string ext_directory = "./ext/";
    app.add_option("--ext-directory", ext_directory,
                   "Path where external libraries are located. DEFAULT: ./ext/");

    std::string data_directory = "./data/";
    app.add_option("--data-directory", data_directory,
                   "Path where additional data like dbSNP and BIOGRID is located. DEFAULT: ./data/");

    int num_threads = 0;
    app.add_option("--num-threads", num_threads,
                   "The number of threads. 0 means all available cores. DEFAULT: use all available cores (" +
                   std::to_string(omp_get_num_procs()) + " on this machine)")->check(CLI::NonNegativeNumber);


    // filters
    bool filter_rsids = false;
    app.add_flag("--filter-rsids", filter_rsids, "Filter out any variants which IDs do not start with 'rs'");

    bool filter_duplicates = false;
    app.add_flag("--filter-duplicates", filter_duplicates, "Apply plink's duplicate filter");

    bool filter_missingness = false;
    app.add_flag("--filter-missing", filter_missingness, "Filters out all genotypes and individuals with > 0.2 missing values.");

    bool filter_sex_discrepancy = false;
    app.add_flag("--filter-sex-discrepancy", filter_sex_discrepancy, "Filters out samples with discrepancy between reported sex and sex according to genome.");

    bool filter_allele_frequency = false;
    app.add_flag("--filter-allele-frequency", filter_allele_frequency, "Apply MAF filter with threshold 0.05 for sets with #individuals < 50,000 and 0.01 else.");

    bool filter_hardy_weinberg = false;
    app.add_flag("--filter-hardy-weinberg", filter_hardy_weinberg, "Filter based on the Hard-Weinberg test.");

    bool filter_hh_missing = false;
    app.add_flag("--filter-hh-missing", filter_hh_missing, "Erases heterzygous haploid and nonnormal Y chromosome genotype calls.");

    bool filter_snp_network = false;
    app.add_flag("--filter-snp-network", filter_snp_network, "Filters out all SNPs that would not be part of the dbSNP/BioGRID based SNP-SNP-interaction network used by NeEDL.");

    bool all_filters = false;
    app.add_flag("--all-filters", all_filters, "Applies all available filters.");


    // Parse the options.
    CLI11_PARSE(app, argc, argv);

    if (num_threads < 0) {
        throw epi::Error("Invalid thread number provided.");
    }
    if (num_threads == 0) {
        num_threads = omp_get_num_procs();
    }
    omp_set_num_threads(num_threads);
    Logger::logLine(
            "OpenMP threads: " + std::to_string(omp_get_num_procs()) + " available, " + std::to_string(num_threads) +
            " selected");

    epi::SequentialJob seq;

    epi::expand_multi_param("--input-format", input_formats, std::string("AUTO"), input_files.size(), false);
    epi::expand_multi_param("--override-phenotype", override_phenotypes, std::string("NO"), input_files.size(), false);


    // convert all input files
    std::vector<std::string> converted_inputs;
    for (size_t i = 0; i < input_files.size(); i++) {
        std::string outfile = output_directory + "input_" + std::to_string(i + 1);
        if (input_formats[i] == "AUTO") {
            // if automatic filetype detection --> check if json file exists
            if (
                    (boost::filesystem::is_regular_file(input_files[i]) &&
                    input_files[i].substr(input_files[i].size() - 5) == ".json")
                    ) {

                input_formats[i] = "JSON";
            } else if (boost::filesystem::is_regular_file(input_files[i] + ".json")) {
                input_formats[i] = "JSON";
                input_files[i] += ".json";
            }
        }

        if (input_formats[i] == "JSON") {
            seq.add(std::make_shared<epi::PlinkConvertJsonToPlink>(input_files[i], outfile, phenotype, num_categories, ext_directory, num_threads));
        } else {
            seq.add(std::make_shared<epi::PlinkInputConverter>(input_formats[i], input_files[i], outfile, ext_directory,
                                                               num_threads));
        }
        if (override_phenotypes[i] != "NO") {
            std::string infile = outfile;
            outfile = output_directory + "input_pheno_" + std::to_string(i + 1);
            seq.add(std::make_shared<epi::PlinkSetPhenotype>(override_phenotypes[i], infile, outfile, ext_directory, num_threads));
            seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(infile, std::vector<std::string>{".bim", ".bed", ".fam", ".nosex" }));
        }
        converted_inputs.push_back(outfile);
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "input file " + std::to_string(i + 1)));
    }

    // combine datasets
    std::string current_input_file = output_directory + "combined";
    seq.add(std::make_shared<epi::PlinkMergeDatasets>(converted_inputs, current_input_file, ext_directory, num_threads));
    seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(current_input_file, phenotype, "input combined"));

    // remove individual datasets
    for (auto & conv_input : converted_inputs) {
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(conv_input, std::vector<std::string>{".bim", ".bed", ".fam", ".nosex" }));
    }


    // filter
    if (all_filters) {
        filter_rsids = filter_duplicates = filter_missingness = filter_sex_discrepancy = filter_allele_frequency = filter_hardy_weinberg = filter_hh_missing = filter_snp_network = true;
    }

    if (filter_rsids) {
        std::string outfile = output_directory + "filtered_rsids";
        seq.add(std::make_shared<epi::PlinkFilterRSIds>(current_input_file, outfile, ext_directory, num_threads));
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "filter_rsids"));
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));
        current_input_file = outfile;
    }

    if (filter_duplicates) {
        std::string outfile = output_directory + "filtered_dups";
        seq.add(std::make_shared<epi::PlinkFilterDuplicates>(current_input_file, outfile, ext_directory, num_threads));
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "filter_duplicates"));
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));
        current_input_file = outfile;
    }

    if (filter_missingness) {
        std::string outfile = output_directory + "filtered_missing";
        seq.add(std::make_shared<epi::PlinkFilterMissing>(current_input_file, outfile, ext_directory, "0.2", num_threads));
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "filter_missingness"));
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));
        current_input_file = outfile;
    }

    if (filter_sex_discrepancy) {
        std::string outfile = output_directory + "filtered_sex";
        seq.add(std::make_shared<epi::PlinkFilterSexDiscrepancy>(current_input_file, outfile, ext_directory, num_threads));
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "filter_sex_discrepancy"));
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));
        current_input_file = outfile;
    }

    if (filter_allele_frequency) {
        std::string outfile = output_directory + "filtered_maf";
        seq.add(std::make_shared<epi::PlinkFilterMAF>(current_input_file, outfile, ext_directory, num_threads));
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "filter_allele_frequency"));
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));
        current_input_file = outfile;
    }

    if (filter_hardy_weinberg) {
        std::string outfile = output_directory + "filtered_hwe";
        seq.add(std::make_shared<epi::PlinkFilterHWE>(current_input_file, outfile, ext_directory, num_threads));
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "filter_hardy_weinberg"));
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));
        current_input_file = outfile;
    }

    if (filter_hh_missing) {
        std::string outfile = output_directory + "filtered_hhmiss";
        seq.add(std::make_shared<epi::PlinkSetHHMissing>(current_input_file, outfile, ext_directory, num_threads));
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "filter_hh_missing"));
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));
        current_input_file = outfile;
    }

    // always do strict missingness filtering as NeEDL does not support any missing genotype data
    {
        std::string outfile = output_directory + "filtered_missingstrict";
        seq.add(std::make_shared<epi::PlinkFilterMissing>(current_input_file, outfile, ext_directory, "0", num_threads));
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "filter_strict_missing"));
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));
        current_input_file = outfile;
    }

    if (filter_snp_network) {
        // create snp network
        seq.add(std::make_shared<epi::InternalBimLoader>(current_input_file + ".bim"));
        seq.add(std::make_shared<epi::DbSNPAnnotator>(data_directory));
        seq.add(std::make_shared<epi::SameAnnotationConnector>());
        seq.add(std::make_shared<epi::BioGridConnector>(data_directory));
        seq.add(std::make_shared<epi::NetworkStatsPrinter>(false));

        std::string outfile = output_directory + "filtered_network";
        seq.add(std::make_shared<epi::PlinkFilterSNPNetwork>(current_input_file, outfile, ext_directory, num_threads));
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "filter_snp_network"));
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));
        current_input_file = outfile;
    }

    // optionally shuffle the phenotypes of all individuals
    if (shuffle_phenotype) {
        std::string outfile = output_directory + "shuffled_pheno";
        seq.add(std::make_shared<epi::PlinkShufflePhenotype>(current_input_file, outfile));
        seq.add(std::make_shared<epi::PlinkCollectDatasetStats>(outfile, phenotype, "shuffle_phenotype"));
        seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));
        current_input_file = outfile;
    }


    // process disease SNPs if provided
    if (!disease_snps_file.empty()) {
        seq.add(std::make_shared<epi::DiseaseSNPReader>(disease_snps_file, current_input_file));
    }


    // create requested output
    if (make_all_formats) {
        make_json = make_bim_bed_fam = make_ped_map = make_tped_tfam = make_vcf = make_macoed = make_linden = true;
    }
    if (!make_json && !make_bim_bed_fam && !make_ped_map && !make_tped_tfam && !make_vcf && !make_macoed && !make_linden) {
        throw epi::Error("No output format specified. Please provide at least one --make-... command.");
    }

    if (make_json) {
        seq.add(std::make_shared<epi::PlinkConvertPlinkBinToJson>(current_input_file, output_directory + "dataset.json", ext_directory, phenotype, num_categories, num_threads));
    }

    if (make_macoed && (!make_all_formats || phenotype == "DICHOTOMOUS")) {
        seq.add(std::make_shared<epi::PlinkConvertPlinkBinToMACOEDInput>(current_input_file, output_directory + "dataset.macoed", ext_directory, phenotype, num_threads));
    }

    if (make_linden && (!make_all_formats || phenotype == "DICHOTOMOUS")) {
        seq.add(std::make_shared<epi::PlinkConvertPlinkBinToLINDENInput>(current_input_file, output_directory + "dataset", ext_directory, phenotype, num_threads));
    }

    output_directory += "dataset";

    if (make_bim_bed_fam) {
        seq.add(std::make_shared<epi::PlinkOutputConverter>("BIM_BED_FAM", current_input_file, output_directory, ext_directory, num_threads));
    }

    if (make_ped_map) {
        seq.add(std::make_shared<epi::PlinkOutputConverter>("PED_MAP", current_input_file, output_directory, ext_directory, num_threads));
    }

    if (make_tped_tfam) {
        seq.add(std::make_shared<epi::PlinkOutputConverter>("TPED_TFAM", current_input_file, output_directory, ext_directory, num_threads));
    }

    if (make_vcf) {
        seq.add(std::make_shared<epi::PlinkOutputConverter>("VCF", current_input_file, output_directory, ext_directory, num_threads));
    }


    // print statistics about filtering and clean up
    seq.add(std::make_shared<epi::PlinkPrintDatasetStats>(phenotype));
    seq.add(std::make_shared<epi::PlinkRemoveTempFiles>(current_input_file, std::vector<std::string>{".bim", ".bed", ".fam" }));

    auto data = std::make_shared<epi::DataModel>(true);
    seq.run(data);

    return 0;
}

