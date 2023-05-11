//
// Created by juli on 01.05.23.
//

#include "PlinkConvertJsonToPlink.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkConvertJsonToPlink::PlinkConvertJsonToPlink(std::string input_path, std::string output_path, std::string phenotype, size_t num_categories, std::string ext_path, int num_threads) {
        this->input_path = input_path;
        this->output_path = output_path;
        this->phenotype = phenotype;
        this->num_categories = num_categories;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
    }

    void PlinkConvertJsonToPlink::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("Convert json format to .tped/.tfam files.");

        long num_snps = 0, num_inds = 0;
        std::vector<unsigned char> genotypes;
        std::vector<std::string> phenotypes;
        struct snp_data {
            std::string name;
            std::string chromosome;
            std::string position;
            std::string allele_from;
            std::string allele_to;
            unsigned long long genotype_sum;
        };
        std::vector<snp_data> snp_data;

        // load json file
        Logger::logLine("Parse json file");
        {
            rapidjson::Document json_doc;
            try {
                // std::ifstream json_file (filename);
                // rapidjson::IStreamWrapper json_file_wrapper(json_file);
                // json_doc.ParseStream(json_file_wrapper);

                FILE *fp = fopen(input_path.c_str(), "r");
                char readBuffer[65536];
                rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
                json_doc.ParseStream(is);

            }
            catch (...) {
                throw Error("The file " + input_path + " cannot be opened.");
            }

            {
                std::string reported_phenotype = toUpperCase(json_doc["model_type"].GetString());
                size_t reported_categories = 2;
                if (reported_phenotype != "QUANTITATIVE") {
                    reported_categories = json_doc["num_categories"].GetUint();
                }
                // fix model type in file but do not allow it to specified wrong when using epiJSON
                if (reported_phenotype == "CATEGORICAL" && reported_categories == 2) {
                    reported_phenotype = "DICHOTOMOUS";
                }
                if (reported_phenotype != phenotype) {
                    throw epi::Error("The reported phenotype of the json file does not match the phenotype provided to epiJSON (json file: " + reported_phenotype + ", epiJSON: " + phenotype + ")");
                }
                if (reported_categories != num_categories) {
                    throw epi::Error("The reported number of categories of the json file does not match the number of categories provided to epiJSON (json file: " + std::to_string(reported_categories) + ", epiJSON: " + std::to_string(num_categories) + ")");
                }
            }

            try {
                num_snps = json_doc["num_snps"].GetInt();
            }
            catch (...) {
                throw Error(
                        "The input file must contain a field \"num_snps\" whose value must be convertible to an int greater 0.");
            }

            if (num_snps <= 0) {
                throw Error(
                        "The input file must contain a field \"num_snps\" whose value must be convertible to an int greater 0.");
            }
            try {
                num_inds = json_doc["num_inds"].GetInt();
            }
            catch (...) {
                throw Error(
                        "The input file must contain a field \"num_inds\" whose value must be convertible to an int greater 0.");
            }
            if (num_inds <= 0) {
                throw Error(
                        "The input file must contain a field \"num_inds\" whose value must be convertible to an int greater 0.");
            }

            // Load genotypes.
            if (!json_doc.HasMember("genotype")) {
                throw Error("The input file must contain a field \"genotype\".");
            }
            SNP snp{0};
            Ind ind{0};
            for (auto &row: json_doc["genotype"].GetArray()) {
                ind = 0;
                for (auto &cell: row.GetArray()) {
                    try {
                        genotypes.emplace_back(cell.GetUint());
                    }
                    catch (...) {
                        throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) +
                                    " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
                    }
                    if (genotypes.back() > 2) {
                        throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) +
                                    " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
                    }
                    ind++;
                }
                if (ind != num_inds) {
                    throw Error(
                            "The actual number of individuals " + std::to_string(ind) + " for SNP " +
                            std::to_string(snp) +
                            " in the \"genotype\" field does not match the number of individuals " +
                            std::to_string(num_inds) + " specified in the \"num_inds\" field.");
                }
                snp++;
            }
            if (snp != num_snps) {
                throw Error(
                        "The actual number of SNPs in the \"genotype\" field  does not match the specified number of SNPs specified in the \"num_snps\" field.");
            }

            // Load phenotypes.
            if (!json_doc.HasMember("phenotype")) {
                throw Error("The input file must contain a field \"phenotype\".");
            }
            ind = 0;
            for (const auto &phenotype: json_doc["phenotype"].GetArray()) {
                std::string phenotype_parsed;
                if (phenotype.IsString()) phenotype_parsed = phenotype.GetString();
                else if (phenotype.IsInt()) phenotype_parsed = std::to_string(phenotype.GetInt());
                else if (phenotype.IsDouble()) phenotype_parsed = std::to_string(phenotype.GetDouble());
                else if (phenotype.IsFloat()) phenotype_parsed = std::to_string(phenotype.GetFloat());
                else throw Error("Phenotype has unknown datatype. Need to be Integer or String.");
                phenotypes.emplace_back(phenotype_parsed);
                ind++;
            }
            if (ind != num_inds) {
                throw Error(
                        "The actual number of individuals in the \"phenotype\" field  does not match the specified number of individuals specified in the \"num_inds\" field.");
            }

            // Load RS IDs.
            if (!json_doc.HasMember("snps")) {
                throw Error("The input file must contain a field \"snps\".");
            }
            snp = 0;
            for (auto &snp_info: json_doc["snps"].GetArray()) {
                int counter = 0;
                struct snp_data snp_d;
                for (auto &cell: snp_info.GetArray()) {

                    if (counter == 0) {
                        std::string rs_id_val = cell.GetString();
                        snp_d.name = cell.GetString();
                    } else if (counter == 1) {
                        snp_d.chromosome = cell.GetString();
                    } else if (counter == 2) {
                        snp_d.position = cell.GetString();
                    } else if (counter == 3) {
                        snp_d.allele_from = cell.GetString();
                    } else if (counter == 4) {
                        snp_d.allele_to = cell.GetString();
                    } else {
                        break;
                    }
                    counter++;
                }
                snp_data.push_back(snp_d);
                snp++;
            }

            if (snp != num_snps) {
                throw Error(
                        "The actual number of SNPs in the \"snps\" field  does not match the specified number of SNPs specified in the \"num_snps\" field.");
            }
        }


        // create tfam file
        Logger::logLine("Create tfam file");
        std::ofstream tfam_file (output_path + ".tfam");

        bool is_dichotomous = phenotype == "DICHOTOMOUS";
        for (size_t i = 0; i < num_inds; ++i) {
            // use index as within family id: according to plink documentation is not allowed to be 0
            tfam_file << "0 " << (i+1) << " 0 0 0 ";
            if (is_dichotomous) {
                tfam_file << std::stoi(phenotypes[i]) + 1;
            } else {
                tfam_file << phenotypes[i];
            }
            tfam_file << '\n';
        }

        tfam_file.close();

        // create tped file
        Logger::logLine("Create tped file");
        std::ofstream tped_file (output_path + ".tped");

        for (size_t i = 0; i < num_snps; ++i) {
            tped_file << snp_data[i].chromosome << ' ' << snp_data[i].name << " 0 " << snp_data[i].position;
            for (size_t j = 0; j < num_inds; ++j) {
                std::string variant_specifier;
                switch(genotypes[i * num_inds + j]) {
                    case 2:
                        variant_specifier = snp_data[i].allele_from + ' ' + snp_data[i].allele_from;
                        break;
                    case 1:
                        variant_specifier = snp_data[i].allele_to + ' ' + snp_data[i].allele_from;
                        break;
                    case 0:
                        variant_specifier = snp_data[i].allele_to + ' ' + snp_data[i].allele_to;
                        break;
                    default:
                        break;
                }
                tped_file << ' ' << variant_specifier;
            }

            tped_file << '\n';
        }

        tped_file.close();

        // create reference file
        Logger::logLine("Create ref file");
        std::ofstream ref_file (output_path + ".ref");

        for (size_t i = 0; i < num_snps; ++i) {
            ref_file << snp_data[i].name << ' ' << snp_data[i].allele_from << ' ' << snp_data[i].allele_to << '\n';
        }
        ref_file.close();

        run_plink_command(ext_path, {
                "--threads",
                std::to_string(num_threads),
                "--tfile",
                output_path,
                "--a1-allele",
                output_path + ".ref",
                "2",
                "1",
                "--noweb",
                "--make-bed",
                "--out",
                output_path
        });

        remove_plink_temp_files(output_path, { ".log", ".hh", ".nosex", ".tfam", ".ref", ".tped" });

        logger.stop();
    }


} // epi