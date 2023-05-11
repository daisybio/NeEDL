//
// Created by juli on 25.04.23.
//

#include "PlinkInputConverter.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkInputConverter::PlinkInputConverter(std::string input_format, std::string input_path, std::string output_path,
                                             std::string ext_path, int num_threads, bool impute_sex) {

        if (input_format == "AUTO") {
          // try to detect the input file type automatically
          if (
                  boost::filesystem::is_regular_file(input_path + ".bim") &&
                  boost::filesystem::is_regular_file(input_path + ".fam") &&
                  boost::filesystem::is_regular_file(input_path + ".bed")
            ) {
              this->input_format = BIM_BED_FAM;
          } else if (
                  boost::filesystem::is_regular_file(input_path + ".ped") &&
                  boost::filesystem::is_regular_file(input_path + ".map")
                  ) {
              this->input_format = PED_MAP;
          } else if (
                  boost::filesystem::is_regular_file(input_path + ".tped") &&
                  boost::filesystem::is_regular_file(input_path + ".tfam")
                  ) {
              this->input_format = TPED_TFAM;
          } else if (
                  boost::filesystem::is_regular_file(input_path) &&
                  input_path.substr(input_path.size() - 4) == ".vcf"
                  ) {
              this->input_format = VCF;
          } else if (
                  boost::filesystem::is_regular_file(input_path + ".vcf")
                  ) {
              input_path += ".vcf";
              this->input_format = VCF;
          } else if (
                  boost::filesystem::is_regular_file(input_path) &&
                  input_path.substr(input_path.size() - 4) == ".bcf"
                  ) {
              this->input_format = BCF2;
          } else if (
                  boost::filesystem::is_regular_file(input_path + ".bcf")
                  ) {
              input_path += ".bcf";
              this->input_format = BCF2;
          } else {
              throw epi::Error("Cannot detect input format automatically for path '" + input_path + "'. Please specify format explicitly.");
          }
        } else if (input_format == "BIM_BED_FAM") {
            this->input_format = BIM_BED_FAM;
        } else if (input_format == "PED_MAP") {
            this->input_format = PED_MAP;
        } else if (input_format == "TPED_TFAM") {
            this->input_format = TPED_TFAM;
        } else if (input_format == "VCF") {
            this->input_format = VCF;
        } else if (input_format == "BCF2") {
            this->input_format = BCF2;
        } else {
            throw epi::Error("Unknown input format " + input_format);
        }

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
        this->impute_sex = impute_sex;
    }

    void PlinkInputConverter::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("convert input format to .bim/.bed/.fam");

        std::string input_arg;
        switch(input_format) {
            case PED_MAP:
                input_arg = "--file";
                break;
            case TPED_TFAM:
                input_arg = "--tfile";
                break;
            case VCF:
                input_arg = "--vcf";
                break;
            case BCF2:
                input_arg = "--bcf";
                break;
            case BIM_BED_FAM:
                input_arg = "--bfile";
                break;
            default:
                input_arg = "";
        }

        std::vector<std::string> arguments = {
                "--threads",
                std::to_string(num_threads),
                input_arg,
                input_path,
                "--noweb",
                "--make-bed",
                "--out",
                output_path
        };

        if (impute_sex) {
            arguments.push_back("--impute-sex");
        }

        run_plink_command(ext_path, arguments);

        remove_plink_temp_files(output_path, { ".log", ".hh", ".nosex" });

        logger.stop();
    }
} // epi