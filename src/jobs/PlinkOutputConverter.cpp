//
// Created by juli on 25.04.23.
//

#include "PlinkOutputConverter.hpp"
#include "../util/TimeLogger.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    PlinkOutputConverter::PlinkOutputConverter(std::string output_format, std::string input_path,
                                               std::string output_path, std::string ext_path, int num_threads) {
        if (output_format == "BIM_BED_FAM") {
            this->output_format = BIM_BED_FAM;
        } else if (output_format == "PED_MAP") {
            this->output_format = PED_MAP;
        } else if (output_format == "TPED_TFAM") {
            this->output_format = TPED_TFAM;
        } else if (output_format == "VCF") {
            this->output_format = VCF;
        } else {
            throw epi::Error("Unknown output format " + output_format);
        }

        this->input_path = input_path;
        this->output_path = output_path;
        this->ext_path = ext_path;
        this->num_threads = num_threads;
    }

    void PlinkOutputConverter::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger ("convert output format");

        std::vector<std::string> output_args;
        switch(output_format) {
            case PED_MAP:
                output_args = { "--recode" };
                break;
            case TPED_TFAM:
                output_args = { "--recode", "transpose" };
                break;
            case VCF:
                output_args = { "--recode", "vcf" };
                break;
            case BIM_BED_FAM:
                output_args = { "--make-bed" };
                break;
            default:
                output_args = {};
        }

        std::vector<std::string> args = {
                "--threads",
                std::to_string(num_threads),
                "--bfile",
                input_path,
                "--noweb"
        };

        args.insert(args.end(), output_args.begin(), output_args.end());

        args.push_back("--out");
        args.push_back(output_path);

        run_plink_command(ext_path, args);

        remove_plink_temp_files(output_path, { ".log", ".hh", ".nosex" });

        logger.stop();

    }
} // epi