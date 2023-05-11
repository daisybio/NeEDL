//
// Created by juli on 17.06.22.
//

#include "MAFFilter.hpp"
#include "../util/FileFinder.hpp"
#include "../util/helper_functions.hpp"
#include "../util/TimeLogger.hpp"


namespace epi {

    MAFFilter::MAFFilter(double filter_cutoff) {
        this->filter_cutoff = filter_cutoff;
        this->has_additional_maf_info = false;
    }

    MAFFilter::MAFFilter(double filter_cutoff, std::string additional_maf_file) {
        this->filter_cutoff = filter_cutoff;
        this->has_additional_maf_info = true;

        FileFinder ff;
        ff.add_ends_with(".csv");
        this->additional_maf_file = find_file_by_ending(std::move(additional_maf_file), ff);
    }

    void MAFFilter::run(std::shared_ptr <DataModel> data) {
        TimeLogger logger("MAF filter");

        if (!data->snpStorage->has_maf_information()) {
            // no instance maf --> check if we have additional maf data
            if (!has_additional_maf_info) {
                throw epi::Error("No MAF data in input file and no additional MAF file provided.");
            } else {
                CSVParser parser_add_maf_information;
                parser_add_maf_information.parse(additional_maf_file, '\t');

                std::vector<double> maf_data;
                for (std::size_t x = 0; x < parser_add_maf_information.num_rows(); x++) {
                    maf_data.emplace_back(std::stod(parser_add_maf_information.cell(x, 0)));
                }

                data->snpStorage->set_maf_information(maf_data);
            }
        }

        // store mma_scores and mark as removed if threshold is not reached
        for (const auto &snp : data->snpStorage->all()) {
            double maf = data->snpStorage->snp_get_maf(snp);
            if (maf >= filter_cutoff) data->snpStorage->set_removed(snp, true);
        }

        logger.stop();
    }

    rapidjson::Value MAFFilter::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("MAFFilter"), doc.GetAllocator());
        obj.AddMember("filter_cutoff", rapidjson::Value().SetDouble(filter_cutoff), doc.GetAllocator());
        obj.AddMember("has_additional_maf_info", rapidjson::Value().SetBool(has_additional_maf_info), doc.GetAllocator());
        if (has_additional_maf_info) {
            obj.AddMember("additional_maf_file", rapidjson::Value().SetString(additional_maf_file.c_str(), additional_maf_file.size(), doc.GetAllocator()), doc.GetAllocator());
        }
        return obj;
    }
} // epi