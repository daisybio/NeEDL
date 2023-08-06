//
// Created by juli on 26.05.22.
//

#include <utility>
#include "InstanceLoader.hpp"
#include "../util/helper_functions.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    InstanceLoader::InstanceLoader(std::string path_, std::string input_format_, std::string phenotype_, size_t num_categories_, std::string covariates_file_) {
        std::string ending = ".csv";
        input_format_ = toUpperCase(input_format_);
        inputFormatStr = input_format_;
        if (input_format_ == "JSON_EPIGEN") {
            inputFormat = options::InputFormat::JSON_EPIGEN;
            ending = ".json";
        } else if (input_format_ == "NEEDL_BIN") {
            inputFormat = options::InputFormat::NEEDL_BIN;
            ending = ".bin";
        } else if (input_format_ == "CSV_SNPS_AS_ROWS_FIRST") {
            inputFormat = options::InputFormat::CSV_SNPS_AS_ROWS_FIRST;
        } else if (input_format_ == "CSV_SNPS_AS_ROWS_LAST") {
            inputFormat = options::InputFormat::CSV_SNPS_AS_ROWS_LAST;
        } else if (input_format_ == "CSV_SNPS_AS_COLUMNS_FIRST") {
            inputFormat = options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST;
        } else if (input_format_ == "CSV_SNPS_AS_COLUMNS_LAST") {
            inputFormat = options::InputFormat::CSV_SNPS_AS_COLUMNS_LAST;
        } else {
            throw epi::Error("Invalid input format.");
        }

        // check path
        FileFinder ff;
        ff.add_ends_with(ending);
        this->path = find_file_by_ending(std::move(path_), ff);

        phenotype_ = toUpperCase(phenotype_);
        phenotypeStr = phenotype_;
        if (phenotype_ == "QUANTITATIVE") phenotype = epi::options::PhenoType::QUANTITATIVE;
        else if (phenotype_ == "CATEGORICAL" || phenotype_ == "DICHOTOMOUS") phenotype = epi::options::PhenoType::CATEGORICAL;
        else {
            throw epi::Error("Invalid phenotype " + phenotype_ + ".");
        }

        if (num_categories_ < 2) {
            throw epi::Error("At least 2 categories needed.");
        }
        num_categories = num_categories_;

        // check covariates file if provided
        if (!covariates_file_.empty()) {
            FileFinder cov_ff;
            cov_ff.add_ends_with(".csv");
            this->covariates_file = find_file_by_ending(covariates_file_, cov_ff);
        }
    }

    void InstanceLoader::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("loading GWAS data");
        Logger::logLine("input file: " + path);

        if (phenotype == options::PhenoType::QUANTITATIVE) {
            auto instance = std::make_shared<Instance<epi::QuantitativePhenoType>>();
            instance->load(inputFormat, path);
            if (!covariates_file.empty()) {
                instance->load_cov(options::InputFormat::CSV_COV, covariates_file);
            }
            auto storage = std::make_shared<SNPStorage_PhenoType<epi::QuantitativePhenoType>>(instance);
            auto snpStorage = std::static_pointer_cast<SNPStorage>(storage);
            data->snpStorage = snpStorage;
            SNPStorage::currentSnpStorage = snpStorage;
        } else if (phenotype == options::PhenoType::CATEGORICAL) {
            auto instance = std::make_shared<Instance<epi::CategoricalPhenoType>>(num_categories);
            instance->load(inputFormat, path);
            if (!covariates_file.empty()) {
                instance->load_cov(options::InputFormat::CSV_COV, covariates_file);
            }
            auto storage = std::make_shared<SNPStorage_PhenoType<epi::CategoricalPhenoType>>(instance);
            auto snpStorage = std::static_pointer_cast<SNPStorage>(storage);
            data->snpStorage = snpStorage;
            SNPStorage::currentSnpStorage = snpStorage;
        }

        logger.stop();
    }

    rapidjson::Value InstanceLoader::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("InstanceLoader"), doc.GetAllocator());
        obj.AddMember("path", rapidjson::Value().SetString(path.c_str(), path.size(), doc.GetAllocator()), doc.GetAllocator());
        obj.AddMember("input_format", rapidjson::Value().SetString(inputFormatStr.c_str(), inputFormatStr.size(), doc.GetAllocator()), doc.GetAllocator());
        obj.AddMember("phenotype", rapidjson::Value().SetString(phenotypeStr.c_str(), phenotypeStr.size(), doc.GetAllocator()), doc.GetAllocator());
        if (phenotype != options::PhenoType::QUANTITATIVE) {
            obj.AddMember("num_categories", rapidjson::Value().SetUint64(num_categories),doc.GetAllocator());
        }
        return obj;
    }

    const std::string InstanceLoader::getInputFilePath() const {
        return path;
    }
} // epi