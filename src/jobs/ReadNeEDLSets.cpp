//
// Created by juli on 21.10.22.
//

#include "ReadNeEDLSets.hpp"
#include "../util/FileFinder.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    ReadNeEDLSets::ReadNeEDLSets(const std::string& path) {
        FileFinder ff;
        ff.add_ends_with(".csv");
        this->path = find_file_by_ending(std::move(path), ff);

        // check if column indices are valid
        std::ifstream file(this->path);
        std::string first_line;
        std::getline(file, first_line);
        auto splits = string_split(first_line, '\t');
        bool found = false;
        for (size_t i = 0; i < splits.size(); i++) {
            if (splits[i] == "RS_IDS")  {
                found = true;
                this->column_index = i;
            }
        }

        if (!found) {
            throw epi::Error("Column RS_IDS not found in SNP set file.");
        }

        file.close();
    }

    void ReadNeEDLSets::run(std::shared_ptr<DataModel> data) {
        CSVParser parser;
        parser.parse(path, '\t', '\"', true, false);

        // search for set property columns
        std::vector<std::pair<size_t, std::string>> property_cols;
        for (size_t i = 0; i < parser.num_columns(); i++) {
            const std::string& colname = parser.cell(0, i);
            if (colname == "RS_IDS") continue;
            if (colname.substr(0, 6) == "RANK (") continue;
            if (options::is_epistasis_model_string(colname)) continue;
            if (colname == "ANNOTATIONS") continue;
            if (colname.substr(0, 16) == "NUM_INDIVIDUALS_") continue;
            if (colname.substr(0, 17) == "FREQ_INDIVIDUALS_") continue;
            if (colname.substr(0, 12) == "INDIVIDUALS_") continue;

            property_cols.emplace_back(i, colname);
        }

        data->snpSetStorage.clear();
        for (size_t i = 1; i < parser.num_rows(); i++) {
            const std::string& snp_string = parser.cell(i, column_index);
            auto rs_ids = string_split(snp_string, ';');
            SNPSet set;
            for (auto & id : rs_ids) {
                auto snp = data->snpStorage->by_name(id);
                set += snp;
            }

            for (auto & col : property_cols) {
                set.set_attribute(col.second, parser.cell(i, col.first));
            }

            data->snpSetStorage.push_back(set);
        }
    }

    rapidjson::Value ReadNeEDLSets::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("ReadNeEDLSets"), doc.GetAllocator());
        obj.AddMember("path", rapidjson::Value().SetString(path.c_str(), path.size(), doc.GetAllocator()), doc.GetAllocator());
        return obj;
    }
} // epi