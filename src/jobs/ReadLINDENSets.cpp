//
// Created by juli on 21.10.22.
//

#include "ReadLINDENSets.hpp"
#include "../util/FileFinder.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    ReadLINDENSets::ReadLINDENSets(const std::string& path, bool ignore_unknown_snps) {
        FileFinder ff;
        ff.add_ends_with(".reciprocalPairs");
        this->path = find_file_by_ending(std::move(path), ff);

        // check if column indices are valid
        std::ifstream file(this->path);
        std::string first_line;
        std::getline(file, first_line);
        auto splits = string_split(first_line, '\t');
        bool found1 = false, found2 = false;
        for (size_t i = 0; i < splits.size(); i++) {
            if (splits[i] == "snp1")  {
                found1 = true;
                this->column1_index = i;
            } else if (splits[i] == "snp2") {
                found2 = true;
                this->column2_index = i;
            }
        }

        if (!found1 || !found2) {
            throw epi::Error("Column snp1 or snp2 not found in SNP set file.");
        }

        file.close();


        this->ignore_unknown_snps = ignore_unknown_snps;
    }

    void ReadLINDENSets::run(std::shared_ptr<DataModel> data) {
        CSVParser parser;
        parser.parse(path, '\t');

        data->snpSetStorage.clear();
        for (size_t i = 1; i < parser.num_rows(); i++) {
            const std::string& snp1_string = parser.cell(i, column1_index);
            const std::string& snp2_string = parser.cell(i, column2_index);
            SNPSet set;

            try {
                set += data->snpStorage->by_name(snp1_string);
            } catch (epi::SNPNotFoundError &err) {
                if (!ignore_unknown_snps) throw err;
            }

            try {
                set += data->snpStorage->by_name(snp2_string);
            } catch (epi::SNPNotFoundError &err) {
                if (!ignore_unknown_snps) throw err;
            }

            if (set.size() == 0) continue;

            set.set_attribute("LINDEN_rank", std::to_string(i));

            data->snpSetStorage.push_back(set);
        }
    }

    rapidjson::Value ReadLINDENSets::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("ReadLINDENSets"), doc.GetAllocator());
        obj.AddMember("path", rapidjson::Value().SetString(path.c_str(), path.size(), doc.GetAllocator()), doc.GetAllocator());
        return obj;
    }
} // epi