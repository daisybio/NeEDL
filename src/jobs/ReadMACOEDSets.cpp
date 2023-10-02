//
// Created by juli on 21.10.22.
//

#include "ReadMACOEDSets.hpp"
#include "../util/FileFinder.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    ReadMACOEDSets::ReadMACOEDSets(const std::string& path, bool ignore_unknown_snps) {
        FileFinder ff;
        ff.add_ends_with(".txt");
        this->path = find_file_by_ending(std::move(path), ff);

       // do not perform a column check here because these files are really fucked up


        this->ignore_unknown_snps = ignore_unknown_snps;
    }

    void ReadMACOEDSets::run(std::shared_ptr<DataModel> data) {
        data->snpSetStorage.clear();

        std::ifstream file(this->path);
        std::string line;
        while(file.good() && !file.eof()) {
            std::getline(file, line);
            std::vector<std::string> splits = string_split(line, '\t');

            if (splits.size() < 3) continue;
            if (splits[0] != "final") continue;

            const std::string& snp1_string = splits[1];
            const std::string& snp2_string = splits[2];
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

            data->snpSetStorage.push_back(set);
        }
    }

    rapidjson::Value ReadMACOEDSets::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("ReadMACOEDSets"), doc.GetAllocator());
        obj.AddMember("path", rapidjson::Value().SetString(path.c_str(), path.size(), doc.GetAllocator()), doc.GetAllocator());
        return obj;
    }
} // epi