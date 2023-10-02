//
// Created by juli on 04.08.22.
//

#include "JointDegreeAnalyser.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    JointDegreeAnalyser::JointDegreeAnalyser(std::string name) {
        this->name = name;
    }

    void JointDegreeAnalyser::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("Creating joint degree matrix");

        if (data->outputDirectory == nullptr) {
            throw epi::Error("Output directory need to be specificed to save the network.");
        }

        if (data->snpNetwork == nullptr) {
            throw epi::Error("No network available to analyze.");
        }

        std::vector<SNP_t> snp_list;
        std::unordered_map<SNP_t, size_t, SNP_t::SNPHash> degree_map;
        std::set<size_t> existing_degrees;
        auto snps_iter = data->snpNetwork->all();
        for (auto snp: snps_iter) {
            snp_list.push_back(snp);
            size_t deg = data->snpNetwork->get_degree(snp);
            degree_map.insert({ snp,  deg});
            existing_degrees.insert(deg);
        }

        std::unordered_map<size_t, size_t> degree_index_map;
        size_t deg_i = 0;
        for (auto deg : existing_degrees) {
            degree_index_map.insert({ deg, deg_i ++ });
        }

        // create matrix
        std::vector<size_t> mat (existing_degrees.size() * existing_degrees.size(), 0);
        std::vector<size_t> tot_vec (existing_degrees.size(), 0);
        for (auto & snp_out : snp_list) {
            size_t out_degree = degree_index_map[degree_map[snp_out]];
            tot_vec[out_degree] ++;

            auto neighbours = data->snpNetwork->get_adjacent_snps(snp_out);
            for (auto & snp_in : neighbours) {
                size_t in_degree = degree_index_map[degree_map[snp_in]];
                mat[in_degree * existing_degrees.size() + out_degree] ++;
            }
        }

        // create file stream
        auto file = data->outputDirectory->get_ofstream(name, ".csv");

        // csv-header
        for (auto &deg : existing_degrees) {
            file << deg << '\t';
        }
        file << "TOTAL\n";

        // matrix body
        for (size_t i = 0; i < existing_degrees.size(); i++) {
            for (size_t j = 0; j < existing_degrees.size(); j++) {
                file << mat[i * existing_degrees.size() + j] << '\t';
            }

            // write total
            file << tot_vec[i] << '\n';
        }

        // finalize stuff
        file.close();
        logger.stop();
    }

    rapidjson::Value JointDegreeAnalyser::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("JointDegreeAnalyser"), doc.GetAllocator());
        obj.AddMember("name", rapidjson::Value().SetString(name.c_str(), name.size(), doc.GetAllocator()), doc.GetAllocator());
        return obj;
    }
} // epi