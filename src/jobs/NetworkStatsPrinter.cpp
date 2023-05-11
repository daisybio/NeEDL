//
// Created by juli on 27.05.22.
//

#include "NetworkStatsPrinter.hpp"

namespace epi {
    void NetworkStatsPrinter::run(std::shared_ptr<DataModel> data) {
        if (data->snpNetwork == nullptr) {
            Logger::logLine("network stats: no network created yet");
        } else {
            auto num_nodes = data->snpNetwork->num_nodes();
            auto num_edges = data->snpNetwork->num_edges();
            Logger::logLine("network stats: ");
            Logger::logLine("    nodes:                "+ std::to_string(num_nodes));
            Logger::logLine("    edges:                "+ std::to_string(num_edges));
            if (num_nodes > 1) Logger::logLine("    % of possible edges:  "+ Logger::to_string(100.0 * double(num_edges) / (double(num_nodes - 1) * double(num_nodes) * .5)) + " %");

            // determine the number of annotations
            std::unordered_set<std::string> annotations;
            auto anno_map = data->snpStorage->get_annotations_map();
            for (auto &snp : data->snpNetwork->get_network_snps()) {
                auto annos = data->snpStorage->snp_get_annotations(snp);
                annotations.insert(annos.begin(), annos.end());
            }

            Logger::logLine("    annotations:          " + std::to_string(annotations.size()));

            if (printTimeExpensiveStats) {
                Logger::logLine("    connected:        " + std::string(data->snpNetwork->is_connected() ? "yes" : "no"));
                Logger::logLine("    diameter:         " + std::to_string(data->snpNetwork->get_network_diameter()));
            }
        }
    }

    NetworkStatsPrinter::NetworkStatsPrinter(bool printTimeExpensiveSats) {
        this->printTimeExpensiveStats = printTimeExpensiveSats;
    }

    rapidjson::Value NetworkStatsPrinter::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("NetworkStatsPrinter"), doc.GetAllocator());
        obj.AddMember("print_time_expensive_stats", rapidjson::Value().SetBool(printTimeExpensiveStats), doc.GetAllocator());
        return obj;
    }
} // epi