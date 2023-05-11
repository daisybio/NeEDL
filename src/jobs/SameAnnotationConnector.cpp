//
// Created by juli on 26.05.22.
//

#include "SameAnnotationConnector.hpp"
#include "../util/TimeLogger.hpp"
#include "../data_model/SNPNetwork.hpp"

namespace epi {
    void SameAnnotationConnector::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("connecting items with same label");

        if (data->snpNetwork == nullptr) {
            data->snpNetwork = std::make_shared<SNPNetwork>();
        }

        std::unordered_set<SNP_t, SNP_t::SNPHash> nodes;
        std::vector<SNPEdge> edges;

        auto annotations_map = data->snpStorage->get_annotations_map();
        for (auto & anno : annotations_map) {
            std::vector<SNP_t> nodes_repr;
            for (auto & node : anno.second) {
                auto no = data->snpStorage->by_instance_id(node);
                if (!data->snpStorage->snp_is_removed(no)) {
                    nodes_repr.push_back(no);
                }
            }

            if (nodes_repr.size() > 1) {
                // add nodes to network
                for (auto & node : nodes_repr) {
                    nodes.insert(node);
                }
                for (size_t i = 0; i < nodes_repr.size(); i++) {
                    for (size_t j = i + 1; j < nodes_repr.size(); j++) {
                        edges.emplace_back( nodes_repr[i], nodes_repr[j] );
                    }
                }
            }
        }

        // insert nodes and edges
        data->snpNetwork->add_nodes(nodes.begin(), nodes.end());
        data->snpNetwork->add_edges(edges.begin(), edges.end(), "SAME_TAG");

        logger.stop();
    }

    rapidjson::Value SameAnnotationConnector::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("SameAnnotationConnector"), doc.GetAllocator());
        return obj;
    }
} // epi