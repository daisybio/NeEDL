//
// Created by juli on 22.06.22.
//

#include "SaveNetwork.hpp"
#include "../util/helper_functions.hpp"
#include "../util/TimeLogger.hpp"
#include <sqlite3.h>

namespace epi {
    SaveNetwork::SaveNetwork(std::string type, std::string name) {
        if (type == "ADJACENCY_MATRIX_JSON") {
           network_type = ADJACENCY_MATRIX_JSON;
        } else if (type == "ADJACENCY_MATRIX_CSV") {
            network_type = ADJACENCY_MATRIX_CSV;
        } else if (type == "ADJACENCY_LIST_JSON") {
            network_type = ADJACENCY_LIST_JSON;
        } else if (type == "NODE_EDGE_LIST") {
            network_type = NODE_EDGE_LISTS;
        } else if (type == "SQLITE") {
            network_type = SQLITE;
        } else {
            throw epi::Error("Unknown network save type " + type);
        }
        this->network_type_str = type;
        this->name = name;
    }

    void SaveNetwork::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("saving network");

        if (data->outputDirectory == nullptr) {
            throw epi::Error("Output directory need to be specificed to save the network.");
        }

        if (data->snpNetwork == nullptr) {
            throw epi::Error("No network available to save.");
        }


        if (network_type == ADJACENCY_MATRIX_JSON) {
            auto adjacency_list = data->snpNetwork->get_adjacency_list();

            std::unordered_map<SNP_t, size_t, SNP_t::SNPHash> snp_index_map;
            for (size_t i = 0; i < adjacency_list.size(); i++) snp_index_map[adjacency_list[i].first] = i;

            auto num_nodes = data->snpNetwork->num_nodes();
            std::vector<bool> adjacency_matrix(num_nodes * num_nodes);
            for (auto & snp_entry : adjacency_list) {
                auto index1 = snp_index_map[snp_entry.first];
                for (auto & adjacent : snp_entry.second) {
                    auto index2 = snp_index_map[adjacent];
                    adjacency_matrix[index1 + index2 * num_nodes] = true;
                    adjacency_matrix[index2 + index1 * num_nodes] = true;
                }
            }

            auto file = data->outputDirectory->get_ofstream(name, ".csv");
            // also save the snp names
            file << "{\"rs_ids\":[";
            for (size_t i = 0; i < num_nodes; i++) {
                if (i > 0) file << ',';
                file << '"' << maskForJson(data->snpStorage->snp_get_name(adjacency_list[i].first)) << '"';
            }
            file << "],\"adjacency_matrix\":[";
            for (size_t i = 0; i < num_nodes; i++) {
                file << '[';
                if (i > 0) file << ',';
                for (size_t j = 0; j < num_nodes; j++) {
                    if (j > 0) file << ',';
                    file << (adjacency_matrix[i + j * num_nodes] ? 1 : 0);
                }
                file << ']';
            }
            file << "]}";
            file.close();
        } else if (network_type == ADJACENCY_MATRIX_CSV) {
            auto adjacency_list = data->snpNetwork->get_adjacency_list();

            std::unordered_map<SNP_t, size_t, SNP_t::SNPHash> snp_index_map;
            for (size_t i = 0; i < adjacency_list.size(); i++) snp_index_map[adjacency_list[i].first] = i;

            auto num_nodes = data->snpNetwork->num_nodes();
            std::vector<bool> adjacency_matrix(num_nodes * num_nodes);
            for (auto & snp_entry : adjacency_list) {
                auto index1 = snp_index_map[snp_entry.first];
                for (auto & adjacent : snp_entry.second) {
                    auto index2 = snp_index_map[adjacent];
                    adjacency_matrix[index1 + index2 * num_nodes] = true;
                    adjacency_matrix[index2 + index1 * num_nodes] = true;
                }
            }

            auto file = data->outputDirectory->get_ofstream(name, ".csv");
            file << "RS_ID\t";
            for (size_t i = 0; i < num_nodes; i++) {
                if (i > 0) file << ',';
                file << data->snpStorage->snp_get_name(adjacency_list[i].first);
            }
            file << '\n';

            for (size_t i = 0; i < num_nodes; i++) {
                file << data->snpStorage->snp_get_name(adjacency_list[i].first);
                for (size_t j = 0; j < num_nodes; j++) {
                    file << '\t' << (adjacency_matrix[i + j * num_nodes] ? '1' : '0');
                }
                file << '\n';
            }
            file.close();
        } else if (network_type == ADJACENCY_LIST_JSON) {
            auto file = data->outputDirectory->get_ofstream(name, ".json");
            auto network_snps = data->snpNetwork->get_network_snps();
            file << '{';
            bool is_first = true;
            for (auto & snp : network_snps) {
                if (is_first) is_first = false;
                else file << ',';
                file << '"' << data->snpStorage->snp_get_name(snp) << "\":[";
                auto adjacent_nodes = data->snpNetwork->get_adjacent_snps(snp);
                for (size_t j = 0; j < adjacent_nodes.size(); j++) {
                    if (j > 0) file << ',';
                     file << '"' << data->snpStorage->snp_get_name(adjacent_nodes[j]) << '"';
                }
                file << ']';
            }
            file << '}';
            file.close();
        } else if (network_type == NODE_EDGE_LISTS) {
            // node list
            auto node_file = data->outputDirectory->get_ofstream(name, "_nodes.csv");
            node_file << "ID (BASE32)\tRS_ID\tAnnotations";

            std::vector<SNP_t> nodes;
            std::set<std::string> attrib_names;
            auto network_snps = data->snpNetwork->get_network_snps();
            for (auto &node: network_snps) {
                nodes.push_back(node);
                auto keys = data->snpStorage->snp_get_variable_attribute_keys(node);
                attrib_names.insert(keys.begin(), keys.end());
            }

            for (auto &attrib: attrib_names) {
                node_file << '\t' << attrib;
            }
            node_file << '\n';

            for (size_t i = 0; i < nodes.size(); i++) {
                node_file << numberToBase32(nodes[i].value) << '\t';
                auto annos = data->snpStorage->snp_get_annotations(nodes[i]);
                node_file << data->snpStorage->snp_get_name(nodes[i]) << '\t';
                bool is_first = true;
                for (auto &anno: annos) {
                    if (is_first) is_first = false;
                    else node_file << ';';
                    node_file << anno;
                }

                for (auto &attrib_name: attrib_names) {
                    node_file << '\t' << data->snpStorage->snp_get_variable_attribute(nodes[i], attrib_name);
                }
                node_file << '\n';
            }
            node_file.close();
            // edge list
            auto edge_file = data->outputDirectory->get_ofstream(name, "_edges.json");

            auto edge_labels = data->snpNetwork->get_all_edge_tags();
            edge_file << "{\"labels\":[";
            bool first = true;
            for (const auto &item: edge_labels) {
                if (first) first = false;
                else edge_file << ',';
                edge_file << '"' << maskForJson(item) << '"';
            }
            edge_file << "],\"edges\":{";

            // std::unordered_set<SNPEdge, SNPEdgeHash> already_printed;

            bool first_edge = true;
            for (auto &node: network_snps) {
                std::string from = numberToBase32(node.value);
                auto adjacent_nodes = data->snpNetwork->get_adjacent_snps(node);
                bool printed_start = false;
                for (auto &adjacent_node: adjacent_nodes) {
                    // SNPEdge curr (node, adjacent_node);
                    // if (already_printed.find(curr) != already_printed.end()) continue;
                    // already_printed.insert(curr);

                    if (!printed_start) {
                        printed_start = true;

                        if (first_edge) first_edge = false;
                        else edge_file << "},";

                        edge_file << '"' << from << "\":{";
                    } else edge_file << ',';


                    std::string to = numberToBase32(adjacent_node.value);
                    auto labels = data->snpNetwork->get_edge_label_ids({node, adjacent_node});
                    edge_file << '"' << to << "\":[";

                    bool first_label = true;
                    for (auto &l: labels) {
                        if (first_label) first_label = false;
                        else edge_file << ',';
                        edge_file << l;
                    }

                    edge_file << "]";
                }
            }
            edge_file << "}}}";
            edge_file.close();
            /*
            // edge list
            auto edge_file = data->outputDirectory->get_ofstream(name, "_edges.csv");
            edge_file << "SNP1\tSNP2\ttag\n";

            std::unordered_set<SNPEdge, SNPEdgeHash> already_printed;

            for (auto & node : network_snps) {
                std::string from = numberToBase32(node.value);
                auto adjacent_nodes = data->snpNetwork->get_adjacent_snps(node);
                for (auto & adjacent_node : adjacent_nodes) {
                    SNPEdge curr (node, adjacent_node);
                    if (already_printed.find(curr) != already_printed.end()) continue;
                    already_printed.insert(curr);
                    std::string to = numberToBase32(adjacent_node.value);
                    auto labels = data->snpNetwork->get_edge_labels({node, adjacent_node });
                    std::string labels_str = "";
                    bool first_label = true;
                    for (auto & l : labels) {
                        if (first_label) first_label = false;
                        else labels_str += ';';
                        labels_str += l;
                    }
                    edge_file << from << '\t' << to << '\t' << labels_str << '\n';
                }
            }
            edge_file.close();
             */
        } else if (network_type == SQLITE) {
            char *zErrMsg = nullptr;
            int rc;

            // open db file
            auto db_file = data->outputDirectory->get_free_filepath(name, ".sqlite3");
            size_t num_threads = 1; // omp_get_max_threads();
            sqlite3 **db = new sqlite3*[num_threads];
            for (size_t thri = 0; thri < num_threads; thri++) {
                rc = sqlite3_open_v2(db_file.c_str(), &db[thri], SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_NOMUTEX, nullptr);
                if (rc) {
                    throw epi::Error("Cannot create sqlite3 database at " + db_file);
                }

                // configure DB
                std::string configuration_stmt = "PRAGMA journal_mode = WAL;"
                                                 "PRAGMA synchronous = FULL;"
                                                 "PRAGMA locking_mode = NORMAL;"
                                                 "PRAGMA cache_size = -1000000";
                rc = sqlite3_exec(db[thri], configuration_stmt.c_str(), nullptr, nullptr, &zErrMsg);
                if (rc) {
                    throw epi::Error("SQLITE-error: " + std::string(zErrMsg));
                }
            }

            // create tables
            std::vector<SNP_t> nodes;
            std::set<std::string> attrib_names;
            auto network_snps = data->snpNetwork->get_network_snps();
            for (auto &node: network_snps) {
                nodes.push_back(node);
                auto keys = data->snpStorage->snp_get_variable_attribute_keys(node);
                attrib_names.insert(keys.begin(), keys.end());
            }

            std::string create_stmt = "CREATE TABLE node_annotations ("
                                      "id           INT         PRIMARY KEY NOT NULL,"
                                      "name         VARCHAR     NOT NULL"
                                      ") WITHOUT ROWID;"
                                      "CREATE TABLE has_annotation ("
                                      "node         INT         NOT NULL,"
                                      "annotation   INT         NOT NULL,"
                                      "FOREIGN KEY(node) REFERENCES nodes(id),"
                                      "FOREIGN KEY(annotation) REFERENCES node_annotation(id),"
                                      "PRIMARY KEY (node, annotation)"
                                      ") WITHOUT ROWID;"
                                      "CREATE TABLE nodes ("
                                      "id           INT         PRIMARY KEY NOT NULL,"
                                      "name         VARCHAR     NOT NULL";
            for (auto & key : attrib_names) {
                create_stmt += ",\n" + encode_sql_column_name(key) + " VARCHAR";
            }
            create_stmt += ") WITHOUT ROWID;"
                           "CREATE TABLE edges ("
                           "node1       INT     NOT NULL,"
                           "node2       INT     NOT NULL,";

            auto edge_labels = data->snpNetwork->get_all_edge_tags();
            for (auto &label : edge_labels) {
                create_stmt += encode_sql_column_name(label) + " BOOLEAN,";
            }

            create_stmt += "FOREIGN KEY(node1) REFERENCES nodes(id),"
                           "FOREIGN KEY(node2) REFERENCES nodes(id),"
                           "PRIMARY KEY (node1, node2)"
                           ") WITHOUT ROWID;";

            rc = sqlite3_exec(db[0], create_stmt.c_str(), nullptr, nullptr, &zErrMsg);
            if (rc) {
                throw epi::Error("SQLITE-error: " + std::string(zErrMsg));
            }

            std::string begin_transaction = "BEGIN";
            rc = sqlite3_exec(db[0], begin_transaction.c_str(), nullptr, nullptr, &zErrMsg);
            if (rc) {
                throw epi::Error("SQLITE-error: " + std::string(zErrMsg));
            }

            // insert nodes
            std::unordered_set<std::string> network_node_annotations;
            std::string node_insert_stmt = "INSERT INTO nodes (id, name ";
            for (auto &key : attrib_names) node_insert_stmt += ", " + encode_sql_column_name(key);
            node_insert_stmt += ") VALUES (?1, ?2";
            for (size_t i = 0; i < attrib_names.size(); i++) node_insert_stmt += ", ?" + std::to_string(i + 3);
            node_insert_stmt += ")";
            sqlite3_stmt *node_insert_stmt_sqlite;
            rc = sqlite3_prepare_v2(db[0], node_insert_stmt.c_str(), -1, &node_insert_stmt_sqlite, nullptr);
            if (rc) {
                throw epi::Error("SQLITE-error when preparing node insert statement");
            }
            for (auto &node : nodes) {
                sqlite3_bind_int(node_insert_stmt_sqlite, 1, node.value);
                auto node_name = data->snpStorage->snp_get_name(node);
                sqlite3_bind_text(node_insert_stmt_sqlite, 2, node_name.c_str(), node_name.size(), SQLITE_STATIC);

                auto annos = data->snpStorage->snp_get_annotations(node);
                network_node_annotations.insert(annos.begin(), annos.end());

                int parameter_num = 3;
                std::vector<std::string> attrib_data;
                for (auto &attrib_name: attrib_names) {
                    auto attrib_value = data->snpStorage->snp_get_variable_attribute(node, attrib_name);
                    if (attrib_value.empty()) {
                        sqlite3_bind_null(node_insert_stmt_sqlite, parameter_num);
                    } else {
                        attrib_data.push_back(attrib_value);
                        sqlite3_bind_text(node_insert_stmt_sqlite, parameter_num, attrib_data.back().c_str(), attrib_data.back().size(), SQLITE_STATIC);
                    }

                    parameter_num ++;
                }
                sqlite3_step(node_insert_stmt_sqlite);
                sqlite3_reset(node_insert_stmt_sqlite);
            }
            sqlite3_finalize(node_insert_stmt_sqlite);

            // insert node annotations
            std::string insert_annos_stmt = "INSERT INTO node_annotations (id, name) VALUES (?1, ?2)";
            sqlite3_stmt *insert_annos_stmt_sqlite;
            rc = sqlite3_prepare_v2(db[0], insert_annos_stmt.c_str(), -1, &insert_annos_stmt_sqlite, nullptr);
            if (rc) {
                throw epi::Error("SQLITE-error when preparing node_annotation insert statement");
            }
            std::string has_anno_stmt = "INSERT INTO has_annotation (node, annotation) VALUES (?1, ?2)";
            sqlite3_stmt *has_anno_stmt_sqlite;
            rc = sqlite3_prepare_v2(db[0], has_anno_stmt.c_str(), -1, &has_anno_stmt_sqlite, nullptr);
            if (rc) {
                throw epi::Error("SQLITE-error when preparing has_annotation insert statement");
            }
            auto all_annotations = data->snpStorage->get_annotations_map();

            size_t anno_id = 0;
            for (auto &anno : network_node_annotations) {
                sqlite3_bind_int(insert_annos_stmt_sqlite, 1, anno_id);
                sqlite3_bind_text(insert_annos_stmt_sqlite, 2, anno.c_str(), anno.size(), SQLITE_STATIC);
                sqlite3_step(insert_annos_stmt_sqlite);

                for (auto & node : all_annotations[anno]) {
                    sqlite3_bind_int(has_anno_stmt_sqlite, 1, node);
                    sqlite3_bind_int(has_anno_stmt_sqlite, 2, anno_id);
                    sqlite3_step(has_anno_stmt_sqlite);
                    sqlite3_reset(has_anno_stmt_sqlite);
                }

                sqlite3_reset(insert_annos_stmt_sqlite);
                anno_id++;
            }
            sqlite3_finalize(insert_annos_stmt_sqlite);
            sqlite3_finalize(has_anno_stmt_sqlite);


            std::string commit_transaction = "COMMIT";
            rc = sqlite3_exec(db[0], commit_transaction.c_str(), nullptr, nullptr, &zErrMsg);
            if (rc) {
                throw epi::Error("SQLITE-error: " + std::string(zErrMsg));
            }

            // insert edges and has_label
            std::string insert_edges_stmt = "INSERT INTO edges (node1, node2";
            for (size_t i = 0; i < edge_labels.size(); i++) insert_edges_stmt += ", " + encode_sql_column_name(edge_labels[i]);
            insert_edges_stmt += ") VALUES (?1, ?2";
            for (size_t i = 0; i < edge_labels.size(); i++) insert_edges_stmt += ", ?" + std::to_string(i + 3);
            insert_edges_stmt += ")";
            auto **insert_edges_stmt_sqlite = new sqlite3_stmt*[num_threads];
            for (size_t thri = 0; thri < num_threads; thri++) {
                rc = sqlite3_exec(db[thri], begin_transaction.c_str(), nullptr, nullptr, &zErrMsg);
                if (rc) {
                    throw epi::Error("SQLITE-error: " + std::string(zErrMsg));
                }
                rc = sqlite3_prepare_v2(db[thri], insert_edges_stmt.c_str(), -1, &insert_edges_stmt_sqlite[thri], nullptr);
                if (rc) {
                    throw epi::Error("SQLITE-error when preparing edges insert statement");
                }
            }

// #pragma omp parallel for default(none) shared(db, network_snps, data, insert_edges_stmt_sqlite, has_label_stmt_sqlite, edge_num, std::cout, commit_transaction, rc, begin_transaction, zErrMsg) private(edge_range_end, curr_edge_num)
            for (auto snp_iter = network_snps.begin(); snp_iter != network_snps.end(); snp_iter++) {
                size_t thri = omp_get_thread_num();

                auto node = *snp_iter;
                auto adjacent_nodes = data->snpNetwork->get_adjacent_edges(node);
                for (auto &adjacent_edge: adjacent_nodes) {
                    if (adjacent_edge.first.bits.snp1 != node.value) continue;
                    // if (node.value > adjacent_node.value) continue;
                    sqlite3_bind_int(insert_edges_stmt_sqlite[thri], 1, adjacent_edge.first.bits.snp1);
                    sqlite3_bind_int(insert_edges_stmt_sqlite[thri], 2, adjacent_edge.first.bits.snp2);

                    std::vector<bool> contain_label (edge_labels.size(), false);

                    // also insert all labels
                    auto labels = adjacent_edge.second;
                    for (auto &label : labels) {
                        contain_label[label] = true;
                    }

                    for (size_t i = 0; i < contain_label.size(); i++)
                        sqlite3_bind_int(insert_edges_stmt_sqlite[thri], i + 3, contain_label[i] ? 1 : 0);

                    sqlite3_step(insert_edges_stmt_sqlite[thri]);
                    sqlite3_reset(insert_edges_stmt_sqlite[thri]);
                }
            }


            for (size_t thri = 0; thri < num_threads; thri++) {
                rc = sqlite3_exec(db[thri], commit_transaction.c_str(), nullptr, nullptr, &zErrMsg);
                if (rc) {
                    throw epi::Error("SQLITE-error: " + std::string(zErrMsg));
                }
                sqlite3_finalize(insert_edges_stmt_sqlite[thri]);

                sqlite3_close_v2(db[thri]);
            }
            delete[] insert_edges_stmt_sqlite;
            delete[] db;
        }

        logger.stop();
    }

    std::string SaveNetwork::encode_sql_column_name(std::string name) {
        std::stringstream res;
        res << '"';
        for (auto c : name) {
            if ((c >= 65 && c <= 90) // A-Z
            || (c >= 97 && c <= 122) // a-z
            || (c >= 48 && c <= 57) // 0-9
            || c == 45) // -
                res << c;
            else res << '_';
        }
        res << '"';
        return res.str();
    }

    rapidjson::Value SaveNetwork::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("SaveNetwork"), doc.GetAllocator());
        obj.AddMember("network_type", rapidjson::Value().SetString(network_type_str.c_str(), network_type_str.size(), doc.GetAllocator()), doc.GetAllocator());
        obj.AddMember("name", rapidjson::Value().SetString(name.c_str(), name.size(), doc.GetAllocator()), doc.GetAllocator());
        return obj;
    }
} // epi