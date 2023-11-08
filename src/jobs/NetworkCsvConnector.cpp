//
// Created by juli on 25.08.22.
//

#include "NetworkCsvConnector.hpp"
#include "../util/FileFinder.hpp"
#include "../util/helper_functions.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    NetworkCsvConnector::NetworkCsvConnector(std::string name, std::string path, size_t column1, size_t column2, bool has_header,
                                             char csv_separator, char col1_separator, char col2_separator)
            : NetworkCsvConnector(name, has_header, csv_separator, col1_separator, col2_separator) {

        FileFinder ff;
        ff.add_ends_with(".csv");
        this->path = find_file_by_ending(std::move(path), ff);

        check_columns(column1, column2);
    }

    NetworkCsvConnector::NetworkCsvConnector(std::string name, std::string path, std::string column1, std::string column2,
                                             char csv_separator, char col1_separator, char col2_separator)
            : NetworkCsvConnector(name, true, csv_separator, col1_separator, col2_separator) {

        FileFinder ff;
        ff.add_ends_with(".csv");
        this->path = find_file_by_ending(std::move(path), ff);

        check_columns(column1, column2);
    }

    NetworkCsvConnector::NetworkCsvConnector(std::string name, bool has_header, char csv_separator, char col1_separator,
                                             char col2_separator) {
        this->has_header = has_header;
        this->csv_separator = csv_separator;
        this->col1_separator = col1_separator;
        this->col2_separator = col2_separator;
        this->name = name;
    }

    void NetworkCsvConnector::check_columns(size_t col1, size_t col2) {
        this->column1 = col1;
        this->column2 = col2;
        // check if column indices are valid
        std::ifstream file(this->path);
        std::string first_line;
        std::getline(file, first_line);
        auto splits = string_split(first_line, csv_separator);
        if (column1 >= splits.size()) {
            throw epi::Error("Invalid column index " + std::to_string(column1) + ". Network csv file only has " +
                             std::to_string(splits.size()) + " columns.");
        }
        if (column2 >= splits.size()) {
            throw epi::Error(
                    "Invalid column index " + std::to_string(column2) + ". Network csv file only has " +
                    std::to_string(splits.size()) + " columns.");
        }

        file.close();
    }

    void NetworkCsvConnector::check_columns(const std::string &col1, const std::string &col2) {
        // check if column indices are valid
        std::ifstream file(this->path);
        std::string first_line;
        std::getline(file, first_line);
        auto splits = string_split(first_line, csv_separator);
        auto col1_pos = std::find(splits.begin(), splits.end(), col1);
        if (col1_pos == splits.end()) {
            throw epi::Error("Invalid column name " + col1 + ".");
        } else this->column1 = col1_pos - splits.begin();

        auto col2_pos = std::find(splits.begin(), splits.end(), col2);
        if (col2_pos == splits.end()) {
            throw epi::Error("Invalid column name " + col2 + ".");
        } else this->column2 = col2_pos - splits.begin();

        file.close();
    }

    void NetworkCsvConnector::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("connecting SNPs with CSV map " + name);

        if (data->snpNetwork == nullptr) {
            data->snpNetwork = std::make_shared<SNPNetwork>();
        }

        // extract mappings
        std::vector<std::pair<SNP_t, std::string>> all_annotations;
        CSVParser parser;
        parser.parse(path, csv_separator);

        std::vector<std::vector<std::string>> col1_splits (parser.num_rows()), col2_splits(parser.num_rows());
        if (col1_separator != -1) {
#pragma omp parallel for default(none) shared(has_header, column1, col1_separator, parser, col1_splits)
            for (size_t i = has_header ? 1 : 0; i < parser.num_rows(); i++) {
                col1_splits[i] = string_split(parser.cell(i, column1), col1_separator);
            }
        }

        if (col2_separator != -1) {
#pragma omp parallel for default(none) shared(has_header, column2, col2_separator, parser, col2_splits)
            for (size_t i = has_header ? 1 : 0; i < parser.num_rows(); i++) {
                col2_splits[i] = string_split(parser.cell(i, column2), col2_separator);
            }
        }

        size_t num_threads = omp_get_max_threads();

        std::vector<std::vector<SNP_t>> nodes (num_threads);
        std::vector<std::vector<SNPEdge>> edges (num_threads);

        std::vector<std::vector<std::string>> col1 (num_threads), col2 (num_threads);
#pragma omp parallel for default(none) shared(col1, col2, parser, col1_splits, col2_splits, data, nodes, edges)
        for (size_t i = has_header ? 1 : 0; i < parser.num_rows(); i++) {
            size_t thr = omp_get_thread_num();
            if (col1_separator == -1) {
                col1[thr].clear();
                col1[thr].push_back(parser.cell(i, column1));
            } else {
                // string_split(snps, parser.cell(i, snp_column), snp_separator);
                col1[thr] = col1_splits[i];
            }
            if (col2_separator == -1) {
                col2[thr].clear();
                col2[thr].push_back(parser.cell(i, column2));
            } else {
                // string_split(annotations, parser.cell(i, annotation_column), annotation_separator);
                col2[thr] = col2_splits[i];
            }

            for (auto & gene_symbol_1 : col1[thr]) {
                for (auto & gene_symbol_2 : col2[thr]) {

                    if (gene_symbol_1 != gene_symbol_2) {

                        std::vector<SNP_t> snps_symbol_1 = data->snpStorage->by_annotation(gene_symbol_1);
                        std::vector<SNP_t> snps_symbol_2 = data->snpStorage->by_annotation(gene_symbol_2);

                        if (!snps_symbol_1.empty() && !snps_symbol_2.empty()) {
                            nodes[thr].insert(nodes[thr].end(), snps_symbol_1.begin(), snps_symbol_1.end());
                            nodes[thr].insert(nodes[thr].end(), snps_symbol_2.begin(), snps_symbol_2.end());

                            for (auto &snp1: snps_symbol_1) {
                                for (auto &snp2: snps_symbol_2) {
                                    edges[thr].emplace_back(snp1, snp2);
                                }
                            }
                        }
                    }
                }
            }

            // save partial result if we already have 1.000.000 edges
            if (edges[thr].size() > 10000000) {
#pragma omp critical
                {
                    data->snpNetwork->add_nodes(nodes[thr].begin(), nodes[thr].end());
                    data->snpNetwork->add_edges(edges[thr].begin(), edges[thr].end(), name);
                }

                edges[thr].clear();
                nodes[thr].clear();
            }
        }

        // insert nodes and edges
        // std::unordered_set<SNP_t, SNP_t::SNPHash> nodes_set;
        // for(auto & x : nodes) nodes_set.insert(x);
        for (size_t i = 0; i < num_threads; i++) {
            data->snpNetwork->add_nodes(nodes[i].begin(), nodes[i].end());
            data->snpNetwork->add_edges(edges[i].begin(), edges[i].end(), name);
        }


        logger.stop();
    }

    std::string NetworkCsvConnector::get_name() const {
        return name;
    }

    rapidjson::Value NetworkCsvConnector::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("NetworkCsvConnector"), doc.GetAllocator());
        obj.AddMember("name", rapidjson::Value().SetString(name.c_str(), name.size(), doc.GetAllocator()), doc.GetAllocator());
        obj.AddMember("path", rapidjson::Value().SetString(path.c_str(), path.size(), doc.GetAllocator()), doc.GetAllocator());
        obj.AddMember("has_header", rapidjson::Value().SetBool(has_header), doc.GetAllocator());
        obj.AddMember("column1", rapidjson::Value().SetUint64(column1), doc.GetAllocator());
        obj.AddMember("column2", rapidjson::Value().SetUint64(column2), doc.GetAllocator());
        std::string csv_sep_str; csv_sep_str += csv_separator;
        obj.AddMember("csv_separator", rapidjson::Value().SetString(csv_sep_str.c_str(), csv_sep_str.size(), doc.GetAllocator()), doc.GetAllocator());

        if (col1_separator == -1) {
            obj.AddMember("col1_separator", rapidjson::Value(), doc.GetAllocator());
        } else {
            std::string col1_sep_str;
            col1_sep_str += col1_separator;
            obj.AddMember("col1_separator",
                          rapidjson::Value().SetString(col1_sep_str.c_str(), col1_sep_str.size(), doc.GetAllocator()),
                          doc.GetAllocator());
        }

        if (col2_separator == -1) {
            obj.AddMember("col2_separator", rapidjson::Value(), doc.GetAllocator());
        } else {
            std::string col2_sep_str;
            col2_sep_str += col2_separator;
            obj.AddMember("col2_separator",
                          rapidjson::Value().SetString(col2_sep_str.c_str(), col2_sep_str.size(), doc.GetAllocator()),
                          doc.GetAllocator());
        }

        return obj;
    }
} // epi