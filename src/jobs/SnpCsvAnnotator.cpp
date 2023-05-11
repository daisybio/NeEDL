//
// Created by juli on 25.08.22.
//

#include <thread>
#include "SnpCsvAnnotator.hpp"
#include "../util/FileFinder.hpp"
#include "../util/helper_functions.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {

    SnpCsvAnnotator::SnpCsvAnnotator(bool has_header, char csv_separator, char snp_separator,
                                     char annotation_separator) {
        this->has_header = has_header;
        this->csv_separator = csv_separator;
        this->snp_separator = snp_separator;
        this->annotation_separator = annotation_separator;
    }

    SnpCsvAnnotator::SnpCsvAnnotator(std::string path, size_t snp_column, size_t annotation_column, bool has_header,
                                     char csv_separator, char snp_separator, char annotation_separator)
            : SnpCsvAnnotator(has_header, csv_separator, snp_separator, annotation_separator) {

        FileFinder ff;
        ff.add_ends_with(".csv");
        this->path = find_file_by_ending(std::move(path), ff);

        check_columns(snp_column, annotation_column);
    }


    SnpCsvAnnotator::SnpCsvAnnotator(std::string path, std::string snp_column, std::string annotation_column,
                                     char csv_separator, char snp_separator, char annotation_separator)
            : SnpCsvAnnotator(true, csv_separator, snp_separator, annotation_separator) {

        FileFinder ff;
        ff.add_ends_with(".csv");
        this->path = find_file_by_ending(std::move(path), ff);

        check_columns(snp_column, annotation_column);
    }

    void SnpCsvAnnotator::check_columns(size_t snp_col, size_t annotation_col) {
        this->snp_column = snp_col;
        this->annotation_column = annotation_col;
        // check if column indices are valid
        std::ifstream file(this->path);
        std::string first_line;
        std::getline(file, first_line);
        auto splits = string_split(first_line, csv_separator);
        if (snp_column >= splits.size()) {
            throw epi::Error("Invalid column index " + std::to_string(snp_column) + ". Annotation csv file only has " +
                             std::to_string(splits.size()) + " columns.");
        }
        if (annotation_column >= splits.size()) {
            throw epi::Error(
                    "Invalid column index " + std::to_string(annotation_column) + ". Annotation csv file only has " +
                    std::to_string(splits.size()) + " columns.");
        }

        file.close();
    }

    void SnpCsvAnnotator::check_columns(const std::string &snp_col, const std::string &annotation_col) {
        // check if column indices are valid
        std::ifstream file(this->path);
        std::string first_line;
        std::getline(file, first_line);
        auto splits = string_split(first_line, csv_separator);
        auto snp_col_pos = std::find(splits.begin(), splits.end(), snp_col);
        if (snp_col_pos == splits.end()) {
            throw epi::Error("Invalid column name " + snp_col + ".");
        } else this->snp_column = snp_col_pos - splits.begin();

        auto anno_col_pos = std::find(splits.begin(), splits.end(), annotation_col);
        if (anno_col_pos == splits.end()) {
            throw epi::Error("Invalid column name " + annotation_col + ".");
        } else this->annotation_column = anno_col_pos - splits.begin();

        file.close();
    }

    void SnpCsvAnnotator::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("annotating SNPs with CSV map");

        // extract mappings
        std::vector<std::pair<SNP_t, std::string>> all_annotations;
        CSVParser parser;
        parser.parse(path, csv_separator);

        std::vector<std::vector<std::string>> snp_splits (parser.num_rows()), annotation_splits(parser.num_rows());
        if (snp_separator != -1) {
#pragma omp parallel for default(none) shared(has_header, snp_column, snp_separator, parser, snp_splits)
            for (size_t i = has_header ? 1 : 0; i < parser.num_rows(); i++) {
                snp_splits[i] = string_split(parser.cell(i, snp_column), snp_separator);
            }
        }

        if (annotation_separator != -1) {
#pragma omp parallel for default(none) shared(has_header, annotation_column, annotation_separator, parser, annotation_splits)
            for (size_t i = has_header ? 1 : 0; i < parser.num_rows(); i++) {
                annotation_splits[i] = string_split(parser.cell(i, annotation_column), annotation_separator);
            }
        }

        // size_t anno_max = 0, anno_sum = 0, anno_count = 0;

        std::vector<std::string> snps, annotations;
        for (size_t i = has_header ? 1 : 0; i < parser.num_rows(); i++) {
            if (snp_separator == -1) {
                snps.clear();
                snps.push_back(parser.cell(i, snp_column));
            } else {
                // string_split(snps, parser.cell(i, snp_column), snp_separator);
                snps = snp_splits[i];
            }
            if (annotation_separator == -1) {
                annotations.clear();
                annotations.push_back(parser.cell(i, annotation_column));
            } else {
                // string_split(annotations, parser.cell(i, annotation_column), annotation_separator);
                annotations = annotation_splits[i];
            }

            for (auto &snp: snps) {
                if (data->snpStorage->contains_name(snp)) {
                    auto snp_t = data->snpStorage->by_name(snp);

                    for (auto &anno: annotations){
                        /*
                        anno_max = std::max(anno_max, anno.size());
                        anno_sum += anno.size();
                        anno_count ++;
                         */
                        all_annotations.emplace_back(snp_t, anno);
                    }
                }
            }
        }

        // std::cout << anno_max << ", " << anno_count << ", " << anno_sum << std::endl;

        // apply mappings
        data->snpStorage->add_SNP_annotations(all_annotations);

        logger.stop();
    }

    rapidjson::Value SnpCsvAnnotator::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("SnpCsvAnnotator"), doc.GetAllocator());
        obj.AddMember("path", rapidjson::Value().SetString(path.c_str(), path.size(), doc.GetAllocator()), doc.GetAllocator());
        obj.AddMember("has_header", rapidjson::Value().SetBool(has_header), doc.GetAllocator());
        obj.AddMember("snp_column", rapidjson::Value().SetUint64(snp_column), doc.GetAllocator());
        obj.AddMember("annotation_column", rapidjson::Value().SetUint64(annotation_column), doc.GetAllocator());

        std::string csv_sep_str; csv_sep_str += csv_separator;
        obj.AddMember("csv_separator", rapidjson::Value().SetString(csv_sep_str.c_str(), csv_sep_str.size(), doc.GetAllocator()), doc.GetAllocator());

        if (snp_separator == -1) {
            obj.AddMember("snp_separator", rapidjson::Value(), doc.GetAllocator());
        } else {
            std::string col1_sep_str;
            col1_sep_str += snp_separator;
            obj.AddMember("snp_separator",
                          rapidjson::Value().SetString(col1_sep_str.c_str(), col1_sep_str.size(), doc.GetAllocator()),
                          doc.GetAllocator());
        }

        if (annotation_separator == -1) {
            obj.AddMember("annotation_separator", rapidjson::Value(), doc.GetAllocator());
        } else {
            std::string col2_sep_str;
            col2_sep_str += annotation_separator;
            obj.AddMember("annotation_separator",
                          rapidjson::Value().SetString(col2_sep_str.c_str(), col2_sep_str.size(), doc.GetAllocator()),
                          doc.GetAllocator());
        }

        return obj;

    }

    SnpCsvAnnotator SnpCsvAnnotator::parse_from_source_string(const std::string &method) {
        auto escapedCharMap = getEscapedCharMap();
        auto splits = epi::string_split(method, '|');
        if (splits.size() < 4) throw epi::Error("Not enough arguments to specify a snp annotation method");

        std::string path = splits[0];
        bool has_header = epi::toUpperCase(splits[1]) == "YES";
        std::string snp_column = splits[2];
        std::string annotation_column = splits[3];

        bool snp_col_cast = false, anno_col_cast = false;
        size_t snp_col_num = try_parse_number(snp_column, snp_col_cast);
        size_t anno_col_num = try_parse_number(annotation_column, anno_col_cast);

        char csv_sep = ';';
        if (splits.size() > 4 && !splits[4].empty()) {
            if (splits[4].size() > 1) {
                if (escapedCharMap.find(splits[4]) != escapedCharMap.end()) {
                    splits[4][0] = escapedCharMap[splits[4]];
                } else throw epi::Error("Separator can only be a single character");
            }
            csv_sep = splits[4][0];
        }

        char snp_sep = ';';
        if (splits.size() > 5 && !splits[5].empty()) {
            if (escapedCharMap.find(splits[5]) != escapedCharMap.end()) {
                splits[5][0] = escapedCharMap[splits[5]];
            } else throw epi::Error("Separator can only be a single character");
            if (!splits[5].empty()) snp_sep = splits[5][0];
        }

        char anno_sep = ';';
        if (splits.size() > 6 && !splits[6].empty()) {
            if (escapedCharMap.find(splits[6]) != escapedCharMap.end()) {
                splits[6][0] = escapedCharMap[splits[6]];
            } else throw epi::Error("Separator can only be a single character");
            if (!splits[6].empty()) anno_sep = splits[6][0];
        }


        if (snp_col_cast && anno_col_cast) {
            return epi::SnpCsvAnnotator(path, snp_col_num, anno_col_num, has_header, csv_sep, snp_sep, anno_sep);
        } else if (has_header) {
            return epi::SnpCsvAnnotator(path, snp_column, annotation_column, csv_sep, snp_sep, anno_sep);
        } else {
            throw epi::Error("Error in SNP annotation source: A csv without header was given but the columns are referenced by column name.");
            return epi::SnpCsvAnnotator("", 0, 0);
        }
    }


} // epi