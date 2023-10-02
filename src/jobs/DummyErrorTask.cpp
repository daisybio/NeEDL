//
// Created by juli on 02.10.23.
//

#include "DummyErrorTask.hpp"

#include "../util/FileFinder.hpp"
#include "../util/helper_functions.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {
    void DummyErrorTask::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger("Run CSVParser to cause the error");

        // extract mappings
        std::vector<std::pair<SNP_t, std::string>> all_annotations;
        CSVParser parser;
        parser.parse(path, csv_separator);


        logger.stop();
    }

    DummyErrorTask::DummyErrorTask(std::string name, std::string path, std::string column1, std::string column2,
                                   char csv_separator, char col1_separator, char col2_separator)
            : DummyErrorTask(name, true, csv_separator, col1_separator, col2_separator) {

        FileFinder ff;
        // ff.add_ends_with(".csv");
        this->path = find_file_by_ending(std::move(path), ff);

        check_columns(column1, column2);
    }

    DummyErrorTask::DummyErrorTask(std::string name, bool has_header, char csv_separator, char col1_separator,
                                   char col2_separator) {
        this->has_header = has_header;
        this->csv_separator = csv_separator;
        this->col1_separator = col1_separator;
        this->col2_separator = col2_separator;
        this->name = name;
    }

    void DummyErrorTask::check_columns(const std::string &col1, const std::string &col2) {
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
} // epi