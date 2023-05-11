//
// Created by juli on 25.08.22.
//

#ifndef GENEPISEEKER_NETWORKCSVCONNECTOR_HPP
#define GENEPISEEKER_NETWORKCSVCONNECTOR_HPP

#include "Job.hpp"

namespace epi {

    class NetworkCsvConnector : public Job {
    public:
        NetworkCsvConnector(std::string name, std::string path, size_t column1, size_t column2, bool has_header = false, char csv_separator = ',', char col1_separator = -1, char col2_separator = -1);
        NetworkCsvConnector(std::string name, std::string path, std::string  column1, std::string column2, char csv_separator = ',', char col1_separator = -1, char col2_separator = -1);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;

        std::string get_name() const;

    protected:
        explicit NetworkCsvConnector(std::string name, bool has_header = false, char csv_separator = ',', char col1_separator = -1, char col2_separator = -1);

        void check_columns(size_t col1, size_t col2);
        void check_columns(const std::string& col1, const std::string& col2);

        std::string path;
        std::string name;
        bool has_header = true;
        size_t column1, column2;
        char csv_separator, col1_separator, col2_separator;
    };

} // epi

#ifdef HEADER_ONLY
#include "NetworkCsvConnector.cpp"
#endif


#endif //GENEPISEEKER_NETWORKCSVCONNECTOR_HPP
