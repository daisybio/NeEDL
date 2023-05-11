//
// Created by juli on 25.08.22.
//

#ifndef GENEPISEEKER_SNPCSVANNOTATOR_HPP
#define GENEPISEEKER_SNPCSVANNOTATOR_HPP

#include "Job.hpp"

namespace epi {

    class SnpCsvAnnotator : public Job {
    public:
        SnpCsvAnnotator(
                std::string path,
                size_t snp_column,
                size_t annotation_column,
                bool has_header = false,
                char csv_separator = ',',
                char snp_separator = -1,
                char annotation_separator = -1
        );
        SnpCsvAnnotator(
                std::string path,
                std::string snp_column,
                std::string annotation_column,
                char csv_separator = ',',
                char snp_separator = -1,
                char annotation_separator = -1
        );

        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;

        static SnpCsvAnnotator parse_from_source_string(const std::string &method);
    protected:
        explicit SnpCsvAnnotator(bool has_header = false, char csv_separator = ',', char snp_separator = -1, char annotation_separator = -1);

        void check_columns(size_t snp_col, size_t annotation_col);
        void check_columns(const std::string& snp_col, const std::string& annotation_col);

        std::string path;
        bool has_header = true;
        size_t snp_column, annotation_column;
        char csv_separator, snp_separator, annotation_separator;
    };

} // epi

#ifdef HEADER_ONLY
#include "SnpCsvAnnotator.cpp"
#endif


#endif //GENEPISEEKER_SNPCSVANNOTATOR_HPP
