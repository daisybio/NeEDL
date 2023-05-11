//
// Created by juli on 22.06.22.
//

#ifndef GENEPISEEKER_SAVENETWORK_HPP
#define GENEPISEEKER_SAVENETWORK_HPP

#include "Job.hpp"

namespace epi {

    class SaveNetwork : public Job {
    public:
        SaveNetwork(std::string type, std::string name);
        void run(std::shared_ptr<DataModel> data) override;
        rapidjson::Value getConfig(rapidjson::Document &doc) override;
    private:
        static std::string encode_sql_column_name(std::string name);

        enum NETWORK_TYPES {
            ADJACENCY_MATRIX_JSON,
            ADJACENCY_MATRIX_CSV,
            ADJACENCY_LIST_JSON,
            NODE_EDGE_LISTS,
            SQLITE
        };

        NETWORK_TYPES network_type;
        std::string network_type_str;
        std::string name = "network";
    };

} // epi

#ifdef HEADER_ONLY
#include "SaveNetwork.cpp"
#endif


#endif //GENEPISEEKER_SAVENETWORK_HPP
