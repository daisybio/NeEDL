//
// Created by juli on 02.10.23.
//

#ifndef NEEDL_DUMMYERRORTASK_HPP
#define NEEDL_DUMMYERRORTASK_HPP

#include "Job.hpp"

namespace epi {

    class DummyErrorTask : public Job {
    public:
        DummyErrorTask(std::string name, std::string path, std::string  column1, std::string column2, char csv_separator = ',', char col1_separator = -1, char col2_separator = -1);
        void run(std::shared_ptr<DataModel> data) override;

    protected:
        explicit DummyErrorTask(std::string name, bool has_header = false, char csv_separator = ',', char col1_separator = -1, char col2_separator = -1);

        void check_columns(const std::string& col1, const std::string& col2);

        std::string path;
        std::string name;
        bool has_header = true;
        size_t column1, column2;
        char csv_separator, col1_separator, col2_separator;
    };

} // epi

#ifdef HEADER_ONLY
#include "DummyErrorTask.cpp"
#endif

#endif //NEEDL_DUMMYERRORTASK_HPP
