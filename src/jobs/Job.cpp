//
// Created by juli on 25.05.22.
//

#include "Job.hpp"
#include "../../ext/rapidjson/filewritestream.h"
#include "../../ext/rapidjson/writer.h"

namespace epi {

    rapidjson::Value Job::getConfig(rapidjson::Document &doc) {
        // return null
        return rapidjson::Value();
    }

    void Job::save_pipeline_config(std::shared_ptr<Job> job, std::string path) {
        rapidjson::Document doc;
        auto config = job->getConfig(doc);
        doc.Swap(config);

        FILE* fp = fopen(path.c_str(), "wb");

        char writeBuffer[65536];
        rapidjson::FileWriteStream os(fp, writeBuffer, sizeof(writeBuffer));

        rapidjson::Writer<rapidjson::FileWriteStream> writer(os);
        doc.Accept(writer);

        fclose(fp);
    }
}