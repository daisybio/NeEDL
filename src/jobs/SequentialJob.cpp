//
// Created by juli on 25.05.22.
//

#include "SequentialJob.hpp"
#include "../util/Logger.hpp"

#include <utility>
#include <string>

namespace epi {
    SequentialJob::SequentialJob(std::vector<std::shared_ptr<Job>> job_list) {
        this->job_list = std::move(job_list);
    }

    void SequentialJob::add(const std::shared_ptr<Job> &job) {
        this->job_list.push_back(job);
    }

    void SequentialJob::add(std::vector<std::shared_ptr<Job>> job_list_) {
        this->job_list.insert(this->job_list.end(), job_list_.begin(), job_list_.end());
    }

    void SequentialJob::run(std::shared_ptr<DataModel> data) {
        size_t i = 0;
        for (auto &job: this->job_list) {
            Logger::logLine(
                    "Sequential pipeline " + std::to_string(++i) + " of " + std::to_string(this->job_list.size()));
            job->run(data);
        }
    }

    rapidjson::Value SequentialJob::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("SequentialJob"), doc.GetAllocator());

        rapidjson::Value job_li(rapidjson::kArrayType);
        for (size_t i = 0; i < job_list.size(); ++i) {
            auto config = job_list[i]->getConfig(doc);
            job_li.PushBack(config, doc.GetAllocator());
        }
        obj.AddMember("jobs", job_li, doc.GetAllocator());

        return obj;
    }
}
