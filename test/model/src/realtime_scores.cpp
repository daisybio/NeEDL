//
// Created by juli on 07.06.21.
//

#include <CLI11.hpp>
#include <uwebsockets/App.h>
#include "../../ext/rapidjson/document.h"

#define HEADER_ONLY
#define ALLOW_MULTIPLE_SNPSTORAGES
#include "../../../src/jobs/InstanceLoader.hpp"
#include "../../../src/util/helper_functions.hpp"

// using namespace epi;

std::string json_to_str(const rapidjson::Value &object) {
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    object.Accept(writer);

    return buffer.GetString();
}

int main(int argc, char **argv) {
    // parse args
    CLI::App cli("Realtime scores - A tool to evaluate SNP-sets as a realtime API service");
    std::vector<std::string> gwas_input_format;
    // cli.add_option("--gwas-input-format", gwas_input_format, "The format of the input GWAS file.")->required();
    std::vector<std::string> gwas_input_phenotype;
    // cli.add_option("--gwas-input-phenotype", gwas_input_phenotype, "The phenotype of the GWAS data.")->required();
    std::vector<std::string> gwas_input_file;
    // cli.add_option("--gwas-input-file", gwas_input_file, "The file containing GWAS data.")->required();
    std::vector<std::string> input_name;
    // cli.add_option("--input-name", input_name, "The name used to reference the dataset through the API.")->required();
    std::vector<unsigned int> num_categories;
    // cli.add_option("--num-categories", num_categories, "Number of threads to be used for score calculation");

    std::string input_dir;
    cli.add_option("--input-dir", input_dir, "A directory containing all datasets.")->required();

    std::string bind_addr = "127.0.0.1";
    cli.add_option("--bind-addr", bind_addr, "Address to bind the server to. Default is 127.0.0.1");
    uint16_t bind_port = 3000;
    cli.add_option("--bind-port", bind_port, "Port to bind the server to. Default is 3000");


    std::vector<std::string> models;
    cli.add_option("--model", models,
                   "Add a model to the list of models that are applied on the SNP sets. If none is specified, all available models will be used.");


    unsigned int num_threads = 2;
    cli.add_option("--num-threads", num_threads,
                   "Number of threads to be used for score calculation. Default: 1. 0 means all available");

    // Parse the options.
    CLI11_PARSE(cli, argc, argv);

    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
    }

    boost::filesystem::path gwas_dir_b(input_dir);
    if (is_directory(gwas_dir_b)) {
        for (auto it = directory_iterator(gwas_dir_b); it != directory_iterator(); ++it) {
            path path = it->path();
            if (is_regular_file(path)) {
                std::string file = path.filename().string();

                auto splits = epi::string_split(file, '.');
                if (splits.size() < 3) {
                    throw epi::Error("File " + file + " does not have the right format.");
                } else {
                    std::string name = splits[0];
                    std::string format = splits[1];
                    std::string phenotype = splits[2];
                    unsigned int ncat = 2;
                    if (splits.size() >= 4) {
                        ncat = std::stoi(splits[3]);
                    }

                    std::cout << "dataset: { name: \"" << name << "\", format: \"" << format << "\", phenotype: \""
                              << phenotype << "\", num_categories: " << ncat << " }" << std::endl;

                    gwas_input_file.push_back(path.c_str());
                    gwas_input_format.push_back(format);
                    gwas_input_phenotype.push_back(phenotype);
                    num_categories.push_back(ncat);
                    input_name.push_back(name);
                }
            }
        }
    }


    // process datasets
    auto num_datasets = gwas_input_file.size();

    // process models
    const std::vector<std::string> all_available_models = epi::options::get_all_epistasis_scores();
    if (models.empty()) {
        // set all available models
        models = {all_available_models.begin(), all_available_models.end()};
    }

    std::sort(models.begin(), models.end());
    models.erase(std::unique(models.begin(), models.end()), models.end());

    std::unordered_map<std::string, epi::options::EpistasisScore> models_map;
    for (auto &m: models) models_map.insert({m, epi::options::epistasis_score_from_string(m)});



    // create instances
    std::unordered_map<std::string, std::shared_ptr<epi::DataModel>> instances;
    for (size_t i = 0; i < num_datasets; i++) {
        instances[input_name[i]] = std::make_shared<epi::DataModel>(false);
        auto il = epi::InstanceLoader(gwas_input_file[i], gwas_input_format[i], gwas_input_phenotype[i],
                                      num_categories[i]);
        il.run(instances[input_name[i]]);
        instances[input_name[i]]->snpStorage->init_model_containers(num_threads);
    }

    std::vector<size_t> thread_no (num_threads);
    for (size_t i = 0; i< num_threads; ++i) thread_no[i] = i;

    std::vector<std::thread *> threads(num_threads);
    std::mutex stderr_mutex;

    std::transform(thread_no.begin(), thread_no.end(), threads.begin(), [&stderr_mutex, &instances, &models_map, &bind_addr, &bind_port](size_t thread_index) {
       return new std::thread([&stderr_mutex, &instances, &models_map, &bind_addr, &bind_port, &thread_index]() {
           uWS::App().get("/needl/score/:dataset/:scores/:snps", [&stderr_mutex, &instances, &models_map, &thread_index](auto *res, auto *req) {
               std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();

               std::string dataset = std::string(req->getParameter(0));
               std::string scores_str = std::string(req->getParameter(1));
               std::string snps_str = std::string(req->getParameter(2));

               rapidjson::Document doc;
               auto alloc = doc.GetAllocator();

               // find dataset
               auto dm_ptr = instances.find(dataset);
               if (dm_ptr == instances.end()) {
                   stderr_mutex.lock();
                   std::cerr << "Unknown dataset: " << dataset << std::endl;
                   stderr_mutex.unlock();
                   rapidjson::Value json(rapidjson::kObjectType);
                   json.AddMember("ok", false, alloc);
                   json.AddMember("reason", "unknown dataset", alloc);
                   res->writeStatus("404 Not Found");
                   res->end(json_to_str(json));
               } else {
                   auto &dm = dm_ptr->second;

                   // find the selected models
                   auto scores_s = epi::string_split(scores_str, ';');
                   std::sort(scores_s.begin(), scores_s.end());
                   scores_s.erase(std::unique(scores_s.begin(), scores_s.end()), scores_s.end());

                   rapidjson::Value unknown_scores(rapidjson::kArrayType);
                   std::vector<epi::options::EpistasisScore> selected_models;
                   for (auto &user_score: scores_s) {
                       auto item = models_map.find(user_score);
                       if (item == models_map.end()) {
                           unknown_scores.PushBack(rapidjson::Value(user_score.c_str(), alloc).Move(), alloc);
                       } else {
                           selected_models.push_back(item->second);
                       }
                   }

                   try {
                       auto snps_s = epi::string_split(snps_str, ';');
                       if (snps_s.size() > 10) {
                           stderr_mutex.lock();
                           std::cerr << "Too many SNPs (#snps = " << snps_s.size() << ")" << std::endl;
                           stderr_mutex.unlock();
                           rapidjson::Value json(rapidjson::kObjectType);
                           json.AddMember("ok", false, alloc);
                           json.AddMember("reason", "too many SNPs (max. 10 allowed)", alloc);

                           if (unknown_scores.Capacity() > 0)
                               json.AddMember("unknown_or_unavailable_scores", unknown_scores.Move(), alloc);

                           res->writeStatus("400 Bad Request");
                           res->end(json_to_str(json));
                       } else {
                           std::vector<epi::SNP_t> snps_t;
                           snps_t.reserve(snps_s.size());
                           rapidjson::Value snps_j(rapidjson::kArrayType);
                           for (auto &s: snps_s) {
                               snps_t.push_back(dm->snpStorage->by_name(s));
                               snps_j.PushBack(rapidjson::Value(s.c_str(), alloc).Move(), alloc);
                           }
                           epi::SNPSet set(snps_t);
                           rapidjson::Value json_scores(rapidjson::kObjectType);
                           for (auto &model: selected_models) {
                               json_scores.AddMember(
                                       rapidjson::Value(epi::options::epistasis_score_to_string(model).c_str(),
                                                        alloc).Move(),
                                       dm->snpStorage->calculate_score(set, model, thread_index),
                                       alloc
                               );
                           }
                           rapidjson::Value json(rapidjson::kObjectType);
                           json.AddMember("ok", true, alloc);
                           json.AddMember("scores", json_scores.Move(), alloc);
                           json.AddMember("SNPs", snps_j.Move(), alloc);

                           if (unknown_scores.Capacity() > 0)
                               json.AddMember("unknown_or_unavailable_scores", unknown_scores.Move(), alloc);

                           std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
                           auto timeMs = std::chrono::duration_cast<std::chrono::microseconds>(
                                   endTime - startTime).count();
                           json.AddMember("time_microseconds", timeMs, alloc);
                           res->end(json_to_str(json));
                       }
                   } catch (epi::SNPNotFoundError &err) {
                       stderr_mutex.lock();
                       std::cerr << "Unknown SNP " << err.get_name() << std::endl;
                       stderr_mutex.unlock();
                       rapidjson::Value json(rapidjson::kObjectType);
                       json.AddMember("ok", false, alloc);
                       json.AddMember("reason",
                                      rapidjson::Value(("unknown SNP " + err.get_name()).c_str(), alloc).Move(),
                                      alloc);
                       json.AddMember("unknown_snp", rapidjson::Value(err.get_name().c_str(), alloc).Move(), alloc);
                       if (unknown_scores.Capacity() > 0)
                           json.AddMember("unknown_or_unavailable_scores", unknown_scores.Move(), alloc);
                       res->writeStatus("404 Not Found");
                       res->end(json_to_str(json));
                   } catch (epi::Error &err) {
                       stderr_mutex.lock();
                       std::cerr << "epi::Error  " << err.what() << std::endl;
                       stderr_mutex.unlock();
                       rapidjson::Value json(rapidjson::kObjectType);
                       json.AddMember("ok", false, alloc);
                       json.AddMember("reason", "internal error", alloc);
                       if (unknown_scores.Capacity() > 0)
                           json.AddMember("unknown_or_unavailable_scores", unknown_scores.Move(), alloc);
                       res->writeStatus("500 Internal Error");
                       res->end(json_to_str(json));
                   }
               }
           }).listen(bind_addr, bind_port, [&stderr_mutex, &bind_addr, &bind_port](auto *listen_socket) {
               stderr_mutex.lock();
               if (listen_socket) {
                   std::cout << "[Thread-" << std::this_thread::get_id() << "]  Listening on " << bind_addr << ':' << bind_port << std::endl;
               } else {
                   std::cerr << "[Thread-" << std::this_thread::get_id() << "]  Failed to listen on " << bind_addr << ':' << bind_port << std::endl;
               }
               stderr_mutex.unlock();
           }).run();
       });
    });

    std::for_each(threads.begin(), threads.end(), [](std::thread *t) {
        t->join();
    });

    return 0;
}

