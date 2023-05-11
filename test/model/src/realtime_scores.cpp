//
// Created by juli on 07.06.21.
//

#include <CLI11.hpp>
#include <drogon/drogon.h>

#define HEADER_ONLY
#define ALLOW_MULTIPLE_SNPSTORAGES
#include "../../../src/jobs/InstanceLoader.hpp"
#include "../../../src/util/helper_functions.hpp"

// using namespace epi;

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


    int num_threads = 2;
    cli.add_option("--num-threads", num_threads,
                   "Number of threads to be used for score calculation. Default: 2. 0 means all available");

    // Parse the options.
    CLI11_PARSE(cli, argc, argv);

    // process num threads
    drogon::app().setThreadNum(num_threads);
    num_threads = drogon::app().getThreadNum();

    boost::filesystem::path gwas_dir_b(input_dir);
    if(is_directory(gwas_dir_b)) {
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

                    LOG_WARN << "dataset: { name: \"" << name << "\", format: \"" << format << "\", phenotype: \"" << phenotype << "\", num_categories: " << ncat << " }";

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

    /*
    // process num_categories
    if (num_categories.empty()) {
        num_categories.push_back(2);
    }
    if (num_categories.size() == 1) {
        for (size_t i = 1; i < num_datasets; i++) num_categories.push_back(num_categories[0]);
    }

    // process input format
    if (gwas_input_format.size() == 1) {
        for (size_t i = 1; i < num_datasets; i++) gwas_input_format.push_back(gwas_input_format[0]);
    }

    // process phenotype
    if (gwas_input_phenotype.size() == 1) {
        for (size_t i = 1; i < num_datasets; i++) gwas_input_phenotype.push_back(gwas_input_phenotype[0]);
    }

    if (
            input_name.size() != num_datasets ||
            num_categories.size() != num_datasets ||
            gwas_input_phenotype.size() != num_datasets ||
            gwas_input_format.size() != num_datasets
       ) {
        throw epi::Error("You did not provide the right amount of arguments of each type.");
    }

     */


    // process models
    const std::vector<std::string> all_available_models =  { "BAYESIAN",
                   "PENETRANCE_NLL", "PENETRANCE_LLH", "PENETRANCE_AIC", "PENETRANCE_BIC",
                   "REGRESSION_NLL", "REGRESSION_LLH", "REGRESSION_AIC", "REGRESSION_BIC",
                   "REGRESSION_NLL-GAIN", "REGRESSION_LLH-GAIN", "REGRESSION_AIC-GAIN", "REGRESSION_BIC-GAIN",
                   "VARIANCE" };
    if (models.empty()) {
        // set all available models
        models = { all_available_models.begin(), all_available_models.end() };
    }

    std::sort(models.begin(), models.end());
    models.erase(std::unique(models.begin(), models.end()), models.end());

    std::unordered_map<std::string, epi::options::EpistasisScore> models_map;
    for (auto & m : models) models_map.insert({ m, epi::options::epistasis_score_from_string(m) });



    // create instances
    std::unordered_map<std::string, std::shared_ptr<epi::DataModel>> instances;
    for (size_t i = 0; i < num_datasets; i++) {
        instances[input_name[i]] = std::make_shared<epi::DataModel>(false);
        auto il = epi::InstanceLoader(gwas_input_file[i], gwas_input_format[i], gwas_input_phenotype[i], num_categories[i]);
        il.run(instances[input_name[i]]);
        instances[input_name[i]]->snpStorage->init_model_containers(num_threads);
    }

    drogon::app().registerHandler(
            "/needl/score/{dataset}/{scores}/{snps}",
            [&instances, &models_map](const drogon::HttpRequestPtr &,
               std::function<void(const drogon::HttpResponsePtr &)> &&callback,
               const std::string &dataset, const std::string &scores_str, const std::string &snps_str) {

                std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();

                // find dataset
                auto dm_ptr = instances.find(dataset);
                if (dm_ptr == instances.end()) {
                    LOG_WARN << "Unknown dataset: " << dataset;
                    Json::Value json;
                    json["ok"] = false;
                    json["reason"] = "unknown dataset";
                    auto resp = drogon::HttpResponse::newHttpJsonResponse(json);
                    resp->setStatusCode(drogon::k404NotFound);
                    callback(resp);
                } else {
                    auto &dm = dm_ptr->second;

                    // find the selected models
                    auto scores_s = epi::string_split(scores_str, ';');
                    std::sort(scores_s.begin(), scores_s.end());
                    scores_s.erase(std::unique(scores_s.begin(), scores_s.end()), scores_s.end());

                    Json::Value unknown_scores;
                    std::vector<epi::options::EpistasisScore> selected_models;
                    for (auto &user_score : scores_s) {
                        auto item = models_map.find(user_score);
                        if (item == models_map.end()) {
                            unknown_scores.append(user_score);
                        } else {
                            selected_models.push_back(item->second);
                        }
                    }

                    try {
                        auto snps_s = epi::string_split(snps_str, ';');
                        if (snps_s.size() > 10) {
                            LOG_WARN << "Too many SNPs (#snps = " << snps_s.size() << ")";
                            Json::Value json;
                            json["ok"] = false;
                            json["reason"] = "too many SNPs (max. 10 allowed)";
                            if (unknown_scores.size() > 0) json["unknown_or_unavailable_scores"] = unknown_scores;
                            auto resp = drogon::HttpResponse::newHttpJsonResponse(json);
                            resp->setStatusCode(drogon::k400BadRequest);
                            callback(resp);
                        } else {
                            std::vector<epi::SNP_t> snps_t;
                            snps_t.reserve(snps_s.size());
                            Json::Value snps_j;
                            for (auto &s: snps_s) {
                                snps_t.push_back(dm->snpStorage->by_name(s));
                                snps_j.append(s);
                            }
                            epi::SNPSet set(snps_t);
                            Json::Value json_scores;
                            for (auto &model: selected_models) {
                                json_scores[epi::options::epistasis_score_to_string(model)] = dm->snpStorage->calculate_score(set, model,
                                                                                  drogon::app().getCurrentThreadIndex());
                            }
                            Json::Value json;
                            json["ok"] = true;
                            json["scores"] = json_scores;
                            json["SNPs"] = snps_j;

                            if (unknown_scores.size() > 0) json["unknown_or_unavailable_scores"] = unknown_scores;

                            std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
                            auto timeMs = std::chrono::duration_cast<std::chrono::microseconds>(
                                    endTime - startTime).count();
                            json["time_microseconds"] = timeMs;
                            callback(drogon::HttpResponse::newHttpJsonResponse(json));
                        }
                    } catch (epi::SNPNotFoundError &err) {
                        LOG_WARN << "Unknown SNP " << err.get_name();
                        Json::Value json;
                        json["ok"] = false;
                        json["reason"] = "unknown SNP " + err.get_name();
                        json["unknown_snp"] = err.get_name();
                        if (unknown_scores.size() > 0) json["unknown_or_unavailable_scores"] = unknown_scores;
                        auto resp = drogon::HttpResponse::newHttpJsonResponse(json);
                        resp->setStatusCode(drogon::k404NotFound);
                        callback(resp);
                    } catch (epi::Error &err) {
                        LOG_WARN << "epi::Error  " << err.what();
                        Json::Value json;
                        json["ok"] = false;
                        if (unknown_scores.size() > 0) json["unknown_or_unavailable_scores"] = unknown_scores;
                        auto resp = drogon::HttpResponse::newHttpJsonResponse(json);
                        resp->setStatusCode(drogon::k500InternalServerError);
                        callback(resp);
                    }
                }
            },
            {drogon::Get});


    LOG_INFO << "Server running on " << bind_addr << ':' << bind_port << " with " << num_threads << " threads";
    drogon::app().addListener(bind_addr, bind_port);
    drogon::app().run();
    return 0;
}

