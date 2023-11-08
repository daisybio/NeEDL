//
// Created by juli on 22.06.22.
//

#include "WriteSets.hpp"
#include "../util/TerminalTable.hpp"
#include "../util/TimeLogger.hpp"
#include <set>
#include <utility>
#include "../util/helper_functions.hpp"

namespace epi {

    WriteSets::WriteSets(std::string rank_model, std::vector<std::string> scores, std::string name, bool individual_snps) {
        this->name = std::move(name);
        this->individual_snps = individual_snps;

        if (!rank_model.empty()) {
            scores.insert(scores.begin(), rank_model);
            has_rank_model = true;
            this->rank_model = options::epistasis_score_from_string(rank_model);
        } else {
            has_rank_model = false;
        }

        std::set<options::EpistasisScore> models;
        std::transform(scores.begin(), scores.end(), std::inserter(models, models.end()), [] (auto & model_str) -> auto {
           return options::epistasis_score_from_string(model_str);
        });
        /*
        if (has_rank_model) {
            models.erase(this->rank_model);
        }
         */

        epi_models = {models.begin(), models.end()};
    }

    void WriteSets::run(std::shared_ptr<DataModel> data) {
        TimeLogger logger(std::string("writing ") + (individual_snps ? "individual SNPs" : "sets") + " (+ calculating additional scores)");

        std::vector<SNPSet> all_sets;
        if (individual_snps) {
            std::unordered_set<SNP_t, SNP_t::SNPHash> all_snps;
            for (auto & set : data->snpSetStorage) {
                all_snps.insert(set.begin(), set.end());
            }
            all_sets.reserve(all_snps.size());
            std::transform(all_snps.begin(), all_snps.end(), std::back_inserter(all_sets), [] (auto & snp) -> auto {
                return SNPSet({ snp });
            });
        } else {
            all_sets = { data->snpSetStorage.begin(), data->snpSetStorage.end() };
        }

        // find out all available attributes
        std::set<std::string> attrib_names;
        for (auto & set : all_sets) {
            auto keys = set.get_attribute_keys();
            attrib_names.insert(keys.begin(), keys.end());
        }

        // create heading
        std::vector<std::string> heading_names;
        if (has_rank_model) heading_names.emplace_back("RANK (" + options::epistasis_score_to_string(rank_model) + ")");
        heading_names.emplace_back("RS_IDS");

        for (auto & model : epi_models) heading_names.push_back(options::epistasis_score_to_string(model));
        heading_names.insert(heading_names.end(), attrib_names.begin(), attrib_names.end());

        heading_names.push_back("ANNOTATIONS");

        // add absolute and percentage column for every category --> how many individuals match this configuration in each category?
        std::vector<size_t> num_individuals_per_phenotype;
        if (data->snpStorage->is_categorical()) {
            num_individuals_per_phenotype = data->snpStorage->get_num_individuals_per_category();
            for (size_t i = 0; i < data->snpStorage->num_categories(); i++) {
                heading_names.push_back("NUM_INDIVIDUALS_" + std::to_string(i));
                heading_names.push_back("FREQ_INDIVIDUALS_" + std::to_string(i));
                heading_names.push_back("INDIVIDUALS_" + std::to_string(i));
            }
        }

        // calculate all scores multi-threaded
#pragma omp parallel for default(none) shared(all_sets)
        for (size_t i = 0; i < all_sets.size(); i++) {
            for(auto & model : epi_models) {
                all_sets[i].calculate_score(model);
            }
        }

        // rank if requested
        if (has_rank_model) {
            bool need_minimizing = data->snpStorage->need_minimizing_score(rank_model);
            std::sort(all_sets.begin(), all_sets.end(), [&, need_minimizing] (auto a, auto b) -> bool {
                double score_a = a.calculate_score(rank_model);
                double score_b = b.calculate_score(rank_model);
                return (need_minimizing && score_a < score_b) || (!need_minimizing && score_a > score_b);
            });
        }


        if (Logger::consoleLoggingActive()) {
            // print table to console
            TerminalTable table(150, 10);
            table.add_row(heading_names);

            // go over all sets
            size_t pos = 1;
            for (auto & set : all_sets) {
                std::vector<std::string> row_vals;

                if (has_rank_model) {
                    row_vals.push_back(std::to_string(pos));
                }

                row_vals.push_back(set.get_snp_string());

                for(auto & model : epi_models) {
                    row_vals.push_back(Logger::to_string(set.calculate_score(model)));
                }

                for (auto & key : attrib_names) {
                    row_vals.push_back(set.get_attribute(key));
                }

                auto annotations = set.get_annotations();
                row_vals.push_back(epi::string_join(annotations, ";"));


                if (data->snpStorage->is_categorical()) {
                    auto inds = data->snpStorage->get_individuals_per_category(set);
                    for (size_t i = 0; i < inds.size(); i++) {
                        row_vals.push_back(std::to_string(inds[i].size()));
                        double ratio = double(inds[i].size()) / double(num_individuals_per_phenotype[i]);
                        row_vals.push_back(Logger::to_string(ratio));
                        row_vals.push_back(epi::string_join(inds[i], ";", std::to_string));

                    }
                }

                table.add_row(row_vals);
                pos ++;
            }

            table.print();
        }

        bool has_outfile = false;
        std::string path;
        if (!ext_outfile_path.empty()) {
            path = ext_outfile_path;
            has_outfile = true;
        } else if (data->outputDirectory != nullptr) {
            // write to file
            path = data->outputDirectory->get_free_filepath(name, ".csv");
            has_outfile = true;
        }

        if (has_outfile) {
            std::ofstream file (path);

            // write header
            bool is_first = true;
            for (auto & heading_name : heading_names) {
                if (is_first) is_first = false;
                else file << '\t';

                file << heading_name;
            }
            file << '\n';

            // write csv body
            size_t pos = 1;
            for (auto & set : all_sets) {
                std::vector<std::string> row_vals;

                if (has_rank_model) {
                    row_vals.push_back(std::to_string(pos));
                }

                row_vals.push_back(set.get_snp_string());

                for(auto & model : epi_models) {
                    row_vals.push_back(Logger::to_string(set.calculate_score(model)));
                }

                for (auto & key : attrib_names) {
                    row_vals.push_back(set.get_attribute(key));
                }


                auto annotations = set.get_annotations();
                row_vals.push_back(epi::string_join(annotations, ";"));

                if (data->snpStorage->is_categorical()) {
                    auto inds = data->snpStorage->get_individuals_per_category(set);
                    for (size_t i = 0; i < inds.size(); i++) {
                        row_vals.push_back(std::to_string(inds[i].size()));
                        double ratio = double(inds[i].size()) / double(num_individuals_per_phenotype[i]);
                        row_vals.push_back(Logger::to_string(ratio));
                        row_vals.push_back(epi::string_join(inds[i], ";", std::to_string));

                    }
                }

                is_first = true;
                for (auto & cell : row_vals) {
                    if (is_first) is_first = false;
                    else file << '\t';

                    file << cell;
                }
                file << '\n';

                pos ++;
            }
        }

        logger.stop();
    }

    void WriteSets::outfile_path(std::string path) {
        ext_outfile_path = path;
    }

    rapidjson::Value WriteSets::getConfig(rapidjson::Document &doc) {
        rapidjson::Value obj(rapidjson::kObjectType);
        obj.AddMember("JOB", rapidjson::Value().SetString("WriteSets"), doc.GetAllocator());
        obj.AddMember("name", rapidjson::Value().SetString(name.c_str(), name.size(), doc.GetAllocator()), doc.GetAllocator());
        if(!ext_outfile_path.empty()) {
            obj.AddMember("ext_outfile_path",
                          rapidjson::Value().SetString(ext_outfile_path.c_str(), ext_outfile_path.size(),
                                                       doc.GetAllocator()), doc.GetAllocator());
        }

        obj.AddMember("individual_snps", rapidjson::Value().SetBool(individual_snps), doc.GetAllocator());
        obj.AddMember("has_rank_model", rapidjson::Value().SetBool(has_rank_model), doc.GetAllocator());
        if (has_rank_model) {
            auto rm_str = options::epistasis_score_to_string(rank_model);
            obj.AddMember("rank_model", rapidjson::Value().SetString(rm_str.c_str(), rm_str.size(), doc.GetAllocator()), doc.GetAllocator());
        }

        rapidjson::Value model_arr(rapidjson::kArrayType);
        for (auto & model : epi_models) {
            auto m_str = options::epistasis_score_to_string(model);
            model_arr.PushBack(rapidjson::Value().SetString(m_str.c_str(), m_str.size(), doc.GetAllocator()), doc.GetAllocator());
        }

        return obj;
    }


} // epi