
#ifndef GENEPISEEKER_SNPSTORAGE_IPP
#define GENEPISEEKER_SNPSTORAGE_IPP

#include "SNPStorage.hpp"

namespace epi {

    template<class PhenoType>
    SNPStorage_PhenoType<PhenoType>::SNPStorage_PhenoType(std::shared_ptr<epi::Instance<PhenoType>> instance) : SNPStorage() {
        // this->score_calculation_max_snp_set_sizes = std::vector<size_t>(omp_get_num_threads(), 0);
        this->instance = instance;

        // create snp_data
        snp_data = std::vector<SNPData>(instance->rs_ids_.size());

        // add names to lookup map
        for (size_t i = 0; i < instance->rs_ids_.size(); i++) {
            const auto &name = instance->rs_ids_[i];
            if (name_map.find(name) != name_map.end()) {
                throw epi::Error("SNP name " + name + " found multiple times.");
            }

            name_map.insert({ name, i });

            if (name.size() >= 2 && name.substr(0, 2) == "rs") {
                std::string name_without_rs_prefix = name.substr(2);
                if (name_map.find(name_without_rs_prefix) != name_map.end()) {
                    throw epi::Error("SNP name " + name_without_rs_prefix + " found multiple times (automatic trimming of 'rs' for RS-IDs).");
                }

                name_map.insert({ name_without_rs_prefix, i });
            }
        }

        init_model_containers(omp_get_max_threads());
    }


    template<class PhenoType>
    void SNPStorage_PhenoType<PhenoType>::init_model_containers(size_t num_threads) {
        bayesian_models = { num_threads, nullptr };
        penetrance_models = { num_threads, nullptr };
        regression_models = { num_threads, nullptr };
        variance_models = { num_threads, nullptr };
        penetrance_model_scores = { num_threads, options::EpistasisScore::PENETRANCE_NLL };
        regression_model_scores = { num_threads, options::EpistasisScore::REGRESSION_NLL };
    }

    template<class PhenoType>
    const std::string &epi::SNPStorage_PhenoType<PhenoType>::snp_get_name(const SNP_t &snp) const {
        return instance->rs_ids_[snp.value];
    }

    template<class PhenoType>
    const std::string &epi::SNPStorage_PhenoType<PhenoType>::snp_get_chromosome(const SNP_t &snp) const {
        if (instance->rs_ids_chromosomes_.empty()) throw epi::Error("Chromosome information needed but not provided.");
        return instance->rs_ids_chromosomes_[snp.value];
    }

    template<class PhenoType>
    double epi::SNPStorage_PhenoType<PhenoType>::snp_get_maf(const SNP_t &snp) const {
        if (instance->rs_ids_maf_.empty()) throw epi::Error("MAF information needed but not provided.");
        return instance->rs_ids_maf_[snp.value];
    }

    template<class PhenoType>
    void SNPStorage_PhenoType<PhenoType>::save(std::string path) {
        instance->save_bin(path);
    }

    template<class PhenoType>
    double SNPStorage_PhenoType<PhenoType>::calculate_score(const SNPSet & set, options::EpistasisScore score, size_t thread_index) {
        std::vector<SNP> set_instance {};
        for (auto &snp : set) set_instance.push_back(snp.value);

        /*
        if (score_calculation_max_snp_set_sizes[omp_get_thread_num()] < set.size()) {
            score_calculation_max_snp_set_sizes[omp_get_thread_num()] = set.size();
            Logger::logLine("Thread " + std::to_string(omp_get_thread_num()) + ": max set size now " + std::to_string(set.size()) + " SNPs.");
        }
        */

        request_score_method(score, thread_index);

        switch(epistasis_model_from_epistasis_score(score)) {
            case options::EpistasisModel::BAYESIAN_MODEL:
                return bayesian_models[thread_index]->evaluate(set_instance);
            case options::EpistasisModel::PENETRANCE_MODEL:
                return penetrance_models[thread_index]->evaluate(set_instance);
            case options::EpistasisModel::REGRESSION_MODEL:
                return regression_models[thread_index]->evaluate(set_instance);
            case options::EpistasisModel::VARIANCE_MODEL:
                return variance_models[thread_index]->evaluate(set_instance);
        }

        return 0;
    }

    template<class PhenoType>
    void SNPStorage_PhenoType<PhenoType>::request_score_method(options::EpistasisScore score, size_t current_thread) {
        switch (epistasis_model_from_epistasis_score(score)) {
            case options::EpistasisModel::BAYESIAN_MODEL:
                if (bayesian_models[current_thread] == nullptr)
                    bayesian_models[current_thread] = std::make_shared<BayesianModel<PhenoType>>(instance.get());
                break;
            case options::EpistasisModel::PENETRANCE_MODEL:
                if (penetrance_models[current_thread] == nullptr)
                    penetrance_models[current_thread] = std::make_shared<PenetranceModel<PhenoType>>(instance.get());
                break;
            case options::EpistasisModel::REGRESSION_MODEL:
                if (regression_models[current_thread] == nullptr)
                    regression_models[current_thread] = std::make_shared<RegressionModel<PhenoType>>(instance.get());
                break;
            case options::EpistasisModel::VARIANCE_MODEL:
                if (variance_models[current_thread] == nullptr)
                    variance_models[current_thread] = std::make_shared<VarianceModel<PhenoType>>(instance.get());
                break;
        }

        // set required options if necessary
        switch(score) {
            case options::EpistasisScore::PENETRANCE_NLL:
                if (penetrance_model_scores[current_thread] != score) {
                    penetrance_models[current_thread]->set_options("--score NLL");
                    penetrance_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::PENETRANCE_LLH:
                if (penetrance_model_scores[current_thread] != score) {
                    penetrance_models[current_thread]->set_options("--score LLH");
                    penetrance_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::PENETRANCE_AIC:
                if (penetrance_model_scores[current_thread] != score) {
                    penetrance_models[current_thread]->set_options("--score AIC");
                    penetrance_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::PENETRANCE_BIC:
                if (penetrance_model_scores[current_thread] != score) {
                    penetrance_models[current_thread]->set_options("--score BIC");
                    penetrance_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::REGRESSION_NLL:
                if (regression_model_scores[current_thread] != score) {
                    regression_models[current_thread]->set_options("--score NLL");
                    regression_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::REGRESSION_LLH:
                if (regression_model_scores[current_thread] != score) {
                    regression_models[current_thread]->set_options("--score LLH");
                    regression_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::REGRESSION_AIC:
                if (regression_model_scores[current_thread] != score) {
                    regression_models[current_thread]->set_options("--score AIC");
                    regression_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::REGRESSION_BIC:
                if (regression_model_scores[current_thread] != score) {
                    regression_models[current_thread]->set_options("--score BIC");
                    regression_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::REGRESSION_NLL_GAIN:
                if (regression_model_scores[current_thread] != score) {
                    regression_models[current_thread]->set_options("--score NLL-GAIN");
                    regression_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::REGRESSION_LLH_GAIN:
                if (regression_model_scores[current_thread] != score) {
                    regression_models[current_thread]->set_options("--score LLH-GAIN");
                    regression_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::REGRESSION_AIC_GAIN:
                if (regression_model_scores[current_thread] != score) {
                    regression_models[current_thread]->set_options("--score AIC-GAIN");
                    regression_model_scores[current_thread] = score;
                }
                break;
            case options::EpistasisScore::REGRESSION_BIC_GAIN:
                if (regression_model_scores[current_thread] != score) {
                    regression_models[current_thread]->set_options("--score BIC-GAIN");
                    regression_model_scores[current_thread] = score;
                }
                break;
            default:
                break;
        }
    }

    template<class PhenoType>
    bool SNPStorage_PhenoType<PhenoType>::has_maf_information() const {
        return !instance->rs_ids_maf_.empty();
    }

    template<class PhenoType>
    void SNPStorage_PhenoType<PhenoType>::set_maf_information(std::vector<double> maf_data) {
        if (maf_data.size() != instance->rs_ids_.size()) {
            throw epi::Error("Cannot insert " + std::to_string(maf_data.size()) + " MAF values into instance with " + std::to_string(instance->rs_ids_.size()) + " SNPs.");
        }
    }

    template<class PhenoType>
    size_t SNPStorage_PhenoType<PhenoType>::num_snps() {
        return instance->rs_ids_.size();
    }

    template<class PhenoType>
    bool SNPStorage_PhenoType<PhenoType>::is_categorical() const {
        return instance->categorical_phenotypes();
    }

    template<class PhenoType>
    bool SNPStorage_PhenoType<PhenoType>::need_minimizing_score(options::EpistasisScore score, size_t thread_index) {
        request_score_method(score, thread_index);
        switch(epistasis_model_from_epistasis_score(score)) {
            case options::EpistasisModel::BAYESIAN_MODEL:
                return bayesian_models[thread_index]->model_sense() == options::ModelSense::MINIMIZE;
            case options::EpistasisModel::PENETRANCE_MODEL:
                return penetrance_models[thread_index]->model_sense() == options::ModelSense::MINIMIZE;
            case options::EpistasisModel::REGRESSION_MODEL:
                return regression_models[thread_index]->model_sense() == options::ModelSense::MINIMIZE;
            case options::EpistasisModel::VARIANCE_MODEL:
                return variance_models[thread_index]->model_sense() == options::ModelSense::MINIMIZE;
        }

        return false;
    }

    template<class PhenoType>
    double SNPStorage_PhenoType<PhenoType>::calculate_monte_carlo_score(const SNPSet &set, options::EpistasisScore score,
                                                                        size_t num_permutations, size_t thread_index) {
        std::vector<SNP> set_instance {};
        for (auto &snp : set) set_instance.push_back(snp.value);

        request_score_method(score, thread_index);

        switch(epistasis_model_from_epistasis_score(score)) {
            case options::EpistasisModel::BAYESIAN_MODEL:
                return bayesian_models[thread_index]->monte_carlo_p_value(set_instance, num_permutations);
            case options::EpistasisModel::PENETRANCE_MODEL:
                return penetrance_models[thread_index]->monte_carlo_p_value(set_instance, num_permutations);
            case options::EpistasisModel::REGRESSION_MODEL:
                return regression_models[thread_index]->monte_carlo_p_value(set_instance, num_permutations);
            case options::EpistasisModel::VARIANCE_MODEL:
                return variance_models[thread_index]->monte_carlo_p_value(set_instance, num_permutations);
        }

        return 0;
    }

    template<class PhenoType>
    std::vector<std::vector<size_t>> SNPStorage_PhenoType<PhenoType>::get_individuals_per_category(const SNPSet & set) const {
        if (!is_categorical()) {
            throw epi::Error("Tried to access individuals by category, even though the phenotype is not categorical.");
        }

        return {};
    }



    template<class PhenoType>
    size_t SNPStorage_PhenoType<PhenoType>::num_categories() const {
        if (!is_categorical()) {
            throw epi::Error("Tried to access number of categories, even though the phenotype is not categorical.");
        }
        return {};
    }

    template<class PhenoType>
    std::vector<size_t> SNPStorage_PhenoType<PhenoType>::get_num_individuals_per_category() const {
        if (!is_categorical()) {
            throw epi::Error("Tried to access number of categories, even though the phenotype is not categorical.");
        }
        return {};
    }


}

#endif // GENEPISEEKER_SNPSTORAGE_IPP