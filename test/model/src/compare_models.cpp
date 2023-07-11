/*******************************************************************************
 *                                                                             *
 *   Copyright (C) 2020 by David B. Blumenthal                                 *
 *                                                                             *
 *   This file is part of NeEDL.                                        *
 *                                                                             *
 *   NeEDL is free software: you can redistribute it and/or modify it   *
 *   under the terms of the GNU General Public License as published by         *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   NeEDL is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with NeEDL. If not, see <http://www.gnu.org/licenses/>.      *
 *                                                                             *
 ******************************************************************************/


/*!
 * @file  compare_models.cpp
 * @brief Command line tool to compare different epistasis models.
 * @details Compile via the option <tt>\--target compare_models</tt>.
 * For more information, execute the compiled binary with the option <tt>\--help</tt>.
 */

#define HEADER_ONLY
#include "../../../src/model/all_models.hpp"
#include "../../../src/util/progress_bar.hpp"
#include <CLI11.hpp>
#include "util.hpp"

int parse_options(int argc, char* argv[], epi::ModelTestCLIOptions & options) {

    // Maps for mapping strings to enums.
    std::map<std::string, epi::options::InputFormat> format_map{{"JSON_EPIGEN", epi::options::InputFormat::JSON_EPIGEN}, {"CSV_SNPS_AS_ROWS_FIRST", epi::options::InputFormat::CSV_SNPS_AS_ROWS_FIRST},
                                                                {"CSV_SNPS_AS_ROWS_LAST", epi::options::InputFormat::CSV_SNPS_AS_ROWS_LAST}, {"CSV_SNPS_AS_COLUMNS_FIRST", epi::options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST},
                                                                {"CSV_SNPS_AS_COLUMNS_LAST", epi::options::InputFormat::CSV_SNPS_AS_COLUMNS_LAST}};
    std::map<std::string, epi::options::PhenoType> pheno_type_map{{"QUANTITATIVE", epi::options::PhenoType::QUANTITATIVE}, {"CATEGORICAL", epi::options::PhenoType::CATEGORICAL}, {"DICHOTOMOUS", epi::options::PhenoType::CATEGORICAL}};
    std::map<std::string, epi::options::EpistasisModel> model_map{{"VARIANCE_MODE", epi::options::EpistasisModel::VARIANCE_MODEL}, {"BAYESIAN_MODEL", epi::options::EpistasisModel::BAYESIAN_MODEL},
                                                                  {"PENETRANCE_MODEL", epi::options::EpistasisModel::PENETRANCE_MODEL}, {"REGRESSION_MODEL", epi::options::EpistasisModel::REGRESSION_MODEL}};

    // Function to check that argument is an integer greater than 1.
    std::function<std::string(std::string &)> check_integer_greater_1 = [] (std::string & arg) -> std::string {
        int arg_as_int;
        bool failed{false};
        try {
            arg_as_int = std::stoi(arg);
        }
        catch (...) {
            failed = true;
        }
        if (not failed) {
            failed = (arg_as_int < 2);
        }
        if (failed) {
            return "arguments must be integers greater than 1";
        }
        return "";
    };

    // Required options.
    CLI::App app("Run this application to compare different epistasis models.");
    app.add_option("--input-directory", options.input_directory, "The directory of the input files.")->required()->check(CLI::ExistingDirectory);
    app.add_option("--input-format", options.input_format, "The format of the input file.")->required()->transform(CLI::CheckedTransformer(format_map, CLI::ignore_case));
    app.add_option("--pheno-type", options.pheno_type, "The type of the phenotypes.")->required()->transform(CLI::CheckedTransformer(pheno_type_map, CLI::ignore_case));

    // Optional options.
    options.num_categories = 2;
    app.add_option("--num-categories", options.num_categories, "The number of categories for categorical phenotypes. DEFAULT: 2.")->check(CLI::Validator(check_integer_greater_1, "GREATER1"));
    options.num_random_solutions = 20;
    app.add_option("--num-random-solutions", options.num_random_solutions, "The number of random solutions used in evaluation. DEFAULT: 20.")->check(CLI::Validator(check_integer_greater_1, "GREATER1"));
    app.add_option("--size-range-random-solutions", options.size_range_random_solutions, "The range for the size of the random solutions. DEFAULT: number of disease SNPs +/- 1.")->expected(2)->check(CLI::Validator(check_integer_greater_1, "GREATER1"));
    app.add_option("--exclude", options.excluded_models, "Models that should be excluded from the comparison. DEFAULT: none.")->transform(CLI::CheckedTransformer(model_map, CLI::ignore_case));
    app.add_option("--input-files", options.input_files, "Input files containing the epistasis instances on which the comparison should be run. DEFAULT: all files contained in the input directory.");
    app.add_option("--disease-snps", options.disease_snps, "The disease SNPs. Must be provided unless specified in input files. DEFAULT: disease SNPs provided in input files.")->check(CLI::NonNegativeNumber);
    options.output_directory = "../res/";
    app.add_option("--output-directory", options.output_directory, "The output directory. DEFAULT: ../res/.")->check(CLI::ExistingDirectory);
    std::random_device rng;
    options.seed = static_cast<std::size_t>(rng());
    app.add_option("--seed", options.seed, "The seed to generate the random solutions. DEFAULT: drawn from the system.")->check(CLI::NonNegativeNumber);
    options.num_threads = 1;
    app.add_option("--num-threads", options.num_threads, "The number of threads. DEFAULT: 1.")->check(CLI::PositiveNumber);

    // Parse the options.
    CLI11_PARSE(app, argc, argv);
    if (options.input_directory.back() != '/') {
        options.input_directory += "/";
    }
    if (options.output_directory.back() != '/') {
        options.output_directory += "/";
    }

    // Get all input files if none are specified.
    if (options.input_files.empty()) {
        std::string suffix(".csv");
        if (options.input_format == epi::options::InputFormat::JSON_EPIGEN) {
            suffix = ".json";
        }
        boost::filesystem::path dir(options.input_directory);
        std::string file_name;
        for (const boost::filesystem::directory_entry & dir_entry : boost::filesystem::directory_iterator(dir)) {
            if (boost::filesystem::is_regular_file(dir_entry.path())) {
                file_name = dir_entry.path().filename().string();
                if (file_name.size() > suffix.size() and file_name.substr(file_name.size() - suffix.size()) == suffix) {
                    options.input_files.emplace_back(file_name);
                }
            }
        }
    }
    return 0;
}

void specify_model_variants(const epi::ModelTestCLIOptions & options, std::vector<epi::ModelTestComparedModel> & models_with_options) {

    if (options.num_threads == 1) {
        std::cout << "-- Specifying the model variants" << std::flush;
    }
    epi::ModelTestComparedModel model_with_options;
    model_with_options.options = "";
    model_with_options.predict = false;

    // Add BayesianModel.
    if (std::find(options.excluded_models.begin(), options.excluded_models.end(), epi::options::EpistasisModel::BAYESIAN_MODEL) == options.excluded_models.end()) {
        model_with_options.model = epi::options::EpistasisModel::BAYESIAN_MODEL;
        models_with_options.emplace_back(model_with_options);
    }

    // Add VarianceModel.
    if (std::find(options.excluded_models.begin(), options.excluded_models.end(), epi::options::EpistasisModel::VARIANCE_MODEL) == options.excluded_models.end()) {
        model_with_options.model = epi::options::EpistasisModel::VARIANCE_MODEL;
        models_with_options.emplace_back(model_with_options);
    }

    // Add variants of PenetranceModel.
    if (std::find(options.excluded_models.begin(), options.excluded_models.end(), epi::options::EpistasisModel::PENETRANCE_MODEL) == options.excluded_models.end()) {
        model_with_options.model = epi::options::EpistasisModel::PENETRANCE_MODEL;
        std::vector<std::string> penetrance_scores{"NLL", "AIC"};
        for (const auto & score : penetrance_scores) {
            model_with_options.options = "--score " + score;
            if (score == "NLL") {
                model_with_options.predict = true;
            }
            else {
                model_with_options.predict = false;
            }
            models_with_options.emplace_back(model_with_options);
        }
    }

    // Add variants of RegressionModel.
    if (std::find(options.excluded_models.begin(), options.excluded_models.end(), epi::options::EpistasisModel::REGRESSION_MODEL) == options.excluded_models.end()) {
        std::vector<std::string> regression_scores{"NLL", "AIC", "NLL-GAIN", "AIC-GAIN"};
        model_with_options.model = epi::options::EpistasisModel::REGRESSION_MODEL;
        for (const auto & score : regression_scores) {
            model_with_options.options = "--score " + score;
            if (score == "NLL") {
                model_with_options.predict = true;
            }
            else {
                model_with_options.predict = false;
            }
            models_with_options.emplace_back(model_with_options);
        }
    }
    if (options.num_threads == 1) {
        std::cout << "\n";
    }
}

template <class PhenoType>
void load_instance_and_construct_models(const std::string & input_file, epi::ModelTestCLIOptions & options, const std::vector<epi::ModelTestComparedModel> & models_with_options,
                                        epi::Instance<PhenoType> & instance, std::vector<epi::EpistasisModel<PhenoType> *> & epistasis_models) {

    // Load the instance.
    if (options.num_threads == 1) {
        std::cout << "-- Loading the instance and constructing the models." << std::flush;
    }
    instance.load(options.input_format, options.input_directory + input_file);
    if (options.disease_snps.empty()) {
        options.disease_snps = instance.disease_snps();
    }
    if (options.size_range_random_solutions.empty()) {
        options.size_range_random_solutions = {options.disease_snps.size() - 1, options.disease_snps.size() + 1};
    }

    // Construct and initialize an epistasis model for all selected model variants.
    for (std::size_t pos{0}; pos < models_with_options.size(); pos++) {
        switch (models_with_options.at(pos).model) {
            case epi::options::EpistasisModel::BAYESIAN_MODEL:
                epistasis_models.emplace_back(new epi::BayesianModel<PhenoType>(&instance));
                break;
            case epi::options::EpistasisModel::PENETRANCE_MODEL:
                epistasis_models.emplace_back(new epi::PenetranceModel<PhenoType>(&instance));
                break;
            case epi::options::EpistasisModel::REGRESSION_MODEL:
                epistasis_models.emplace_back(new epi::RegressionModel<PhenoType>(&instance));
                break;
            case epi::options::EpistasisModel::VARIANCE_MODEL:
                epistasis_models.emplace_back(new epi::VarianceModel<PhenoType>(&instance));
                break;
            default:
                throw epi::Error("Invalid epistasis model.");
                break;
        }
        epistasis_models.back()->set_options(models_with_options.at(pos).options);
        epistasis_models.back()->initialize();
    }
    if (options.num_threads == 1) {
        std::cout << "\n";
    }
}

template <class PhenoType>
void compute_scores(const epi::ModelTestCLIOptions & options, std::vector<epi::ModelTestComparedModel> & models_with_options,
                    epi::Instance<PhenoType> & instance, std::vector<epi::EpistasisModel<PhenoType> *> & epistasis_models) {

    // Compute scores for disease SNPs.
    epi::ProgressBar progress_bar(models_with_options.size());
    if (options.num_threads == 1) {
        std::cout << "\r-- Computing scores for disease SNPS: " << progress_bar << std::flush;
    }
    for (std::size_t pos{0}; pos < models_with_options.size(); pos++) {
        models_with_options.at(pos).score_disease_snps = epistasis_models.at(pos)->evaluate_track_time(options.disease_snps);
        models_with_options.at(pos).runtimes.emplace_back(epistasis_models.at(pos)->evaluation_time());
        if (options.num_threads == 1) {
            progress_bar.increment();
            std::cout << "\r-- Computing scores for disease SNPS: " << progress_bar << std::flush;
        }
    }
    if (options.num_threads == 1) {
        std::cout << "\n";
    }

    // Compute scores for random solutions.
    std::vector<epi::SNP> all_snps;
    for (epi::SNP snp{0}; snp < instance.num_snps(); snp++) {
        all_snps.emplace_back(snp);
    }
    std::vector<epi::SNP> random_solution;
    std::mt19937 urng;
    urng.seed(options.seed);
    std::uniform_int_distribution<std::size_t> uni(options.size_range_random_solutions.at(0), options.size_range_random_solutions.at(1));
    if (options.num_threads == 1) {
        progress_bar.reset(options.num_random_solutions * models_with_options.size());
        std::cout << "\r-- Computing scores for random solutions: " << progress_bar << std::flush;
    }
    for (std::size_t sol_id{0}; sol_id < options.num_random_solutions; sol_id++) {
        std::shuffle(all_snps.begin(), all_snps.end(), urng);
        random_solution.clear();
        for (std::size_t snp_id{0}; snp_id < uni(urng); snp_id++) {
            random_solution.push_back(all_snps.at(snp_id));
        }
        for (std::size_t pos{0}; pos < models_with_options.size(); pos++) {
            models_with_options.at(pos).scores_random_solutions.emplace_back(epistasis_models.at(pos)->evaluate_track_time(random_solution));
            models_with_options.at(pos).runtimes.emplace_back(epistasis_models.at(pos)->evaluation_time());
            if (options.num_threads == 1) {
                progress_bar.increment();
                std::cout << "\r-- Computing scores for random solutions: " << progress_bar << std::flush;
            }
        }
    }
    if (options.num_threads == 1) {
        std::cout << "\n";
    }
}

double predictions_to_score(const epi::ModelTestCLIOptions & options, const epi::Instance<epi::QuantitativePhenoType> & instance, std::vector<epi::QuantitativePhenoType> predictions) {
    double mean_absolute_error{0.0};
    for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
        mean_absolute_error += std::fabs(instance.phenotype(ind) - predictions.at(ind));
    }
    return mean_absolute_error / static_cast<double>(instance.num_inds());
}

double predictions_to_score(const epi::ModelTestCLIOptions & options, const epi::Instance<epi::CategoricalPhenoType> & instance, std::vector<epi::CategoricalPhenoType> predictions) {
    std::vector<std::size_t> true_positives(options.num_categories, 0.0);
    std::vector<std::size_t> class_counts(options.num_categories, 0);
    for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
        if (instance.phenotype(ind) == predictions.at(ind)) {
            true_positives.at(instance.phenotype(ind))++;
        }
        class_counts.at(instance.phenotype(ind))++;
    }
    double balanced_accuracy{0.0};
    std::size_t num_empty_categories{0};
    for (std::size_t class_id{0}; class_id < options.num_categories; class_id++) {
        if (class_counts.at(class_id) > 0) {
            balanced_accuracy += static_cast<double>(true_positives.at(class_id)) / static_cast<double>(class_counts.at(class_id));
        }
        else {
            num_empty_categories++;
        }
    }
    return balanced_accuracy / static_cast<double>(options.num_categories - num_empty_categories);
}

template <class PhenoType>
void test_predictions(const std::string & input_file, const epi::ModelTestCLIOptions & options, std::vector<epi::ModelTestComparedModel> & models_with_options,
                      epi::Instance<PhenoType> & instance, std::vector<epi::EpistasisModel<PhenoType> *> & epistasis_models) {

    for (auto & model_with_options : models_with_options) {
        model_with_options.mean_prediction_time = 0.0;
        model_with_options.prediction_score = 0.0;
    }
    std::size_t num_predict_models{0};
    for (const auto & model_with_options : models_with_options) {
        if (model_with_options.predict) {
            num_predict_models++;
        }
    }
    epi::ProgressBar progress_bar(10 * num_predict_models);
    if (options.num_threads == 1) {
        std::cout << "\r-- Testing predictions: " << progress_bar << std::flush;
    }
    for (std::size_t fold_id{0}; fold_id < 5; fold_id++) {
        instance.load(options.input_format, options.input_directory + input_file, 5, fold_id, epi::options::DataPurpose::TRAINING);
        for (std::size_t pos{0}; pos < models_with_options.size(); pos++) {
            if (not models_with_options.at(pos).predict) {
                continue;
            }
            epistasis_models.at(pos)->evaluate(options.disease_snps);
            progress_bar.increment();
            if (options.num_threads == 1) {
                std::cout << "\r-- Testing predictions: " << progress_bar << std::flush;
            }
        }
        instance.load(options.input_format,  options.input_directory + input_file, 5, fold_id, epi::options::DataPurpose::VALIDATION);
        for (std::size_t pos{0}; pos < models_with_options.size(); pos++) {
            std::vector<PhenoType> predictions;
            double sum_prediction_time{0.0};
            if (not models_with_options.at(pos).predict) {
                continue;
            }
            for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
                predictions.emplace_back(epistasis_models.at(pos)->predict_track_time(ind));
                sum_prediction_time += epistasis_models.at(pos)->prediction_time();
            }
            models_with_options.at(pos).mean_prediction_time += sum_prediction_time / static_cast<double>(instance.num_inds());
            models_with_options.at(pos).prediction_score += predictions_to_score(options, instance, predictions);
            if (options.num_threads == 1) {
                progress_bar.increment();
                if (options.num_threads == 1) {
                    std::cout << "\r-- Testing predictions: " << progress_bar << std::flush;
                }
            }
        }
    }
    for (auto & model_with_options : models_with_options) {
        model_with_options.mean_prediction_time /= 5.0;
        model_with_options.prediction_score /= 5.0;
    }
    if (options.num_threads == 1) {
        std::cout << "\n";
    }

}

template <class PhenoType>
void aggregate_results(const epi::ModelTestCLIOptions & options, const std::vector<epi::EpistasisModel<PhenoType> *> & epistasis_models, std::vector<epi::ModelTestComparedModel> & models_with_options, Eigen::MatrixXd & correlation_matrix) {

    // Compute p-values and mean runtimes for all tested models, and populate correlation matrix.
    epi::ProgressBar progress_bar(models_with_options.size());
    if (options.num_threads == 1) {
        std::cout << "\r-- Aggregating the results: " << progress_bar << std::flush;
    }
    epi::options::TestDirection direction;
    Eigen::MatrixXd random_scores_matrix(models_with_options.size(), models_with_options.at(0).scores_random_solutions.size());
    for (std::size_t pos{0}; pos < models_with_options.size(); pos++) {

        // Compute mean runtime.
        models_with_options.at(pos).mean_runtime = 0;
        for (auto runtime : models_with_options.at(pos).runtimes) {
            models_with_options.at(pos).mean_runtime += runtime;
        }
        models_with_options.at(pos).mean_runtime /= static_cast<double>(models_with_options.at(pos).runtimes.size());

        // Populate random scores matrix.
        for (std::size_t score_id{0}; score_id < models_with_options.at(pos).scores_random_solutions.size(); score_id++) {
            random_scores_matrix(pos, score_id) = models_with_options.at(pos).scores_random_solutions.at(score_id);
        }

        // Compute p-value of one-sample t-test.
        if (epistasis_models.at(pos)->model_sense() == epi::options::ModelSense::MINIMIZE) {
            direction = epi::options::TestDirection::UPPER_TAILED;
        }
        else {
            direction = epi::options::TestDirection::LOWER_TAILED;
        }
        try {
            models_with_options.at(pos).p_value_from_score = epi::misc::one_sample_t_test(models_with_options.at(pos).scores_random_solutions, models_with_options.at(pos).score_disease_snps, direction);
        }
        catch (...) {
            std::stringstream model_name;
            model_name << models_with_options.at(pos).model;
            throw epi::Error("Error in one_sample_t_test for " + model_name.str() + "[" + models_with_options.at(pos).options + "]");
        }
        if (options.num_threads == 1) {
            progress_bar.increment();
            std::cout << "\r-- Aggregating the results: " << progress_bar << std::flush;
        }
    }
    if (options.num_threads == 1) {
        std::cout << "\n";
    }

    // Compute correlations.
    Eigen::MatrixXd product_matrix = random_scores_matrix * random_scores_matrix.transpose();
    for (std::size_t pos_1{0}; pos_1 < models_with_options.size(); pos_1++) {
        for (std::size_t pos_2{0}; pos_2 < models_with_options.size(); pos_2++) {
            correlation_matrix(pos_1, pos_2) = product_matrix(pos_1, pos_2) / (std::sqrt(product_matrix(pos_1, pos_1)) * std::sqrt(product_matrix(pos_2, pos_2)));
        }
    }

    // Compute mean correlations.
    for (std::size_t pos_1{0}; pos_1 < models_with_options.size(); pos_1++) {
        models_with_options.at(pos_1).mean_correlation = 0;
        for (std::size_t pos_2{0}; pos_2 < models_with_options.size(); pos_2++) {
            if (pos_2 != pos_1) {
                models_with_options.at(pos_1).mean_correlation += correlation_matrix(pos_1, pos_2);
            }
        }
        models_with_options.at(pos_1).mean_correlation /= static_cast<double>(models_with_options.size() - 1);
    }
}

void write_results(const std::string & input_file, const epi::ModelTestCLIOptions & options, const std::vector<epi::ModelTestComparedModel> & models_with_options, const Eigen::MatrixXd & correlation_matrix) {

    // Contruct the prefix of the result files.
    std::string prefix(options.output_directory);
    std::string suffix(".csv");
    if (options.input_format == epi::options::InputFormat::JSON_EPIGEN) {
        suffix = ".json";
    }
    prefix += input_file.substr(0, input_file.size() - suffix.size());


    // Write aggregated results.
    std::ofstream results((prefix + "_AGG.csv").c_str());
    results << "model,evaluation-time,neg-log-p-value,mean-correlation,prediction-time,";
    if (options.pheno_type == epi::options::PhenoType::CATEGORICAL) {
        results << "balanced-accuracy\n";
    }
    else {
        results << "mean-absolute-error\n";
    }
    for (const auto & model : models_with_options) {
        results << model;
    }
    results.close();

    // Write correlation matrix.
    results.open((prefix + "_CORR.csv").c_str(), std::ios::out);
    for (const auto & model_with_options : models_with_options) {
        results << "," << model_with_options.model << "[" << model_with_options.options << "]";
    }
    results << "\n";
    for (std::size_t pos_1{0}; pos_1 < models_with_options.size(); pos_1++) {
        results << models_with_options.at(pos_1).model << "[" << models_with_options.at(pos_1).options << "]";
        for (std::size_t pos_2{0}; pos_2 < models_with_options.size(); pos_2++) {
            results << "," << correlation_matrix(pos_1, pos_2);
        }
        results << "\n";
    }
    results.close();

    // Print information about output files.
    if (options.num_threads == 1) {
        std::cout << "-- Result files: ";
        std::cout << prefix + "_AGG.csv" << ", ";
        std::cout << prefix + "_CORR.csv" << "\n";
    }
}

template <class PhenoType>
void compare_models(const std::string & input_file, epi::ModelTestCLIOptions & options, epi::Instance<PhenoType> & instance) {

    // Specify all model variants that are included in the comparison.
    std::vector<epi::ModelTestComparedModel> models_with_options;
    specify_model_variants(options, models_with_options);

    // Load epistasis instance and construct models.
    std::vector<epi::EpistasisModel<PhenoType> *> epistasis_models;
    load_instance_and_construct_models(input_file, options, models_with_options, instance, epistasis_models);

    // Compare the selected models.
    compute_scores(options, models_with_options, instance, epistasis_models);
    test_predictions(input_file, options, models_with_options, instance, epistasis_models);

    // Aggregate the results.
    Eigen::MatrixXd correlation_matrix(models_with_options.size(), models_with_options.size());
    aggregate_results(options, epistasis_models, models_with_options, correlation_matrix);

    // Free memory.
    for (auto model : epistasis_models) {
        delete model;
    }

    // Write results.
    write_results(input_file, options, models_with_options, correlation_matrix);
}

int main(int argc, char* argv[]) {

    std::cout << "\n";
    std::cout << "**************************************************\n";
    std::cout << "                        NeEDL                   \n";
    std::cout << "             Comparing Epistasis Models           \n";
    std::cout << "**************************************************\n";

    // Parse the input options.
    epi::ModelTestCLIOptions options;
    parse_options(argc, argv, options);


#ifdef _OPENMP
    omp_set_num_threads(options.num_threads);
#pragma omp parallel for if (options.num_threads > 1) schedule(dynamic)
#endif
    for (std::size_t pos = 0; pos < options.input_files.size(); pos++) {

#pragma omp critical
        {
            std::cout << "-- Processing instance " << pos + 1 << " of " << options.input_files.size() << ".\n";
        }

        epi::ModelTestCLIOptions private_options(options);

        // Run the comparison.
        std::string input_file(private_options.input_files.at(pos));
        if (options.pheno_type == epi::options::PhenoType::CATEGORICAL) {
            epi::Instance<epi::CategoricalPhenoType> instance(options.num_categories);
            compare_models(input_file, private_options, instance);
        }
        else {
            epi::Instance<epi::QuantitativePhenoType> instance;
            compare_models(input_file, private_options, instance);
        }
    }

    // Return 0 if the program terminated as expected.
    std::cout << "**************************************************\n";
    return 0;
}
