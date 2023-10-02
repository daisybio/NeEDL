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
 * @file  regression_model.cpp
 * @brief Unit tests for epi::RegressionModel class.
 * @details Compile via the option <tt>\--target regression_model</tt> or <tt>\--target unit</tt>.
 */

#define HEADER_ONLY
#include "../../../src/model/regression_model.hpp"
#include <catch.hpp>

#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

#include "../../../src/jobs/InstanceLoader.hpp"
#include "../../../src/jobs/DummyErrorTask.hpp"


TEST_CASE("Quantitative Instance") {
    epi::Instance<epi::QuantitativePhenoType> instance;
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/quantitative/2_disease_snps/exponential/0_1_ASW.json"));
    epi::RegressionModel<epi::QuantitativePhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{4, 30, 43};
    std::vector<epi::SNP> disease_snp_set_shuffled{30, 4, 43};
    std::vector<std::string> scores{"LLH", "NLL", "AIC", "BIC", "LLH-GAIN", "NLL-GAIN", "AIC-GAIN", "BIC-GAIN"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--max-itrs 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));

        // from max
        // REQUIRE(model.evaluate(disease_snp_set) > 0);

        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<double> predictions;
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_NOTHROW(predictions.emplace_back(model.predict(ind)));
        }
        CHECK_NOTHROW(model.save("../init/quantitative_model.ini"));
        CHECK_THAT(model.monte_carlo_p_value(disease_snp_set, 100, 0),
                   Catch::WithinAbs(model.monte_carlo_p_value(disease_snp_set_shuffled, 100, 0), 0.0001));
        CHECK_NOTHROW(model.load("../init/quantitative_model.ini"));
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_THAT(predictions.at(ind), Catch::WithinAbs(model.predict(ind), 0.0001));
        }
    }
}

TEST_CASE("Categorical Instance 1") {
    epi::Instance<epi::CategoricalPhenoType> instance(2);
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json"));
    epi::RegressionModel<epi::CategoricalPhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{69, 33};
    std::vector<epi::SNP> disease_snp_set_shuffled{33, 69};

    // std::vector<std::string> scores{"LLH", "NLL", "AIC", "BIC", "LLH-GAIN", "NLL-GAIN", "AIC-GAIN", "BIC-GAIN"};
    std::vector<std::string> scores{"NLL-GAIN", "AIC-GAIN", "BIC-GAIN"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--max-itrs 5 --score " + score)));
        double check1 = model.evaluate(disease_snp_set);
        double check2 = model.evaluate(disease_snp_set_shuffled);
        CHECK_THAT(check1, Catch::WithinAbs(check2, 0.0001));

        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<double> predictions;
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_NOTHROW(predictions.emplace_back(model.predict(ind)));
        }
        CHECK_NOTHROW(model.save("../init/categorical_model.ini"));
        CHECK_THAT(model.monte_carlo_p_value(disease_snp_set, 100, 0), Catch::WithinAbs(model.monte_carlo_p_value(disease_snp_set_shuffled, 100, 0), 0.0001));
        CHECK_NOTHROW(model.load("../init/categorical_model.ini"));
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK(predictions.at(ind) == model.predict(ind));
        }
    }
}

TEST_CASE("Categorical Instance 2") {
    epi::Instance<epi::CategoricalPhenoType> instance(2);
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST, "../../../data/MACOED/NME/00.1600.0.antesnp100.csv"));
    epi::RegressionModel<epi::CategoricalPhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{0, 1};
    std::vector<epi::SNP> disease_snp_set_shuffled{1, 0};
    std::vector<std::string> scores{"LLH", "NLL", "AIC", "BIC", "LLH-GAIN", "NLL-GAIN", "AIC-GAIN", "BIC-GAIN"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--max-itrs 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));

        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<double> predictions;
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_NOTHROW(predictions.emplace_back(model.predict(ind)));
        }
        CHECK_NOTHROW(model.save("../init/categorical_model.ini"));
        CHECK_THAT(model.monte_carlo_p_value(disease_snp_set, 100, 0), Catch::WithinAbs(model.monte_carlo_p_value(disease_snp_set_shuffled, 100, 0), 0.0001));
        CHECK_NOTHROW(model.load("../init/categorical_model.ini"));
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK(predictions.at(ind) == model.predict(ind));
        }
    }
}

TEST_CASE("Load Categorical Instance 1 - JSOPN_EPIGEN and CSV_COV") {
    epi::Instance<epi::CategoricalPhenoType> instance(2);
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN,
                                  "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::RegressionModel<epi::CategoricalPhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
}

TEST_CASE("Quantitative Instance - Cov Scores") {
    epi::Instance<epi::QuantitativePhenoType> instance;
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/quantitative/2_disease_snps/exponential/0_1_ASW.json"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::RegressionModel<epi::QuantitativePhenoType> model(&instance);
    model.cov_activate();
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{4, 30, 43};
    std::vector<epi::SNP> disease_snp_set_shuffled{30, 4, 43};
    std::vector<std::string> scores{"CLG-L-LC", "CLG-Q-QC", "CLG-Q-LC"};

    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--max-itrs 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));


        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<double> predictions;
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_NOTHROW(predictions.emplace_back(model.predict(ind)));
        }

        Eigen::MatrixXd model1 = model.get_interaction_model_();
        std::string model1_str = model.get_score();

        CHECK_NOTHROW(model.save("../init/quantitative_model.ini"));
        CHECK_THAT(model.monte_carlo_p_value(disease_snp_set, 100, 0),
                   Catch::WithinAbs(model.monte_carlo_p_value(disease_snp_set_shuffled, 100, 0), 0.0001));
        CHECK_NOTHROW(model.load("../init/quantitative_model.ini"));

        CHECK_THAT(model1_str, Catch::Matchers::Equals(model.get_score()));

        Eigen::MatrixXd model2 = model.get_interaction_model_();
        std::vector<double> expected_vec(model1.data(), model1.data() + model1.size());
        std::vector<double> actual_vec(model2.data(), model2.data() + model2.size());
        CHECK_THAT(expected_vec, Catch::Approx(actual_vec).epsilon(0.1));

        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_THAT(predictions.at(ind), Catch::WithinAbs(model.predict(ind), 0.0001));
        }
    }
}

TEST_CASE("Categorical Instance 1 - Cov Scores") {
    epi::Instance<epi::CategoricalPhenoType> instance(2);
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::RegressionModel<epi::CategoricalPhenoType> model(&instance);
    model.cov_activate();
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{69, 33};
    std::vector<epi::SNP> disease_snp_set_shuffled{33, 69};

    std::vector<std::string> scores{"CLG-Q-QC", "CLG-L-LC", "CLG-Q-LC"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--max-itrs 5 --score " + score)));
        double check1 = model.evaluate(disease_snp_set);
        double check2 = model.evaluate(disease_snp_set_shuffled);
        CHECK_THAT(check1, Catch::WithinAbs(check2, 0.0001));

        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<double> predictions;
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_NOTHROW(predictions.emplace_back(model.predict(ind)));
        }
        Eigen::MatrixXd model1 = model.get_interaction_model_();
        std::string model1_str = model.get_score();
        CHECK_NOTHROW(model.save("../init/categorical_model.ini"));
        CHECK_THAT(model.monte_carlo_p_value(disease_snp_set, 100, 0), Catch::WithinAbs(model.monte_carlo_p_value(disease_snp_set_shuffled, 100, 0), 0.0001));
        CHECK_NOTHROW(model.load("../init/categorical_model.ini"));
        CHECK_THAT(model1_str, Catch::Matchers::Equals(model.get_score()));

        Eigen::MatrixXd model2 = model.get_interaction_model_();
        std::vector<double> expected_vec(model1.data(), model1.data() + model1.size());
        std::vector<double> actual_vec(model2.data(), model2.data() + model2.size());
        CHECK_THAT(expected_vec, Catch::Approx(actual_vec).epsilon(0.1));
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK(predictions.at(ind) == model.predict(ind));
        }
    }
}

TEST_CASE("Categorical Instance 2 - Cov Scores") {
    epi::Instance<epi::CategoricalPhenoType> instance(2);
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST, "../../../data/MACOED/NME/00.1600.0.antesnp100.csv"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::RegressionModel<epi::CategoricalPhenoType> model(&instance);
    model.cov_activate();
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{0, 1};
    std::vector<epi::SNP> disease_snp_set_shuffled{1, 0};
    std::vector<std::string> scores{"CLG-Q-QC", "CLG-L-LC", "CLG-Q-LC"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--max-itrs 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));

        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<double> predictions;
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_NOTHROW(predictions.emplace_back(model.predict(ind)));
        }
        CHECK_NOTHROW(model.save("../init/categorical_model.ini"));
        CHECK_THAT(model.monte_carlo_p_value(disease_snp_set, 100, 0), Catch::WithinAbs(model.monte_carlo_p_value(disease_snp_set_shuffled, 100, 0), 0.0001
        ));
        CHECK_NOTHROW(model.load("../init/categorical_model.ini"));
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK(predictions.at(ind) == model.predict(ind));
        }
    }
}

TEST_CASE("Quantitative Instance - Cov activate/deactivate") {
    epi::Instance<epi::QuantitativePhenoType> instance;
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/quantitative/2_disease_snps/exponential/0_1_ASW.json"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::RegressionModel<epi::QuantitativePhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{4, 30, 43};
    std::vector<epi::SNP> disease_snp_set_shuffled{30, 4, 43};
    std::vector<std::string> scores{"CLG-Q-QC", "CLG-L-LC", "CLG-Q-LC", "LLH", "NLL", "AIC", "BIC", "LLH-GAIN", "NLL-GAIN", "AIC-GAIN", "BIC-GAIN"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--max-itrs 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));


        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<double> predictions;
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_NOTHROW(predictions.emplace_back(model.predict(ind)));
        }

        Eigen::MatrixXd model1 = model.get_interaction_model_();
        std::string model1_str = model.get_score();

        CHECK_NOTHROW(model.save("../init/quantitative_model.ini"));
        CHECK_THAT(model.monte_carlo_p_value(disease_snp_set, 100, 0),
                   Catch::WithinAbs(model.monte_carlo_p_value(disease_snp_set_shuffled, 100, 0), 0.0001));
        CHECK_NOTHROW(model.load("../init/quantitative_model.ini"));

        CHECK_THAT(model1_str, Catch::Matchers::Equals(model.get_score()));

        Eigen::MatrixXd model2 = model.get_interaction_model_();
        std::vector<double> expected_vec(model1.data(), model1.data() + model1.size());
        std::vector<double> actual_vec(model2.data(), model2.data() + model2.size());
        CHECK_THAT(expected_vec, Catch::Approx(actual_vec).epsilon(0.1));

        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_THAT(predictions.at(ind), Catch::WithinAbs(model.predict(ind), 0.0001));
        }
    }
}

TEST_CASE("Categorical Instance 1 - Cov activate/deactivate") {
    epi::Instance<epi::CategoricalPhenoType> instance(2);
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::RegressionModel<epi::CategoricalPhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{69, 33};
    std::vector<epi::SNP> disease_snp_set_shuffled{33, 69};

    std::vector<std::string> scores{"CLG-Q-QC", "CLG-L-LC", "CLG-Q-LC", "LLH", "NLL", "AIC", "BIC", "LLH-GAIN", "NLL-GAIN", "AIC-GAIN", "BIC-GAIN"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--max-itrs 5 --score " + score)));
        double check1 = model.evaluate(disease_snp_set);
        double check2 = model.evaluate(disease_snp_set_shuffled);
        CHECK_THAT(check1, Catch::WithinAbs(check2, 0.0001));

        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<double> predictions;
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_NOTHROW(predictions.emplace_back(model.predict(ind)));
        }
        Eigen::MatrixXd model1 = model.get_interaction_model_();
        std::string model1_str = model.get_score();
        CHECK_NOTHROW(model.save("../init/categorical_model.ini"));
        CHECK_THAT(model.monte_carlo_p_value(disease_snp_set, 100, 0), Catch::WithinAbs(model.monte_carlo_p_value(disease_snp_set_shuffled, 100, 0), 0.0001));
        CHECK_NOTHROW(model.load("../init/categorical_model.ini"));
        CHECK_THAT(model1_str, Catch::Matchers::Equals(model.get_score()));

        Eigen::MatrixXd model2 = model.get_interaction_model_();
        std::vector<double> expected_vec(model1.data(), model1.data() + model1.size());
        std::vector<double> actual_vec(model2.data(), model2.data() + model2.size());
        CHECK_THAT(expected_vec, Catch::Approx(actual_vec).epsilon(0.1));
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK(predictions.at(ind) == model.predict(ind));
        }
    }
}

TEST_CASE("Categorical Instance 1 - No NaN values") {
    epi::Instance<epi::CategoricalPhenoType> instance(2);
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::RegressionModel<epi::CategoricalPhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);

    std::vector<std::string> problematic_rs_ids = { "rs11589207", "rs10927631" };
    std::vector<epi::SNP> problematic_snps(2);
    // ugly way to find RS-IDs in dataset
    for (size_t i = 0; i < problematic_rs_ids.size(); ++i) {
        bool found = false;
        for (size_t j = 0; j < instance.rs_ids_.size(); ++j) {
            if (problematic_rs_ids[i] == instance.rs_ids_[j]) {
                problematic_snps[i] = j;
                found = true;
            }
        }
        CHECK(found);
    }

    std::vector<std::string> scores{"CLG-Q-QC", "CLG-L-LC", "CLG-Q-LC", "LLH", "NLL", "AIC", "BIC" };
    for(size_t i = 0; i < 1; ++i) {
        for (const auto &score: scores) {
            REQUIRE_NOTHROW(model.set_options(std::string("--score " + score)));
            REQUIRE_NOTHROW(model.initialize());
            REQUIRE_NOTHROW(model.cov_activate());

            double value = model.evaluate(problematic_snps);
            CHECK(!std::isnan(value));
            std::cout << value << std::endl;
        }
    }
}

TEST_CASE("Categorical Instance 1 - NeEDL with CSVParser") {
    // load the instance via InstanceLoader Job
    omp_set_num_threads(1);
    auto data = std::make_shared<epi::DataModel>(true);

    epi::InstanceLoader(
                                           "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json",
                                           "JSON_EPIGEN",
                                           "DICHOTOMOUS",
                                           2,
                                           "../../../data/COV_TEST/EpiGEN_RND_COV.csv"
                                   ).run(data);

    // this is the error-causing code (src/jobs/DummyErrorTask)
    // in the DummyErrorTask class I reduced the problematic code essentially to a CSVParser call after which the scores are nan
    epi::DummyErrorTask(
            "BIOGRID",
            "../../../data/BIOGRID/BIOGRID-ORGANISM-Homo_sapiens-3.5.182.tab2.txt",
            "Official Symbol Interactor A",
            "Official Symbol Interactor B",
            '\t',
            -1,
            -1).run(data);


    // the error does not occur when adding the CSVParser code directly to the test here


    // construct SNP set
    epi::SNPSet set {{data->snpStorage->by_name("rs11589207"),data->snpStorage->by_name("rs10927631")}};
    data->snpSetStorage.push_back(set);

    // test scores
    // auto all_scores = epi::options::get_all_epistasis_scores(true);
    std::vector<std::string> all_scores = { "REGRESSION_COV_NLL" };
    for (const auto & score : all_scores) {
        auto score_rep = epi::options::epistasis_score_from_string(score);
        double score_value = set.calculate_score(score_rep);

        CHECK(!std::isnan(score_value));
        std::cout << score_value << std::endl;
    }

    data->snpStorage.reset();
}

// this test is just here to show that the above issue really comes from the call to epi::DummyErrorTask
TEST_CASE("Categorical Instance 1 - NeEDL without CSVParser") {
    // load the instance via InstanceLoader Job
    omp_set_num_threads(1);
    auto data = std::make_shared<epi::DataModel>(true);

    epi::InstanceLoader(
            "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json",
            "JSON_EPIGEN",
            "DICHOTOMOUS",
            2,
            "../../../data/COV_TEST/EpiGEN_RND_COV.csv"
    ).run(data);

    // construct SNP set
    epi::SNPSet set {{data->snpStorage->by_name("rs11589207"),data->snpStorage->by_name("rs10927631")}};
    data->snpSetStorage.push_back(set);

    // test scores
    // auto all_scores = epi::options::get_all_epistasis_scores(true);
    std::vector<std::string> all_scores = { "REGRESSION_COV_NLL" };
    for (const auto & score : all_scores) {
        auto score_rep = epi::options::epistasis_score_from_string(score);
        double score_value = set.calculate_score(score_rep);

        CHECK(!std::isnan(score_value));
        std::cout << score_value << std::endl;
    }

    data->snpStorage.reset();
}
