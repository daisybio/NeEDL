/*******************************************************************************
 *                                                                             *
 *   Copyright (C) 2020 by David B. Blumenthal                                 *
 *                                                                             *
 *   This file is part of GenEpiSeeker.                                        *
 *                                                                             *
 *   GenEpiSeeker is free software: you can redistribute it and/or modify it   *
 *   under the terms of the GNU General Public License as published by         *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   GenEpiSeeker is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with GenEpiSeeker. If not, see <http://www.gnu.org/licenses/>.      *
 *                                                                             *
 ******************************************************************************/


/*!
 * @file  penetrance_model.cpp
 * @brief Unit tests for epi::PenetranceModel class.
 * @details Compile via the option <tt>\--target penetrance_model</tt> or <tt>\--target unit</tt>.
 */

#define HEADER_ONLY
#include "../../../src/model/penetrance_model.hpp"
#include <catch.hpp>

TEST_CASE("Quantitative Instance") {
    epi::Instance<epi::QuantitativePhenoType> instance;
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/quantitative/2_disease_snps/exponential/0_1_ASW.json"));
    epi::PenetranceModel<epi::QuantitativePhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{4, 30, 43};
    std::vector<epi::SNP> disease_snp_set_shuffled{30, 4, 43};
    std::vector<std::string> scores{"LLH", "NLL", "AIC", "BIC"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--min-cell-size 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));
        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<double> predictions;
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_NOTHROW(predictions.emplace_back(model.predict(ind)));
        }
        CHECK_NOTHROW(model.save("../init/quantitative_model.ini"));
        CHECK_THAT(model.monte_carlo_p_value(disease_snp_set, 100, 0), Catch::WithinAbs(model.monte_carlo_p_value(disease_snp_set_shuffled, 100, 0), 0.0001));
        CHECK_NOTHROW(model.load("../init/quantitative_model.ini"));
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_THAT(predictions.at(ind), Catch::WithinAbs(model.predict(ind), 0.0001));
        }
    }
}

TEST_CASE("Categorical Instance 1") {
    epi::Instance<epi::CategoricalPhenoType> instance(2);
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json"));
    epi::PenetranceModel<epi::CategoricalPhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{55, 63};
    std::vector<epi::SNP> disease_snp_set_shuffled{63, 55};
    std::vector<std::string> scores{"LLH", "NLL", "AIC", "BIC"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--min-cell-size 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));
        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<std::size_t> predictions;
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
    epi::PenetranceModel<epi::CategoricalPhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{0, 1};
    std::vector<epi::SNP> disease_snp_set_shuffled{1, 0};
    std::vector<std::string> scores{"LLH", "NLL", "AIC", "BIC"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--min-cell-size 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));
        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<std::size_t> predictions;
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

TEST_CASE("LOAD COV") {
    epi::Instance<epi::QuantitativePhenoType> instance;
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/quantitative/2_disease_snps/exponential/0_1_ASW.json"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::PenetranceModel<epi::QuantitativePhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
}

TEST_CASE("Quantitative Instance - Cov") {
    epi::Instance<epi::QuantitativePhenoType> instance;
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/quantitative/2_disease_snps/exponential/0_1_ASW.json"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::PenetranceModel<epi::QuantitativePhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    REQUIRE_NOTHROW(model.cov_activate());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{4, 30, 43};
    std::vector<epi::SNP> disease_snp_set_shuffled{30, 4, 43};
    std::vector<std::string> scores{"LLH", "NLL", "AIC", "BIC"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--min-cell-size 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));
        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<double> predictions;
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_NOTHROW(predictions.emplace_back(model.predict(ind)));
        }
        CHECK_NOTHROW(model.save("../init/quantitative_model.ini"));
        CHECK_THAT(model.monte_carlo_p_value(disease_snp_set, 100, 0), Catch::WithinAbs(model.monte_carlo_p_value(disease_snp_set_shuffled, 100, 0), 0.0001));
        CHECK_NOTHROW(model.load("../init/quantitative_model.ini"));
        for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
            CHECK_THAT(predictions.at(ind), Catch::WithinAbs(model.predict(ind), 0.0001));
        }
    }
}

TEST_CASE("Categorical Instance 1 - Cov") {
    epi::Instance<epi::CategoricalPhenoType> instance(2);
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::PenetranceModel<epi::CategoricalPhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    REQUIRE_NOTHROW(model.cov_activate());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{55, 63};
    std::vector<epi::SNP> disease_snp_set_shuffled{63, 55};
    std::vector<std::string> scores{"LLH", "NLL", "AIC", "BIC"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--min-cell-size 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));
        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<std::size_t> predictions;
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

TEST_CASE("Categorical Instance 2 - Cov") {
    epi::Instance<epi::CategoricalPhenoType> instance(2);
    REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST, "../../../data/MACOED/NME/00.1600.0.antesnp100.csv"));
    REQUIRE_NOTHROW(instance.load_cov(epi::options::InputFormat::CSV_COV, "../../../data/COV_TEST/EpiGEN_RND_COV.csv"));
    epi::PenetranceModel<epi::CategoricalPhenoType> model(&instance);
    REQUIRE_NOTHROW(model.initialize());
    REQUIRE_NOTHROW(model.cov_activate());
    CHECK(model.is_predictive());
    CHECK(model.model_sense() == epi::options::ModelSense::MINIMIZE);
    std::vector<epi::SNP> disease_snp_set{0, 1};
    std::vector<epi::SNP> disease_snp_set_shuffled{1, 0};
    std::vector<std::string> scores{"LLH", "NLL", "AIC", "BIC"};
    for (const auto & score : scores) {
        REQUIRE_NOTHROW(model.set_options(std::string("--min-cell-size 5 --score " + score)));
        CHECK_THAT(model.evaluate(disease_snp_set), Catch::WithinAbs(model.evaluate(disease_snp_set_shuffled), 0.0001));
        CHECK_NOTHROW(model.evaluate_track_time(disease_snp_set));
        std::vector<std::size_t> predictions;
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


