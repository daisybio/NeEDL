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
 * @file  instance.cpp
 * @brief Unit tests for epi::Instance class.
 * @details Compile via the option <tt>\--target instance</tt> or <tt>\--target unit</tt>.
 */

#define HEADER_ONLY
#include "../../../src/model/instance.hpp"
#include <catch.hpp>

TEST_CASE("Quantitative Instance") {
	epi::Instance<epi::QuantitativePhenoType> instance;
	REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/quantitative/2_disease_snps/exponential/0_1_ASW.json"));
	std::size_t total_num_inds{instance.num_inds()};
	std::size_t min_fold_size{total_num_inds / 3};
	std::size_t sum_num_validation_inds{0};
	std::vector<std::size_t> accepted_fold_sizes{min_fold_size, min_fold_size + 1};
	for (std::size_t fold_id{0}; fold_id < 3; fold_id++) {
		REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/quantitative/2_disease_snps/exponential/0_1_ASW.json", 3, fold_id, epi::options::DataPurpose::TRAINING));
		std::size_t num_training_inds{instance.num_inds()};
		REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/quantitative/2_disease_snps/exponential/0_1_ASW.json", 3, fold_id, epi::options::DataPurpose::VALIDATION));
		std::size_t num_validation_inds{instance.num_inds()};
		sum_num_validation_inds += num_validation_inds;
		REQUIRE(num_training_inds + num_validation_inds == total_num_inds);
		REQUIRE_THAT(accepted_fold_sizes, Catch::VectorContains(num_validation_inds));
	}
	REQUIRE(sum_num_validation_inds == total_num_inds);
	for (epi::SNP snp{0}; snp < instance.num_snps(); snp++) {
		for (auto geno = instance.genotypes_of_all_inds_begin(snp); geno != instance.genotypes_of_all_inds_end(snp); geno++) {
			REQUIRE_NOTHROW(static_cast<unsigned>(*geno));
		}
	}
	for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
		for (auto geno = instance.genotypes_at_all_snps_begin(ind); geno != instance.genotypes_at_all_snps_end(ind); geno++) {
			REQUIRE_NOTHROW(static_cast<unsigned>(*geno));
		}
	}
	std::vector<epi::SNP> snp_set{2, 23, 74, 66, 88};
	for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {

		epi::QuantitativePhenoType pheno;
		REQUIRE_NOTHROW(pheno = instance.phenotype(ind));

		std::vector<epi::GenoType> genotype;
		instance.genotype_at_snp_set(snp_set, ind, genotype);
		std::size_t genotype_id{instance.genotype_at_snp_set(snp_set, ind)};

		std::vector<epi::GenoType> reconstructed_genotype;
		epi::misc::id_to_genotype(genotype_id, snp_set.size(), reconstructed_genotype);
		REQUIRE_THAT(genotype, Catch::Equals(reconstructed_genotype));

		std::vector<epi::Ind> inds_with_genotype_at_snp_set;
		instance.inds_with_genotype_at_snp_set(snp_set, genotype, inds_with_genotype_at_snp_set);
		REQUIRE_THAT(inds_with_genotype_at_snp_set, Catch::VectorContains(ind));

		instance.inds_with_genotype_at_snp_set(snp_set, genotype_id, inds_with_genotype_at_snp_set);
		REQUIRE_THAT(inds_with_genotype_at_snp_set, Catch::VectorContains(ind));
	}
}

TEST_CASE("Categorical Instance") {
	epi::Instance<epi::CategoricalPhenoType> instance(2);
	REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json"));
	std::size_t total_num_inds{instance.num_inds()};
	std::size_t min_fold_size{total_num_inds / 3};
	std::size_t sum_num_validation_inds{0};
	std::vector<std::size_t> accepted_fold_sizes{min_fold_size, min_fold_size + 1};
	for (std::size_t fold_id{0}; fold_id < 3; fold_id++) {
		REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json", 3, fold_id, epi::options::DataPurpose::TRAINING));
		std::size_t num_training_inds{instance.num_inds()};
		REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::JSON_EPIGEN, "../../../data/EpiGEN/dichotomous/2_disease_snps/exponential/1_1_ASW.json", 3, fold_id, epi::options::DataPurpose::VALIDATION));
		std::size_t num_validation_inds{instance.num_inds()};
		sum_num_validation_inds += num_validation_inds;
		REQUIRE(num_training_inds + num_validation_inds == total_num_inds);
		REQUIRE_THAT(accepted_fold_sizes, Catch::VectorContains(num_validation_inds));
	}
	REQUIRE(sum_num_validation_inds == total_num_inds);
	for (epi::SNP snp{0}; snp < instance.num_snps(); snp++) {
		for (auto geno = instance.genotypes_of_all_inds_begin(snp); geno != instance.genotypes_of_all_inds_end(snp); geno++) {
			REQUIRE_NOTHROW(static_cast<unsigned>(*geno));
		}
	}
	for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {
		for (auto geno = instance.genotypes_at_all_snps_begin(ind); geno != instance.genotypes_at_all_snps_end(ind); geno++) {
			REQUIRE_NOTHROW(static_cast<unsigned>(*geno));
		}
	}
	std::vector<epi::SNP> snp_set{2, 23, 74, 66, 88};
	for (epi::Ind ind{0}; ind < instance.num_inds(); ind++) {

		epi::QuantitativePhenoType pheno;
		REQUIRE_NOTHROW(pheno = instance.phenotype(ind));

		std::vector<epi::GenoType> genotype;
		instance.genotype_at_snp_set(snp_set, ind, genotype);
		std::size_t genotype_id{instance.genotype_at_snp_set(snp_set, ind)};

		std::vector<epi::GenoType> reconstructed_genotype;
		epi::misc::id_to_genotype(genotype_id, snp_set.size(), reconstructed_genotype);
		REQUIRE_THAT(genotype, Catch::Equals(reconstructed_genotype));

		std::vector<epi::Ind> inds_with_genotype_at_snp_set;
		instance.inds_with_genotype_at_snp_set(snp_set, genotype, inds_with_genotype_at_snp_set);
		REQUIRE_THAT(inds_with_genotype_at_snp_set, Catch::VectorContains(ind));

		instance.inds_with_genotype_at_snp_set(snp_set, genotype_id, inds_with_genotype_at_snp_set);
		REQUIRE_THAT(inds_with_genotype_at_snp_set, Catch::VectorContains(ind));
	}
}

TEST_CASE("Categorical CSV Instances") {
	epi::Instance<epi::CategoricalPhenoType> instance(2);
	REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST, "../../../data/MACOED/NME/00.1600.0.antesnp100.csv"));
	std::size_t total_num_inds{instance.num_inds()};
	std::size_t min_fold_size{total_num_inds / 3};
	std::size_t sum_num_validation_inds{0};
	std::vector<std::size_t> accepted_fold_sizes{min_fold_size, min_fold_size + 1};
	for (std::size_t fold_id{0}; fold_id < 3; fold_id++) {
		REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST, "../../../data/MACOED/NME/00.1600.0.antesnp100.csv", 3, fold_id, epi::options::DataPurpose::TRAINING));
		std::size_t num_training_inds{instance.num_inds()};
		REQUIRE_NOTHROW(instance.load(epi::options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST, "../../../data/MACOED/NME/00.1600.0.antesnp100.csv", 3, fold_id, epi::options::DataPurpose::VALIDATION));
		std::size_t num_validation_inds{instance.num_inds()};
		sum_num_validation_inds += num_validation_inds;
		REQUIRE(num_training_inds + num_validation_inds == total_num_inds);
		REQUIRE_THAT(accepted_fold_sizes, Catch::VectorContains(num_validation_inds));
	}
	REQUIRE(sum_num_validation_inds == total_num_inds);
}
