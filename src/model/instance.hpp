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
 * @file  instance.hpp
 * @brief Declaration of epi::Instance.
 */

#ifndef SRC_MODEL_INSTANCE_HPP_
#define SRC_MODEL_INSTANCE_HPP_

#include "../util/types.hpp"
#include "../util/csv_parser.hpp"
#include "../util/misc.hpp"

namespace epi {

/*!
 * @brief Contains the epistasis instance.
 * @tparam PhenoType The type of the phenotypes. Use epi::CategoricalPhenoType for categorical phenotypes and epi::QuantitativePhenoType for quantitative phenotypes.
 * @details All objectives and algorithms used for epistasis detection operate on top of this class.
 */
    template<class PhenoType>
    class Instance {

    public:

        /*!
         * @brief Constant iterator for iterating over the genotypes of all individuals at a given SNP.
         */
        typedef std::vector<GenoType>::const_iterator const_ind_iterator;

        /*!
     * @brief Constructs and empty instance.
     * @param[in] num_categories The number of categories for categorical phenotypes. Has no effect if @p PhenoType is set to epi::QuantititativePhenoType.
     */
        Instance(std::size_t num_categories = 2);

        /*!
         * @brief Constant iterator for iterating over the genotype of an individual at all SNPs.
         */
        class const_snp_iterator {

        public:

            typedef std::vector<GenoType>::const_iterator::iterator_category iterator_category;

            typedef std::vector<GenoType>::const_iterator::difference_type difference_type;

            typedef std::vector<GenoType>::const_iterator::value_type value_type;

            typedef std::vector<GenoType>::const_iterator::pointer pointer;

            typedef std::vector<GenoType>::const_iterator::reference reference;

            const_snp_iterator(std::vector<GenoType>::const_iterator itr, std::size_t offset, std::size_t pos);

            const_snp_iterator operator++();

            const_snp_iterator operator++(int);

            const value_type & operator*();

            bool operator==(const const_snp_iterator & rhs);

            bool operator!=(const const_snp_iterator & rhs);

        private:

            std::vector<GenoType>::const_iterator itr_;

            std::size_t offset_;

            std::size_t pos_;
        };




        /*!
         * @brief Loads the instance.
         * @param[in] input_format Format of the input file.
         * @param[in] filename Path to the input file.
         * - If @p input_format is set to epi::options::InputFormat::CSV_SNPS_AS_ROWS_FIRST, the input must be a comma-separated file whose rows
         *   represent the SNPs and whose columns represent the individuals. The first column must contain the reference SNP IDs (e.g. rs206437),
         *   and the last row must contain the phenotypes of the individuals. All other cells must contain the number of minor alleles of the
         *   individual corresponding to the current column at the SNP corresponding to the current row (either 0, 1, or 2).
         * - If @p input_format is set to epi::options::InputFormat::CSV_SNPS_AS_ROWS_LAST, the input must be a comma-separated file whose rows
         *   represent the SNPs and whose columns represent the individuals. The last column must contain the reference SNP IDs (e.g. rs206437),
         *   and the last row must contain the phenotypes of the individuals. All other cells must contain the number of minor alleles of the
         *   individual corresponding to the current column at the SNP corresponding to the current row (either 0, 1, or 2).
         * - If @p input_format is set to epi::options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST, the input must be a comma-separated file whose columns
         *   represent the SNPs and whose rows represent the individuals. The first row must contain the reference SNP IDs (e.g. rs206437),
         *   and the last column must contain the phenotypes of the individuals. All other cells must contain the number of minor alleles of the
         *   individual corresponding to the current row at the SNP corresponding to the current column (either 0, 1, or 2).
         * - If @p input_format is set to epi::options::InputFormat::CSV_SNPS_AS_COLUMNS_LAST, the input must be a comma-separated file whose columns
         *   represent the SNPs and whose rows represent the individuals. The last row must contain the reference SNP IDs (e.g. rs206437),
         *   and the last column must contain the phenotypes of the individuals. All other cells must contain the number of minor alleles of the
         *   individual corresponding to the current row at the SNP corresponding to the current column (either 0, 1, or 2).
         * - If @p input_format is set to epi::options::InputFormat::JSON_EPIGEN, the input file must be in the JSON format generated by EpiGEN.
         *   For more details, consult the documentation of EpiGEN: https://github.com/baumbachlab/epigen.
         * @param[in] num_folds If set to a value greater than 1, this parameter specifies the number of folds to be used for cross-validation.
         * @param[in] fold_id If @p num_folds is set to a value greater than 1, this parameter specifies which fold of the data should be loaded.
         * @param[in] cv_data Specifies whether the loaded data should be used for training of for validation.
         */
        void load(options::InputFormat input_format, const std::string & filename, std::size_t num_folds = 1, std::size_t fold_id = 0, options::DataPurpose cv_data = options::DataPurpose::TRAINING);

        void load_cov(options::InputFormat input_format, const std::string &filename);
        /*!
         * @brief Returns the number of SNPs in the instance.
         * @return The number of SNPs.
         */
        std::size_t num_snps() const;

        /*!
         * @brief Returns the number of individuals in the instance.
         * @return The number of individuals.
         */
        std::size_t num_inds() const;

        /*!
         * @brief Returns the number of categories for categorical phenotypes.
         * @return The number of categories for categorical phenotypes.
         * @note Since the number of categories is meaningful only for categorical phenotypes,
         * this member function is implemented only if the template parameter @p PhenoType equals epi::CategoricalPhenoType.
         */
        std::size_t num_categories() const;

        /*!
         * @brief Returns the genotype of an individual at a given SNP.
         * @param[in] snp The SNP for which the genotype should be returned.
         * @param[in] ind The individual whose genotype should be returned.
         * @return The genotype of individual @p ind at SNP @p snp.
         */
        GenoType genotype_at_snp(SNP snp, Ind ind) const;

        /*!
         * @brief Provides access to the genotype of an individual at a given SNP set.
         * @param[in] snp_set The SNP set for which the genotype should be returned.
         * @param[in] ind The individual whose genotype should be returned.
         * @param[out] genotype The genotype of individual @p ind at the SNP set @p snp_set.
         */
        void genotype_at_snp_set(const std::vector<SNP> & snp_set, Ind ind, std::vector<GenoType> & genotype) const;

        /*!
         * @brief Returns the genotype (as integer ID) of an individual at a given SNP set.
         * @param[in] snp_set The SNP set for which the genotype should be returned.
         * @param[in] ind The individual whose genotype should be returned.
         * @return The genotype (as integer ID) of individual @p ind at the SNP set @p snp_set.
         */
        std::size_t genotype_at_snp_set(const std::vector<SNP> & snp_set, Ind ind) const;

        /*!
         * @brief Collects all individuals with a given genotype at a given SNP set.
         * @param[in] snp_set The reference SNP set.
         * @param[in] genotype The genotype of the selected individuals at the SNP set @p snp_set.
         * @param[out] inds All individuals whose genotype at the SNP set @p snp_set equal @p genotype.
         * @note The parameters @p snp_set and @p genotype must be of the same size.
         */
        void inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, const std::vector<GenoType> & genotype, std::vector<Ind> & inds) const;

        /*!
         * @brief Collects all individuals with a given genotype at a given SNP set.
         * @param[in] snp_set The reference SNP set.
         * @param[in] genotype_id The genotype (as integer ID) of the selected individuals at the SNP set @p snp_set.
         * @param[out] inds All individuals whose genotype at the SNP set @p snp_set equal @p genotype.
         * @note The parameters @p snp_set and @p genotype must be of the same size.
         */
        void inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, std::size_t genotype_id, std::vector<Ind> & inds) const;

        /*!
         * @brief Counts the individuals with a given genotype at a given SNP set.
         * @param[in] snp_set The reference SNP set.
         * @param[in] genotype The target genotype at the SNP set @p snp_set.
         * @return The number of individuals whose genotype at the SNP set @p snp_set equals @p genotype.
         * @note The parameters @p snp_set and @p genotype must be of the same size.
         */
        std::size_t num_inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, const std::vector<GenoType> & genotype) const;

        /*!
         * @brief Counts the individuals with a given genotype at a given SNP set.
         * @param[in] snp_set The reference SNP set.
         * @param[in] genotype_id The target genotype (as integer ID) at the SNP set @p snp_set.
         * @return The number of individuals whose genotype at the SNP set @p snp_set equals @p genotype.
         * @note The parameters @p snp_set and @p genotype must be of the same size.
         */
        std::size_t num_inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, std::size_t genotype_id) const;

        /*!
         * @brief Returns the phenotype of an individual.
         * @param[in] ind The individuals whose phenotype should be returned.
         * @return The phenotype of individual @p ind.
         */
        PhenoType phenotype(Ind ind) const;


        /*!
         * @brief Provides read-only access to the genotypes of all individuals at a given SNP.
         * @param[in] snp The SNP for which the genotypes should be provided.
         * @return An constant iterator pointing to the beginning of the genotypes of interest.
         */
        const_ind_iterator genotypes_of_all_inds_begin(SNP snp) const;

        /*!
         * @brief Provides read-only access to the genotypes of all individuals at a given SNP.
         * @param[in] snp The SNP for which the genotypes should be provided.
         * @return An constant iterator pointing to the end of the genotypes of interest.
         */
        const_ind_iterator genotypes_of_all_inds_end(SNP snp) const;

        /*!
         * @brief Provides read-only access to the genotypes of a given individuals at all SNPs.
         * @param[in] ind The individual for which the genotypes should be provided.
         * @return An constant iterator pointing to the beginning of the genotypes of interest.
         */
        const_snp_iterator genotypes_at_all_snps_begin(Ind ind) const;

        /*!
         * @brief Provides read-only access to the genotypes of a given individuals at all SNPs.
         * @param[in] ind The individual for which the genotypes should be provided.
         * @return An constant iterator pointing to the end of the genotypes of interest.
         */
        const_snp_iterator genotypes_at_all_snps_end(Ind ind) const;

        /*!
         * @brief Shuffles the phenotypes.
         */
        void shuffle_phenotypes();

        /*!
         * @brief Restores the original phenotypes.
         */
        void restore_phenotypes();

        /*!
         * @brief Sets the seed used for shuffling the phenotypes.
         * @param[in] seed The selected seed.
         */
        void set_seed(std::size_t seed);

        /*!
         * @brief Checks whether the phenotypes are quantitative.
         * @return Boolean true if @p PhenoType equals epi::QuantitativePhenoType.
         */
        bool quantitative_phenotypes() const;

        /*!
         * @brief Checks whether the phenotypes are categorical.
         * @return Boolean true if @p PhenoType equals epi::CategoricalPhenoType.
         */
        bool categorical_phenotypes() const;

        /*!
         * @brief Provides access to disease related SNPs (if available).
         * @return Constant reference to vector containing the disease related SNP.
         */
        const std::vector<SNP> & disease_snps() const;

        /*!
         * @brief Sets the disease SNPs.
         * @param[in] disease_snps The vector of disease SNPs.
         */
        void set_disease_snps(const std::vector<SNP> & disease_snps);

        /*!
         * @brief Returns textual SNP descriptor.
         * @param[in] snp The SNP whose descriptor should be returned.
         * @return Returns the descriptor of SNP @p snp provided in the input file.
         */
        std::string snp_descriptor(SNP snp) const;

        /*!
         * @brief stores the instance in a compact binary format that can be loaded more quickly than a json file
         * @param filename path to the output file
         */
        void save_bin(const std::string & filename);

        std::vector<std::string> rs_ids_;
        std::vector<std::string> rs_ids_chromosomes_;
        std::vector<double> rs_ids_maf_;

        size_t num_covs() const;

        Eigen::MatrixXd get_covariates() const;

        Eigen::VectorXd get_covariates_at_ind(Ind ind) const;

    private:

        std::size_t num_categories_;

        std::size_t num_snps_;

        std::size_t num_inds_;

        std::vector<GenoType> genotypes_;

        std::vector<PhenoType> phenotypes_;

        std::vector<PhenoType> original_phenotypes_;

        std::vector<SNP> disease_snps_;

        std::mt19937 urng_;

        Eigen::MatrixXd covariates_;

        void load_csv_(const std::string & filename, bool snps_as_rows, bool snp_info_comes_first, std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data);

        void load_json_(const std::string & filename, std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data);

        void load_bin_(const std::string & filename, std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data);

        PhenoType parse_phenotype_(const std::string & field, Ind ind) const;

        std::size_t construct_folds_(std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data, std::vector<bool> & skip);

        void load_csv_cov_(const std::string &filename);
    };

}

#include "instance.ipp"

#ifdef HEADER_ONLY
#include "instance.cpp"
#endif


#endif /* SRC_MODEL_INSTANCE_HPP_ */
