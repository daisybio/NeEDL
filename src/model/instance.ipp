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
 * @file  instance.ipp
 * @brief Definition of epi::Instance.
 */

//#ifndef SRC_MODEL_INSTANCE_IPP_
// #define SRC_MODEL_INSTANCE_IPP_

#include "../../ext/rapidjson/document.h"
#include "../../ext/rapidjson/istreamwrapper.h"
#include "../../ext/rapidjson/filereadstream.h"
#include "instance.hpp"

namespace epi {

    template<class PhenoType>
    Instance<PhenoType>::
    const_snp_iterator::
    const_snp_iterator(std::vector<GenoType>::const_iterator itr, std::size_t offset, std::size_t pos):
            itr_(itr),
            offset_{offset},
            pos_{pos} {}

    template<class PhenoType>
    typename Instance<PhenoType>::const_snp_iterator
    Instance<PhenoType>::
    const_snp_iterator::
    operator++() {
        itr_ += offset_;
        pos_ += offset_;
        return *this;
    }

    template<class PhenoType>
    typename  Instance<PhenoType>::const_snp_iterator
    Instance<PhenoType>::
    const_snp_iterator::
    operator++(int) {
        const_snp_iterator temp(*this);
        itr_ += offset_;
        pos_ += offset_;
        return temp;
    }

    template<class PhenoType>
    const typename Instance<PhenoType>::const_snp_iterator::value_type &
    Instance<PhenoType>::
    const_snp_iterator::
    operator*() {
        return *itr_;
    }

    template<class PhenoType>
    bool
    Instance<PhenoType>::
    const_snp_iterator::
    operator==(const const_snp_iterator & rhs) {
        return pos_ == rhs.pos_;
    }

    template<class PhenoType>
    bool
    Instance<PhenoType>::
    const_snp_iterator::
    operator!=(const const_snp_iterator & rhs) {
        return pos_ != rhs.pos_;
    }

    template<class PhenoType>
    Instance<PhenoType>::
    Instance(std::size_t num_categories):
            rs_ids_(),
            rs_ids_chromosomes_(),
            num_categories_{num_categories},
            num_snps_{undefined_uint()},
            num_inds_{undefined_uint()},
            genotypes_(),
            phenotypes_(),
            original_phenotypes_(),
            disease_snps_(),
            urng_(),
            covariates_(),
            has_cov_(false) {}

    template<class PhenoType>
    void
    Instance<PhenoType>::
    load(options::InputFormat input_format, const std::string & filename, std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data) {
        genotypes_.clear();
        phenotypes_.clear();
        original_phenotypes_.clear();
        rs_ids_.clear();
        rs_ids_chromosomes_.clear();
        switch (input_format) {
            case options::InputFormat::CSV_SNPS_AS_ROWS_FIRST:
                load_csv_(filename, true, true, num_folds, fold_id, cv_data);
                break;
            case options::InputFormat::CSV_SNPS_AS_ROWS_LAST:
                load_csv_(filename, true, false, num_folds, fold_id, cv_data);
                break;
            case options::InputFormat::CSV_SNPS_AS_COLUMNS_FIRST:
                load_csv_(filename, false, true, num_folds, fold_id, cv_data);
                break;
            case options::InputFormat::CSV_SNPS_AS_COLUMNS_LAST:
                load_csv_(filename, false, false, num_folds, fold_id, cv_data);
                break;
            case options::InputFormat::JSON_EPIGEN:
                load_json_(filename, num_folds, fold_id, cv_data);
                break;
            case options::InputFormat::NEEDL_BIN:
                load_bin_(filename, num_folds, fold_id, cv_data);
                break;
            default:
                throw Error("Unsupported input format.");
        }
        original_phenotypes_ = phenotypes_;
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    load_cov(options::InputFormat input_format, const std::string & filename) {
        covariates_ = Eigen::MatrixXd();
        switch (input_format) {
            case options::InputFormat::CSV_COV:
                load_csv_cov_(filename);
                break;
            default:
                throw Error("Unsupported input format.");
        }
        has_cov_ = true;
    }

    template<class PhenoType>
    std::size_t
    Instance<PhenoType>::
    num_snps() const {
        return num_snps_;
    }

    template<class PhenoType>
    std::size_t
    Instance<PhenoType>::
    num_inds() const {
        return num_inds_;
    }

    template<class PhenoType>
    bool
    Instance<PhenoType>::
    has_cov() {
        return has_cov_;
    }

    template<class PhenoType>
    std::size_t
    Instance<PhenoType>::
    num_covs() const {
        return covariates_.cols();
    }

    template<class PhenoType>
    GenoType
    Instance<PhenoType>::
    genotype_at_snp(SNP snp, Ind ind) const {
        return genotypes_.at(snp * num_inds_ + ind);
    }

    template<class PhenoType>
    Eigen::MatrixXd
    Instance<PhenoType>::
    get_covariates() const {
        return covariates_;
    }

    template<class PhenoType>
    Eigen::VectorXd
    Instance<PhenoType>::
    get_covariates_at_ind(Ind ind) const {
        return covariates_.row(ind);
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    genotype_at_snp_set(const std::vector<SNP> & snp_set, Ind ind, std::vector<GenoType> & genotype) const {
        genotype.clear();
        for (SNP snp : snp_set) {
            genotype.emplace_back(genotype_at_snp(snp, ind));
        }
    }

    template<class PhenoType>
    std::size_t
    Instance<PhenoType>::
    genotype_at_snp_set(const std::vector<SNP> & snp_set, Ind ind) const {
        std::size_t genotype_id{0};
        std::size_t exponent{snp_set.size() - 1};
        for (SNP snp : snp_set) {
            genotype_id += genotype_at_snp(snp, ind) * misc::uint_pow(3, exponent--);
        }
        return genotype_id;
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, const std::vector<GenoType> & genotype, std::vector<Ind> & inds) const {
        inds.clear();
        for (Ind ind{0}; ind < num_inds_; ind++) {
            bool match{true};
            for (std::size_t pos{0}; pos < snp_set.size(); pos++) {
                if (genotype_at_snp(snp_set.at(pos), ind) != genotype.at(pos)) {
                    match = false;
                    break;
                }
            }
            if (match) {
                inds.emplace_back(ind);
            }
        }
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, std::size_t genotype_id, std::vector<Ind> & inds) const {
        std::vector<GenoType> genotype;
        misc::id_to_genotype(genotype_id, snp_set.size(), genotype);
        return inds_with_genotype_at_snp_set(snp_set, genotype, inds);
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::inds_with_nonzero_genotype_at_snp_set(const std::vector<SNP> &snp_set, std::vector<Ind> & inds) const {
        inds.clear();

        for (Ind ind{0}; ind < num_inds_; ind++) {
            bool match{true};
            for (std::size_t pos{0}; pos < snp_set.size(); pos++) {
                if (genotype_at_snp(snp_set.at(pos), ind) == 0) {
                    match = false;
                    break;
                }
            }
            if (match) {
                inds.emplace_back(ind);
            }
        }
    }

    template<class PhenoType>
    std::size_t
    Instance<PhenoType>::
    num_inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, const std::vector<GenoType> & genotype) const {
        std::size_t num_matches{0};
        for (Ind ind{0}; ind < num_inds_; ind++) {
            bool match{true};
            for (std::size_t pos{0}; pos < snp_set.size(); pos++) {
                if (genotype_at_snp(snp_set.at(pos), ind) != genotype.at(pos)) {
                    match = false;
                    break;
                }
            }
            if (match) {
                num_matches++;
            }
        }
        return num_matches;
    }

    template<class PhenoType>
    std::size_t
    Instance<PhenoType>::
    num_inds_with_genotype_at_snp_set(const std::vector<SNP> & snp_set, std::size_t genotype_id) const {
        std::vector<GenoType> genotype;
        misc::id_to_genotype(genotype_id, snp_set.size(), genotype);
        return num_inds_with_genotype_at_snp_set(snp_set, genotype);
    }

    template<class PhenoType>
    PhenoType
    Instance<PhenoType>::
    phenotype(Ind ind) const {
        return phenotypes_.at(ind);
    }

    template<class PhenoType>
    std::vector<PhenoType> Instance<PhenoType>::all_phenotypes() const {
        return { phenotypes_.begin(), phenotypes_.end() };
    }

    template<class PhenoType>
    typename Instance<PhenoType>::const_ind_iterator
    Instance<PhenoType>::
    genotypes_of_all_inds_begin(SNP snp) const {
        return (genotypes_.cbegin() + (snp * num_inds_));
    }

    template<class PhenoType>
    typename Instance<PhenoType>::const_ind_iterator
    Instance<PhenoType>::
    genotypes_of_all_inds_end(SNP snp) const {
        return (genotypes_.cbegin() + ((snp + 1) * num_inds_));
    }

    template<class PhenoType>
    typename Instance<PhenoType>::const_snp_iterator
    Instance<PhenoType>::
    genotypes_at_all_snps_begin(Ind ind) const {
        return const_snp_iterator(genotypes_.cbegin() + ind, num_inds_, ind);
    }

    template<class PhenoType>
    typename Instance<PhenoType>::const_snp_iterator
    Instance<PhenoType>::
    genotypes_at_all_snps_end(Ind ind) const {
        return const_snp_iterator(genotypes_.cbegin() + ind, num_inds_, ind + num_inds_ * num_snps_);
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    shuffle_phenotypes() {
        std::shuffle(phenotypes_.begin(), phenotypes_.end(), urng_);
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    restore_phenotypes() {
        phenotypes_ = original_phenotypes_;
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    set_seed(std::size_t seed) {
        urng_.seed(seed);
    }

    template<class PhenoType>
    bool
    Instance<PhenoType>::
    categorical_phenotypes() const {
        return not quantitative_phenotypes();
    }

    template<class PhenoType>
    const std::vector<SNP> &
    Instance<PhenoType>::
    disease_snps() const {
        if (disease_snps_.empty()) {
            throw Error("No disease SNP set available.");
        }
        return disease_snps_;
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    set_disease_snps(const std::vector<SNP> & disease_snps) {
        disease_snps_ = disease_snps;
        for (SNP snp : disease_snps_) {
            if (snp >= num_snps_) {
                throw Error("The selected disease SNP " + std::to_string(snp) + " does not exist.");
            }
        }
        if (disease_snps_.empty()) {
            throw Error("The selected disease SNP set is empty.");
        }
        std::set<double> test_unique(disease_snps_.begin(), disease_snps_.end());
        if (test_unique.size() != disease_snps_.size()) {
            throw Error("The selected disease SNP set contains duplicates.");
        }
    }

    template<class PhenoType>
    std::string
    Instance<PhenoType>::
    snp_descriptor(SNP snp) const {
        if (snp >= num_snps_) {
            throw Error("The selected SNP " + std::to_string(snp) + " does not exist.");
        }
        return rs_ids_.at(snp);
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    load_csv_(const std::string & filename, bool snps_as_rows, bool snp_info_comes_first, std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data) {
        CSVParser csv_parser;
        csv_parser.parse(filename);
        num_snps_ = snps_as_rows ? csv_parser.num_rows() - 1 : csv_parser.num_columns() - 1;
        num_inds_ = snps_as_rows ? csv_parser.num_columns() - 1 : csv_parser.num_rows() - 1;
        std::string field;
        std::string index;

        // Determine which individuals should be skipped.
        std::vector<bool> skip;
        std::size_t num_skipped_inds{construct_folds_(num_folds, fold_id, cv_data, skip)};

        // Determine indices of first and of last individual as well as of SNP information.
        Ind first_ind{0};
        Ind last_ind{num_inds_ - 1};
        std::size_t index_snp_info{num_inds_};
        std::size_t skip_ind_shift{0};
        if (snp_info_comes_first) {
            skip_ind_shift = 1;
            first_ind = 1;
            last_ind = num_inds_;
            index_snp_info = 0;
        }

        // Load genotypes.
        for (SNP snp{0}; snp < num_snps_; snp++) {
            for (Ind ind{first_ind}; ind <= last_ind + num_skipped_inds; ind++) {
                if (not skip.at(ind - skip_ind_shift)) {
                    field = snps_as_rows ? csv_parser.cell(snp, ind) : csv_parser.cell(ind, snp);
                    try {
                        genotypes_.emplace_back(static_cast<GenoType>(std::stoul(field)));
                    }
                    catch (...) {
                        throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) + " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
                    }
                    if (genotypes_.back() > 2) {
                        throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) + " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
                    }
                }
            }
        }

        // Load phenotypes.
        for (Ind ind{first_ind}; ind <= last_ind + num_skipped_inds; ind++) {
            if (not skip.at(ind - skip_ind_shift)) {
                field = snps_as_rows ? csv_parser.cell(num_snps_, ind) : csv_parser.cell(ind, num_snps_);
                phenotypes_.emplace_back(parse_phenotype_(field, ind));
            }
        }

        // Load RS IDs.
        for (SNP snp{0}; snp < num_snps_; snp++) {
            field = snps_as_rows ? csv_parser.cell(snp, index_snp_info) : csv_parser.cell(index_snp_info, snp);
            rs_ids_.push_back(field);
        }
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    load_json_(const std::string & filename, std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data) {
        rapidjson::Document json_doc;
        try {
            // std::ifstream json_file (filename);
            // rapidjson::IStreamWrapper json_file_wrapper(json_file);
            // json_doc.ParseStream(json_file_wrapper);

            FILE* fp = fopen(filename.c_str(), "r");
            char readBuffer[65536];
            rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
            json_doc.ParseStream(is);

        }
        catch (...) {
            throw Error("The file " + filename + " cannot be opened.");
        }
        try {
            num_snps_ = json_doc["num_snps"].GetInt();
        }
        catch (...) {
            throw Error("The input file must contain a field \"num_snps\" whose value must be convertible to an int greater 0.");
        }
        if (num_snps_ <= 0) {
            throw Error("The input file must contain a field \"num_snps\" whose value must be convertible to an int greater 0.");
        }
        try {
            num_inds_ = json_doc["num_inds"].GetInt();
        }
        catch (...) {
            throw Error("The input file must contain a field \"num_inds\" whose value must be convertible to an int greater 0.");
        }
        if (num_inds_ <= 0) {
            throw Error("The input file must contain a field \"num_inds\" whose value must be convertible to an int greater 0.");
        }

        // Determine which individuals should be skipped.
        std::vector<bool> skip;
        std::size_t num_skipped_inds{construct_folds_(num_folds, fold_id, cv_data, skip)};

        // Load genotypes.
        if (!json_doc.HasMember("genotype")) {
            throw Error("The input file must contain a field \"genotype\".");
        }
        SNP snp{0};
        Ind ind{0};
        for (auto & row : json_doc["genotype"].GetArray()) {
            ind = 0;
            for (auto & cell : row.GetArray()) {
                if (not skip.at(ind)) {

                    try {
                        genotypes_.emplace_back(cell.GetUint());
                    }
                    catch (...) {
                        throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) + " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
                    }
                    if (genotypes_.back() > 2) {
                        throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) + " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
                    }
                }
                ind++;
            }
            if (ind != num_inds_ + num_skipped_inds) {
                throw Error("The actual number of individuals " + std::to_string(ind) + " for SNP " + std::to_string(snp) + " in the \"genotype\" field does not match the number of individuals " + std::to_string(num_inds_ + num_skipped_inds) + " specified in the \"num_inds\" field.");
            }
            snp++;
        }
        if (snp != num_snps_) {
            throw Error("The actual number of SNPs in the \"genotype\" field  does not match the specified number of SNPs specified in the \"num_snps\" field.");
        }

        // Load phenotypes.
        if (!json_doc.HasMember("phenotype")) {
            throw Error("The input file must contain a field \"phenotype\".");
        }
        ind = 0;
        for (const auto & phenotype : json_doc["phenotype"].GetArray()) {
            if (not skip.at(ind)) {
                std::string phenotype_parsed;
                if (phenotype.IsString()) phenotype_parsed = phenotype.GetString();
                else if (phenotype.IsInt()) phenotype_parsed = std::to_string(phenotype.GetInt());
                else if (phenotype.IsDouble()) phenotype_parsed = std::to_string(phenotype.GetDouble());
                else if (phenotype.IsFloat()) phenotype_parsed = std::to_string(phenotype.GetFloat());
                else throw Error("Phenotype has unknown datatype. Need to be Integer or String.");
                phenotypes_.emplace_back(parse_phenotype_(phenotype_parsed, ind));
            }
            ind++;
        }
        if (ind != num_inds_ + num_skipped_inds) {
            throw Error("The actual number of individuals in the \"phenotype\" field  does not match the specified number of individuals specified in the \"num_inds\" field.");
        }

        // Load RS IDs.
        if (!json_doc.HasMember("snps")) {
            throw Error("The input file must contain a field \"snps\".");
        }
        snp = 0;
        for (auto & snp_info : json_doc["snps"].GetArray()) {
            int counter = 0;
            for (auto & cell : snp_info.GetArray()) {

                if(counter==0)
                {
                    std::string rs_id_val = cell.GetString();
                    rs_ids_.push_back(rs_id_val);
                }
                else if(counter==1)
                {
                    rs_ids_chromosomes_.emplace_back(cell.GetString());
                }
                else
                {
                    break;
                }
                counter++;
            }
            snp++;
        }
        for (const auto & snp_info : json_doc["mafs"].GetArray())
        {
            double maf = snp_info.GetDouble();
            rs_ids_maf_.emplace_back(maf);
        }
        if (snp != num_snps_) {
            throw Error("The actual number of SNPs in the \"snps\" field  does not match the specified number of SNPs specified in the \"num_snps\" field.");
        }

        // Load disease SNPs (if available).
        if (!json_doc.HasMember("disease_snps")) {
            return;
        }
        disease_snps_.clear();
        for (auto & disease_snp : json_doc["disease_snps"].GetArray()) {
            disease_snps_.emplace_back(disease_snp.GetUint());
            if (disease_snps_.back() >= num_snps_) {
                throw Error("The selected disease SNP " + std::to_string(disease_snps_.back()) + " does not exist");
            }
        }
        std::set<double> test_unique(disease_snps_.begin(), disease_snps_.end());
        if (test_unique.size() != disease_snps_.size()) {
            throw Error("The selected disease SNP set contains duplicates.");
        }
    }



    template<class PhenoType>
    std::size_t
    Instance<PhenoType>::
    construct_folds_(std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data, std::vector<bool> & skip) {
        if (num_folds > num_inds_) {
            throw Error("The specified number of folds exceeds the number of individuals in the input data.");
        }
        if (fold_id >= num_folds) {
            throw Error("The specified fold ID has to be smaller than the specified number of folds.");
        }
        bool skip_most_inds{cv_data == options::DataPurpose::TRAINING or num_folds == 1 ? false : true};
        skip = std::vector<bool>(num_inds_, skip_most_inds);
        if (num_folds == 1) {
            return 0;
        }
        std::size_t size_fold{num_inds_ / num_folds};
        std::size_t remainder{num_inds_ % num_folds};
        Ind ind_in_fold{fold_id * size_fold + std::min(fold_id, remainder)};
        Ind last_ind_in_fold{(fold_id + 1) * size_fold - 1 + std::min(fold_id + 1, remainder)};
        std::size_t num_inds_in_fold{0};
        while (ind_in_fold <= last_ind_in_fold) {
            skip[ind_in_fold++] = not skip_most_inds;
            num_inds_in_fold++;
        }
        std::size_t num_skipped_inds{skip_most_inds ? num_inds_ - num_inds_in_fold : num_inds_in_fold};
        if (skip_most_inds) {
            num_inds_ = num_inds_in_fold;
        }
        else {
            num_inds_ -= num_inds_in_fold;
        }
        return num_skipped_inds;
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    save_bin(const std::string & filename) {
        std::ofstream file;
        file.open(filename, std::ofstream::out | std::ofstream::binary);

        if (!file.is_open() || !file.good()) {
            throw Error("The file " + filename + " cannot be opened.");
        }

        long val = num_snps_;
        file.write(reinterpret_cast<char*>(&val), sizeof(val));

        val = num_inds_;
        file.write(reinterpret_cast<char*>(&val), sizeof(val));


        // save genotypes.
        size_t geno_pos = 0;
        char *genotypes_buffer = new char[num_inds_];
        for (size_t snp = 0; snp < num_snps_; snp++) {
            for (size_t ind = 0; ind < num_inds_; ind++) {
                genotypes_buffer[ind] = genotypes_[geno_pos++];
            }
            file.write(genotypes_buffer, sizeof(char) * num_inds_);
        }
        delete[] genotypes_buffer;

        // save phenotypes.
        char dataType = quantitative_phenotypes() ? 1 : 2;
        file.write(&dataType, sizeof(dataType));
        if (dataType == 1) { // double
            double *phenotype_buffer = new double[num_inds_];
            for (size_t ind = 0; ind < num_inds_; ind++) {
                phenotype_buffer[ind] = phenotypes_[ind];
            }
            file.write(reinterpret_cast<char*>(phenotype_buffer), sizeof(double) * num_inds_);
            delete[] phenotype_buffer;
        } else if (dataType == 2) {
            unsigned long *phenotype_buffer = new unsigned long[num_inds_];
            for (size_t ind = 0; ind < num_inds_; ind++) {
                phenotype_buffer[ind] = phenotypes_[ind];
            }
            file.write(reinterpret_cast<char*>(phenotype_buffer), sizeof(unsigned long) * num_inds_);
            delete[] phenotype_buffer;
        }

        char hasChromosomeInfo = rs_ids_chromosomes_.empty() ? 0 : 1;
        file.write(&hasChromosomeInfo, sizeof(hasChromosomeInfo));

        char hasMafInfo = rs_ids_maf_.empty() ? 0 : 1;
        file.write(&hasMafInfo, sizeof(hasMafInfo));

        // save RS IDs.
        for (size_t snp = 0; snp < num_snps_; snp++) {
            // name
            const std::string& name = rs_ids_[snp];
            unsigned char name_length = name.size();
            file.write(reinterpret_cast<char*>(&name_length), sizeof(name_length));
            file.write(name.data(), sizeof(char) * name_length);

            // chromosome
            if (!rs_ids_chromosomes_.empty()) {
                const std::string &chromosome = rs_ids_chromosomes_[snp];
                unsigned char chromosome_length = chromosome.size();
                file.write(reinterpret_cast<char *>(&chromosome_length), sizeof(chromosome_length));
                file.write(chromosome.data(), sizeof(char) * chromosome_length);
            }

            // maf
            if (!rs_ids_maf_.empty()) {
                double maf = rs_ids_maf_[snp];
                file.write(reinterpret_cast<char *>(&maf), sizeof(maf));
            }
        }

        // save disease SNPs (if available).
        unsigned long num_disease_snps = disease_snps_.size();
        file.write(reinterpret_cast<char*>(&num_disease_snps), sizeof(num_disease_snps));
        if (num_disease_snps == 0){
            return;
        }

        unsigned long *disease_snps_buffer = new unsigned long[num_disease_snps];

        for (size_t dis_snp = 0; dis_snp < num_disease_snps; dis_snp++) {
            disease_snps_buffer[dis_snp] = disease_snps_[dis_snp];
        }
        file.write(reinterpret_cast<char*>(&disease_snps_buffer), sizeof(unsigned long) * num_disease_snps);
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    load_bin_(const std::string & filename, std::size_t num_folds, std::size_t fold_id, options::DataPurpose cv_data) {
        std::ifstream file;
        file.open(filename, std::ifstream::in | std::ifstream::binary);

        if (!file.is_open() || !file.good()) {
            throw Error("The file " + filename + " cannot be opened.");
        }
        try {
            long val = 0;
            file.read(reinterpret_cast<char*>(&val), sizeof(val));
            num_snps_ = val;
        }
        catch (...) {
            throw Error("The input file must contain a field \"num_snps\" whose value must be convertible to an int greater 0.");
        }
        if (num_snps_ <= 0) {
            throw Error("The input file must contain a field \"num_snps\" whose value must be convertible to an int greater 0.");
        }
        try {
            long val = 0;
            file.read(reinterpret_cast<char*>(&val), sizeof(val));
            num_inds_ = val;
        }
        catch (...) {
            throw Error("The input file must contain a field \"num_inds\" whose value must be convertible to an int greater 0.");
        }
        if (num_inds_ <= 0) {
            throw Error("The input file must contain a field \"num_inds\" whose value must be convertible to an int greater 0.");
        }

        // Determine which individuals should be skipped.
        std::vector<bool> skip;
        std::size_t num_skipped_inds{construct_folds_(num_folds, fold_id, cv_data, skip)};

        // Load genotypes.
        char *genotypes_buffer = new char[num_inds_];
        for (size_t snp = 0; snp < num_snps_; snp++) {
            file.read(genotypes_buffer, sizeof(char) * num_inds_);

            for (size_t ind = 0; ind < num_inds_; ind++) {
                if(not skip.at(ind)) {
                    genotypes_.push_back(genotypes_buffer[ind]);
                    if (genotypes_.back() > 2) {
                        throw Error("The input file contains invalid genotype for individual " + std::to_string(ind) +
                                    " at SNP " + std::to_string(snp) + ". Expected: 0, 1, or 2.");
                    }
                }
            }
        }
        delete[] genotypes_buffer;

        // Load phenotypes.
        char dataType = 0;
        file.read(&dataType, sizeof(dataType));
        if (dataType == 1) { // double
            double *phenotype_buffer = new double[num_inds_];
            file.read(reinterpret_cast<char*>(phenotype_buffer), sizeof(double) * num_inds_);
            for (size_t ind = 0; ind < num_inds_; ind++) {
                if (not skip.at(ind)) {
                    phenotypes_.emplace_back(parse_phenotype_(std::to_string(phenotype_buffer[ind]), ind));
                }
            }
            delete[] phenotype_buffer;
        } else if (dataType == 2) {
            unsigned long *phenotype_buffer = new unsigned long[num_inds_];
            file.read(reinterpret_cast<char*>(phenotype_buffer), sizeof(unsigned long) * num_inds_);
            for (size_t ind = 0; ind < num_inds_; ind++) {
                if (not skip.at(ind)) {
                    phenotypes_.emplace_back(parse_phenotype_(std::to_string(phenotype_buffer[ind]), ind));
                }
            }
            delete[] phenotype_buffer;
        }

        char hasChromosomeInfo = 0;
        file.read(&hasChromosomeInfo, sizeof(hasChromosomeInfo));

        char hasMafInfo = 0;
        file.read(&hasMafInfo, sizeof(hasMafInfo));

        // Load RS IDs.
        for (size_t snp = 0; snp < num_snps_; snp++) {
            // name
            unsigned char name_length = 0;
            file.read(reinterpret_cast<char*>(&name_length), sizeof(name_length));

            char *name_buffer = new char[name_length];
            file.read(name_buffer, sizeof(char) * name_length);
            std::string name (name_buffer, name_length);

            rs_ids_.push_back(name);

            // chromosome
            if(hasChromosomeInfo == 1) {
                unsigned char chromosome_length = 0;
                file.read(reinterpret_cast<char *>(&chromosome_length), sizeof(chromosome_length));

                char *chromosome_buffer = new char[chromosome_length];
                file.read(chromosome_buffer, sizeof(char) * chromosome_length);
                std::string chromosome(chromosome_buffer, chromosome_length);

                rs_ids_chromosomes_.push_back(chromosome);
            }


            // maf
            if (hasMafInfo == 1) {
                double maf = 0;
                file.read(reinterpret_cast<char *>(&maf), sizeof(maf));
                rs_ids_maf_.push_back(maf);
            }
        }

        // Load disease SNPs (if available).
        unsigned long num_disease_snps = 0;
        file.read(reinterpret_cast<char*>(&num_disease_snps), sizeof(num_disease_snps));
        if (num_disease_snps == 0){
            return;
        }

        disease_snps_.clear();
        unsigned long *disease_snps_buffer = new unsigned long[num_disease_snps];
        file.read(reinterpret_cast<char*>(&disease_snps_buffer), sizeof(unsigned long) * num_disease_snps);
        for (size_t dis_snp = 0; dis_snp < num_disease_snps; dis_snp++) {
            disease_snps_.push_back(disease_snps_buffer[dis_snp]);
        }
        std::set<SNP> test_unique(disease_snps_.begin(), disease_snps_.end());
        if (test_unique.size() != disease_snps_.size()) {
            throw Error("The selected disease SNP set contains duplicates.");
        }
    }

    template<class PhenoType>
    void
    Instance<PhenoType>::
    load_csv_cov_(const std::string & filename) {
        std::ifstream csv_file(filename);
        if (!csv_file) {
            throw Error("The file " + filename + " cannot be opened.");
        }

        std::string line;
        std::vector<std::string> headers;
        std::getline(csv_file, line);
        std::istringstream iss(line);
        std::string header;
        std::getline(iss, header, ','); // read the first header separately
        headers.push_back(header);
        while (std::getline(iss, header, ',')) {
            headers.push_back(header);
        }

        covariates_.resize(0, headers.size()-1);

        while (std::getline(csv_file, line)) {
            std::istringstream iss(line);
            std::vector<double> covariates_row(headers.size() - 1);
            for (std::size_t i = 0; i < headers.size() - 1; ++i) {
                std::string value_str;
                std::getline(iss, value_str, ','); // specify comma as delimiter
                double value;
                try {
                    value = std::stod(value_str);
                }
                catch (...) {
                    throw Error("The input file contains an invalid value: " + value_str);
                }
                covariates_row[i] = value;
            }
            covariates_.conservativeResize(covariates_.rows() + 1, Eigen::NoChange);
            covariates_.row(covariates_.rows() - 1) = Eigen::Map<Eigen::RowVectorXd>(covariates_row.data(), covariates_row.size());
        }

    }
}



//#endif /* SRC_MODEL_INSTANCE_IPP_ */
