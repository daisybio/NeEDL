//
// Created by juli on 25.05.22.
//

#include "SNPStorage.hpp"

#include <utility>
#include "SNP_t.hpp"
#include "../util/types.hpp"
#include "../util/helper_functions.hpp"

namespace epi {
    const std::string SNPStorage::EMPTY_STR = "";
    bool SNPStorage::any_snp_storage_created = false;
    std::shared_ptr<SNPStorage> SNPStorage::currentSnpStorage = nullptr;

    void SNPStorage::clear() {
        snp_data.clear();
        name_map.clear();
    }

    SNP_t SNPStorage::by_instance_id(size_t instance_id) {
        if (instance_id >= num_snps()) {
            throw epi::SNPNotFoundError("SNP with instance id '" + std::to_string(instance_id) + "' not found", std::to_string(instance_id));
        }
        return SNP_t(instance_id);
    }

    SNP_t SNPStorage::by_name(const std::string& name) {
        auto ptr = name_map.find(name);
        if (ptr == name_map.end()) {
            throw epi::SNPNotFoundError("SNP with name '" + name + "' not found", name);
        }
        return by_instance_id(ptr->second);
    }

    bool SNPStorage::contains_name(const std::string& name) {
        auto ptr = name_map.find(name);
        return ptr != name_map.end();
    }

    void SNPStorage::add_SNP_annotations(const std::vector<std::pair<SNP_t, std::string>>& annotations) {
        for (auto & anno : annotations) {
            snp_data[anno.first.value].annotations.push_back(anno.second);
            auto map_item = annotations_map.find(anno.second);
            if (map_item == annotations_map.end()) {
                std::vector<size_t> v;
                v.push_back(anno.first.value);
                annotations_map.insert({anno.second, v });
            } else {
                map_item->second.push_back(anno.first.value);
            }
        }
    }

    std::vector<SNP_t> SNPStorage::by_annotation(const std::string &name, bool include_removed) {
        auto map_item = annotations_map.find(name);
        if (map_item == annotations_map.end()) {
            return {};
        } else {
            std::vector<SNP_t> items;
            items.reserve(map_item->second.size());
            for (auto & x : map_item->second) {
                if (include_removed || !snp_data[x].is_removed) items.push_back(SNP_t(x));
            }
            return items;
        }
    }

    std::unordered_map<std::string, std::vector<size_t>> SNPStorage::get_annotations_map() {
        return annotations_map;
    }

    void SNPStorage::set_removed(SNP_t snp, bool is_removed) {
        snp_data[snp.value].is_removed = is_removed;
    }

    SNPStorage::all_iterable SNPStorage::all() {
        return SNPStorage::all_iterable(this);
    }

    double SNPStorage::get_snp_mma(SNP_t snp) {
        return snp_data[snp.value].mma_score;
    }

    void SNPStorage::set_snp_mma(SNP_t snp, double mma) {
        snp_data[snp.value].mma_score = mma;
    }

    bool SNPStorage::snp_is_removed(SNP_t snp) {
        return snp_data[snp.value].is_removed;
    }

    /*
    SNP SNPStorage::snp_get_instance_id(const SNP_t &snp) const {
        return 0;
    }

    const std::string &SNPStorage::snp_get_name(const SNP_t &snp) const {
        return EMPTY_STR;
    }

    const std::string &SNPStorage::snp_get_chromosome(const SNP_t &snp) const {
        return EMPTY_STR;
    }

    double SNPStorage::snp_get_maf(const SNP_t &snp) const {
        return 0;
    }

    void SNPStorage::save(std::string path) {

    }
     */
    SNPStorage::all_iterator::all_iterator(size_t _max, size_t _pos) {
        max = _max;
        val = std::min(_pos, max);
    }

    SNPStorage::all_iterator &SNPStorage::all_iterator::operator++() {
        val = std::min(val + 1, max);
        return *this;
    }

    SNPStorage::all_iterator SNPStorage::all_iterator::operator++(int) {
        all_iterator retval = *this;
        ++(*this);
        return retval;
    }

    bool SNPStorage::all_iterator::operator!=(SNPStorage::all_iterator other) const {
        return val != other.val;
    }

    bool SNPStorage::all_iterator::operator==(SNPStorage::all_iterator other) const  {
        return val == other.val;
    }

    SNP_t SNPStorage::all_iterator::operator*() const {
        return SNP_t(val);
    }

    SNPStorage::all_iterator& SNPStorage::all_iterator::operator+=(long n) {
        val = std::min(val + n, max);
        return *this;
    }

    SNPStorage::all_iterator& SNPStorage::all_iterator::operator-=(long n) {
        val = std::min(val - n, max);
        return *this;
    }


    SNPStorage::all_iterator operator+(SNPStorage::all_iterator iter, long n) {
        iter += n;
        return iter;
    }

    SNPStorage::all_iterator operator-(SNPStorage::all_iterator iter, long n) {
        iter -= n;
        return iter;
    }


    long operator-(SNPStorage::all_iterator left, SNPStorage::all_iterator right) {
       return left.val - right.val;
    }

    std::vector<std::string> SNPStorage::snp_get_annotations(const SNP_t& snp) {
        return snp_data[snp.value].annotations;
    }

    void SNPStorage::snp_set_variable_attribute(const SNP_t &snp, const std::string& key, const std::string& value) {
        auto res = snp_data[snp.value].variable_attributes.insert({ key, value });
        if (!res.second) res.first->second = value;
    }

    std::vector<std::string> SNPStorage::snp_get_variable_attribute_keys(const SNP_t &snp) {
        const auto &attributes = snp_data[snp.value].variable_attributes;
        std::vector<std::string> keys;
        for (auto & item : attributes) {
            keys.push_back(item.first);
        }
        return keys;

    }

    std::string SNPStorage::snp_get_variable_attribute(const SNP_t &snp, const std::string& key) {
        const auto &attributes = snp_data[snp.value].variable_attributes;
        auto item = attributes.find(key);
        if (item == attributes.end()) {
            return "";
        } else {
            return item->second;
        }
    }

    SNPStorage::SNPStorage() {
#ifndef ALLOW_MULTIPLE_SNPSTORAGES
        if (any_snp_storage_created) {
            throw epi::Error("Tried to create more than one instances of SNPStorage. This is currently not supported.");
        }
#endif
        any_snp_storage_created = true;
    }

    SNPStorage::~SNPStorage() {
        any_snp_storage_created = false;
    }

    void
    SNPStorage::snp_set_or_add_variable_attribute(const SNP_t &snp, const std::string &key, const std::string &value,
                                                  const char sep) {
        auto item = snp_data[snp.value].variable_attributes.find(key);
        if (item == snp_data[snp.value].variable_attributes.end())
            snp_data[snp.value].variable_attributes.insert({ key, value });
        else {
            auto splits = string_split(item->second, sep);
            splits.push_back(value);
            std::sort(splits.begin(), splits.end());
            auto last = std::unique(splits.begin(), splits.end());
            splits.erase(last, splits.end());
            item->second = "";
            bool is_first = true;
            for (const auto &s: splits) {
                if (is_first) is_first = false;
                else item->second += sep;

                item->second += s;
            }
        }
    }

    SNP SNPStorage::snp_get_instance_id(const SNP_t &snp) const {
        return snp.value;
    }

    SNP_t SNPStorage::all_iterator::operator[](long n) {
        return SNP_t(n);
    }

    SNPStorage::all_iterable::all_iterable(SNPStorage *_owner) {
        owner = _owner;
    }

    SNPStorage::all_iterator SNPStorage::all_iterable::begin() {
        return {owner->snp_data.size(), 0};
    }

    SNPStorage::all_iterator SNPStorage::all_iterable::end() {
        return {owner->snp_data.size(), size_t(-1)};
    }
    template<>
    std::vector<std::vector<size_t>> SNPStorage_PhenoType<CategoricalPhenoType>::get_individuals_per_category(
            const SNPSet &set) const {

        std::vector<SNP> set_instance {};
        for (auto &snp : set) {
            set_instance.push_back(snp.value);
        }

        std::vector<epi::Ind> individuals;
        instance->inds_with_nonzero_genotype_at_snp_set(set_instance, individuals);

        std::vector<std::vector<size_t>> count_map (instance->num_categories());

        for (auto & ind : individuals) {
            auto phenotype = instance->phenotype(ind);
            count_map[phenotype].push_back(ind);
        }

        return count_map;
    }

    template<>
    std::vector<size_t> SNPStorage_PhenoType<CategoricalPhenoType>::get_num_individuals_per_category() const {
        auto all_phenotypes = instance->all_phenotypes();

        std::vector<size_t> count_map (instance->num_categories(), 0);

        for (auto & phenotype : all_phenotypes) {
            count_map[phenotype] ++;
        }

        return count_map;
    }


    template<>
    size_t SNPStorage_PhenoType<CategoricalPhenoType>::num_categories() const {
        return instance->num_categories();
    }

    SNPStorage_WithoutGeno::SNPStorage_WithoutGeno(std::string plink_file_path) {
        // read bim file
        CSVParser bim_parser;
        bim_parser.parse(plink_file_path, '\t');


        // create snp_data
        snp_data = std::vector<SNPData>(bim_parser.num_rows());

        // add names to lookup map
        for (size_t i = 0; i < bim_parser.num_rows(); i++) {
            const auto &name = bim_parser.cell(i, 1);
            if (name_map.find(name) != name_map.end()) {
                throw epi::Error("SNP name " + name + " found multiple times.");
            }

            name_map.insert({ name, i });
            snp_names.push_back(name);
            snp_chromosomes.push_back(bim_parser.cell(i, 0));

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

    void SNPStorage_WithoutGeno::init_model_containers(size_t num_threads) {

    }

    const std::string &SNPStorage_WithoutGeno::snp_get_name(const SNP_t &snp) const {
        return snp_names[snp.value];
    }

    const std::string &SNPStorage_WithoutGeno::snp_get_chromosome(const SNP_t &snp) const {
       return snp_chromosomes[snp.value];
    }

    double SNPStorage_WithoutGeno::snp_get_maf(const SNP_t &snp) const {
        throw epi::Error("MAF information not available in BIM-only SNPStorage.");
    }

    bool SNPStorage_WithoutGeno::has_maf_information() const {
        return false;
    }

    void SNPStorage_WithoutGeno::set_maf_information(std::vector<double> maf_data) {
        throw epi::Error("MAF information cannot be set for BIM-only SNPStorage.");
    }

    void SNPStorage_WithoutGeno::save(std::string path) {
        throw epi::Error("BIM-only SNPStorage cannot be saved.");
    }

    double
    SNPStorage_WithoutGeno::calculate_score(const SNPSet &set, options::EpistasisScore model, size_t thread_index) {
        throw epi::Error("Score calculation not available in BIM-only SNPStorage.");
    }

    double SNPStorage_WithoutGeno::calculate_monte_carlo_score(const SNPSet &set, options::EpistasisScore,
                                                               size_t num_permutations, size_t thread_index) {
        throw epi::Error("Monte-carlo p-value calculation not available in BIM-only SNPStorage.");
    }

    void SNPStorage_WithoutGeno::request_score_method(options::EpistasisScore model, size_t thread_index) {

    }

    bool SNPStorage_WithoutGeno::need_minimizing_score(options::EpistasisScore model, size_t thread_index) {
        return false;
    }

    bool SNPStorage_WithoutGeno::is_categorical() const {
        return false;
    }

    size_t SNPStorage_WithoutGeno::num_snps() {
        return snp_names.size();
    }

    size_t SNPStorage_WithoutGeno::num_categories() const {
        return 0;
    }

    std::vector<std::vector<size_t>> SNPStorage_WithoutGeno::get_individuals_per_category(const SNPSet &set) const {
        throw epi::Error("Individuals per category not available in BIM-only SNPStorage.");
        return {};
    }

    std::vector<size_t> SNPStorage_WithoutGeno::get_num_individuals_per_category() const {
        throw epi::Error("Individuals per category not available in BIM-only SNPStorage.");
        return {};
    }

} // epi