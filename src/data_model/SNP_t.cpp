//
// Created by juli on 25.05.22.
//

#include "SNP_t.hpp"

#include <utility>

namespace epi {

    SNP_t::SNP_t(const SNP_t& snp) : value(snp.value ) { }

    SNP_t::SNP_t() : value(INVALID_SNP) { }

    /*
    void SNP::set_maf(double maf_) {
        data->maf = maf_;
    }

    size_t SNP::get_instance_id() const {
        return data->instance_id;
    }
     */

    SNP_t::SNP_t(SNP val) : value(val) { }

    bool SNP_t::operator<(const SNP_t &rhs) const {
        return this->value < rhs.value;
    }

    bool SNP_t::operator==(const SNP_t &rhs) const {
        return this->value == rhs.value;
    }

    SNP_t& SNP_t::operator=(const SNP_t &rhs) {
        if (this == &rhs) return *this;
        this->value = rhs.value;
        return *this;
    }

    bool SNP_t::isValid() const {
        return value != INVALID_SNP;
    }

    /*
    const std::string &SNP::get_name() const {
        return this->data->name;
    }

    const std::string &SNP::get_chromosome() const {
        return this->data->chromosome;
    }

    double SNP::get_maf() const {
        return this->data->maf;
    }
     */

    size_t SNP_t::SNPHash::operator()(const SNP_t &t) const {
        return std::hash<size_t>{}(t.value);
    }

    bool SNPEdge::operator==(const SNPEdge &rhs) const {
        return val == rhs.val;
    }

    SNPEdge::SNPEdge(const SNP_t& first_, const SNP_t& second_) {
        bits = { std::min(first_.value, second_.value), std::max(first_.value, second_.value) };
    }

    bool SNPEdge::operator<(const SNPEdge &rhs) const {
        if (bits.snp1 == rhs.bits.snp1) return bits.snp2 < rhs.bits.snp2;
        return bits.snp1 < rhs.bits.snp1;
    }

    SNPEdge::SNPEdge() {
        bits = { 0, 0 };
    }



    size_t SNPEdgeHash::operator()(const SNPEdge &t) const {
        return std::hash<uint_fast64_t>{}(t.val);
    }


    SNPSet::SNPSet() {
        snps = std::make_shared<std::vector<SNP_t>>();
    }

    SNPSet::SNPSet(const std::vector<SNP_t> &snps) {
        if (snps.size() > MAXIMUM_SNP_SET_SIZE) {
            throw epi::Error("Tried to create a snp-set with size > " + std::to_string(MAXIMUM_SNP_SET_SIZE));
        }
        this->snps = std::make_shared<std::vector<SNP_t>>(snps.begin(), snps.end());
        std::sort(this->snps->begin(), this->snps->end());
    }

    SNPSet::SNPSet(const SNP_t& snp) {
        snps = std::make_shared<std::vector<SNP_t>>();
        snps->push_back(snp);
    }

    std::vector<SNP_t>::const_iterator SNPSet::begin() const {
        return snps->cbegin();
    }

    std::vector<SNP_t>::const_iterator SNPSet::end() const {
        return snps->cend();
    }

    void SNPSet::set_attribute(const std::string& key, const std::string& value) {
        if (attributes == nullptr) attributes = std::make_shared<std::map<std::string, std::string>>();
        auto res = attributes->insert({ key, value });
        if (!res.second) res.first->second = value;
    }

    std::string SNPSet::get_attribute(const std::string& key) const  {
        if (attributes == nullptr) return "";
       auto item = attributes->find(key);
       if (item == attributes->end()) {
           return "";
       } else {
           return item->second;
       }
    }

    std::map<std::string, std::string> SNPSet::get_attributes() const {
        if (attributes == nullptr) return {};
        return { *attributes };
    }

    std::vector<std::string> SNPSet::get_attribute_keys() const {
        if (attributes == nullptr) return {};

        std::vector<std::string> keys;
        for (auto & item : *attributes) {
            keys.push_back(item.first);
        }
        return keys;
    }

    bool SNPSet::operator<(const SNPSet &rhs) const {
        if (snps->size() == rhs.snps->size()) {
            for (size_t i = 0; i < snps->size(); i++) {
                if (snps->at(i) == rhs.snps->at(i)) continue;
                return snps->at(i) < rhs.snps->at(i);
            }
            return false;
        } else return snps->size() == rhs.snps->size();
    }

    bool SNPSet::operator==(const SNPSet &rhs) const {
        if (snps->size() != rhs.snps->size()) return false;

        for (size_t i = 0; i < snps->size(); i++) {
            if (!(snps->at(i) == rhs.snps->at(i))) return false;
        }
        return true;
    }

    size_t SNPSet::size() const {
        return snps->size();
    }

    SNP_t SNPSet::operator[](size_t index) const {
        return snps->at(index);
    }

    SNPSet operator+(const SNPSet &set, const SNP_t &snp) {
        SNPSet new_set = set.clone();
        new_set += snp;
        return new_set;
    }

    SNPSet operator-(const SNPSet &set, const SNP_t &snp) {
        SNPSet new_set = set.clone();
        new_set -= snp;
        return new_set;
    }

    SNPSet &SNPSet::operator+=(const SNP_t &snp) {
        if (std::find(snps->begin(), snps->end(), snp) == snps->end()) {
            snps->insert(snps->end(), snp);

            if (snps->size() > MAXIMUM_SNP_SET_SIZE) {
                throw epi::Error("Tried to create a snp-set with size > " + std::to_string(MAXIMUM_SNP_SET_SIZE));
            }

            std::sort(snps->begin(), snps->end());
            scores = nullptr;
            scores_calculated = 0;
        }

        return *this;
    }

    SNPSet &SNPSet::operator-=(const SNP_t &snp) {
        auto it = std::find(snps->begin(), snps->end(), snp);
        if (it != snps->end()) {
            snps->erase(it);
            scores = nullptr;
            scores_calculated = 0;
        }

        return *this;
    }

    SNPSet &SNPSet::operator+=(const SNPSet &other) {
        size_t size_before = snps->size();
        snps->insert(snps->end(), other.snps->begin(), other.snps->end());
        std::sort(snps->begin(), snps->end());
        snps->erase( unique( snps->begin(), snps->end() ), snps->end() );

        if (snps->size() != size_before) {
            if (snps->size() > MAXIMUM_SNP_SET_SIZE) {
                throw epi::Error("Tried to create a snp-set with size > " + std::to_string(MAXIMUM_SNP_SET_SIZE));
            }
            scores = nullptr;
            scores_calculated = 0;
        }

        return *this;
    }

    SNPSet SNPSet::clone() const {
        SNPSet new_set;

        new_set.snps = std::make_shared<std::vector<SNP_t>>(snps->begin(), snps->end());
        if (attributes != nullptr) new_set.attributes = std::make_shared<std::map<std::string, std::string>>(attributes->begin(), attributes->end());
        new_set.scores_calculated = scores_calculated;
        if (scores != nullptr) new_set.scores = std::make_shared<std::vector<double>>(scores->begin(), scores->end());

        return new_set;
    }

    SNPSet operator+(const SNPSet &set1, const SNPSet &set2) {
        SNPSet new_set = set1.clone();
        new_set += set2;
        return new_set;
    }

    void SNPSet::clear_attributes() {
        if (attributes != nullptr) attributes->clear();
    }

    const std::vector<SNP_t> SNPSet::vec() const {
        return *snps;
    }

    SNP_t SNPSet::at(size_t pos) {
        return snps->at(pos);
    }

    SNPSet::operator const std::vector<SNP_t> &() const {
        return *snps;
    }


    size_t SNPSetHash::operator()(const SNPSet &t) const {
        size_t hash = 0;
        auto hasher = SNP_t::SNPHash{};
        for (auto & snp : t) {
            hash ^= hasher(snp);
        }
        return hash;
    }
} // epi