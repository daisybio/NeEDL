//
// Created by juli on 25.05.22.
//

#ifndef GENEPISEEKER_SNP_T_HPP
#define GENEPISEEKER_SNP_T_HPP

#include <igraph.h>
#include <algorithm>
#include <string>
#include "../util/types.hpp"

#define MAXIMUM_SNP_SET_SIZE 10

namespace epi {

    class SNPStorage;
    class SNPData;
    class SNPNetwork;
    class LDTester;
    class SaveNetwork;

    union SNPEdge;

    class SNP_t {
    public:

        static const SNP INVALID_SNP = -1;

        class SNPHash {
        public:
            // id is returned as hash function
            size_t operator()(const SNP_t& t) const;
        };

        SNP_t(const SNP_t& snp);
        SNP_t(); // initialize to an empty/invalid SNP
        bool isValid() const;

        /*
        void set_maf(double maf_);
        size_t get_instance_id() const;
        const std::string& get_name() const;
        const std::string& get_chromosome() const;
        double get_maf() const;
         */

        bool operator<(const SNP_t& rhs) const;
        bool operator==(const SNP_t& rhs) const;
        SNP_t& operator=(const SNP_t&rhs);

        SNP value;

    private:

        friend class SNPNetwork;
        friend class SNPNetwork_base;
        friend class SNPStorage;
        friend union SNPEdge;
        friend class LDTester;
        friend class SaveNetwork;

        explicit SNP_t(SNP val);
    };

    union SNPEdge {
        struct {
            uint_fast32_t snp1 : 32;
            uint_fast32_t snp2 : 32;
        } bits{};
        uint_fast64_t val;

        SNPEdge(const SNP_t& first_, const SNP_t& last_);
        /**
         * constructs an invalid edge (the edge 0 <--> 0)
         */
        SNPEdge();
        bool operator==(const SNPEdge& rhs) const;
        bool operator<(const SNPEdge& rhs) const;
        std::string get_snp_string() const;
    };
/*
    class SNPEdge : public std::pair<SNP, SNP> {
    public:
        SNPEdge(const SNP& first_, const SNP& last_);
        bool operator==(const SNPEdge& rhs) const;
        bool operator<(const SNPEdge& rhs) const;
    };
    */

    class SNPEdgeHash {
    public:
        // id is returned as hash function
        size_t operator()(const SNPEdge& t) const;
    };

    class SNPSet {
    public:
        SNPSet();
        SNPSet(const std::vector<SNP_t> &snps);
        SNPSet(const SNP_t& snp);

        std::vector<SNP_t>::const_iterator begin() const;
        std::vector<SNP_t>::const_iterator end() const;

        void set_attribute(const std::string& key, const std::string& value);
        std::string get_attribute(const std::string& key) const;
        std::map<std::string, std::string> get_attributes() const;
        std::vector<std::string> get_attribute_keys() const;
        void clear_attributes();
        std::string get_snp_string() const;
        std::vector<std::string> get_annotations();

        double calculate_score(options::EpistasisScore model);

        SNP_t at(size_t pos);


        bool operator<(const SNPSet& rhs) const;
        bool operator==(const SNPSet& rhs) const;

        explicit operator const std::vector<SNP_t>& () const;
        const std::vector<SNP_t> vec() const;

        // attributes are copied/merged always --> scores are only kept if set did not change
        SNPSet& operator+=(const SNP_t &snp);
        SNPSet& operator+=(const SNPSet &other);
        SNPSet& operator-=(const SNP_t &snp);
        friend SNPSet operator+(const SNPSet &set, const SNP_t &snp);
        friend SNPSet operator+(const SNPSet &set1, const SNPSet &set2);
        friend SNPSet operator-(const SNPSet &set, const SNP_t &snp);
        SNPSet clone() const;


        size_t size() const;
        SNP_t operator[](size_t index) const;

    protected:
        std::shared_ptr<std::vector<SNP_t>> snps = nullptr;
        std::shared_ptr<std::map<std::string,std::string>> attributes = nullptr;
        std::shared_ptr<std::vector<double>> scores = nullptr;
        uint_fast32_t scores_calculated = 0;
    };

    class SNPSetHash {
    public:
        // id is returned as hash function
        size_t operator()(const SNPSet& t) const;
    };

} // epi

#ifdef HEADER_ONLY
#include "SNP_t.cpp"
#endif

#endif //GENEPISEEKER_SNP_T_HPP
