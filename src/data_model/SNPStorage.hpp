//
// Created by juli on 25.05.22.
//

#ifndef GENEPISEEKER_SNPSTORAGE_HPP
#define GENEPISEEKER_SNPSTORAGE_HPP

#include <string>
#include <igraph.h>
#include <vector>
#include <unordered_map>
#include "SNP_t.hpp"
#include "../model/instance.hpp"
#include "../model/all_models.hpp"

namespace epi {
    class SNPData {
    public:
        SNPData() = default;
        virtual ~SNPData() {};
    private:
        // attributes
        std::vector<std::string> annotations {};
        std::map<std::string, std::string> variable_attributes {};
        double mma_score = 0;

        // this is used for filtering
        bool is_removed = false;

        // give SNPNetwork access to private arguments
        friend class SNP_t;
        friend class SNPStorage;
        friend union SNPEdge;
    };

    class SNPStorage {
    public:

    class all_iterator : public std::iterator<std::random_access_iterator_tag, SNP_t, long, const SNP_t *, SNP_t> {
    size_t val = 0;
    size_t max = 0;
    public:
        all_iterator(size_t max, size_t pos);
        all_iterator& operator++();
        all_iterator operator++(int);

        all_iterator& operator+=(long n);
        friend all_iterator operator+(all_iterator iter, long n);
        all_iterator& operator-=(long n);
        friend all_iterator operator-(all_iterator iter, long n);
        friend long operator-(all_iterator left, all_iterator right);
        SNP_t operator[](long n);

        bool operator==(all_iterator other) const;
        bool operator!=(all_iterator other) const;

        SNP_t operator*() const;
    };

    class all_iterable {
    public:
        explicit all_iterable(SNPStorage *_owner);
        all_iterator begin();
        all_iterator end();
    private:
        SNPStorage *owner;
    };

        void clear();

        all_iterable all();

        SNP_t by_instance_id(SNP instance_id);
        std::vector<SNP_t> by_annotation(const std::string& name, bool include_removed = false);
        SNP_t by_name(const std::string& name);
        void set_removed(SNP_t snp, bool is_removed = true);
        void set_snp_mma(SNP_t snp, double mma);
        double get_snp_mma(SNP_t snp);
        bool snp_is_removed(SNP_t snp);

        bool contains_name(const std::string& name);

        void add_SNP_annotations(const std::vector<std::pair<SNP_t, std::string>>& annotations);
        std::unordered_map<std::string, std::vector<size_t>> get_annotations_map();
        std::vector<std::string> snp_get_annotations(const SNP_t& snp);

        void snp_set_variable_attribute(const SNP_t& snp, const std::string& key, const std::string& value);
        std::vector<std::string> snp_get_variable_attribute_keys(const SNP_t& snp);
        std::string snp_get_variable_attribute(const SNP_t& snp, const std::string& key);
        void snp_set_or_add_variable_attribute(const SNP_t& snp, const std::string& key, const std::string& value, const char sep);


        // functions that need access to the instance
        virtual void init_model_containers(size_t num_threads) = 0;
        virtual SNP snp_get_instance_id(const SNP_t& snp) const;
        virtual const std::string& snp_get_name(const SNP_t& snp) const = 0;
        virtual const std::string& snp_get_chromosome(const SNP_t& snp) const = 0;
        virtual double snp_get_maf(const SNP_t& snp) const = 0;
        virtual bool has_maf_information() const = 0;
        virtual void set_maf_information(std::vector<double> maf_data) = 0;
        virtual void save(std::string path) = 0;
        inline double calculate_score(const SNPSet & set, options::EpistasisScore score) {
            return calculate_score(set, score, omp_get_thread_num());
        };
        virtual double calculate_score(const SNPSet & set, options::EpistasisScore, size_t thread_index) = 0;
        inline double calculate_monte_carlo_score(const SNPSet & set, options::EpistasisScore score, size_t num_permutations) {
            return calculate_monte_carlo_score(set, score, num_permutations, omp_get_thread_num());
        };
        virtual double calculate_monte_carlo_score(const SNPSet & set, options::EpistasisScore score, size_t num_permutations, size_t thread_index) = 0;
        virtual void request_score_method(options::EpistasisScore score, size_t thread_index) = 0;
        inline bool need_minimizing_score(options::EpistasisScore score) {
            return need_minimizing_score(score, omp_get_thread_num());
        };
        virtual bool need_minimizing_score(options::EpistasisScore, size_t thread_index) = 0;
        virtual size_t num_snps() = 0;
        virtual bool is_categorical() const = 0;
        virtual size_t num_categories() const = 0;
        virtual std::vector<std::vector<size_t>> get_individuals_per_category(const SNPSet & set) const = 0;
        virtual std::vector<size_t> get_num_individuals_per_category() const = 0;
        virtual void shuffle_phenotypes() = 0;


        SNPStorage();
        virtual ~SNPStorage();

        // global access to SNPStorage
        static std::shared_ptr<SNPStorage> currentSnpStorage;
    protected:
        static const std::string EMPTY_STR;

        static bool any_snp_storage_created;

        std::vector<SNPData> snp_data{};
        std::unordered_map<std::string, SNP> name_map;
        std::unordered_map<std::string, std::vector<SNP>> annotations_map;

    };

    /*!
     * this class is only used for the preprocessing pipeline and allows to create a SNP storage from a BIM file without genotype information
     */
    class SNPStorage_WithoutGeno : public SNPStorage {
    public:
        explicit SNPStorage_WithoutGeno(std::string plink_file_path);

        void init_model_containers(size_t num_threads) override;
        const std::string& snp_get_name(const SNP_t& snp) const override;
        const std::string& snp_get_chromosome(const SNP_t& snp) const override;
        double snp_get_maf(const SNP_t& snp) const override;
        bool has_maf_information() const override;
        void set_maf_information(std::vector<double> maf_data) override;
        void save(std::string path) override;
        double calculate_score(const SNPSet & set, options::EpistasisScore model, size_t thread_index) override;
        double calculate_monte_carlo_score(const SNPSet & set, options::EpistasisScore, size_t num_permutations, size_t thread_index) override;
        void request_score_method(options::EpistasisScore model, size_t thread_index) override;
        bool need_minimizing_score(options::EpistasisScore model, size_t thread_index) override;
        bool is_categorical() const override;
        size_t num_snps() override;
        size_t num_categories() const override;
        std::vector<std::vector<size_t>> get_individuals_per_category(const SNPSet & set) const override;
        std::vector<size_t> get_num_individuals_per_category() const override;
        void shuffle_phenotypes() override;

    private:
        std::vector<std::string> snp_names;
        std::vector<std::string> snp_chromosomes;

    };

    template<class PhenoType>
    class SNPStorage_PhenoType : public SNPStorage {
    public:
        explicit SNPStorage_PhenoType(std::shared_ptr<epi::Instance<PhenoType>> instance);

        void init_model_containers(size_t num_threads) override;


        const std::string& snp_get_name(const SNP_t& snp) const override;
        const std::string& snp_get_chromosome(const SNP_t& snp) const override;
        double snp_get_maf(const SNP_t& snp) const override;
        bool has_maf_information() const override;
        void set_maf_information(std::vector<double> maf_data) override;
        void save(std::string path) override;
        double calculate_score(const SNPSet & set, options::EpistasisScore model, size_t thread_index) override;
        double calculate_monte_carlo_score(const SNPSet & set, options::EpistasisScore, size_t num_permutations, size_t thread_index) override;
        void request_score_method(options::EpistasisScore model, size_t thread_index) override;
        bool need_minimizing_score(options::EpistasisScore model, size_t thread_index) override;
        bool is_categorical() const override;
        size_t num_snps() override;
        size_t num_categories() const override;
        std::vector<std::vector<size_t>> get_individuals_per_category(const SNPSet & set) const override;
        std::vector<size_t> get_num_individuals_per_category() const override;
        void shuffle_phenotypes() override;

    private:
        std::shared_ptr<epi::Instance<PhenoType>> instance;

        std::vector<std::shared_ptr<BayesianModel<PhenoType>>> bayesian_models;
        std::vector<std::shared_ptr<PenetranceModel<PhenoType>>> penetrance_models;
        std::vector<options::EpistasisScore> penetrance_model_scores;
        std::vector<std::shared_ptr<RegressionModel<PhenoType>>> regression_models;
        std::vector<options::EpistasisScore> regression_model_scores;
        std::vector<std::shared_ptr<VarianceModel<PhenoType>>> variance_models;
        // std::vector<size_t> score_calculation_nmax_snp_set_sizes;
    };




} // epi

#include "SNPStorage.ipp"

#ifdef HEADER_ONLY
#include "SNP_t_ext.cpp"
#include "SNPStorage.cpp"
#endif


#endif //GENEPISEEKER_SNPSTORAGE_HPP
