//
// Created by eve on 15.08.21.
//

#ifndef GENEPISEEKER_JSON_FORMAT_HPP
#define GENEPISEEKER_JSON_FORMAT_HPP

#import <string>
#import <vector>
#import <unordered_set>
#import <boost/algorithm/string.hpp>
#import <boost/algorithm/string/replace.hpp>
#include <fstream>
#include <iostream>

using namespace std;


class JSON_format {
public:
    // Attributes related to the input file
    string json_path_;
    string json_filename_;

    // Actual class attributes to store json data
    int num_snps_;
    int num_inds_;
    string model_type_;
    int num_categories_;
    vector<vector<string>> genotype_;
    vector<string> phenotype_;
    vector<vector<string>> snps_;    // [0]snp_id   [1]chromosome   [2]genomic_position   [3]wildtype   [4]variant
    vector<int> disease_snps_;
    vector<double> mafs_;

public:
    // constructor declaration
    JSON_format();
    JSON_format(string json_filepath);

    void add_phenotypes(vector<string> phenotype);
    void add_json_entry(int index, vector<string> genotype, vector<string> snp, bool disease_snp);
    unordered_set<string> get_snp_list();
    string get_snp_string(string sep);
    string get_genotype_string_samplewise(int index, string sep);
    string get_genotype_string_snpwise(int index, string sep);
    string get_phenotype_string(string sep, string case_code, string control_code);

    void read_json(string json_filepath);
    void write_json(string out_file);

    int find(string snp);
    bool is_disease_snp(string snp);
    void calculate_mafs();

    void report_object_dimension();
    void report_object_probes();
    void write_snp_annotation(string out_dir);

private:
    // internal util methods for reading in the json file
    vector<string> tokenize(string line);
    string clean_token(string token, int mode);
    vector<string> split_at_string(string s, string delimiter);
    void assign_tokens(vector<string> tokens);

    // internal util methods for representing json data as string
    string genotype_to_json_string(int index);
    string phenotypes_to_json_string();
    string snp_to_json_string(int index);
    string disease_snps_to_json_string();
    string mafs_to_json_string();
};


#include "JSON_format.ipp"

#endif //GENEPISEEKER_JSON_FORMAT_HPP
