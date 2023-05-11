// combine_JSON_files.cpp
// Created by eve on 17.05.21.
//


#import <string>
#import <vector>
#import <unordered_set>
#include <iostream>
#include "JSON_format.hpp"

using namespace std;



class JSON_combiner{

private:
    string json_case_filepath_;
    string json_control_filepath_;
    string out_filepath_;

    JSON_format json_case = JSON_format();
    JSON_format json_control = JSON_format();
    JSON_format json_out = JSON_format();


    void print_help();
    const string currentDateTime();

    vector<string> get_intersecting_snpIDs();


public:
    JSON_combiner();
    JSON_combiner(string json_case, string json_control, string out_file);

    int handle_arguments(int argc, char* argv[]);
    string combine_JSON_files();
};



/*
int main(int argc, char* argv[]) {  //string json_1_filepath, string json_2_filepath, string out_file
    JSON_combiner jc = JSON_combiner();
    int success = jc.handle_arguments(argc, argv);
    if(success == 0)
        return 0;

    jc.combine_JSON_files();

    return 0;
}
 */



//########################################################################
//                              Constructors
//########################################################################

JSON_combiner::JSON_combiner() {
    cout << "\n************************************************************" << endl;
    cout << "                   Combine two JSON files                   " << endl;
    cout << "************************************************************\n" << endl << endl;
}

JSON_combiner::JSON_combiner(string json_1, string json_2, string out_file) {
    cout << "\n************************************************************" << endl;
    cout << "                   Combine two JSON files                   " << endl;
    cout << "************************************************************\n" << endl << endl;

    if(json_1 != "")
        json_case_filepath_ = json_1;
    else throw runtime_error("WARNING: filepath to case JSON not set.");

    if(json_2 != "")
        json_control_filepath_ = json_2;
    else throw runtime_error("WARNING: filepath to control JSON not set.");

    if(out_file != "")
        out_filepath_ = out_file;
    else throw runtime_error("WARNING: output filepath, where combined JSON should be written to, not set.");
}



//########################################################################
//
//                      Main functionality methode
//
//########################################################################

string JSON_combiner::combine_JSON_files() {
    vector<string> split;
    boost::split(split, out_filepath_, boost::is_any_of("/"));
    string out_dir = split[0];
    for(int i = 1; i < split.size()-1; i++)
        out_dir += "/" + split[i] ;


    // -----------------------------------------------------------------------
    // Read first JSON file
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tRead first json file " << json_case_filepath_ << endl;
    json_case.read_json(json_case_filepath_);
    json_case.report_object_dimension();
    json_case.report_object_probes();
    json_case.write_snp_annotation(out_dir);


    // -----------------------------------------------------------------------
    // Read second JSON file
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tRead second json file " << json_control_filepath_ << endl;
    json_control.read_json(json_control_filepath_);
    json_control.report_object_dimension();
    json_control.report_object_probes();
    json_control.write_snp_annotation(out_dir);


    // -----------------------------------------------------------------------
    // Determine intersecting SNP-IDs
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tDetermine intersecting SNP-IDs" << endl;
    vector<string> intersecting_ids = get_intersecting_snpIDs();



    // -----------------------------------------------------------------------
    // Build combined JSON object
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tBuild combined JSON object" << endl;
    JSON_format json_out = JSON_format();
    json_out.num_snps_ = intersecting_ids.size();
    json_out.num_inds_ = json_case.num_inds_ + json_control.num_inds_;
    if(json_case.model_type_ == json_control.model_type_)
        json_out.model_type_ = json_case.model_type_;
    else json_out.model_type_ = "ambitious";
    json_out.num_categories_ = max(json_case.num_categories_, json_control.num_categories_);
    cout << "\t\t" << currentDateTime() << "\tSuccessfully set general JSON annotation" << endl;

    json_out.phenotype_ = json_case.phenotype_;
    json_out.add_phenotypes(json_control.phenotype_);
    cout << "\t\t" << currentDateTime() << "\tSuccessfully set phenotypes" << endl;

    for(int i = 0; i < intersecting_ids.size(); i++){
        int json_1_index = json_case.find(intersecting_ids[i]);
        json_out.add_json_entry(i, json_case.genotype_[json_1_index], json_case.snps_[json_1_index], json_case.is_disease_snp(intersecting_ids[i]));
        int json_2_index = json_control.find(intersecting_ids[i]);
        json_out.add_json_entry(i, json_control.genotype_[json_2_index],json_control.snps_[json_2_index],json_control.is_disease_snp(intersecting_ids[i]));
    }
    cout << "\t\t" << currentDateTime() << "\tSuccessfully set genotypes, snp annotations and disease snps" << endl;

    json_out.calculate_mafs();
    cout << "\t\t" << currentDateTime() << "\tSuccessfully calculated mafs" << endl;
    json_out.report_object_dimension();


    // -----------------------------------------------------------------------
    // Write JSON
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tWrite JSON to " << out_filepath_ << endl;
    json_out.write_json(out_filepath_);


    // -----------------------------------------------------------------------
    // Wrap up
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tFinished" << endl;

    return(out_filepath_);
}



//########################################################################
//                     Basic administration functions
//########################################################################

int JSON_combiner::handle_arguments(int argc, char* argv[]) {
    if(argc == 1 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h"){
        print_help();
        return 0;
    }
    bool json_case_set = false;
    bool json_control_set = false;
    bool outfile_set = false;

    for(int i = 1; i < argc; i=i+2){
        string att_tag = std::string(argv[i]);
        if(i+1 >= argc){
            cout << "ATTENTION! Incorrect number of attributes detected" << endl << endl;
            print_help();
            return 0;
        }
        string att = std::string(argv[i+1]);
        if(att.substr(0, 2) == "--"){
            cout << "ATTENTION! Missing argument value detected for " << att_tag << endl << endl;
            print_help();
            return 0;
        }

        if(att_tag == "--json_case"){
            json_case_filepath_ = att;
            json_case_set = true;
        }
        else if(att_tag == "--json_control"){
            json_control_filepath_ = att;
            json_control_set = true;
        }
        else if(att_tag == "--out_file"){
            out_filepath_ = att;
            outfile_set = true;
        }
        else{
            cout << "ATTENTION! Not supported attribute detected: " << att_tag << endl;
            print_help();
            return 0;
        }
    }

    if(!json_case_set || !json_control_set){
        cout << "ATTENTION! Please set two mandatory JSON file paths" << endl << endl;
        print_help();
        return 0;
    }
    if(!outfile_set){
        cout << "ATTENTION! Please set mandatory output file path" << endl << endl;
        print_help();
        return 0;
    }
    return 0;
}

void JSON_combiner::print_help() {
    cout << "Usage message" << endl << "-------------" << endl << endl;
    cout << "Mandatory arguments:" << endl;
    cout << "--json_case\tPath to one of the JSON files that should be combined." << endl;
    cout << "--json_control\tPath to the other JSON file that should be combined." << endl;
    cout << "--out_file\t\tAbsolute path with filename of the resulting combined JSON." << endl << endl;
    cout << "Optional arguments:" << endl;
    cout << "--help, -h\tShow help message and exit." << endl << endl;
}

const string JSON_combiner::currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "[%Y/%m/%d-%X]", &tstruct);

    return buf;
}


vector<string> JSON_combiner::get_intersecting_snpIDs() {
    unordered_set<string> ids_case = json_case.get_snp_list();
    unordered_set<string> ids_control = json_control.get_snp_list();
    vector<string> intersecting_ids;

    for(string snp : ids_case){
        if(ids_control.find(snp) == ids_control.end())
            continue;
        intersecting_ids.push_back(snp);
    }

    cout << "\t\t" << currentDateTime() << "\tSuccessfully determined intersecting SNP-IDs" << endl;
    cout << "\t\t - Number of found SNP-IDs: " << intersecting_ids.size() << endl;

    return intersecting_ids;
}