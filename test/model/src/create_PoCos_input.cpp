//
// create_PoCos_input.cpp
// Created by eve on 18.05.21.
//

#include "JSON_format.ipp"
#include <fstream>
#include <iostream>

using namespace std;



string json_filepath_;
string out_dir_;

JSON_format json = JSON_format();



int handle_arguments(int argc, char* argv[]);

void print_help();

void write_genotype_file();

void write_phenotype_file();

void write_snp_names();

const string currentDateTime();



int main(int argc, char* argv[]) {
    cout << "\n************************************************************" << endl;
    cout << "                     Create PoCos input                     " << endl;
    cout << "************************************************************\n" << endl;

    int success = handle_arguments(argc, argv);
    if(success == 0)
        return 0;

    cout << currentDateTime() << "\tRead JSON file from " << json_filepath_ << endl;
    json.read_json(json_filepath_);
    json.report_object_probes();

    cout << currentDateTime() << "\tBegin creating PoCos input files to " << out_dir_ << endl;
    write_phenotype_file();
    cout << currentDateTime() << "\t\tPhenotype file successfully written" << endl;
    write_snp_names();
    cout << currentDateTime() << "\t\tSNP names file successfully written" << endl;
    write_genotype_file();
    cout << currentDateTime() << "\t\tGenotype file successfully written" << endl;
}



int handle_arguments(int argc, char* argv[]){
    if(argc == 1 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h"){
        print_help();
        return 0;
    }
    bool json_set = false;
    bool out_set = false;

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

        if(att_tag == "--json_file"){
            json_filepath_ = att;
            json_set = true;
        }
        else if(att_tag == "--out_dir"){
            out_dir_ = att;
            out_set = true;
        }
        else{
            cout << "ATTENTION! Not supported attribute detected: " << att_tag << endl;
            print_help();
            return 0;
        }
    }

    if(!json_set){
        cout << "ATTENTION! Please set mandatory JSON file path" << endl << endl;
        print_help();
        return 0;
    }
    if(!out_set){
        cout << "ATTENTION! Please set mandatory output directory" << endl << endl;
        print_help();
        return 0;
    }
    return 1;
}

void print_help(){
    cout << "Usage message" << endl << "-------------" << endl << endl;
    cout << "Mandatory arguments:" << endl;
    cout << "--json_file\tPath to the JSON file." << endl;
    cout << "--out_dir\t\tDirectory where output should be written to." << endl << endl;
    cout << "Optional arguments:" << endl;
    cout << "--help, -h\tShow help message and exit." << endl << endl;
}



void write_genotype_file() {
    ofstream $out;
    $out.open(out_dir_ + "/" + json.json_filename_ + "_genotype.txt");
    for(int i = 0; i < json.genotype_.size(); i++){
        $out << json.get_genotype_string_snpwise(i, "\t") << "\n";
    }
    $out.close();
}

void write_phenotype_file() {
    ofstream $out;
    $out.open(out_dir_ + "/" + json.json_filename_ + "_phenotype.txt");
    $out << json.get_phenotype_string("\n", "1", "2");
    $out.close();
}

void write_snp_names() {
    ofstream $out;
    $out.open(out_dir_ + "/" + json.json_filename_ + "_snpNames.txt");
    $out << json.get_snp_string("\n");
    $out.close();
}


const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "[%Y/%m/%d-%X]", &tstruct);

    return buf;
}