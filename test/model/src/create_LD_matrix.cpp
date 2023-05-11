// create_LD_matrix.cpp
// Created by eve on 04.12.2020.
//

#include <iostream>
#include <string>
#include <stdlib.h>
#include <chrono>
#include "JSON_format.hpp"

using namespace std;



class LD_creator{

private:
    string working_dir_ = "/home/eve/GenEpiSeeker/NeEDL";
    string path_to_plink_ = working_dir_ + "/ext/plink/plink_linux_x86_64_20201019/plink";

    string infile_filepath_;
    string output_;
    string ld_filepath_;

    JSON_format json = JSON_format();
    string ped_filepath_;
    string map_filepath_;
    vector<string> files_to_remove_;


    void remove_files(bool remove);
    const string currentDateTime();

    void write_ped_and_map_files();
    void recode_bfile_to_file();

    void run_plink();


public:
    LD_creator();
    LD_creator(string infile, string out_dir, string working_dir);

    string create_LD_matrix_from_bfile(bool remove_redundant_files);
    string create_LD_matrix_from_json(bool remove_redundant_files);
};



//########################################################################
//                              Constructors
//########################################################################

LD_creator::LD_creator() {
    cout << "\n************************************************************" << endl;
    cout << "                      Create LD matrix                      " << endl;
    cout << "************************************************************\n" << endl;
}

LD_creator::LD_creator(string infile, string output, string working_dir) {
    cout << "\n************************************************************" << endl;
    cout << "                      Create LD matrix                      " << endl;
    cout << "************************************************************\n" << endl;

    if(infile != "")
        infile_filepath_ = infile;
    else throw runtime_error("WARNING: filepath to input file not set.");

    if(output != "")
        output_ = output;
    else throw runtime_error("WARNING: path to output directory not set.");

    if(working_dir != "")
        working_dir_ = working_dir;
}



//########################################################################
//
//                      Main functionality methode
//
//########################################################################

string LD_creator::create_LD_matrix_from_bfile(bool remove_redundant_files) {

    // -----------------------------------------------------------------------
    // Write .ped and .map files to use with PLINK for linkage equilibrium
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tRecode BFILEs to .ped and .map files" << endl;
    recode_bfile_to_file();


    // -----------------------------------------------------------------------
    // Run PLINK and manage output files in run_plink()
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tRun PLINK and create .ld LD matrix file" << endl;
    run_plink();


    // -----------------------------------------------------------------------
    // Cleaning and wrap up
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tFinished" << endl;
    remove_files(remove_redundant_files);

    return output_;
}

string LD_creator::create_LD_matrix_from_json(bool remove_redundant_files) {

    // -----------------------------------------------------------------------
    // Read JSON
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tRead JSON file from " << infile_filepath_ << endl;
    json.read_json(infile_filepath_);


    // -----------------------------------------------------------------------
    // Write .ped and .map files to use with PLINK for linkage equilibrium
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tBegin creating .ped and .map file to " << output_ << endl;
    write_ped_and_map_files();


    // -----------------------------------------------------------------------
    // Run PLINK and manage output files in run_plink()
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tRun PLINK and create .ld LD matrix file" << endl;
    run_plink();


    // -----------------------------------------------------------------------
    // Cleaning and wrap up
    // -----------------------------------------------------------------------
    cout << currentDateTime() << "\tFinished" << endl;
    remove_files(remove_redundant_files);

    return output_ + "/" + json.json_filename_;
}



//########################################################################
//                     Basic administration functions
//########################################################################

void LD_creator::remove_files(bool remove) {
    files_to_remove_.push_back(ped_filepath_);
    files_to_remove_.push_back(map_filepath_);
    if(remove){
        for(string file : files_to_remove_)
            system(("rm " + file).c_str());
    }
    else{
        cout << endl << "Files to remove:" << endl;
        for(string file : files_to_remove_)
            cout << "\t- " << file << endl;
    }
}

const string LD_creator::currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "[%Y/%m/%d-%X]", &tstruct);

    return buf;
}



//########################################################################
//                  Convert JSON to PLINK input files
//########################################################################

void LD_creator::write_ped_and_map_files() {
    ped_filepath_ = output_ + "/" + json.json_filename_ + "_ld_creator.ped";
    map_filepath_ = output_ + "/" + json.json_filename_ + "_ld_creator.map";
    files_to_remove_.push_back(ped_filepath_);
    files_to_remove_.push_back(map_filepath_);

    // -----------------------------------------------------------------------
    // Write .ped file to output directory
    // -----------------------------------------------------------------------
    cout << "\t\t" << currentDateTime() << "\tWrite " << ped_filepath_ << endl;
    ofstream $out_ped(ped_filepath_);
    if($out_ped.is_open()){
        for(int i = 0; i < json.num_inds_; i++){
            // tab or space delimited
            // ped columns: [0]family_id    [1]individual_id    [2]paternal_id  [3]maternal_id  [4]sex  [5] phenotype [6-n]genotypes; haploid with space delim
            string out = "FAM" + to_string(i+1) + " " + to_string(i+1) + " 0 0 0  " + to_string(stoi(json.phenotype_[i])+1);

            for(int j = 0; j < json.num_snps_; j++) {
                if (json.genotype_[j][i] == "0") {
                    out = out + "  " + json.snps_[j][4] + " " + json.snps_[j][4];
                }
                else if (json.genotype_[j][i] == "1") {
                    out = out + "  " + json.snps_[j][3] + " " + json.snps_[j][4];
                }
                else if (json.genotype_[j][i] == "2") {
                    out = out + "  " + json.snps_[j][3] + " " + json.snps_[j][3];
                }
                else cout << currentDateTime() << "\tWARNING\tGenotype detected which is not in [0, 1, 2]: " + json.genotype_[j][i] << endl;
            }

            $out_ped << out << endl;
        }
        $out_ped.close();
        $out_ped.clear();
    }
    else throw runtime_error("Unable to write to PED file " + ped_filepath_);


    // -----------------------------------------------------------------------
    // Write .map file to output directory
    // -----------------------------------------------------------------------
    cout << "\t\t" << currentDateTime() << "\tWrite " << map_filepath_ << endl;
    ofstream $out_map(map_filepath_);
    if($out_map.is_open()){
        // annot columns: [0]rs_id [1]chromosome [2]basepair_position  [3]wildtype  [4]variant
        for(int i = 0; i < json.num_snps_; i++){
            // map columns: chromosome  rs_id   genetic_distance    basepair_position
            $out_map<< json.snps_[i][1] + " " + json.snps_[i][0] + " 0 " + json.snps_[i][2] << endl;
        }
        $out_map.close();
        $out_map.clear();
    }
    else throw runtime_error("Unable to write to MAP file " + map_filepath_);
    ld_filepath_ = output_ + "/" + json.json_filename_;
}

void LD_creator::recode_bfile_to_file() {
    system((path_to_plink_ +
            " --noweb" +
            " --bfile " + infile_filepath_ +
            " --recode" +
            " --out " + output_).c_str());
    ped_filepath_ = output_ + ".ped";
    map_filepath_ = output_ + ".map";
    ld_filepath_ = output_;
}



//########################################################################
//                          Create LD matrix
//########################################################################

void LD_creator::run_plink() {
    system((path_to_plink_ +
            " --ped " + ped_filepath_ +
            " --map " + map_filepath_ +
            " --out " + ld_filepath_ +
            " --allow-no-sex --noweb --r2 --matrix").c_str());

    files_to_remove_.push_back(ld_filepath_ + ".log");
    files_to_remove_.push_back(ld_filepath_ + ".nosex");
}