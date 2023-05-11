// JSON_format.cpp
// Created by eve on 17.05.21.
//

#ifndef GENEPISEEKER_JSON_FORMAT_IPP
#define GENEPISEEKER_JSON_FORMAT_IPP

#import "JSON_format.hpp"



//########################################################################
//                              Constructors
//########################################################################

JSON_format::JSON_format() {
    json_path_ = "";
    json_filename_ = "";
    num_snps_ = 0;
    num_inds_ = 0;
    model_type_ = "";
    num_categories_ = 0;
}

JSON_format::JSON_format(string json_filepath) {
    json_path_ = "";
    json_filename_ = "";
    num_snps_ = 0;
    num_inds_ = 0;
    model_type_ = "";
    num_categories_ = 0;
    read_json(json_filepath);
}



//########################################################################
//                      Get, set and add methods
//########################################################################

void JSON_format::add_phenotypes(vector<string> phenotype) {
    phenotype_.insert(phenotype_.end(), phenotype.begin(), phenotype.end());
//    for(string pheno : phenotype)
//        phenotype_.push_back(pheno);
}

void JSON_format::add_json_entry(int index, vector<string> genotype, vector<string> snp, bool disease_snp) {
    if(genotype_.size() == index){
        vector<string> temp;
        genotype_.push_back(temp);
    }
    for(string geno : genotype)
        genotype_[index].push_back(geno);
//    genotype_[index].insert(genotype_[index].end(), genotype.begin(), genotype.end());
    if(snps_.size() == index)
        snps_.push_back(snp);
    else if(snps_[index] != snp){
        if(snps_[index][0] != snp[0])
            cout << "WARNING: mismatching SNP-IDs " << snps_[index][0] << "\t" << snp[0] << endl;
        else if(snps_[index][1] != snp[1])
            cout << "WARNING: mismatching chromosome annotation " << snps_[index][1] << "\t" << snp[1] << endl;
        else if(snps_[index][2] != snp[2])
            cout << "WARNING: mismatching region annotation " << snps_[index][2] << "\t" << snp[2] << endl;
        else if(snps_[index][3] != snp[3] && snps_[index][4] == snp[4]){
            if(snps_[index][3] == "0")
                snps_[index] = snp;
        }
        else if(snps_[index][4] != snp[4]){
            if(snps_[index][4] == "0")
                snps_[index] = snp;
        }
        else{
            cout << "WARNING: mismatching SNP annotation!" << endl;
            cout << "\tObject SNP annot: " << snps_[index][0] << "\t" << snps_[index][1] << "\t" << snps_[index][2] << "\t" << snps_[index][3] << "\t" << snps_[index][4] << endl;
            cout << "\tTo add SNP annot: " << snp[0] << "\t" << snp[1] << "\t" << snp[2] << "\t" << snp[3] << "\t" << snp[4] << endl;
        }
    }
    if(disease_snp)
        disease_snps_.push_back(index);
}

unordered_set<string> JSON_format::get_snp_list() {
    unordered_set<string> snp_ids;
    for(int i = 0; i < snps_.size(); i++)
        snp_ids.insert(snps_[i][0]);
    return snp_ids;
}

string JSON_format::get_snp_string(string sep) {
    string out = snps_[0][0];
    for(int i = 1; i < snps_.size(); i++){
        out += sep + snps_[i][0];
    }
    return out;
}

string JSON_format::get_genotype_string_samplewise(int index, string sep) {
    string out = genotype_[0][index];
    for(int j = 1; j < genotype_.size(); j++){
        out += sep + genotype_[j][index];
    }
    return out;
}

string JSON_format::get_genotype_string_snpwise(int index, string sep) {
    string out = to_string(stoi(genotype_[index][0]));
    for(int j = 1; j < genotype_[index].size(); j++){
        out += sep + to_string(stoi(genotype_[index][j]));
    }
    return out;
}

string JSON_format::get_phenotype_string(string sep, string case_code = "1", string control_code = "0") {
    string out;
    if(phenotype_[0] == "1")
        out += case_code;
    else if(phenotype_[0] == "0")
        out += control_code;
    else cout << "WARNING: Ambiguous phenotype found for index 0 -> " << phenotype_[0] << endl;

    for(int i = 1; i < phenotype_.size(); i++){
        if(phenotype_[i] == "1")
            out += sep + case_code;
        else if(phenotype_[i] == "0")
            out += sep + control_code;
        else cout << "WARNING: Ambiguous phenotype found for index " << i << " -> " << phenotype_[i] << endl;
    }
    return out;
}



//########################################################################
//                    Methods to read JSON file
//########################################################################

void JSON_format::read_json(string json_filepath) {
    json_path_ = json_filepath;

    vector<string> tokens;
    boost::split(tokens, json_path_, boost::is_any_of("/"));
    json_filename_ = tokens[tokens.size()-1];
    boost::replace_all(json_filename_, ".json", "");

    string line;
    string json;
    ifstream $infile(json_path_);

    if($infile.is_open()){
        while(getline($infile, line)){
            json += line;
        }
        $infile.close();
        $infile.clear();
    }
    else throw runtime_error("Unable to read JSON file " + json_path_);

    assign_tokens(tokenize(json));
}


// -----------------------------------------------------------------------
// Util methods to read JSON file
// -----------------------------------------------------------------------

vector<string> JSON_format::tokenize(string line) {
    vector<string> tokens = split_at_string(line, ": [[");
    line = "";
    vector<string> cleaned;
    vector<string> temp_tokens;

    //process first portion of json string
    boost::split(temp_tokens, tokens[0], boost::is_any_of("\""));
    for(string token : temp_tokens){
        if(token == "{" | token == ": " | token == ", " | token == ""){
            continue;
        }
        cleaned.push_back(clean_token(token, 0));
    }

    //process second portion of json string containing matrix structure
    boost::split(temp_tokens, tokens[1], boost::is_any_of(("\"")));
    for(string token : temp_tokens){
        if(token == ""){
            continue;
        }
        cleaned.push_back(clean_token(token, 1));
    }

    //process last portion of json string containing snp annotation
    temp_tokens = split_at_string(tokens[2], "]], ");
    cleaned.push_back(temp_tokens[0]);
    boost::split(temp_tokens, temp_tokens[1], boost::is_any_of("\""));
    for(string token : temp_tokens){
        if(token == ""){
            continue;
        }
        cleaned.push_back(clean_token(token, 2));
    }

    return cleaned;
}

string JSON_format::clean_token(string token, int mode) {
    if(mode == 0){
        boost::replace_all(token, ",", "");
    }

    else if(mode == 1){
        boost::replace_all(token, ": [", "");
        boost::replace_all(token, "], ", "");
    }

    else if(mode == 2){
        boost::replace_all(token, "}", "");
        boost::replace_all(token, "[", "");
        boost::replace_all(token, "],", "");
        boost::replace_all(token, "]", "");
    }

    boost::replace_all(token, ":", "");
    boost::replace_all(token, " ", "");

    return token;
}

vector<string> JSON_format::split_at_string(string s, string delimiter) {
    vector<string> tokens;
    boost::replace_all(s, delimiter, "%");
    boost::split(tokens, s, boost::is_any_of("%"));
    return tokens;
}

void JSON_format::assign_tokens(vector<string> tokens) {
    for(int i = 0; i < tokens.size(); i = i+2){
        if(tokens[i] == "num_snps"){
            num_snps_ = stoi(tokens[i+1]);
        }

        else if(tokens[i] == "num_inds"){
            num_inds_ = stoi(tokens[i+1]);
        }

        else if(tokens[i] == "model_type"){
            model_type_ = tokens[i+1];
        }

        else if(tokens[i] == "num_categories"){
            num_categories_ = stoi(tokens[i+1]);
        }

        else if(tokens[i] == "genotype"){
            vector<string> snp_genos;
            boost::split(snp_genos, tokens[i+1], boost::is_any_of("[["));
            for(string snp_geno : snp_genos){
                if(snp_geno == ": " || snp_geno == ""){
                    continue;
                }

                vector<string> values;
                string snp_geno_cleaned = clean_token(snp_geno, 2);
                boost::split(values, snp_geno_cleaned, boost::is_any_of(","));
                genotype_.push_back(values);
            }
        }

        else if(tokens[i] == "phenotype"){
            vector<string> values;
            boost::split(values, tokens[i+1], boost::is_any_of(","));
            for(string vals : values){
                if(vals == ""){
                    continue;
                }
                phenotype_.push_back(vals);
            }
        }

        else if(tokens[i] == "snps"){
            vector<string> snp_annots = split_at_string(tokens[i+1], "], [");
            for(string snp_annot : snp_annots){
                vector<string> values;
                boost::split(values, snp_annot, boost::is_any_of("\""));
                vector<string> annots;

                for(string val : values){
                    if(val == "" | val == ", "){
                        continue;
                    }

                    //if(val.size() > 3 && val.substr(0, 3) == "chr"){
                    //    val = val.substr(3, val.size());
                    //}
                    annots.push_back(val);
                }
                snps_.push_back(annots);
            }
        }

        else if(tokens[i] == "disease_snps"){
            vector<string> snp_dis;
            boost::split(snp_dis, tokens[i+1], boost::is_any_of(","));
            for(string dis : snp_dis){
                if(dis == "")
                    continue;
                disease_snps_.push_back(stoi(dis));
            }
        }

        else if(tokens[i] == "mafs"){
            vector<string> maf;
            boost::split(maf, tokens[i+1], boost::is_any_of(","));
            for(string m : maf){
                mafs_.push_back(stod(m));
            }
        }

        else cout << "\tWARNING\tjson field " << tokens[i] << " not included in this reader. Value excluded: " << tokens[i+1] << endl;
    }
}



//########################################################################
//             Methods to write JSON object to output file
//########################################################################

void JSON_format::write_json(string out_file) {
    string out_filepath = out_file;
    if(out_file.substr(out_file.size()-5, out_file.size()) != ".json")
        out_filepath += "_" + to_string(num_snps_) + ".json";

    ofstream $out(out_filepath);
    if($out.is_open()){
        $out << "{" <<
             "\"num_snps\": " << num_snps_ << ", " <<
             "\"num_inds\": " << num_inds_ << ", " <<
             "\"model_type\": \"" << model_type_ << "\", " <<
             "\"num_categories\": " << num_categories_ << ", " <<
             "\"genotype\": [" << genotype_to_json_string(0);
        for(int i = 1; i < genotype_.size(); i++)
            $out << "," << genotype_to_json_string(i);
        $out << "], " <<
             "\"phenotype\": " << phenotypes_to_json_string() << ", " <<
             "\"snps\": [" << snp_to_json_string(0);
        for(int i = 1; i < snps_.size(); i++)
            $out << ", " << snp_to_json_string(i);
        $out << "], " <<
             "\"disease_snps\": " << disease_snps_to_json_string() << ", " <<
             "\"mafs\": " << mafs_to_json_string() << "}";

        $out.close();
        $out.clear();
    }
    else throw runtime_error("Unable to write to JSON file " + out_filepath);
}


// -----------------------------------------------------------------------
// Util methods to write JSON object
// -----------------------------------------------------------------------

string JSON_format::genotype_to_json_string(int index) {
    vector<string> snp_genotype = genotype_[index];
    if (snp_genotype.size() == 0)
        return "[]";

    string geno = "[" + snp_genotype[0];
    for (int i = 1; i < snp_genotype.size(); i++)
        geno += "," + snp_genotype[i];
    return geno + "]";
}

string JSON_format::phenotypes_to_json_string() {
    if (phenotype_.size() == 0)
        return "[]";

    string pheno = "[" + phenotype_[0];
    for(int i = 1; i < phenotype_.size(); i++)
        pheno += ", " + phenotype_[i];
    return pheno + "]";
}

string JSON_format::snp_to_json_string(int index) {
    vector<string> snp = snps_[index];
    if (snp.size() == 0)
        return "[]";

    string s = "[\"" + snp[0] + "\"";
    for(int i = 1; i < snp.size(); i++)
        s += ", \"" + snp[i] + "\"";
    return s + "]";
}

string JSON_format::disease_snps_to_json_string() {
    if (disease_snps_.size() == 0)
        return "[]";

    string disease = "[" + to_string(disease_snps_[0]);
    for(int i = 1; i < disease_snps_.size(); i++)
        disease += ", " + to_string(disease_snps_[i]);
    return disease + "]";
}

string JSON_format::mafs_to_json_string() {
    if (mafs_.size() == 0)
        return "[]";

    string maf = "[" + to_string(static_cast<int>(mafs_[0] * 10000.0)/1000.0);
    for(int i = 1; i < mafs_.size(); i++)
        maf += ", " + to_string(static_cast<int>(mafs_[i] * 10000.0)/1000.0);
    return maf + "]";
}



//########################################################################
//                           Util methods
//########################################################################

int JSON_format::find(string snp) {
    for(int i = 0; i < snps_.size(); i++){
        if(snps_[i][0] == snp)
            return i;
    }
    return -1;
}

bool JSON_format::is_disease_snp(string snp) {
    if(disease_snps_.size() == 0)
        return false;
    for(int i : disease_snps_)
        if(snps_[i-1][0] == snp)
            return true;
    return false;
}

void JSON_format::calculate_mafs() {
    for(int i = 0; i < genotype_.size(); i++){
        int geno_sum = 0;
        for(int j = 0; j < genotype_[i].size(); j++){
            string placebo = genotype_[i][j];
            if(placebo == "NA"){
                cout << "OLLA! NA in sight, captain" << endl;
                continue;
            }
            geno_sum += stoi(genotype_[i][j]);
        }
        mafs_.push_back((geno_sum * 1.0)/(2.0 * num_inds_));
    }
}

void JSON_format::report_object_dimension() {
    cout << "\n************************************************************" << endl;
    cout << "                 Object dimension report                  " << endl;
    cout << "************************************************************\n" << endl;
    cout << "Welcome to the JSON_format object dimension report of: " << endl;
    cout << json_path_ << endl << endl;

    cout << " Annotated dimensions" << endl;
    cout << " - Number of SNPs: " << num_snps_ << endl;
    cout << " - Number of individuals: " << num_inds_ << endl << endl;
    cout << " Dimensions of genotype matrix" << endl;
    cout << " - Rows: " << genotype_.size() << endl;
    cout << " - Columns: " << genotype_[0].size() << endl << endl;
    cout << " Dimension of object lists" << endl;
    cout << " - Number of phenotypes: " << phenotype_.size() << endl;
    cout << " - Number of SNPs: " << snps_.size() << endl;
    cout << " - Number of disease SNPs: " << disease_snps_.size() << endl;
    cout << " - Number of mafs: " << mafs_.size() << endl << endl;
}

void JSON_format::report_object_probes() {
    cout << "Genotype probe:\t" << genotype_[0][0] << "\t" << genotype_[0][1] << "\t" << genotype_[1][0] << "\t" << genotype_[1][1] << endl;
    cout << "Phenotype probe:\t" << phenotype_[0] << "\t" << phenotype_[1] << "\t" << phenotype_[2] << "\t" << phenotype_[3] << endl;
    cout << "SNP annotation probe:" << endl;
    cout << "\t\t\t- " << snp_to_json_string(0) << endl;
    cout << "\t\t\t- " << snp_to_json_string(1) << endl;
    cout << "\t\t\t- " << snp_to_json_string(2) << endl;
    cout << "\t\t\t- " << snp_to_json_string(3) << endl << endl;
}

void JSON_format::write_snp_annotation(string out_dir) {
    string annot_filepath = out_dir + "/" + json_filename_ + "_SNPannot.txt";
    ofstream $out(annot_filepath);
    if($out.is_open()){
        $out << "snp_id\tchromosome\tgenomic_position\twildtype\tvariant" << endl;
        for(int i = 0; i < snps_.size(); i++){
            $out << snps_[i][0];
            for(int j = 1; j < snps_[0].size(); j ++)
                $out << "\t" << snps_[i][j];
            $out << endl;
        }
        $out.close();
        $out.clear();
    }
    else throw runtime_error("Unable to write to SNP annotation file " + annot_filepath);
}


#endif //GENEPISEEKER_JSON_FORMAT_IPP