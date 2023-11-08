//
// Created by juli on 08.11.23.
//

#include "PlinkShufflePhenotype.hpp"
#include "../util/TimeLogger.hpp"

namespace epi {

PlinkShufflePhenotype::PlinkShufflePhenotype(std::string input_path, std::string output_path) {
    this->input_path = input_path;
    this->output_path = output_path;
}


void PlinkShufflePhenotype::run(std::shared_ptr<DataModel> data) {
    TimeLogger logger("shufffle phenotype");

    Logger::logLine("Read fam file");
    CSVParser ind_parser;
    ind_parser.parse(input_path + ".fam", ' ');
    if (ind_parser.num_columns() < 6) ind_parser.parse(input_path + ".fam", '\t');

    // create a list of indices and shuffle them -> output file will get the phenotypes in that order
    std::vector<size_t> index_list;
    index_list.reserve(ind_parser.num_rows());
    for(size_t i = 0; i < ind_parser.num_rows();++i) {
        index_list.push_back(i);
    }
    std::shuffle(index_list.begin(), index_list.end(), data->random_device[omp_get_thread_num()]);
    

    std::ofstream pheno_file(output_path + ".fam");

    for (size_t i = 0; i < ind_parser.num_rows(); ++i) {
        for (size_t col = 0; col < 5; ++col) {
            pheno_file << ind_parser.cell(i, col) << '\t';
        }
        pheno_file << ind_parser.cell(index_list[i], 5) << '\n';
    }
    pheno_file.close();

    // copy bim and bed file
    std::filesystem::copy_file(input_path + ".bed", output_path + ".bed");
    std::filesystem::copy_file(input_path + ".bim", output_path + ".bim");

    
    logger.stop();
}
}
