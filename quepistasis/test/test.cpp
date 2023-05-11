#include <iostream>
#include <vector>
#include <iomanip>
#include <random>
#include <tuple>
#include <typeinfo>
#include "snps_optimization.h"
#include "python_wrapper.h"

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}


int main(int argc, char** argv) {

    signal(SIGSEGV, handler);   // install our handler

    int N_SNPS = 4;
    double* stat_info = new double[N_SNPS*N_SNPS]{
        0, 1, 2, 3,
        0, 0, 4, 5,
        0, 0, 0, 6,
        0, 0, 0, 0
    };
    double* bio_info = new double[N_SNPS*N_SNPS]{
        0, 1, 2, 3,
        0, 0, 4, 5,
        0, 0, 0, 6,
        0, 0, 0, 0
    };
    matrix stat_matrix(stat_info, N_SNPS);
    matrix bio_matrix(bio_info, N_SNPS);

    int N_CLIQUE = 2;
    int K = 3;

    snps_qubo_matrix qubo(N_CLIQUE, N_SNPS);

    std::cout << "QUBO generation:\n";
    qubo.fill_n_max_weighted_k_clique_qubo(stat_matrix, bio_matrix, K, 0.5, 10.0, 0.456, 0.789);

    std::cout << "\n\nPrettyprint QUBO:\n";
    qubo.prettyprint_qubo();

    std::cout << "\n\nPrettyprint Ising:\n";
    qubo.prettyprint_ising();


    // with Simulated Annealing on CPU
    std::cout << "\n\n=================================================\n";
    std::cout << "=============== SIM ANNEALING ===================\n";
    std::cout << "=================================================\n";
    std::cout << "Solve with seed 1234 and random initial state\n";
    uint64_t seed = 1234;
    char* initial_state = nullptr;
    const char* logFilePathSimAnn = "./simannlog";
    std::vector<std::vector<int>> sa_snps_set_list = qubo.solve_cpu_simulated_annealing(10, 1000, seed, initial_state, logFilePathSimAnn);
    std::cout << "\n\nPrint result\n";
    for(std::vector<int> snps_set : sa_snps_set_list) {
        std::cout << "SNPs set: ";
        for(int snp : snps_set) {
            std::cout << snp << " ";
        }
        std::cout << std::endl;
    }


    // with quantum computer QAOA on AZURE
    std::cout << "\n\n=================================================\n";
    std::cout << "================= AZURE QC ======================\n";
    std::cout << "=================================================\n";
    const char* logFilePathQC = "./qclog";
    std::vector<std::vector<int>> qaoa_snps_set_list = qubo.solve_azure_QAOA(
        "none", "", "", "", "", "",
        "ADAM", 1, 1, 1, 0, logFilePathQC);
    std::cout << "\n\nPrint result\n";
    if(qaoa_snps_set_list.empty()) {
        std::cout << "Empty SNP set" << std::endl;
    }
    for(std::vector<int> snps_set : qaoa_snps_set_list) {
        std::cout << "SNPs set: ";
        for(int snp : snps_set) {
            std::cout << snp << " ";
        }
        std::cout << std::endl;
    }


    // with several optimizers on Microsoft Azure QIO platform
    // std::cout << "\n\n=================================================\n";
    // std::cout << "================= AZURE QIO ====================\n";
    // std::cout << "=================================================\n";
    // const char* logFilePathAzureQIO = "./azurelog";
    // std::vector<std::vector<int>> qio_snps_set_list = qubo.solve_azure_optimizer(
    //     "xxxxxxx-TOKEN", "azurequantum", "quantumincudstudent", "westeurope", 1, 12345, logFilePathAzureQIO);
    // std::cout << "\n\nPrint result\n";
    // for(std::vector<int> snps_set : qio_snps_set_list) {
    //     std::cout << "SNPs set: ";
    //     for(int snp : snps_set) {
    //         std::cout << snp << " ";
    //     }
    //     std::cout << std::endl;
    // }


    // with quantum annealing
    // std::cout << "\n\n=================================================\n";
    // std::cout << "================= DWAVE QA ======================\n";
    // std::cout << "=================================================\n";
    // const char* logFilePathQA = "./dwavelog";
    // std::vector<std::vector<int>> qa_snps_set_list = qubo.solve_dwave_quantum_annealing(
    //     "CINE-8359d556e22e5c9132042bc427bb0d0063bc1199", 100, 4, 10, 50, 10, 50, 0.5, logFilePathQA);
    // std::cout << "\n\nPrint result\n";
    // for(std::vector<int> snps_set : qa_snps_set_list) {
    //     std::cout << "SNPs set: ";
    //     for(int snp : snps_set) {
    //         std::cout << snp << " ";
    //     }
    //     std::cout << std::endl;
    // }


    // with ATOS simulated quantum annealing
    // std::cout << "\n\n=================================================\n";
    // std::cout << "================= ATOS SQA ======================\n";
    // std::cout << "=================================================\n";
    // const char* logFilePathAtos = "./atoslog";
    // std::vector<std::vector<int>> atos_snps_set_list = qubo.solve_atos_simulated_quantum_annealing(5, 1024, logFilePathAtos);
    // std::cout << "\n\nPrint result\n";
    // if(atos_snps_set_list.empty()) {
    //     std::cout << "Empty SNP set" << std::endl;
    // }
    // for(std::vector<int> snps_set : atos_snps_set_list) {
    //     std::cout << "SNPs set: ";
    //     for(int snp : snps_set) {
    //         std::cout << snp << " ";
    //     }
    //     std::cout << std::endl;
    // }

    std::cout << "\n\n=================================================\n";
    std::cout << "================= CLOSING =======================\n";
    std::cout << "=================================================\n";
    std::cout << "Before closing Python environment...\n" << std::flush;
    PythonWrapper::get_instance().close();

    return 0;
}
