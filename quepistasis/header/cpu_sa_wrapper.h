#ifndef CPU_SA_WRAPPER_H__
#define CPU_SA_WRAPPER_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <tuple>
#include <cmath>

double get_energy_qubo(std::vector<char> solution, double* Q, int N);

std::vector<double> linear_space(double start, double end, int n);

std::vector<double> geometric_space(double start, double end, int n);

std::vector<double> calculate_beta_schedule(std::vector<double> h, std::vector<int> coupler_start, std::vector<int> coupler_end, std::vector<double> J, int num_betas, char mode);

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> convert_qubo_to_ising(double* Q, double offset, int N);

struct sa_return {
public:
        int num_variables;
        int num_samples;
        char* states;
        double* energies;
        sa_return(int num_variables, int num_samples);
        ~sa_return();
};

void simulated_annealing_ising(
        std::vector<double> h, 
        std::vector<int> coupler_start, 
        std::vector<int> coupler_end,
        std::vector<double> J, 
        int num_samples,
        int num_sweeps,
        uint64_t seed,
        sa_return* ret);

#endif