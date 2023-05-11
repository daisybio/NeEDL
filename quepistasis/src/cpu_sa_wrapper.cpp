#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <cmath>
#include "../header/cpu_sa.h"
#include "../header/cpu_sa_wrapper.h"


std::vector<double> linear_space(double start, double end, int n) {
    std::vector<double> v(n);
    double step = (end - start) / (n - 1);
    for(int i = 0; i < n; i++) {
        v[i] = start + step * i;
    }
    return v;
}


std::vector<double> geometric_space(double start, double end, int n) {
    std::vector<double> geom_schedule = linear_space(std::log10(start), std::log10(end), n);
    for(int i = 0; i < n; i++) {
        geom_schedule[i] = std::pow(10, geom_schedule[i]);
    }
    geom_schedule[0] = start;
    geom_schedule[n-1] = end;
    return geom_schedule;
}


std::vector<double> calculate_beta_schedule(std::vector<double> h, std::vector<int> coupler_start, std::vector<int> coupler_end, std::vector<double> J, int num_betas, char mode)
{
    std::vector<double> h_abs(h.size());
    std::vector<double> J_abs(J.size());

    double (*d_abs)(double) = &std::abs;
    std::transform(h.begin(), h.end(), h_abs.begin(), d_abs);
    std::transform(J.begin(), J.end(), J_abs.begin(), d_abs);

    // calculate min energy step
    double min_h_abs = *std::min_element(h_abs.begin(), h_abs.end());
    double min_J_abs = *std::min_element(J_abs.begin(), J_abs.end());
    double min_delta_energy = std::min(min_h_abs, min_J_abs);

    // calculate max energy step
    for(size_t i = 0; i < coupler_start.size(); i++) {
        int k1 = coupler_start[i];
        int k2 = coupler_end[i];
        h_abs[k1] += J_abs[i];
        h_abs[k2] += J_abs[i];
    }
    double max_delta_energy = *std::max_element(h_abs.begin(), h_abs.end());

    // calculate beta range
    double hot_beta = std::log(2) / max_delta_energy;
    double cold_beta = std::log(100) / min_delta_energy;

    // if mode is linear then return a linear schedule
    if(mode == 'L' || mode == 'l') {
        return linear_space(hot_beta, cold_beta, num_betas);
    }
    
    // otherwise return geometric schedule
    return geometric_space(hot_beta, cold_beta, num_betas);
}


std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> convert_qubo_to_ising(double* Q, double Q_offset, int N) {
    
    double* h_vec = new double[N]{0};
    double* J_vec = new double[N*N]{0};
    double bias;
    double linear_offset = 0, quadratic_offset = 0;

    for(int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            bias = Q[N*i + j];
            if(i == j) {
                h_vec[i] += 0.5 * bias;
                linear_offset += bias;
            } else if(bias != 0) {
                J_vec[i*N+j] = 0.25 * bias;
                h_vec[i] += 0.25 * bias;
                h_vec[j] += 0.25 * bias;
                quadratic_offset += bias;
            }
        }
    }

    double offset = Q_offset + 0.5 * linear_offset + 0.25 * quadratic_offset;

    // fill linear vector
    std::vector<double> h(N);
    h.assign(h_vec, h_vec + N);

    // fill quadratic cost vector
    std::vector<int> coupler_start;
    std::vector<int> coupler_end; 
    std::vector<double> J;

    for(int i = 0; i < N*N; i++) {
        if(J_vec[i] != 0.0) {
            int a = i / N;
            int b = i % N;
            coupler_start.push_back(a);
            coupler_end.push_back(b);
            J.push_back(J_vec[i]);
        }
    }

    return std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double>(h, coupler_start, coupler_end, J, offset);
}


double get_energy_qubo(std::vector<char> solution, double* Q, int N) 
{
    double cost = 0.0;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            cost += solution[i] * Q[i*N+j] * solution[j];
        }
    }
    return cost;
}

/*
std::tuple<std::vector<char>, double> simulated_annealing_qubo(double* Q, double offset, int N, int num_samples, int num_sweeps, uint64_t seed) 
{
    // convert QUBO formulation to Ising formulation
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, int> t = convert_qubo_to_ising(Q, offset, N);
    
    std::vector<double> h = std::get<0>(t); 
    std::vector<int> coupler_start = std::get<1>(t);
    std::vector<int> coupler_end = std::get<2>(t);
    std::vector<double> J = std::get<3>(t);
    int offset = std::get<4>(t);

    std::cout << "H:";
    for(int i = 0; i < h.size(); i++) {
        std::cout << h[i] << " ";
    }
    std::cout << "\n";
    std::cout << "J:";
    for(int i = 0; i < coupler_start.size(); i++) {
        std::cout << "(" << coupler_start[i] << "," << coupler_end[i] << ") " << J[i];
        std::cout << "\n";
    }
    std::cout << "\n";
    std::cout << "offset: " << offset << "\n"; 

    // run simulated annealing
    sa_return ret = simulated_annealing_ising(h, coupler_start, coupler_end, J, num_samples, num_sweeps, seed);
    
    // print ising state vector
    std::cout << "ISING states:\n";
    for(int i = 0; i < num_samples; i++) {
        for(int j = 0; j < N; j++) {
            std::cout << (int)ret.states[i * N + j] << " ";
        }
        std::cout << std::endl;
    }

    // check the solution having better score, then transform it back to QUBO formulation
    int best_solution = 0;
    for(int i = 1; i < ret.num_samples; i++) {
        if(ret.energies[best_solution] > ret.energies[i]) {
            best_solution = i;
        }
    }
    std::cout << "Best solution is: " << best_solution << " (N=" << N << ")" << std::endl;

    // assign the best solution
    std::vector<char> solution(N);
    std::cout << "Copy indices from " << best_solution * N << " to " << ((best_solution+1) * N) << std::endl;
    for(int j = 0; j < N; j++) {
        solution[j] = ret.states[best_solution * N + j];
    }
    std::cout << "Resulting in: ";
    for(int i = 0; i < N; i++) {
        std::cout << (int)(solution[i]) << " ";
    }
    std::cout << std::endl;

    // convert states from Ising (+1, -1) to QUBO (0, 1)
    for(int i = 0; i < solution.size(); i++) {
        solution[i] = (solution[i] + 1) / 2;
    }

    // get energy for the solution
    double energy = get_energy_qubo(solution, Q, N);

    return std::tuple<std::vector<char>, double>(solution, energy);
}*/


sa_return::sa_return(int num_variables, int num_samples) {
    this->num_variables = num_variables;
    this->num_samples = num_samples;
    this->states = new char[num_samples * num_variables];
    this->energies = new double[num_samples];
}


sa_return::~sa_return() {
    delete[] states;
    delete[] energies;
}


bool interrupt_callback(void * const interrupt_function) {
    std::cerr << "Called interrupt callback" << std::endl;
    return 0;
}


void simulated_annealing_ising(
        std::vector<double> h, 
        std::vector<int> coupler_start, 
        std::vector<int> coupler_end,
        std::vector<double> J, 
        int num_samples,
        int num_sweeps,
        uint64_t seed,
        sa_return* ret) 
{
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<uint8_t> dis;

    // Calculate number of variables
    const int num_variables = h.size();

    // Setup schedule of simulated annealing
    int num_sweeps_per_beta = (num_sweeps / 1000) | 1;
    int num_betas = num_sweeps / num_sweeps_per_beta + 1;

    // Calculate beta schedule according to the given formulation
    std::vector<double> beta_schedule = calculate_beta_schedule(h, coupler_start, coupler_end, J, num_betas, 'G');

    // Call simulated annealing, return number of actual samples
    // The function fills the vectors in ret.states and ret.energies
    int actual_n_samples = general_simulated_annealing(
        ret->states, ret->energies, num_samples, h, coupler_start, coupler_end, J,
        num_sweeps_per_beta, beta_schedule, seed,
        interrupt_callback, NULL);

    // overwrite num_samples with the actual value
    ret->num_samples = actual_n_samples;
}