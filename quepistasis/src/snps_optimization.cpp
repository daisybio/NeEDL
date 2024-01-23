#include <tuple>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <ctime>
#include "../header/python_wrapper.h"
#include "../header/cpu_sa_wrapper.h"
#include "../header/snps_optimization.h"


matrix::matrix(double* m, int n) {
    this->n = n;
    this->mat = m;
}


matrix::matrix(int n) {
    this->n = n;
    this->mat = new double[n*n]{0};
}


matrix::~matrix() {
    delete[] this->mat;
}


double& matrix::get(int i, int j) {
    return mat[i * n + j];
}


double& matrix::operator()(int i, int j) {
    return mat[i * n + j];
}


double matrix::get_const(int i, int j) const {
    return mat[i * n + j];
}


double matrix::operator()(int i, int j) const {
    return get_const(i, j);
}

int matrix::get_n() const {
    return n;
}


double* matrix::get_matrix() {
    return mat;
}

matrix linear_combination(double nu, const matrix& a, const matrix& b){

    int N = a.get_n();
    matrix c(N);

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            c(i, j) = nu * a(i, j) + (1 - nu) * b(i, j);
        }
    }

    return c;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vals)
{
    for (const int& val : vals) 
        os << val << " ";
    return os;
}

template <>
std::ostream& operator<<(std::ostream& os, const std::vector<double>& vals)
{
    for (const double& val : vals) 
        os << std::fixed << std::setprecision(3) << val << " ";
    return os;
}


snps_qubo_matrix::snps_qubo_matrix(int n_cliques, int n_snps) : matrix(n_cliques * n_snps) {
    this->n_cliques = n_cliques;
    this->n_snps = n_snps;
    this->offset = 0;
}


snps_qubo_matrix::~snps_qubo_matrix() { }


double& snps_qubo_matrix::access(int first_clique, int first_variable, int second_clique, int second_variable) {
    return get(first_clique * n_snps + first_variable, second_clique * n_snps + second_variable);
}


double snps_qubo_matrix::access(int first_clique, int first_variable, int second_clique, int second_variable) const {
    return get_const(first_clique * n_snps + first_variable, second_clique * n_snps + second_variable);
}

void snps_qubo_matrix::fill_n_max_weighted_k_clique_qubo(const matrix& stat_corr, const matrix& bio_corr, int K, double nu, double lambda0, double lambda1, double lambda2) {

    matrix ss_corr = linear_combination(nu, stat_corr, bio_corr);
    this->fill_n_max_weighted_k_clique_qubo(ss_corr, K, lambda0, lambda1, lambda2);
}

void snps_qubo_matrix::fill_n_max_weighted_k_clique_qubo(const matrix& ss_corr, int K, double lambda0, double lambda1, double lambda2) {

    const int N = ss_corr.get_n();
    // Hide print due to big matrices
    // for(int i = 0; i < N; i++) {
    //     for(int j = 0; j < N; j++) {
    //         std::cout << ss_corr(i, j) << " ";
    //     }
    //     std::cout << "\n";
    // }

    int N_CLIQUES = this->n_cliques;
    int N_SNPS = ss_corr.get_n();
    if (N_SNPS != this->n_snps) {
        std::ostringstream os;
        os << "The correlation matrix has " << N_SNPS << " SNPs while you have previously declared " << this->n_snps << " SNPs";
        throw std::invalid_argument(os.str());
    }

    this->offset = 0;

    // For each of the N cliques to be found add the penality terms enforcing the size of 
    // the clique (to be close to 'k') and prefer the cliques with larger weights
    for(int l = 0; l < N_CLIQUES; l++) {

        // CONSTRAINT: enforcing K size of snps set
        // lambda_0 * (x_1 + x_2 + ... - K)^2
        // = lambda_0 * (x_1^2 + x_2^ + ... + (-K)^2 + 2 x_1 x_2 + 2 x_1 ... - 2 * K * x_1 - 2 * K * x_2 - ...)
        // = x_i^2 -> lambda_0 (1 - 2 * K)
        //   x_i x_j -> lambda_0 * (2)
        for(int i = 0; i < N_SNPS; i++) {
            this->access(l, i, l, i) += lambda0 * (1 - 2 * K);
            for(int j = i+1; j < N_SNPS; j++) {
                this->access(l, i, l, j) += 2 * lambda0;
            }
        }
        this->offset += lambda0 * K * K;

        // REWARD TERM: reward snps set with higher weights
        // - lambda_1 \sum_{i,j} w_{i,j} * x_i * x_j
        // = x_i x_j -> - lambda_1 * w_{i,j}
        //   x_i -> no correlation between the same snp
        for(int i = 0; i < N_SNPS; i++) {
            for(int j = i + 1; j < N_SNPS; j++) {
                this->access(l, i, l, j) += - lambda1 * ss_corr(i, j);
            }
        }
    }

    // Add dissimilarity term 
    // which is a penalty that grows if two cliques are similar to each other
    for(int l = 0; l < N_CLIQUES; l++) {
        for(int m = l+1; m < N_CLIQUES; m++) {
            for(int i = 0; i < N_SNPS; i++) {
                this->access(l, i, m, i) += lambda2;
            }
        }
    }
}


std::tuple<double*, double, int> snps_qubo_matrix::get_qubo() {
    return std::tuple<double*, double, int>(mat, offset, n);
}


std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> snps_qubo_matrix::get_ising() {
    return convert_qubo_to_ising(this->mat, this->offset, this->get_n());
}


void snps_qubo_matrix::prettyprint_qubo() {
    std::cout << "QUBO:\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            std::cout << std::showpos << std::fixed << std::setprecision(3) << get_const(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Offset: " << this->offset << std::endl;
}


void snps_qubo_matrix::prettyprint_ising() {
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> ising = get_ising();

    std::vector<double> h = std::get<0>(ising);
    std::cout << "h value: " << h << std::endl;
    std::vector<int> coupler_start = std::get<1>(ising);
    std::cout << "J s_idx: " << coupler_start << std::endl;
    std::vector<int> coupler_end = std::get<2>(ising);
    std::cout << "J e_idx: " << coupler_end << std::endl;
    std::vector<double> J = std::get<3>(ising);
    std::cout << "J value: " << J << std::endl;
    double offset = std::get<4>(ising);
    std::cout << "offset: " << offset << std::endl;
}


// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================
// =================================================================================================================


std::vector<std::vector<int>> snps_qubo_matrix::solve_cpu_simulated_annealing(int num_samples, int num_sweeps, uint64_t seed,
    char * initial_state, const char* save_path) {

    // convert QUBO into Ising formulation
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> ising = get_ising();
    std::vector<double> h = std::get<0>(ising);
    std::vector<int> coupler_start = std::get<1>(ising);
    std::vector<int> coupler_end = std::get<2>(ising); 
    std::vector<double> J = std::get<3>(ising);
    double offset = std::get<4>(ising);

    // allocate space for the simulated annealing to save the results
    int num_variables = h.size();
    sa_return* ret = new sa_return(num_variables, num_samples);

    // create or assign the initial state
    if(initial_state == nullptr) {
        // start random number generator
        std::mt19937_64 gen(seed + 1); // pick any reproducible seed different from 'seed' itself
        std::uniform_int_distribution<uint8_t> dis;
        // create random initial state (spin variables are either -1 or +1)
        for(int i = 0; i < num_samples * num_variables; i++) {
            ret->states[i] = (dis(gen) % 2) * 2 - 1;
        }
    } else {
        // assign pre-existent initial state
        std::copy(initial_state, initial_state + num_samples * num_variables, ret->states);
    }

    // get time
    // std::time_t rawtime;
    // std::tm * timeinfo;
    // char buffer[80];
    // std::time(&rawtime);
    // timeinfo = std::localtime(&rawtime);
    // strftime (buffer, 80, "%Y%m%d_%H%M%S", timeinfo);

    // format filename with date
    std::stringstream ss_path;
    ss_path << save_path << ".json"; // save path already includes the filename << "_" << buffer << ".json";

    // open file
    std::ofstream of;
    of.open(ss_path.str());
    of << "{";

    // debug print
    // OLD PRINT std::cout << "Initial state (" << num_samples * num_variables << " elements) = {";
    // OLD PRINT for(int i = 0; i < num_samples * num_variables; i++) { std::cout << int(ret->states[i]) << ", "; }
    // OLD PRINT std::cout << "}" << std::endl;

    of << "\"elements\": " << num_samples * num_variables << ", ";
    of << "\"initial_state\": [";
    for(int i = 0; i < num_samples * num_variables; i++) {
        of << int(ret->states[i]) << ", ";
    }
    of << "], " << std::endl;

    // OLD PRINT std::cout << "Initial state (" << num_samples * num_variables << " elements) = {";
    // OLD PRINT for(int i = 0; i < num_samples * num_variables; i++) { std::cout << int((ret->states[i]) + 1)/2 << ", "; }
    // OLD PRINT std::cout << "}" << std::endl;
    
    // call simulated annealing
    simulated_annealing_ising(h, coupler_start, coupler_end, J, num_samples, num_sweeps, seed, ret);

    // debug purposes: add offset to energy (useless in practice)
    for(int i = 0; i < ret->num_samples; i++) { ret->energies[i] += offset; }

    // debug print
    // OLD PRINT std::cout << "Ising energies: ";
    // OLD PRINT for(int i = 0; i < ret->num_samples; i++) { std::cout << ret->energies[i] << " "; } std::cout << std::endl;
    of << "\"ising_energies\": [";
    for(int i = 0; i < num_samples; i++) {
        of << ret->energies[i] << ", ";
    }
    of << "], " << std::endl;

    // add debug prints of the qubo matrix
    of << "\"h\": [";
    for(unsigned int i = 0; i < h.size(); i++) {
        of << h[i] << ", ";
    }
    of << "], " << std::endl;
    of << "\"J\": [";
    for(unsigned int i = 0; i < J.size(); i++) {
        of << J[i] << ", ";
    }
    of << "], " << std::endl;
    of << "\"coupler_start\": [";
    for(unsigned int i = 0; i < coupler_start.size(); i++) {
        of << coupler_start[i] << ", ";
    }
    of << "], " << std::endl;
    of << "\"coupler_end\": [";
    for(unsigned int i = 0; i < coupler_end.size(); i++) {
        of << coupler_end[i] << ", ";
    }
    of << "], " << std::endl;

    
    // find the solution having the smallest energy
    int best_solution_index = std::min_element(ret->energies, ret->energies + ret->num_samples) - ret->energies;
    of << "\"best_solution_index\": " << best_solution_index << ", ";
    char* best_solution = ret->states + best_solution_index*ret->num_variables; // n_clique * n_snps sized char array

    // convert all the solutions to the list of "SNPs set" format
    of << "\"solutions\": [";
    for(int j = 0; j < ret->num_samples; j++) {
        char* current_solution = ret->states + j * ret->num_variables;
        of << "[";
        for(int i = 0; i < n_cliques; i++) { // iterate for each clique of the solution
            std::vector<int> snp_set; // allocate space for the current clique
            of << "[";
            for(int j = 0; j < n_snps; j++) {
                if(current_solution[i * n_snps + j] == 1) { // each variable is -1 if SNP not chosen, +1 if SNP chosen
                    snp_set.push_back(j);
                    of << j << ", ";
                }
            }
            of << "], ";
        }
        of << "], ";
    }
    of << "],";

    // redo everything for the best solution
    std::vector<std::vector<int>> snp_set_list;
    of << "\"best_snp_set_list\": [";
    for(int i = 0; i < n_cliques; i++) { // iterate for each clique of the solution
        std::vector<int> snp_set; // allocate space for the current clique
        of << "[";
        for(int j = 0; j < n_snps; j++) {
            if(best_solution[i * n_snps + j] == 1) { // each variable is -1 if SNP not chosen, +1 if SNP chosen
                snp_set.push_back(j);
                of << j << ", ";
            }
        }
        of << "], ";
        snp_set_list.push_back(snp_set);
    }
    of << "]";

    // close log file
    of << "}";
    of.close();

    // free memory
    delete ret;

    return snp_set_list;
}


std::vector<std::vector<int>> snps_qubo_matrix::get_snp_set_list_from_solution_vector(
    std::vector<int> best_solution, int n_cliques, int n_snps) {

    // Check if there was an error
    if(best_solution.empty()) {
        std::vector<std::vector<int>> empty;
        return empty;
    }

    // Convert the solution to the list of "SNPs set" format
    std::vector<std::vector<int>> snp_set_list;
    for(int i = 0; i < n_cliques; i++) { // iterate for each clique of the solution
        std::vector<int> snp_set; // allocate space for the current clique
        for(int j = 0; j < n_snps; j++) {
            if(best_solution[i * n_snps + j] == 1) { // each variable is -1 if SNP not chosen, +1 if SNP chosen
                snp_set.push_back(j);
            }
        }
        snp_set_list.push_back(snp_set);
    }

    return snp_set_list;
}

std::vector<std::vector<int>> snps_qubo_matrix::solve_cpu_parallel_tempering(
    int num_chains, int num_steps, const char* save_path) {

    // Convert QUBO into Ising formulation
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> ising = get_ising();
    std::vector<double> h = std::get<0>(ising);
    std::vector<int> coupler_start = std::get<1>(ising);
    std::vector<int> coupler_end = std::get<2>(ising); 
    std::vector<double> J = std::get<3>(ising);
    double offset = std::get<4>(ising);

    // Call python
    std::vector<int> best_solution = PythonWrapper::get_instance().run_parallel_tempering(
        h, J, coupler_start, coupler_end, num_chains, num_steps, save_path);

    std::vector<std::vector<int>> snp_set_list = get_snp_set_list_from_solution_vector(best_solution, n_cliques, n_snps);

    return snp_set_list;
}

std::vector<std::vector<int>> snps_qubo_matrix::solve_dwave_quantum_annealing(
    const char* token,
    int num_reads, int solver_idx, double fw_annealing_ramp_time, double fw_annealing_pause_time, 
    double rev_annealing_ramp_time, double rev_annealing_pause_time, double rev_annealing_s_target, 
    const char* save_path) {

    // Convert QUBO into Ising formulation
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> ising = get_ising();
    std::vector<double> h = std::get<0>(ising);
    std::vector<int> coupler_start = std::get<1>(ising);
    std::vector<int> coupler_end = std::get<2>(ising); 
    std::vector<double> J = std::get<3>(ising);
    double offset = std::get<4>(ising);

    // Call python
    std::vector<int> best_solution = PythonWrapper::get_instance().run_quantum_annealer(
        token, h, J, coupler_start, coupler_end,
        num_reads, solver_idx, 
        fw_annealing_ramp_time, fw_annealing_pause_time, 
        rev_annealing_ramp_time, rev_annealing_pause_time, rev_annealing_s_target, 
        save_path);

    std::vector<std::vector<int>> snp_set_list = get_snp_set_list_from_solution_vector(best_solution, n_cliques, n_snps);

    return snp_set_list;
}


std::vector<std::vector<int>> snps_qubo_matrix::solve_azure_optimizer(
    const char* azure_subscription_id, const char* azure_resource_group, const char* azure_name, const char* azure_location,
    int timeout_seconds, int seed, const char* save_path) {

    // Convert QUBO into Ising formulation
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> ising = get_ising();
    std::vector<double> h = std::get<0>(ising);
    std::vector<int> coupler_start = std::get<1>(ising);
    std::vector<int> coupler_end = std::get<2>(ising);
    std::vector<double> J = std::get<3>(ising);
    double offset = std::get<4>(ising);

    // Call python
    std::vector<int> best_solution = PythonWrapper::get_instance().run_azure_optimizers(
        h, J, coupler_start, coupler_end,
        azure_subscription_id, azure_resource_group, azure_name, azure_location,
        timeout_seconds, seed, save_path);

    std::vector<std::vector<int>> snp_set_list = get_snp_set_list_from_solution_vector(best_solution, n_cliques, n_snps);

    return snp_set_list;
}


std::vector<std::vector<int>> snps_qubo_matrix::solve_azure_QAOA(const char* vendor,
        const char* azure_subscription_id, const char* azure_resource_group, const char* azure_name, const char* azure_location,
        const char* azure_backend,
        const char* optimizer, int maxiter, int reps, int n_shots, int is_recursive_qaoa, const char* save_path) {

    // Convert QUBO into Ising formulation
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> ising = get_ising();
    std::vector<double> h = std::get<0>(ising);
    std::vector<int> coupler_start = std::get<1>(ising);
    std::vector<int> coupler_end = std::get<2>(ising);
    std::vector<double> J = std::get<3>(ising);
    double offset = std::get<4>(ising);

    // Call python
    std::vector<int> best_solution = PythonWrapper::get_instance().run_qaoa(
        vendor, azure_subscription_id, azure_resource_group, azure_name, azure_location, azure_backend,
        h, J, coupler_start, coupler_end,
        optimizer, maxiter, reps, n_shots, is_recursive_qaoa, save_path);

    std::vector<std::vector<int>> snp_set_list = get_snp_set_list_from_solution_vector(best_solution, n_cliques, n_snps);

    return snp_set_list;
}


std::vector<std::vector<int>> snps_qubo_matrix::solve_atos_simulated_quantum_annealing(int clique_size, int shots, const char* save_path) {

    // Convert QUBO into Ising formulation
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> ising = get_ising();
    std::vector<double> h = std::get<0>(ising);
    std::vector<int> coupler_start = std::get<1>(ising);
    std::vector<int> coupler_end = std::get<2>(ising);
    std::vector<double> J = std::get<3>(ising);
    double offset = std::get<4>(ising);

    // Call Python
    std::vector<int> best_solution = PythonWrapper::get_instance().run_sqa(
        h, J, coupler_start, coupler_end, clique_size, shots, save_path);

    std::vector<std::vector<int>> snp_set_list = get_snp_set_list_from_solution_vector(best_solution, n_cliques, n_snps);

    return snp_set_list;
}
