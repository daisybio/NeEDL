#ifndef SNPS_QUBO_H__
#define SNPS_QUBO_H__

#include <tuple>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <stdexcept>
#include "cpu_sa_wrapper.h"

// replace it with boost ublas triangular matrix
class matrix {
protected:

    double* mat;
    int n;

public:

    matrix(double* m, int n);

    matrix(int n);

    ~matrix();

    double& get(int i, int j);

    double& operator()(int i, int j);

    double get_const(int i, int j) const;

    double operator()(int i, int j) const;

    int get_n() const;

    double* get_matrix();

};

/**
 * @brief Calculate linear combination of two matrices a, b
 * @param nu: Coefficient of the linear combination
 * @param a: First square matrix
 * @param b: Second square matrix
 * @return matrix corresponding to nu * a + (1-nu) * b
 */
matrix linear_combination(double nu, const matrix& a, const matrix& b);


/**
 * @brief Class in charge of constructing the QUBO formulation, converting it to Ising form and run it on any
 * type of simulated and real quantum hardware
 */
class snps_qubo_matrix : public matrix {

    int n_cliques, n_snps;
    double offset;

    /**
     * @brief Construct the SNP set list from the solution of the Ising problem in vector format
     * @param best_solution: vector of spin solutions
     * @param n_cliques: number of cliques of the solution
     * @param n_snps: number of SNPs of the solution
     * @param return: list of SNPs set
     */
    std::vector<std::vector<int>> get_snp_set_list_from_solution_vector(
        std::vector<int> best_solution, int n_cliques, int n_snps);

public:


    /**
     * @brief Create and allocate space for the QUBO formulation of the problem having the wanted number of output cliques and the 
     * required number of SNPs
     * @param n_cliques number of output cliques
     * @param n_snps number of input SNPs
     */
    snps_qubo_matrix(int n_cliques, int n_snps);


    /**
     * @brief Destructor, deallocation of matrix is within the super class
     */
    ~snps_qubo_matrix();


    /**
     * @brief Write access of the QUBO component representing the term x[1st clique, 1st variable] * x[2nd clique, 2nd variable]
     * @param first_clique clique of the first term
     * @param first_variable variable of the first term
     * @param second_clique clique of the second term
     * @param second_variable variable of the second term
     * @return reference access
     */
    double& access(int first_clique, int first_variable, int second_clique, int second_variable);


    /**
     * @brief Read access of the QUBO component representing the term x[1st clique, 1st variable] * x[2nd clique, 2nd variable]
     * @param first_clique clique of the first term
     * @param first_variable variable of the first term
     * @param second_clique clique of the second term
     * @param second_variable variable of the second term
     * @return copied value of the matrix
     */
    double access(int first_clique, int first_variable, int second_clique, int second_variable) const;


    /**
     * @brief Create the QUBO matrix corresponding to the "N max-weighted k cliques"
     * @param ss_corr: get correlation between i-th and j-th snp
     * @param N_CLIQUES: how many snps set has to be found within the current subproblem
     * @param K: recommended size for the SNPs set that has to be found
     * @param lambda0: hyperparameter setting the strenght of the SNPs size constraint
     * @param lambda1: hyperparameter setting the strength of the SNPs correlation sum penalty term
     * @param lambda2: hyperparameter setting the strength of the SNPs sets dissimilarity term
     */
    void fill_n_max_weighted_k_clique_qubo(const matrix& ss_corr, int K, double lambda0, double lambda1, double lambda2);


    /**
     * @brief Create the QUBO matrix corresponding to the "N max-weighted k cliques"
     * @param stat_corr: get correlation between i-th and j-th snp
     * @param bio_corr:
     * @param N_CLIQUES: how many snps set has to be found within the current subproblem
     * @param K: recommended size for the SNPs set that has to be found
     * @param nu:
     * @param lambda0: hyperparameter setting the strenght of the SNPs size constraint
     * @param lambda1: hyperparameter setting the strength of the SNPs correlation sum penalty term
     * @param lambda2: hyperparameter setting the strength of the SNPs sets dissimilarity term
     */
    void fill_n_max_weighted_k_clique_qubo(const matrix& stat_corr, const matrix& bio_corr, int K, double nu, double lambda0, double lambda1, double lambda2);


    /**
     * @brief Get the object representing the QUBO formulation, consisting of the (double*) qubo matrix Q, the (double) offset and the number of variables (int) N
     * @return QUBO formulation as QUBO matrix (flattened into a 1D vector of size n^2), the offset, the dimensionality of the matrix n
     */
    std::tuple<double*, double, int> get_qubo();


    /**
     * @brief Convert the QUBO formulation into the Ising representation (diagonal h, off-diagonal elements splitted in start index, end index, values and offset)
     * @return QUBO formulation as Ising in terms of h vector, J start index vector, J end index vector, J vector, offset
     */
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, std::vector<double>, double> get_ising();


    /**
     * @brief Pretty-print the QUBO matrix and offset
     */
    void prettyprint_qubo();


    /**
     * @brief Pretty-print the Ising equivalent of the QUBO matrix
     */
    void prettyprint_ising();

    // =================================================================================================================
    // =================================================================================================================
    // =================================================================================================================
    // =================================================================================================================
    // =================================================================================================================

    /**
     * @brief Solve the current QUBO matrix with simulated annealing
     * @param num_samples: number of times the annealing process is independently repeated starting from a different point
     * @param num_sweeps: number of steps of the computation
     * @param seed: random seed
     * @param initial_state: fix that to nullptr in production, otherwise assign a num_var * num_sample char array having values in {-1, +1} for reproducibility purposes
     * @param save_path Path where to save the file including the first part of the file name
     * @return A vector of 'n_cliques' elements, each element represents a suggested SNPs set
     */
    std::vector<std::vector<int>> solve_cpu_simulated_annealing(int num_samples, int num_sweeps, uint64_t seed,
        char* initial_state = nullptr, const char* save_path = nullptr);


    /**
     * @brief Solve the current QUBO matrix with quantum annealing (call Python then D-Wave devices). In case of errors the return value is the empty vector.
     * @param token: D-Wave platform token
     * @param num_reads Number of times the annealing process is repeated with a different, random initial state
     * @param solver_idx Use 4 to use latest D-Wave machine 'Advantage-6.1'
     * @param fw_annealing_ramp_time For FORWARD annealing, it is the time taken to go from 0% to 50% of the annealing process. The same value indicates the time the annealing takes to go from 50% to 100%
     * @param fw_annealing_pause_time For FORWARD annealing, it is the time that the annealing stops at 50%
     * @param rev_annealing_ramp_time For REVERSE annealing, it is the time taken to go from 100% to S_TARGET% of the annealing process. The same value indicates the time the annealing takes to go from S_TARGET% to 100%
     * @param rev_annealing_pause_time For REVERSE annealing, it is the time that the annealing stops at S_TARGET%
     * @param rev_annealing_s_target Indicated the point where the annealing stop running in reverse.
     * @param save_path Path where to save the file including the first part of the file name
     * @return A vector of 'n_cliques' elements, each element represents a suggested SNPs set OR the empty vector in case of errors. 
     */
    std::vector<std::vector<int>> solve_dwave_quantum_annealing(const char* token, int num_reads, int solver_idx,
        double fw_annealing_ramp_time, double fw_annealing_pause_time, double rev_annealing_ramp_time,
        double rev_annealing_pause_time, double rev_annealing_s_target, const char* save_path);


    /**
     * @brief Solve the current QUBO matrix with the classical QUBO optimizers given in Azure QIO environment
     * @param azure_subscription_id:
     * @param azure_resource_group:
     * @param azure_name:
     * @param azure_location:
     * @param timeout_seconds: number of seconds for the computation to complete
     * @param seed: random seed
     * @param save_path: path where to save the file including the first part of the file name
     */
    std::vector<std::vector<int>> solve_azure_optimizer(
        const char* azure_subscription_id, const char* azure_resource_group, const char* azure_name, const char* azure_location,
        int timeout_seconds, int seed, const char* save_path);


    /**
     * @brief Solve the current QUBO matrix with QAOA on simulator or real quantum hw through Azure
     * @param vendor: 'none' to run simulation on actual hardware, 'azure' to use azure subscription
     * @param azure_subscription_id:
     * @param azure_resource_group:
     * @param azure_name:
     * @param azure_location:
     * @param azure_backend: One of the backend available for your account in Azure Quantum or empty string
     * @param optimizer: one between 'ADAM', 'COBYLA', 'NELDER_MEAD' or 'SLSQP'
     * @param maxiter: Number of iteration of the optimization algorithm
     * @param reps: Number of repetitions of the QAOA ansatz
     * @param n_shots: Number of sample estimating the expectation value
     * @param is_recursive_qaoa: Use recursive QAOA algorithm
     * @param save_path: Path where to save the file including the first part of the file name
     * @return A vector of 'n_cliques' elements, each element represents a suggested SNPs set OR the empty vector in case of errors.
     */
    std::vector<std::vector<int>> solve_azure_QAOA(const char* vendor,
        const char* azure_subscription_id, const char* azure_resource_group, const char* azure_name, const char* azure_location,
        const char* azure_backend,
        const char* optimizer, int maxiter, int reps, int n_shots, int is_recursive_qaoa, const char* save_path);


    /**
     * @brief Solve the current QUBO matrix with Simulated Quantum Annealing using ATOS QLM hardware
     * @param clique_size: size of the wanted clique
     * @param shots: number of shots
     * @param save_path: Path where to save the file including the first part of the file name
     * @return A vector of 'n_cliques' elements, each element represents a suggested SNPs set OR the empty vector in case of errors.
     */
    std::vector<std::vector<int>> solve_atos_simulated_quantum_annealing(int clique_size, int shots, const char* save_path);


};

#endif