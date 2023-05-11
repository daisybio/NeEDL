#ifndef PYTHON_WRAPPER_H_
#define PYTHON_WRAPPER_H_

#include <Python.h>
#include <vector>
#include <iostream>

/**
 * @brief Class that wraps the C-Python API instructions
 */
class PythonWrapper {

    /**
     * Object representing the module used to interact with the quantum annealer
     */
    PyObject *pModuleQuantumAnnealer = nullptr;

    /**
     * Object representing the module used to interact with the quantum computers (IBM Azure ecc)
     */
    PyObject *pModuleAzureOptimizer = nullptr;

    /**
     * Object representing the module used to interact with the quantum computers (IBM Azure ecc)
     */
    PyObject *pModuleQAOA = nullptr;

    /**
     * Object representing the function to call "run_quantum_annealer" function in Python
     */
    PyObject *pFuncQuantumAnnealer = nullptr;

    /**
     * Object representing the function to call "run_quantum_annealer" function in Python
     */
    PyObject *pFuncAzureOptimizer = nullptr;

    /**
     * Object representing the function to call "run_qaoa" function in Python
     */
    PyObject *pFuncQAOA = nullptr;

    /**
     * @brief Load a Python file as module
     * @param filepath: path relative to the working directory or absolute
     * @return Python object representing the module
     */
    PyObject* get_module(const char* filepath);

    /**
     * @brief Load a Python function
     * @param functionName: name of the function
     * @param pModule: Python object pointing to the module containing the wanted function
     * @return Python object representing the function
     */
    PyObject* get_function(const char* functionName, PyObject *pModule);

    /**
     * @brief Call a Python function with the given arguments
     * @param pFunc: Python object pointing to the function
     * @param pArgs: Python object pointing to the arguments
     * @return Python object representing the result
     */
    PyObject* call_function(PyObject *pFunc, PyObject *pArgs);

    /**
     * @brief Create a Python object representing a list of floats starting from a vector of C doubles. In the context of Ising optimization, h is the vector of linear coefficients.
     * @param h: vector of double precision floating points elements
     * @param garbage: vector to store the intermediate objects to be deleted
     * @return Python object representing h
     */
    PyObject* create_h(std::vector<double> h, std::vector<PyObject*>& garbage);

    /**
     * @brief Create a Python object representing a dictionary {[int,int]->float}. In the context of Ising optimization, J is the vector of quadratic coefficients.
     * @param J: quadratic coefficients
     * @param start: first index of the coefficient
     * @param end: second index of the coefficient
     * @param garbage: vector to store the intermediate objects to be deleted
     * @return Python object representing J, start, end
     */
    PyObject* create_j(std::vector<double> J, std::vector<int> start, std::vector<int> end, std::vector<PyObject*>& garbage);

    /**
     * Transform a Python object representing the solution of the optimization problem (list of ints) into its C counterpart
     * @param pResult: Python object result of the optimization process
     * @return C vector of integers, solution of the optimization problem
     */
    std::vector<int> unpack_result(PyObject* pResult);

    /**
     * @brief Private constructor. Does nothing, call init to start the environment
     */
    PythonWrapper();

    /**
     * Don't implement copy constructor
     */
    PythonWrapper(PythonWrapper const&);

    /**
     * Don't implement copy operator
     */
    void operator=(PythonWrapper const&);

public:

    /**
     * @brief Initialize the Python environment
     * @return 0 if everything's ok, -1 otherwise
     */
    static PythonWrapper& get_instance();

    /**
     * @brief Initialize the Python environment
     * @return 0 if everything's ok, -1 otherwise
     */
    int init();

    /**
     * @brief Deallocate all the variables and finalize the Python environment
     */
    void close();

    /**
     * @brief Run the quantum annealer method in Python
     * @param token: D-Wave platform token
     * @param h: Ising formulation h coefficients (which in Python results in the object np.nparray[float] of shape (N,))
     * @param j: Ising formulation J coefficients (which - together with start and end params - in Python results in the object dict[tuple(int,int) -> float])
     * @param start: vector Ising formulation J starting index
     * @param end: vector Ising formulation J starting index
     * @param num_reads: Number of samples
     * @param solver_idx: D-Wave solver. Use always 4 for Advantage 6.1
     * @param fw_annealing_ramp_time: For FORWARD annealing, it is the time taken to go from 0% to 50% of the annealing process. The same value indicates the time the annealing takes to go from 50% to 100%
     * @param fw_annealing_pause_time: For FORWARD annealing, it is the time that the annealing stops at 50%
     * @param rev_annealing_ramp_time: For REVERSE annealing, it is the time taken to go from 100% to S_TARGET% of the annealing process. The same value indicates the time the annealing takes to go from S_TARGET% to 100%
     * @param rev_annealing_pause_time: For REVERSE annealing, it is the time that the annealing stops at S_TARGET%
     * @param rev_annealing_s_target: Indicated the point where the annealing stop running in reverse.
     * @param save_path: Path where to save the file including the first part of the file name
     * @return solution to the Ising problem or the empty list in case of errors.
     */
    std::vector<int> run_quantum_annealer(
        const char* token, std::vector<double> h, std::vector<double> J, std::vector<int> start, std::vector<int> end,
        int num_reads, int solver_idx, double fw_annealing_ramp_time, double fw_annealing_pause_time, 
        double rev_annealing_ramp_time, double rev_annealing_pause_time, double rev_annealing_s_target, 
        const char* save_path);

    /**
     * @brief Run the optimizers available of Azure QIO
     * @param h: Ising formulation h coefficients (which in Python results in the object np.nparray[float] of shape (N,))
     * @param j: Ising formulation J coefficients (which - together with start and end params - in Python results in the object dict[tuple(int,int) -> float])
     * @param start: vector Ising formulation J starting index
     * @param end: vector Ising formulation J starting index
     * @param azure_subscription_id:
     * @param azure_resource_group:
     * @param azure_name:
     * @param azure_location:
     * @param timeout_seconds: number of seconds for the computation to complete
     * @param seed: random seed
     * @param save_path: path where to save the file including the first part of the file name
     * @return solution to the Ising problem or the empty list in case of errors.
     */
    std::vector<int> run_azure_optimizers(
        std::vector<double> h, std::vector<double> J, std::vector<int> start, std::vector<int> end,
        const char* azure_subscription_id, const char* azure_resource_group, const char* azure_name, const char* azure_location,
        int timeout_seconds, int seed, const char* save_path);


    /**
     * @brief Run the QAOA method in Python
     * @param vendor: 'none' to run simulation on actual hardware, 'azure' to use azure subscription
     * @param azure_subscription_id:
     * @param azure_resource_group:
     * @param azure_name:
     * @param azure_location:
     * @param azure_backend: One of the backend available for your account in Azure Quantum or empty string
     * @param h: Ising formulation h coefficients (which in Python results in the object np.nparray[float] of shape (N,))
     * @param j: Ising formulation J coefficients (which - together with start and end params - in Python results in the object dict[tuple(int,int) -> float])
     * @param start: vector Ising formulation J starting index
     * @param end: vector Ising formulation J starting index
     * @param optimizer: one between 'ADAM', 'COBYLA', 'NELDER_MEAD' or 'SLSQP'
     * @param maxiter: Number of iteration of the optimization algorithm
     * @param reps: Number of repetitions of the QAOA ansatz
     * @param n_shots: Number of sample estimating the expectation value
     * @param is_recursive_qaoa: 0 if you want to use plain QAOA, != 0 if you want to use Recursive-QAOA
     * @param save_path: Path where to save the file including the first part of the file name
     * @return solution to the Ising problem or the empty list in case of errors.
     */
    std::vector<int> run_qaoa(
        const char* vendor,
        const char* azure_subscription_id, const char* azure_resource_group, const char* azure_name, const char* azure_location,
        const char* azure_backend, std::vector<double> h, std::vector<double> J, std::vector<int> start, std::vector<int> end,
        const char* optimizer, int maxiter, int reps, int n_shots, int is_recursive_qaoa, const char* save_path);


    /**
     * @brief Run the SQA method in Python
     * @param h: Ising formulation h coefficients (which in Python results in the object np.nparray[float] of shape (N,))
     * @param j: Ising formulation J coefficients (which - together with start and end params - in Python results in the object dict[tuple(int,int) -> float])
     * @param start: vector Ising formulation J starting index
     * @param end: vector Ising formulation J starting index
     * @param clique_size: size of the wanted clique
     * @param shots: number of shots
     * @param save_path: Path where to save the file including the first part of the file name
     * @return solution to the Ising problem or the empty list in case of errors.
     */
    std::vector<int> run_sqa(
        std::vector<double> h, std::vector<double> J, std::vector<int> start, std::vector<int> end,
        int clique_size, int shots, const char* save_path);
};

#endif