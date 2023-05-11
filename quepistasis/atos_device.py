import numpy as np
import networkx as nx
import itertools

from qat.core import Variable
from qat.opt import KClique, Ising
from qat.sqa import SQAQPU
from qat.sqa.sqa_qpu import integer_to_spins
from qat.plugins import ScipyMinimizePlugin
from qat.qpus import get_default_qpu

from utils import Trace, get_ising_cost


def get_best_sqa_parameters(h, J, clique_size):
    """
    Interpret h, j as a (non-complete) graph, then through KClique get the best parameters for SQA optimization
    :param h: diagonal coefficients
    :param J: off diagonal coefficients
    :param clique_size: K of k-clique
    :return: best parameters dictionary
    """
    graph = nx.Graph()
    vertices = np.arange(len(h))
    graph.add_nodes_from(vertices)
    avg_weight = np.average([weight for _, weight in J.items()])
    edges = [edge for edge, weight in J.items() if weight > avg_weight]
    graph.add_edges_from(edges)
    k_clique_problem = KClique(graph, clique_size, A=clique_size + 1, B=1)
    return k_clique_problem.get_best_parameters()


def run_atos_simulated_quantum_annealing(h, J, clique_size, shots, save_path):
    """
    Run Simulated Quantum Annealing algorithm on QLM device
    :param h: diagonal Ising coefficients
    :param J: off-diagonal Ising coefficients
    :param clique_size: clique size of the QUBO formulation, used to estimate the SQA parameters
    :param shots: number of shots in the SQA
    :param save_path: path to save the logs
    :return: array of spins
    """
    trace = Trace(print=True, path=save_path)
    trace.add('input', None, 'clique_size', clique_size)
    trace.add('input', None, 'shots', shots)
    trace.add('input', None, 'save_path', save_path)

    # get number of qubits (corresponding also to the variables)
    n = len(h)

    # estimate best SQA parameters
    best_parameters = get_best_sqa_parameters(h, J, clique_size)
    n_monte_carlo_updates = best_parameters['n_monte_carlo_updates']
    n_trotters = best_parameters['n_trotters']
    nsteps_per_gamma = best_parameters['nsteps_per_gamma']
    gamma_max = best_parameters['gamma_max']
    gamma_min = best_parameters['gamma_min']
    temp_max = best_parameters['temp_max']
    temp_min = best_parameters['temp_min']
    n_steps = best_parameters['n_steps']  # n_monte_carlo_updates / (n_trotters * nsteps_per_gamma)
    trace.add('sqa_parameters', None, 'n_monte_carlo_updates', n_monte_carlo_updates)
    trace.add('sqa_parameters', None, 'n_trotters', n_trotters)
    trace.add('sqa_parameters', None, 'nsteps_per_gamma', nsteps_per_gamma)
    trace.add('sqa_parameters', None, 'n_steps', n_steps)
    trace.add('sqa_parameters', None, 'gamma_max', gamma_max)
    trace.add('sqa_parameters', None, 'gamma_min', gamma_min)
    trace.add('sqa_parameters', None, 'temp_max', temp_max)
    trace.add('sqa_parameters', None, 'temp_min', temp_min)

    # 1. construct QUBO
    # SQA maximizes the cost function, thus we need to negate the coefficients
    J_matrix = np.zeros(shape=(n, n))
    for (start, end), coeff in J.items():
        J_matrix[start][end] += -coeff / 2
        J_matrix[end][start] += -coeff / 2
    ising_problem = Ising(J=J_matrix, h=-1 * h)

    # 2. Create a temperature and a gamma schedule
    tmax = 1.0
    t = Variable("t", float)
    temp_t = temp_min * (t / tmax) + temp_max * (1 - t / tmax)
    gamma_t = gamma_min * (t / tmax) + gamma_max * (1 - t / tmax)

    # 3. Create a job and send it to a QPU
    problem_job = ising_problem.to_job(gamma_t=gamma_t, tmax=tmax, nbshots=1)
    best_cost, best_spins = None, None
    for i in range(shots):
        sqa_qpu = SQAQPU(temp_t=temp_t, n_steps=n_steps, n_trotters=n_trotters, seed=shots)
        problem_result = sqa_qpu.submit(problem_job)
        spins = integer_to_spins(problem_result.raw_data[0].state.int, n)
        cost = get_ising_cost(h, J, spins)
        trace.add('candidate_solutions', str(i), 'spins', str(spins.tolist()))
        trace.add('candidate_solutions', str(i), 'cost', cost)
        if best_cost is None or cost < best_cost:
            best_cost = cost
            best_spins = spins
    trace.add('solution', None, 'spins', str(best_spins.tolist()))
    trace.add('solution', None, 'cost', best_cost)
    trace.save()
    return best_spins


def run_atos_qaoa(h, J, depth, seed, save_path):
    """
    Run noiseless QAOA algorithm on QLM device
    :param h: diagonal Ising coefficients
    :param J: off-diagonal Ising coefficients
    :param depth: number of layers of QAOA Ansatz
    :param seed: seed of the optimizer
    :return:
    """

    trace = Trace(print=True, path=save_path)
    trace.add('input', None, 'depth', depth)
    trace.add('input', None, 'seed', seed)
    trace.add('input', None, 'save_path', save_path)

    # randomness is managed by numpy, we cannot seed in the optimizer
    np.random.seed(seed)

    # get number of qubits (corresponding also to the variables)
    n = len(h)

    # SQA maximizes the cost function, thus i need to negate the coefficients
    J_matrix = np.zeros(shape=(n, n))
    for (start, end), coeff in J.items():
        J_matrix[start][end] += coeff / 2
        J_matrix[end][start] += coeff / 2
    ising = Ising(J=J_matrix, h=h)

    # setup training + qpu
    optimizer = ScipyMinimizePlugin(method="COBYLA", tol=1e-5, options={"maxiter": 200})
    trace.add('optimizer', None, 'method', "COBYLA")
    trace.add('optimizer', None, 'tol', 1e-5)
    trace.add('optimizer', None, 'maxiter', 200)
    qpu = get_default_qpu()
    stack = optimizer | qpu

    # run optimization
    job = ising.qaoa_ansatz(depth)
    trace.add('circuit', None, 'depth', "xxx")
    trace.add('circuit', None, 'size', "xxx")
    trace.add('circuit', None, 'single_qubit_gate_number', "xxx")
    trace.add('circuit', None, 'two_qubit_gate_number', "xxx")
    result = stack.submit(job)

    # get value of the parameters and get the most probable states
    sol_job = job(**eval(result.meta_data["parameter_map"]))  # Binding the variables
    sampling_job = sol_job.circuit.to_job()  # Rerunning in 'SAMPLE' mode to get the most probable states
    sol_res = qpu.submit(sampling_job)

    # get the cost of any of the most probable states
    best_cost, best_spins = None, None
    i = 0
    for sample in sol_res:
        if sample.probability > 0.05:
            spins = integer_to_spins(sample.state.int, n)
            cost = get_ising_cost(h, J, spins)
            trace.add('candidate_solutions', str(i), 'spins', str(spins.tolist()))
            trace.add('candidate_solutions', str(i), 'cost', cost)
            i = i + 1
            if best_cost is None or cost < best_cost:
                best_cost = cost
                best_spins = spins

    trace.add('solution', None, 'spins', str(best_spins.tolist()))
    trace.add('solution', None, 'cost', best_cost)
    trace.save()
    return best_spins


def test_run_atos_simulated_quantum_annealing():
    """Run a basic Ising model to check if the library works correctly"""
    # h = [0, 0]
    # J = {(0, 0): 0, (0, 1): 1, (1, 0): 0, (1, 1): 0.5}
    N = 3
    h = np.array([-12.0, 0.5, 0.3])
    J = {(0, 1): 0.2, (0, 2): 0.5, (1, 2): 1.2}
    for variable_assignments in itertools.product([-1, 1], repeat=N):
        print(f"vars={str(variable_assignments):20s} cost={get_ising_cost(h, J, variable_assignments):3.3f}")
    clique_size = 5
    n_shots = 10
    result = run_atos_simulated_quantum_annealing(h, J, clique_size, n_shots)
    print(f"Winning assignment = {result} with cost {get_ising_cost(h, J, result)}")
    return result


def test_run_atos_qaoa():
    """Run a basic Ising model to check if the library works correctly"""
    # h = [0, 0]
    # J = {(0, 0): 0, (0, 1): 1, (1, 0): 0, (1, 1): 0.5}
    N = 3
    h = np.array([-12.0, 0.5, 0.3])
    J = {(0, 1): 0.2, (0, 2): 0.5, (1, 2): 1.2}
    for variable_assignments in itertools.product([-1, 1], repeat=N):
        print(f"vars={str(variable_assignments):20s} cost={get_ising_cost(h, J, variable_assignments):3.3f}")
    result = run_atos_qaoa(h, J, depth=10, seed=12345)
    print(f"Winning assignment = {result} with cost {get_ising_cost(h, J, result)}")
    return result

