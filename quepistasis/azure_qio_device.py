import sys
print("azure_qio_device.py: old sys.path is ", sys.path)
low_prio_path = [p for p in sys.path if "dist-packages" in p]
high_prio_path = [p for p in sys.path if "dist-packages" not in p]
sys.path = low_prio_path + high_prio_path
print("azure_qio_device.py: new sys.path is ", sys.path)

import json
import numpy as np
import itertools
import time
from utils import Trace
from typing import List
from azure.quantum import Workspace
from azure.quantum.optimization import Problem, ProblemType, Term
from azure.quantum.optimization import ParallelTempering, SimulatedAnnealing, PopulationAnnealing, Tabu, SubstochasticMonteCarlo, QuantumMonteCarlo


# convert Ising to AZURE formulation
def create_problem(h, J) -> Problem:
    terms: List[Term] = []
    for i, h_coeff in enumerate(h):
        terms.append(Term(c=h_coeff, indices=[i]))
    for ((start, end), j_coeff) in J.items():
        if start != end: # constant term already added
            terms.append(Term(c=j_coeff, indices=[start, end]))
    # Return an Ising-type problem
    return Problem(name="Epistasis", problem_type=ProblemType.ising, terms=terms)


def run_azure_qio(h, J, azure_subscription_id, azure_resource_group, azure_name, azure_location,
                                    timeout_seconds, seed, save_path):
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
    trace.add('input', None, 'timeout_seconds', timeout_seconds)
    trace.add('input', None, 'seed', seed)
    trace.add('input', None, 'save_path', save_path)

    # get number of qubits (corresponding also to the variables)
    n = len(h)

    # Copy the settings for your workspace below
    workspace = Workspace(
        subscription_id=azure_subscription_id,
        resource_group=azure_resource_group,
        name=azure_name,
        location=azure_location
    )

    # Instantiate solvers
    pt_solver = ParallelTempering(workspace, timeout=timeout_seconds)
    sa_solver = SimulatedAnnealing(workspace, timeout=timeout_seconds, seed=seed)
    pa_solver = PopulationAnnealing(workspace, timeout=timeout_seconds, seed=seed)
    ta_solver = Tabu(workspace, timeout=timeout_seconds, seed=seed)
    smc_solver = SubstochasticMonteCarlo(workspace, timeout=timeout_seconds, seed=seed)
    qmc_solver = QuantumMonteCarlo(workspace, sweeps=2, trotter_number=10, restarts=72,
                                   beta_start=0.1, transverse_field_start=10, transverse_field_stop=0.1, seed=seed)
    solver_list = [pt_solver] # , sa_solver, pa_solver, ta_solver, smc_solver, qmc_solver]
    solver_names = ["ParallelTempering"] # , "SimulatedAnnealing", "PopulationAnnealing", "Tabu", "SubstochasticMonteCarlo", "QuantumMonteCarlo"]

    # create AZURE problem
    problem = create_problem(h, J)

    # solve problem
    solution_list = []
    cost_list = []

    for solver, solver_name in zip(solver_list, solver_names):
        start = time.time()
        result = solver.optimize(problem)
        time_elapsed = time.time() - start
        best_solution = result['configuration']
        best_cost = result['cost']
        trace.add('solution', solver_name, 'time_elapsed', str(time_elapsed))
        trace.add('solution', solver_name, 'best_cost', best_cost)
        trace.add('solution', solver_name, 'best_solution', str(best_solution))
        trace.add('solution', solver_name, 'parameters', str(result['parameters']))
        solution_list.append(best_solution)
        cost_list.append(best_cost)

    # find best solution
    best_index = np.argmin(np.array(cost_list))
    best_cost = cost_list[best_index]
    best_solution = solution_list[best_index]
    trace.add('best_solution', None, 'solver', solver_names[best_index])
    trace.add('best_solution', None, 'solution', str(best_solution))
    trace.add('best_solution', None, 'cost', best_cost)

    # rewrite solution from dict to list
    best_spins = [best_solution[str(i)] for i in range(n)]

    trace.save()
    return best_spins


def test_run_azure_qio():
    """Run a basic Ising model to check if the library works correctly"""
    # h = [0, 0]
    # J = {(0, 0): 0, (0, 1): 1, (1, 0): 0, (1, 1): 0.5}
    N = 4
    h = np.zeros(shape=(N,))
    J = {(i, j): 0.1 for i in range(N) for j in range(i+1, N)}

    password_azure = json.load(open("password_azure.json"))
    azure_subscription_id = password_azure["azure_subscription_id"]
    azure_resource_group = password_azure["azure_resource_group"]
    azure_name = password_azure["azure_name"]
    azure_location = password_azure["azure_location"]
    timeout_seconds = 1
    seed = 12345
    save_path = "."

    spins = run_azure_qio(h, J, azure_subscription_id, azure_resource_group, azure_name, azure_location, timeout_seconds, seed, save_path)
    print("Best spins: ")
    print(spins)
