import sys
print("azure_device.py: old sys.path is ", sys.path)
low_prio_path = [p for p in sys.path if "dist-packages" in p]
high_prio_path = [p for p in sys.path if "dist-packages" not in p]
sys.path = low_prio_path + high_prio_path
print("azure_device.py: new sys.path is ", sys.path)

from azure.quantum.qiskit import AzureQuantumProvider
from azure.quantum import Workspace
import qiskit
import json
from qiskit.circuit import Parameter
from qiskit.utils import QuantumInstance
from qiskit.opflow import PauliSumOp
from qiskit.algorithms import QAOA
from qiskit_optimization.algorithms import MinimumEigenOptimizer, RecursiveMinimumEigenOptimizer
from qiskit_optimization.problems import QuadraticProgram
from qiskit.algorithms.optimizers import ADAM, COBYLA, NELDER_MEAD, SLSQP

import numpy as np
import itertools
from utils import Trace, get_ising_cost


def create_pauli_term(i, j, N):
    """
    Create a Pauli ZZ term at the given indices i, j of size N
    :param i:
    :param j:
    :param N:
    :return:
    """
    return ''.join(['Z' if k in [i, j] else 'I' for k in range(N)])


def convert_ising_to_operator(h, J):
    """

    :param h:
    :param J:
    :return:
    """
    N = len(h)
    pauli_list = []
    for i, h_coeff in enumerate(h):
        pauli_term = (create_pauli_term(i, i, N), h_coeff)
        pauli_list.append(pauli_term)
    for ((start, end), j_coeff) in J.items():
        pauli_term = (create_pauli_term(start, end, N), j_coeff)
        pauli_list.append(pauli_term)
    pauli_operator = PauliSumOp.from_list(pauli_list)
    return pauli_operator, pauli_list


def run_qaoa(vendor, azure_subscription_id, azure_resource_group, azure_name, azure_location, azure_backend, h, J,
             optimizer, maxiter, reps, n_shots, is_recursive_qaoa, save_path):
    """
    Run the quantum annealer
    :param vendor: Hardware vendor, one between 'none', 'azure'
    :param azure_subscription_id: Azure Quantum subscription id
    :param azure_location: Azure Quantum Zone
    :param azure_backend: QPU name, your subscription must have access to such hardware otherwise an error occurs
    :param h: Ising formulation h coefficients (np.nparray[float] of shape (N,))
    :param j: Ising formulation J coefficients (dict[tuple(int,int) -> float])
    :param optimizer: one between 'ADAM', 'COBYLA', 'NELDER_MEAD' or 'SLSQP'
    :param maxiter: Number of iteration of the optimization algorithm
    :param reps: Number of repetitions of the QAOA ansatz
    :param n_shots: Number of sample estimating the expectation value
    :param is_recursive_qaoa: 0 if you want to use plain QAOA, != 0 if you want to use Recursive-QAOA
    :param save_path: Path where to save the file including the first part of the file name
    :return solution to the qubo problem (list[int] of size N) or [] in case of errors.
    """

    trace = Trace(print=True, path=save_path)

    trace.add('input', None, 'vendor', vendor)
    trace.add('input', 'azure', 'azure_subscription_id', 'hidden')
    trace.add('input', 'azure', 'azure_resource_group', azure_resource_group)
    trace.add('input', 'azure', 'azure_name', azure_name)
    trace.add('input', 'azure', 'azure_location', azure_location)
    trace.add('input', 'azure', 'azure_backend', azure_backend)
    trace.add('input', None, 'optimizer', optimizer)
    trace.add('input', None, 'maxiter', maxiter)
    trace.add('input', None, 'reps', reps)
    trace.add('input', None, 'n_shots', n_shots)
    trace.add('input', None, 'is_recursive_qaoa', is_recursive_qaoa)

    # get QPU or simulator device
    if vendor == 'none':
        trace.add('configuration', None, 'device', 'qasm_simulator')
        backend = qiskit.Aer.get_backend('qasm_simulator')
    elif vendor == 'azure':
        workspace = Workspace(
            subscription_id=azure_subscription_id,
            resource_group=azure_resource_group,
            name=azure_name,
            location=azure_location
        )
        provider = AzureQuantumProvider(workspace=workspace)
        trace.add('input', 'azure', 'azure_available_backends', [backend.name() for backend in provider.backends()])
        trace.add('configuration', None, 'device', azure_backend)
        backend = provider.get_backend(azure_backend)

    # convert Ising formulation to Pauli operators
    pauli_operator, pauli_list = convert_ising_to_operator(h, J)
    trace.add('conversion', None, 'pauli_list', pauli_list)
    problem = QuadraticProgram('EPISTASIS_QAOA')
    problem.from_ising(pauli_operator)

    # calculate how many calls will be to the quantum circuit for the QAOA call
    qc_calls_per_qaoa = (1 + 2 * reps) * maxiter + 1  # for each iter 1 global + 1 per parameter; finally a last call
    trace.add('stats', None, 'qc_call_per_qaoa', qc_calls_per_qaoa)
    qc_call_per_rec_qaoa = qc_calls_per_qaoa * (len(h) - 1)  # repeated calls for the given number of variables
    trace.add('stats', None, 'qc_call_per_rec_qaoa', qc_call_per_rec_qaoa)

    # choose optimizer
    optimizer_dict = {'ADAM': ADAM(maxiter=maxiter), 'COBYLA': COBYLA(maxiter=maxiter), 'NELDER_MEAD': NELDER_MEAD(maxiter=maxiter), 'SLSQP': SLSQP(maxiter=maxiter)}
    the_optimizer = optimizer_dict[optimizer]

    # create QAOA callback
    def qaoa_callback(evaluation_count: int, parameters: np.ndarray, eval_mean: float, eval_std: float):
        trace.add('qaoa_trace', f"solver_call_rec_{qaoa_callback.qaoa_solver_calls // qaoa_callback.qc_calls_per_qaoa}_iter_{qaoa_callback.qaoa_solver_calls % qaoa_callback.qc_calls_per_qaoa}", 'eval_count', evaluation_count)
        trace.add('qaoa_trace', f"solver_call_rec_{qaoa_callback.qaoa_solver_calls // qaoa_callback.qc_calls_per_qaoa}_iter_{qaoa_callback.qaoa_solver_calls % qaoa_callback.qc_calls_per_qaoa}", 'eval_parameters', parameters)
        trace.add('qaoa_trace', f"solver_call_rec_{qaoa_callback.qaoa_solver_calls // qaoa_callback.qc_calls_per_qaoa}_iter_{qaoa_callback.qaoa_solver_calls % qaoa_callback.qc_calls_per_qaoa}", 'eval_mean', eval_mean)
        trace.add('qaoa_trace', f"solver_call_rec_{qaoa_callback.qaoa_solver_calls // qaoa_callback.qc_calls_per_qaoa}_iter_{qaoa_callback.qaoa_solver_calls % qaoa_callback.qc_calls_per_qaoa}", 'eval_std', eval_std)
        qaoa_callback.qaoa_solver_calls += 1

    qaoa_callback.qaoa_solver_calls = 0
    qaoa_callback.qc_calls_per_qaoa = qc_calls_per_qaoa

    # create QAOA instance
    quantum_instance = QuantumInstance(backend, shots=n_shots)
    qaoa = QAOA(optimizer=the_optimizer, reps=reps, quantum_instance=quantum_instance, callback=qaoa_callback)
    if is_recursive_qaoa == 0:
        internal_optimizer = None
        optimizer = MinimumEigenOptimizer(qaoa)
    else:
        internal_optimizer = MinimumEigenOptimizer(qaoa)
        optimizer = RecursiveMinimumEigenOptimizer(internal_optimizer)

    # save stats about QAOA circuit length size ecc. on current platform
    params = [Parameter(f'x{i}') for i in range(2*reps)]
    qaoa_qc = qaoa.construct_circuit(params, pauli_operator)[0]
    qaoa_qc_current = qiskit.transpile(qaoa_qc, backend=backend)
    trace.add('stats', 'current_device', 'size', qaoa_qc_current.size())
    trace.add('stats', 'current_device', 'depth', qaoa_qc_current.depth())
    trace.add('stats', 'current_device', 'num_1qubits_op', np.sum([instr.num_qubits == 1 for instr, _, _ in qaoa_qc_current._data]))
    trace.add('stats', 'current_device', 'num_2qubits_op', np.sum([instr.num_qubits == 2 for instr, _, _ in qaoa_qc_current._data]))
    trace.add('stats', 'current_device', 'circuit_qasm', qaoa_qc_current.assign_parameters({params[i]: i for i in range(len(params))}).qasm())

    # run QAOA
    result = optimizer.solve(problem)
    trace.add('output', 'raw', 'result_fval', result.fval)
    trace.add('output', 'raw', 'result_x', result.x)

    # rewrite results in the correct format
    result_variables = [-1 if x == 1 else 1 for x in result.x[::-1]]
    trace.add('output', 'postprocessed', 'result_variables', result_variables)
    trace.add('output', 'postprocessed', 'result_cost', get_ising_cost(h, J, result_variables))

    trace.save()
    return result_variables


def test_run_quantum_universal_computer():
    """Run a basic Ising model to check if the library works correctly"""
    # h = [0, 0]
    # J = {(0, 0): 0, (0, 1): 1, (1, 0): 0, (1, 1): 0.5}
    password_azure = json.load(open("password_azure.json"))
    azure_subscription_id = password_azure["azure_subscription_id"]
    azure_resource_group = password_azure["azure_resource_group"]
    azure_name = password_azure["azure_name"]
    azure_location = password_azure["azure_location"]

    azure_backend = "quantinuum.sim.h1-1e"

    N = 3
    h = np.array([-12.0, 0.5, 0.3])
    J = {(0, 1): 0.2, (0, 2): 0.5, (1, 2): 1.2}
    for vars in itertools.product([-1, 1], repeat=N):
        print(f"vars={str(vars):20s} cost={get_ising_cost(h, J, vars):3.3f}")
    result = run_qaoa('azure', azure_subscription_id, azure_resource_group, azure_name, azure_location, azure_backend,
                      h, J, 'ADAM', 1, 1, 1, 0, 'azurelog')
    return result
