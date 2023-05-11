from dwave.cloud import Client
from dwave.system import DWaveSampler, FixedEmbeddingComposite
from greedy import SteepestDescentSolver
import numpy as np
import datetime
import json
import minorminer
from utils import Trace


def get_rev_annealing_schedule(s_target=0.0, ramp_time=10, pause_time=10):
    """Create schedule for reverse annealing.
    :param s_target: Indicated the point where the annealing stop running in reverse.
    :param ramp_time: It is the time taken to go from 100% to S_TARGET% of the annealing process. The same value indicates the time the annealing takes to go from S_TARGET% to 100%
    :param pause_time: It is the time that the annealing stops at S_TARGET%
    :return annealing (list[[int,int]])
    """
    return [[0, 1.0], [ramp_time, s_target], [ramp_time+pause_time, s_target], [2*ramp_time+pause_time, 1.0]]


def get_annealing_schedule(ramp_time=10, pause_time=0):
    """Create schedule for forward annealing.
    :param ramp_time: It is the time taken to go from 0% to 50% of the annealing process. The same value indicates the time the annealing takes to go from 50% to 100%
    :param pause_time: It is the time that the annealing stops at 50%
    :return annealing (list[[int,int]])
    """
    return [[0, 0.0], [ramp_time, 0.50], [ramp_time+pause_time, 0.50], [2*ramp_time+pause_time, 1.0]]


def run_quantum_annealer(token, h, J, num_reads, solver_idx,
        fw_annealing_ramp_time, fw_annealing_pause_time, 
        rev_annealing_ramp_time, rev_annealing_pause_time, rev_annealing_s_target, 
        save_path):
    """
    Run the quantum annealer
    :param h: Ising formulation h coefficients (np.nparray[float] of shape (N,))
    :param j: Ising formulation J coefficients (dict[tuple(int,int) -> float])
    :param num_reads: number of sampling of the annealer usign each time a different random initial state (https://docs.ocean.dwavesys.com/projects/neal/en/latest/reference/generated/neal.sampler.SimulatedAnnealingSampler.sample.html)
    :param solver_idx: D-Wave solver. Use always 4.
        0 to use 'hybrid_binary_quadratic_model_version2'
        1 to use 'hybrid_discrete_quadratic_model_version1'
        2 to use 'hybrid_constrained_quadratic_model_version1p'
        3 to use 'Advantage_system4.1'
        4 to use 'Advantage_system6.1'
        5 to use 'DW_2000Q_6'
        6 to use 'DW_2000Q_VFYC_6'
    :param fw_annealing_ramp_time: For FORWARD annealing, it is the time taken to go from 0% to 50% of the annealing process. The same value indicates the time the annealing takes to go from 50% to 100%
    :param fw_annealing_pause_time: For FORWARD annealing, it is the time that the annealing stops at 50%
    :param rev_annealing_ramp_time: For REVERSE annealing, it is the time taken to go from 100% to S_TARGET% of the annealing process. The same value indicates the time the annealing takes to go from S_TARGET% to 100%
    :param rev_annealing_pause_time: For REVERSE annealing, it is the time that the annealing stops at S_TARGET%
    :param rev_annealing_s_target: Indicated the point where the annealing stop running in reverse.
    :param save_path: Path where to save the file including the first part of the file name
    :return solution to the qubo problem (list[int] of size N) or [] in case of errors.
    """

    trace = Trace(print=True, path=save_path)
    
    SOLVER_NAMES = [  # see available solvers using client.get_solvers()
        'hybrid_binary_quadratic_model_version2', 
        'hybrid_discrete_quadratic_model_version1', 
        'hybrid_constrained_quadratic_model_version1p',
        'Advantage_system4.1',
        'Advantage_system6.1',
        'DW_2000Q_6',
        'DW_2000Q_VFYC_6']
    
    # trace the input except for the qubo formulation itself which is too large
    trace.add('input', None, 'num_reads', num_reads)
    trace.add('input', None, 'solver_idx', solver_idx)
    trace.add('input', None, 'fw_annealing_ramp_time', fw_annealing_ramp_time)
    trace.add('input', None, 'fw_annealing_pause_time', fw_annealing_pause_time)
    trace.add('input', None, 'rev_annealing_ramp_time', rev_annealing_ramp_time)
    trace.add('input', None, 'rev_annealing_pause_time', rev_annealing_pause_time)
    trace.add('input', None, 'rev_annealing_pause_time', rev_annealing_s_target)
    trace.add('input', None, 'save_path', save_path)
    if h:
        trace.add('input_stats', None, 'h_min', np.min(h))
        trace.add('input_stats', None, 'h_max', np.max(h))
        trace.add('input_stats', None, 'h_mean', np.mean(h))
        trace.add('input_stats', None, 'h_std', np.std(h))
    else:
        trace.add('input_stats', None, 'h_min',  'h empty')
        trace.add('input_stats', None, 'h_max',  'h empty')
        trace.add('input_stats', None, 'h_mean', 'h empty')
        trace.add('input_stats', None, 'h_std',  'h empty')
    j_list = list(J.values())
    if j_list:
        trace.add('input_stats', None, 'J_min', np.min(j_list))
        trace.add('input_stats', None, 'J_max', np.max(j_list))
        trace.add('input_stats', None, 'J_mean', np.mean(j_list))
        trace.add('input_stats', None, 'J_std', np.std(j_list))
    else:
        trace.add('input_stats', None, 'J_min',  'J empty')
        trace.add('input_stats', None, 'J_max',  'J empty')
        trace.add('input_stats', None, 'J_mean', 'J empty')
        trace.add('input_stats', None, 'J_std',  'J empty')

    client = Client.from_config(token=token)
        
    # get the solver from the client object
    solver = client.get_solver(name=SOLVER_NAMES[solver_idx])
    trace.add('configuration', None, 'solver', solver.name)
    trace.add('configuration', None, 'token', 'hidden')

    # create sampler and print its properties
    qpu = DWaveSampler(solver=solver.name, token=token)
    trace.add('qpu_properties', None, 'chip_id', qpu.properties['chip_id'])
    trace.add('qpu_properties', None, 'topology', qpu.properties['topology'])
    trace.add('qpu_properties', None, 'h_range', qpu.properties['h_range'])
    trace.add('qpu_properties', None, 'j_range', qpu.properties['j_range'])
    trace.add('qpu_properties_report', None, 'h_dac_error', [[-5,0.004], [5, 0.0005]])
    trace.add('qpu_properties_report', None, 'J_dac_error', [[-1,0.00225], [1, 0.00025]])

    # print timing properties of the programming cycle and anneal-read cycle
    trace.add('qpu_timing (prog cycle)', None, 'programming_thermalization_range', qpu.properties['programming_thermalization_range'])
    trace.add('qpu_timing (prog cycle)', None, 'default_programming_thermalization', qpu.properties['default_programming_thermalization'])
    trace.add('qpu_timing (anneal-read cycle)', None, 'annealing_time_range', qpu.properties['annealing_time_range'])
    trace.add('qpu_timing (anneal-read cycle)', None, 'default_annealing_time', qpu.properties['default_annealing_time'])
    trace.add('qpu_timing (anneal-read cycle)', None, 'readout_thermalization_range', qpu.properties['readout_thermalization_range'])
    trace.add('qpu_timing (anneal-read cycle)', None, 'default_readout_thermalization', qpu.properties['default_readout_thermalization'])

    # create annealing schedules, both forward and reverse
    fw_anneal_schedule = get_annealing_schedule(fw_annealing_ramp_time, fw_annealing_pause_time)
    rev_anneal_schedule = get_rev_annealing_schedule(rev_annealing_s_target, rev_annealing_ramp_time, rev_annealing_pause_time)
    trace.add('qpu_annealing', 'forward_annealing', 'annealing', fw_anneal_schedule)
    trace.add('qpu_annealing', 'forward_annealing', 'ramp_time', fw_annealing_ramp_time)
    trace.add('qpu_annealing', 'forward_annealing', 'pause_time', fw_annealing_pause_time)
    trace.add('qpu_annealing', 'reverse_annealing', 'annealing', rev_anneal_schedule)
    trace.add('qpu_annealing', 'reverse_annealing', 's_target', rev_annealing_s_target)
    trace.add('qpu_annealing', 'reverse_annealing', 'ramp_time', rev_annealing_ramp_time)
    trace.add('qpu_annealing', 'reverse_annealing', 'pause_time', rev_annealing_pause_time)
    
    # Properties that are not useful for our task
    # * flux_biases
    # * flux_drift_compensation
    # * h_gain_schedule
    # * num_spin_reversal_transforms
    # * programming_thermalization
    # * readout_thermalization
    # * reduce_intersample_correlation

    # find embedding
    start_embedding_time = datetime.datetime.now()
    embedding = minorminer.find_embedding(J.keys(), qpu.properties['couplers'], timeout=5, tries=5, verbose=1)
    end_embedding_time = datetime.datetime.now()
    trace.add('embedding', None, 'time_ms', end_embedding_time - start_embedding_time)

    # check embedding
    if not bool(embedding):
        client.close()
        trace.add('embedding', None, 'exists', False)
        trace.add('error', None, 'error', 'Cannot find an appropriate embedding')
        trace.save()
        return []
    else:
        trace.add('embedding', None, 'exists', True)

    chains = list(embedding.values())
    max_chain_length = np.max([len(chain) for chain in chains])
    trace.add('embedding', None, 'embedding', embedding)
    trace.add('embedding', None, 'max_chain_length', max_chain_length)

    # run forward annealing
    sampler = FixedEmbeddingComposite(qpu, embedding=embedding)
    sampleset = sampler.sample_ising(h, J, num_reads=num_reads, anneal_schedule=fw_anneal_schedule, return_embedding=True)  # run QPU
    best_fw_index = np.argmin(sampleset.record.energy)  # get best energy and samples before preprocessing
    best_fw_energy = sampleset.record.energy[best_fw_index]
    best_fw_sample = sampleset.record.sample[best_fw_index]
    trace.add('fw_annealing', None, 'best_energy', best_fw_energy)
    trace.add('fw_annealing', None, 'best_sample', best_fw_sample)
    for k, v in sampleset.info['timing'].items():
        trace.add('fw_annealing', 'timing', k, v)
    for k, v in sampleset.info['embedding_context'].items():
        if k != 'embedding':
            trace.add('fw_annealing', 'embedding_context', k, v)
    
    # run postprocessing of forward annealing
    sampleset = SteepestDescentSolver().sample_ising(h, J, initial_states=sampleset)  # run postprocessing
    best_fw_index = np.argmin(sampleset.record.energy)  # get best energy and samples before preprocessing
    best_fw_energy = sampleset.record.energy[best_fw_index]
    best_fw_sample = sampleset.record.sample[best_fw_index]
    trace.add('fw_annealing_SteepestDescentSolver', None, 'best_energy', best_fw_energy)
    trace.add('fw_annealing_SteepestDescentSolver', None, 'best_sample', best_fw_sample)
    
    # setup solution of forward annealing as initial state of reverse annealing
    rev_anneal_initial_state = dict(zip(sampleset.variables, best_fw_sample))

    # run QPU with reverse annealing schedule
    sampleset = sampler.sample_ising(h, J, num_reads=num_reads, 
        anneal_schedule=rev_anneal_schedule, 
        initial_state=rev_anneal_initial_state,
        reinitialize_state=True,
        return_embedding=True)  
    best_rev_index = np.argmin(sampleset.record.energy)
    best_rev_energy = sampleset.record.energy[best_rev_index]
    best_rev_sample = sampleset.record.sample[best_rev_index]
    trace.add('rev_annealing', None, 'best_energy', best_rev_energy)
    trace.add('rev_annealing', None, 'best_sample', best_rev_sample)
    for k, v in sampleset.info["timing"].items():
        trace.add('rev_annealing', 'timing', k, v)
    for k, v in sampleset.info["embedding_context"].items():
        if k != 'embedding':
            trace.add('rev_annealing', 'embedding_context', k, v)
    
    # run postprocessing of reverse annealing
    sampleset = SteepestDescentSolver().sample_ising(h, J, initial_states=sampleset)  # run postprocessing
    best_rev_index = np.argmin(sampleset.record.energy)  # get best energy and samples before preprocessing
    best_rev_energy = sampleset.record.energy[best_rev_index]
    best_rev_sample = sampleset.record.sample[best_rev_index]
    trace.add('rev_annealing_SteepestDescentSolver', None, 'best_energy', best_rev_energy)
    trace.add('rev_annealing_SteepestDescentSolver', None, 'best_sample', best_rev_sample)

    client.close()
    trace.save()

    if best_fw_energy > best_fw_energy:
        return best_fw_sample if type(best_fw_sample) == 'list' else best_fw_sample.tolist()
    else:
        return best_rev_sample if type(best_rev_sample) == 'list' else best_rev_sample.tolist()


def test_run_quantum_annealer(the_token):
    """Run a basic Ising model to check if the library works correctly"""
    # h = [0, 0]
    # J = {(0, 0): 0, (0, 1): 1, (1, 0): 0, (1, 1): 0.5}
    N = 4
    h = np.zeros(shape=(N,))
    J = {(i, j): 0.1 for i in range(N) for j in range(N)}
    best_sample = run_quantum_annealer(h, J, the_token, num_reads=100, solver_idx=4,
        fw_annealing_ramp_time=1, fw_annealing_pause_time=5,
        rev_annealing_ramp_time=1, rev_annealing_pause_time=5, rev_annealing_s_target=0.5,
        save_path='dwavelog')
