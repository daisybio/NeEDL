import numpy as np
from pyqubo import Array
import networkx as nx
import neal
from dimod.utilities import qubo_to_ising
from dimod.sampleset import as_samples, infer_vartype, SampleSet

# define the SNP-SNP correlation matrix (assume both statistical and biological information is included)

N_SNPS = 4
graph_matrix = np.array([
    [0, 1, 2, 3],
    [0, 0, 4, 5],
    [0, 0, 0, 6],
    [0, 0, 0, 0]
])
assert graph_matrix.shape == (N_SNPS, N_SNPS)

# define target of the optimization

N_CLIQUES = 2
K = 3

# define the hyper-parameters

LAMBDA0 = 10.0
LAMBDA1 = 0.456
LAMBDA2 = 0.789

# define SA parameters
N_READS = 10
SEED = 1234

# run the test
graph = nx.from_numpy_matrix(graph_matrix, create_using=nx.DiGraph)

def print_expression_as_qubo(expression, N):
    qubo_dict, _ = expression.compile().to_qubo()
    for i in range(N):
        for j in range(N):
            w = qubo_dict.get((f"x[{i}]", f"x[{j}]"), 0)
            print(f"{w: 3.3f} ", end="")
        print()

def construct_n_max_weighted_k_clique(graph, n, k, lambda0, lambda1, lambda2):
    """
    Formulate the problem as WEIGHTED MAX-K-CLIQUE
    The formulation for WEIGHTED MAX-K-CLIQUE is:
    𝑄=∑_(ℓ=1)^𝑁[𝜆_0 (∑_𝑖 𝑥_𝑖ℓ −𝐾)^2-𝜆_1 ∑_((𝑖,𝑗)∉𝐸)〖𝑤_(𝑖ℓ,𝑗ℓ)⋅𝑥_𝑖ℓ 𝑥_𝑗ℓ 〗]+𝑑𝑖𝑠𝑠𝑖𝑚 𝑡𝑒𝑟𝑚
    𝑑𝑖𝑠𝑠𝑖𝑚 𝑡𝑒𝑟𝑚 = 𝜆_2 ∑_(ℓ=1)^𝑁 ∑_(𝑚=ℓ+1)^𝑁 ∑_𝑖^(|𝑉|) ∑_𝑗^(|𝑉|) (𝑥_𝑖ℓ⋅𝑥_𝑗𝑚 )
    :param graph: Graph of SNPs
    :param n: number of cliques
    :param k: size of the clique
    :param lambda0 weight of the first term
    :param lambda1 weight of the second term
    :return: PyQUBO formulation of the max-clique problem
    """
    v = len(graph.nodes)
    x = Array.create('x', shape=(n*v,), vartype='BINARY')
    expression = 0

    # for each of the n cliques
    for i in range(n):
        # first add constraint
        x_row = list(x)[i*v:(i+1)*v]
        expression += lambda0 * (sum(xi for xi in x_row) - k)**2
        # then add reward term for large weight
        for edge in list(graph.edges):
            expression += -lambda1 * graph.edges[edge]['weight'] * x[(i*v) + edge[0]] * x[(i*v) + edge[1]]

    # finally, add dissimilarity term
    for l in range(n):
        for m in range(l+1, n):
            for i in range(v):
                expression += lambda2 * x[(l * v) + i] * x[(m * v) + i]
    return expression


print("QUBO generation:")
expression = construct_n_max_weighted_k_clique(graph, N_CLIQUES, K, LAMBDA0, LAMBDA1, LAMBDA2)
print_expression_as_qubo(expression, N_CLIQUES * N_SNPS)
model = expression.compile()
qubo_bqm, offset = model.to_qubo()
print("Offset:", offset)


print("\n\n\nQUBO to Ising:")
binary_quadratic_model = model.to_bqm()
h, J, offset = qubo_to_ising(qubo_bqm, offset)

h = list(h.items())
get_number_between_brackets = lambda s: int(s.split('[')[-1].split(']')[0])
h.sort(key=lambda item: get_number_between_brackets(item[0]))

print("h index: ", end="")
for hi in h:
    print(f"{hi[0]:5s}", end=" ")
print("\nh value: ", end="")
for hi in h:
    print(f"{hi[1]:4.3f}", end=" ")
print()

J = list(J.items())
J.sort()

print("J s_idx: ", end="")
for ji in J:
    print(f"{ji[0][0]:5s}", end=" ")
print("\nJ e_idx: ", end="")
for ji in J:
    print(f"{ji[0][1]:5s}", end=" ")
print("\nJ value: ", end="")
for ji in J:
    print(f"{ji[1]:4.3f}", end=" ")
print()

print(f"offset: {offset:3.5f}")


model = expression.compile()
bqm = model.to_bqm()
sa = neal.SimulatedAnnealingSampler()

print("\nVariables: ", list(bqm.variables))
variables = list(bqm.variables)
variables.sort(key=get_number_between_brackets)

initial_state = [+1, +1, +0, +0, +0, +1, +0, +1, +0, +0, +1, +1, +0, +0, +1, +1, +1, +0, +1, +1, +1, +1, +0, +0, +1, +0, +1, +0, +1, +1, +1, +0, +0, +1, +1, +1, +0, +1, +0, +0, +0, +1, +0, +0, +0, +1, +1, +0, +0, +0, +1, +0, +1, +1, +1, +1, +1, +0, +0, +1, +0, +1, +0, +0, +0, +0, +0, +1, +1, +1, +1, +0, +0, +0, +0, +1, +1, +1, +1, +1]
initial_state_wrap = []
for i in range(N_READS):
    init_state_i = {}
    for j in range(N_SNPS * N_CLIQUES):
        variable = variables[j] # 'x[k]'
        k = get_number_between_brackets(variable)
        init_state_i[variable] = initial_state[i * N_SNPS * N_CLIQUES + k]
    initial_state_wrap.append(init_state_i)


# initial_state = np.array(initial_state).reshape((N_READS, 8))
print("initial_state:", initial_state_wrap)
initial_states_array, initial_states_variables = as_samples(initial_state_wrap, copy=True)
print("initial_states_array:", initial_states_array)
print("initial_states_variables:", initial_states_variables)
print(bqm.variables ^ initial_states_variables)

sampleset = sa.sample(bqm, num_reads=N_READS, seed=SEED, initial_states=initial_state_wrap)
decoded_samples = model.decode_sampleset(sampleset)

print("\n\nIsing energies: ", [ds.energy for ds in decoded_samples])

best_sample = min(decoded_samples, key=lambda x: x.energy)

