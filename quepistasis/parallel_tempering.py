import numpy as np
from abc import ABC, abstractmethod
from utils import Trace


class ParallelTempering(ABC):

    def __init__(self, num_chains, num_steps, betas):
        """
        Initialize the ParallelTempering instance.

        Parameters:
        - num_chains (int): Number of parallel chains.
        - num_steps (int): Number of MCMC steps for each chain.
        - betas (array-like): Temperatures for each chain.
        """
        self.num_chains = num_chains
        self.num_steps = num_steps
        self.betas = betas
        self.chains = None
        self.log_likes = None

    @abstractmethod
    def initial_probability(self):
        """
        Abstract method to provide the initial probability configuration for each chain.
        """
        pass

    @abstractmethod
    def log_like(self, x):
        """
        Abstract method to provide the log likelihood function for a given configuration.

        Parameters:
        - x (array-like): Configuration for which to calculate the log likelihood.

        Returns:
        - float: Log likelihood.
        """
        pass

    @abstractmethod
    def log_prior(self, x):
        """
        Abstract method to provide the log prior function for a given configuration.

        Parameters:
        - x (array-like): Configuration for which to calculate the log prior.

        Returns:
        - float: Log prior.
        """
        pass

    def metropolis_hastings_proposal(self, x, beta):
        """
        Metropolis-Hastings proposal for binary variables.

        Parameters:
        - x (array-like): Current configuration.
        - beta (float): Inverse temperature for the current chain.

        Returns:
        - array-like: Proposed configuration.
        """
        proposal = np.copy(x)
        flip_index = np.random.randint(len(x))
        proposal[flip_index] *= -1
        return proposal

    def acceptance_probability(self, old_prob, new_prob, beta):
        """
        Calculate the acceptance probability for the Metropolis-Hastings step.

        Parameters:
        - old_prob (float): Log probability of the current configuration.
        - new_prob (float): Log probability of the proposed configuration.
        - beta (float): Inverse temperature for the current chain.

        Returns:
        - float: Acceptance probability.
        """
        return min(1, np.exp(beta * (old_prob - new_prob)))

    def run(self):
        """
        Run the Parallel Tempering MCMC algorithm.

        Returns:
        - tuple: Chains and log likelihoods.
        """
        ndim = len(self.initial_probability())
        self.chains = np.zeros((self.num_chains, self.num_steps, ndim))
        self.log_likes = np.zeros((self.num_chains, self.num_steps))

        for i in range(self.num_chains):
            self.chains[i, 0] = self.initial_probability()
            self.log_likes[i, 0] = self.log_like(self.chains[i, 0]) + self.log_prior(self.chains[i, 0])

        for step in range(1, self.num_steps):
            for i in range(self.num_chains):
                current_chain = self.chains[i, step - 1]
                current_log_prob = self.log_like(current_chain) + self.log_prior(current_chain)

                # Perform Metropolis-Hastings proposal
                proposed_chain = self.metropolis_hastings_proposal(current_chain, self.betas[i])

                # Calculate log likelihood and log prior for proposed chain
                proposed_log_prob = self.log_like(proposed_chain) + self.log_prior(proposed_chain)

                # Accept or reject the proposed chain
                if np.log(np.random.uniform()) < self.acceptance_probability(current_log_prob, proposed_log_prob, self.betas[i]):
                    self.chains[i, step] = proposed_chain
                    self.log_likes[i, step] = proposed_log_prob
                else:
                    self.chains[i, step] = current_chain
                    self.log_likes[i, step] = current_log_prob

            # Exchange states between adjacent chains
            for i in range(self.num_chains - 1):
                diff_prob = self.betas[i + 1] * (self.log_likes[i, step] - self.log_likes[i + 1, step])
                if np.log(np.random.uniform()) < diff_prob:
                    # Swap states between chains i and i+1
                    self.chains[i, step], self.chains[i + 1, step] = self.chains[i + 1, step], self.chains[i, step]
                    self.log_likes[i, step], self.log_likes[i + 1, step] = self.log_likes[i + 1, step], self.log_likes[i, step]

        return self.chains, self.log_likes



class IsingParallelTempering(ParallelTempering):

    def __init__(self, h, J, num_chains, num_steps, betas):
        """
        Initialize the IsingParallelTempering instance.

        Parameters:
        - h (numpy array): External magnetic field.
        - J (dictionary): Coupling strengths between spins (interaction matrix).
        - num_chains (int): Number of parallel chains.
        - num_steps (int): Number of MCMC steps for each chain.
        - betas (array-like): Temperatures for each chain.
        """
        super().__init__(num_chains, num_steps, betas)
        self.h = h
        self.J = J
        self.ising_upper = np.sum(np.abs(self.h)) + np.sum(np.abs(np.array(list(J.values()))))
        self.ndim = len(h)

    def initial_probability(self):
        """
        Generate the initial spin configuration.

        Returns:
        - array-like: Initial spin configuration (+1 or -1).
        """
        return np.random.choice([-1, 1], size=self.ndim)
    
    @staticmethod
    def ising_energy(x, h, J):
        energy = np.sum(h * x)
        for (i, j), value in J.items():
            energy += value * x[i] * x[j]
        return energy

    def log_like(self, x):
        """
        Calculate the log likelihood for a given spin configuration.

        Parameters:
        - x (array-like): Spin configuration for which to calculate the log likelihood.

        Returns:
        - float: Log likelihood
        """
        energy = self.ising_energy(x, self.h, self.J)
        return -1 * (self.ising_upper - energy)

    def log_prior(self, x):
        """
        Calculate the log prior for a given spin configuration.

        Parameters:
        - x (array-like): Spin configuration for which to calculate the log prior.

        Returns:
        - float: Log prior (always 0.0 in this example).
        """
        return 0.0
    
    def get_solution(self):
        """
        Return the solution of the optimization problem. 
        
        Returns:
        - np.ndarray: spin configuration of the best solution.
        - float: Log likelihood of the best solution.
        - float: Energy of the best solution.
        """

        best_chain_index, best_step_index = np.unravel_index(np.argmin(self.log_likes), self.log_likes.shape)
        best_configuration = self.chains[best_chain_index, best_step_index]
        best_log_likelihood = self.log_likes[best_chain_index, best_step_index]
        best_energy = self.ising_energy(best_configuration, self.h, self.J)
        return best_configuration, best_log_likelihood, best_energy


def run_parallel_tempering(h, J, num_chains, num_steps, save_path):

    # trace the input except for the qubo formulation itself which is too large
    trace = Trace(print=True, path=save_path)
    trace.add('input', None, 'num_chains', num_chains)
    trace.add('input', None, 'num_steps', num_steps)
    trace.add('input', None, 'save_path', save_path)
    betas = np.geomspace(1, 1e-2, num_chains)
    ising_pt = IsingParallelTempering(h, J, num_chains, num_steps, betas)
    ising_pt.run()
    spins, loglike, energy = ising_pt.get_solution()
    trace.add('solution', None, 'spins', spins)
    trace.add('solution', None, 'loglike', loglike)
    trace.add('solution', None, 'energy', energy)
    return spins.tolist()

def test_parallel_tempering():
    """
    Run a basic Ising model to check if the code works correctly.
    """
    h = np.array([5, 10, -20])
    J = {(0, 1): 1, (1, 2): 2}

    num_chains = 4
    num_steps = 1000
    betas = np.geomspace(1, 1e-2, num_chains)

    ising_pt = IsingParallelTempering(h, J, num_chains, num_steps, betas)
    ising_pt.run()
    spins, loglike, energy = ising_pt.get_solution()
    print(f"{spins=} {loglike=} {energy=}")
