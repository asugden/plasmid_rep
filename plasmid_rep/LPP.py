from datetime import datetime
import logging
import numpy as np
import os.path
from typing import Any, Dict
import yaml


def crp(n: int, alpha: float):
    """Implementation of the Chinese Restaurant Process

    Alpha of 0.1 leads a cluster of 10 to split up into 1.3 clusters of 8.7 plasmids
    Alpha of 10 leads a cluster of 10 to split up into 7.0 clusters of 1.5 plasmids
    Alpha of >= 100 leads to complete breakup of clusters

    Args:
        n (int): the size of a group to be dissociated/cluster
        alpha (float): the dissociation parameter

    Returns:
        list: a list of "tables"/clusters of plasmids after dissociation

    """
    # Return guaranteed results if possible
    if n == 0: 
        return []
    if n == 1: 
        return [1]
    if alpha >= 100: 
        return [1]*n

    # Initialize random comparisons and return vector
    comp = np.random.random(n - 1)
    out = np.zeros(n, dtype=np.uint16)
    out[0] = 1
    mxt = 0

    # Iterate over all ns greater than 0
    for i in range(1, n):
        prob = _crp_occupied_probability(out[0], i + 1, alpha) # 1-indexed
        table = 0
        while prob < comp[i - 1] and table < mxt + 1:
            prob += _crp_occupied_probability(out[table], i+1, alpha)
            table += 1

        out[table] += 1
        if table > mxt: 
            mxt = table

    return out[0:mxt+1].tolist()


def _crp_occupied_probability(t: int, n: int, alpha: float):
    """Probability of joining a table of size t if it is occupied, given alpha

    Args:
        t (int): the probability of joining a table of size t
        n (int): the size of a group to be dissociated
        alpha (float): aggregation parameter

    Returns:
        float: probability of joining table of size t
    """
    return float(t)/(n - 1 + alpha)


class LatentPlasmidPopulation():
    # Call with population(), simulate(), save()
    def __init__(self):
        # If you decrease the right-tail parameters, you may have to increase this value
        self.maximum_cluster_n = 2000
        # The number of cells to be simulated and sampled to
        self.sample_to_n = None

        # S-phase probability of duplication and probability of equal partitioning
        # These are 0.88 and 0.84, respectively, for EBV, 0.92 for KSHV
        self.s_phase_duplication_prob = 0.92

        # Repulsion-attraction parameter
        # In Nanbo, Sugden, and Sugden 2006, we set the probability of equal partitioning equal to
        # 0.84. At the time, we only simulated this effect during mitosis. Now, to include both
        # attraction and repulsion, we have renamed and inverted the parameter. It is (1 - equal
        # partitioning). For EBV, this should be set to 0.12 (if s-vs-g1 is set to 0, otherwise
        # it should be set to 0.24. For KSHV, it should be set to 1, 0.77 for the fusion protein.

        # Repulsion               Attraction
        # Equal partitioning      Clustering
        # 0 ------------------------------ 1
        self.plasmid_repulsion_attraction = 0.24

        # Right-tail parameters

        # Probability that the average cell will divide each generation and per-plasmid selective
        # disadvantage for cell replication = a - bx - cx^2 where a is average-cell-replication-
        # prob, b is signal-selective-disadvantage-on-cell-replication and c is signal-selective-
        # disadvantage-on-cell-replication-squared. Increase maximum-cluster-n if these values are
        # decreased.
        self.average_cell_replication_prob = 1.0
        self.signal_selective_disadvantage_on_cell_replication = 0
        self.signal_selective_disadvantage_on_cell_replication_squared = 0.00000625

        # Left-tail parameters

        # If clusters "attract", they are clustered until broken up. Cluster dissociation is set
        # by two parameters: the first denotes the probability of some event jostling the cluster
        # and the second is the parameter of the Chinese Restaurant Process that describes the
        # breakup.
        # Alpha of 0.1 leads a cluster of 10 to split up into 1.3 clusters of 8.7 plasmids
        # Alpha of 10 leads a cluster of 10 to split up into 7.0 clusters of 1.5 plasmids
        # Alpha of >= 100 leads to complete breakup of clusters

        # The values from Condor are 0.8, 0.5 after 200,000 compute hours
        # Set to 1.0, 101 to turn off clustering
        self.cluster_jostling = 1.0
        self.cluster_crp_alpha = 101

        # Just in case the probabilities of jostling should be different during the end of s-phase
        # and the beginning of G1, we can set the ratio of breakup probability between them
        # Fusion protein should have this set to 1.0
        self.cluster_jostling_s_vs_g1 = 1.0

        # Positive selection affecting time of replication
        # We will moder this as a sigmoid, and rather than tracking the replication time of every individual
        # cell, we will model replication as a probabilistic process. This also injects more realistic
        # randomness into the model.
        # In particular, the model uses a sigmoid that begins at 0,0: 2(e^(ax)/(1 + e^(ax))) - 1
        # The value 0 removes this selection
        self.positive_selection_coefficient = 0

        # Negative selection specific to clusters
        # We will model negative selection as a decrease in replication probability proportional
        # to the distance from the endogenous chromatin. This will be modeled as exponential decay to use
        # the minimum number of parameters. A cluster will be modeled as a chain. That plasmid
        # directly connected to the chromosome will have the replication probability set by s_phase_duplication_prob
        # As the distance of the chain grows, the probability will decrease as a random number less than
        # s_phase_duplication_prob*e^(-neg_coefficient*x)
        # negative_cluster_selection_coefficient should be positive-- it will be turned negative elsewhere
        self.negative_cluster_selection_coefficient = 0

        # Parameter variables for saving
        self._parameters = {}

        # -------------------------------------------------------------------------------------- #
        # Variables to be used:
        self.n = None
        self.cells = None
        self.zeros = None
        self.n_plasmids_per_cell = None
        self.n_signals_per_cell = None
        self.selection_passed_partner = None

        self.cluster_formation_rates = []
        self.cluster_dissociation_rates = []

        self._total_parents = 0
        self._sample_just_run = False
        self._config = None

    def load_config(self):
        """Load in a config.yaml file describing default parameters"""
        if self._config is None:
            lib_base, _ = os.path.split(os.path.abspath(__file__))
            cfg_path = os.path.join(lib_base, 'config.yaml')
            with open(cfg_path) as fo:
                self._config = yaml.load(fo)
        
        return self._config

    def set_virus(self, virus: str):
        """Set the parmaeters to match those of a virus

        Args:
            virus (str): name of the virus

        """
        assert(virus in ['kshv', 'ebv'])

        self.load_config()
        self.update_parameters(self._config[virus])

    def update_parameters(self, pars: Dict[str, Any]):
        """Update simulation parameters

        Args:
            pars (dict): parameters to update
        
        """
        available_parameters = ['s_phase_duplication_prob', 
                                'plasmid_repulsion_attraction', 
                                'average_cell_replication_prob',
                                'signal_selective_disadvantage_on_cell_replication',
                                'signal_selective_disadvantage_on_cell_replication_squared',
                                'cluster_jostling',
                                'cluster_crp_alpha',
                                'cluster_jostling_s_vs_g1',
                                'positive_selection_coefficient',
                                'negative_cluster_selection_coefficient']
        for key in pars:
            if pars[key] is not None:
                if key not in available_parameters:
                    raise ValueError(f'Variable {key} not in parameters that can be set')
                setattr(self, key, pars[key])

    def population(self, n: int, mu: float = None, sd: float = None, lambda_: float = None):
        """Generate a normally distributed or poisson distributed initial population of cells containing viral particles

        Args:
            n (int): the number of cells to sample
            mu (float): the mean number of plasmids per cell, if using normal distribution
            sd (float): the standard deviation of plasmids per cell, if using normal distribution
            lambda_ (float): the multiplicity of infection per cell, if using poisson distribution

        Returns:
            np.ndarray: the percentage of cells with N plasmids where N is the position in the array

        """
        assert((mu is not None and sd is not None) or lambda_ is not None)

        # Set the number of cells, n
        self.sample_to_n = n
        self.n = self.sample_to_n
        self._sample_just_run = False

        # Create the cell matrix (rows are cells, members of rows are cluster values) and set the counts per line = 0
        self.cells = np.zeros((self.n*2, self.maximum_cluster_n), dtype=np.uint16)
        self.n_plasmids_per_cell = np.zeros(self.n*2, dtype=np.uint16)
        self.n_signals_per_cell = np.zeros(self.n*2, dtype=np.uint16)

        # SELECTION-specific parameters, set to -1 for the burn in stage
        self.selection_passed_partner = np.zeros(self.n*2, dtype=np.int32) - 1

        # Set the number of cells with 0 plasmids to 0. They are simulated separately
        self.zeros = 0

        # Set the number of clusters (of size 1) per cell for each cell and add to cells
        # Set the number of clusters (of size 1) per cell for each cell and add to cells
        if mu is not None and sd is not None:
            plasnums = np.around(np.random.randn(self.n)*sd + mu).astype(int)
            self._total_parents = self.n
        else:
            assert(lambda_ > 0)
            plasnums = []
            while len(plasnums) < self.n:
                plasmids_round = np.random.poisson(lambda_, self.n)
                plasnums = np.concatenate([plasnums, plasmids_round[plasmids_round > 0]])
                self._total_parents += self.n

            self._total_parents -= len(plasnums) - self.n
            plasnums = plasnums[:self.n].astype(int)
            self.zeros = 1.0 - float(len(plasnums))/self._total_parents

        self.n_plasmids_per_cell[:self.n] = plasnums
        self.n_signals_per_cell[:self.n] = np.copy(plasnums)

        for i in range(self.n): 
            self.cells[i, :plasnums[i]] = 1

        return self._current_signals_per_cell()

    def simulate_burn_in(self,
                         generations: int,
                         zeros_after_burn_in: float = 0.0):
        """Run the simulation without saving any variables.

        Args:
            burn_in_generations (int, optional): The number of generations to "burn-in" 
                the simulation, i.e. to come to equilibrium, fter which zeroes will be rest 
                and the true simulationw ill be run. Defaults to 0.
            zeros_after_burn_in (float, optional): the number of cells with zero plasmids
                after burning in. Defaults to 0.0. If you are comparing to real populations 
                and the selection wasn't complete at time 0, set this value to a value 
                0 < zeros-after-burn-in < 1. For normal simulations, 0, for wild-type LANA1, 
                0.1675478927, for ebv, 0.1756, for fusion protein 0.2244897959

        Returns:
            bool: the success of the burn in simulation

        """
        assert(self.n is not None)
        # done is used to check if there was an error during part of the generation
        done = True

        # Run the simulation to come to equilibrium for burn-in-generations 
        logging.info('Burning in...')
        for i in range(generations):
            logging.info('\r\tGen %i/%i', i+1, generations)
            if i > 0: self.sample()

            # Run S & M phases, make sure that they pass successfully
            # S & M phases are combined to account for shared cell fate across mitosis
            done = self.s_m_phase()
            if not done: return False

            # Run G1 phase, make sure that it has passed successfully
            done = self.g1phase()
            if not done: return False

            # Sample, set self.zeros to 0 (no cells with 0 plasmids), and save population
            self.sample()

        # Reset zeros because we assume that the population was under selection before
        # the simulation was begun
        self.zeros = zeros_after_burn_in
        self._sample_just_run = False

        self._parameters['burn_in_generations'] = generations
        self._parameters['zeros_after_burn_in'] = zeros_after_burn_in
        return True

    def simulate(self, generations: int):
        """Run simulation of plasmid replication

        Args:
            generations (int): the number of generations to simulate.

        Returns:
            bool: the success of the simulation

        """
        assert(self.n is not None)
        
        # done is used to check if there was an error during part of the generation
        done = True

        # Run the real simulation for generations, sampling each generation
        logging.info('Simulating...')
        for i in range(generations):
            logging.info('\r\tGen %i/%i', i+1, generations)
            if i > 0:
                self.sample()

            # Run S & M phases
            done = self.s_m_phase()
            if not done: 
                return False

            # Run G1 phase
            done = self.g1phase()
            if not done: 
                return False

        self._parameters['generations'] = generations
        return done

    def save(self, path: str, additional_information: Dict[str, Any]) -> str:
        """Save the output of a run and its associated parameters

        Args:
            path (str): the path into which to save the run
            additional_information (Dict[str, Any]): additional keys and values to be saved
        
        Returns:
            str: the date stamp of the file
        
        """
        assert(self.cells is not None)
        date_stamp = datetime.now().strftime('%y%m%d-%H%M%S')
        for key in vars(self).keys():
            if key not in ['cells', 'n_plasmids_per_cell', 'n_signals_per_cell', 'selection_passed_partner', 
                           'cluster_formation_rates', 'cluster_dissociation_rates'] and key[0] != '_':
                self._parameters[key] = getattr(self, key)
                if type(self._parameters[key]).__module__ == np.__name__:
                    self._parameters[key] = self._parameters[key].item()

        self._parameters['plasmids_per_cell'] = self.get_plasmids_per_cell().tolist()
        self._parameters['plasmids_per_signal'] = self.get_plasmids_per_signal().tolist()
        self._parameters['signals_per_cell'] = self.get_signals_per_cell().tolist()
        self._parameters['date'] = date_stamp

        for key in additional_information:
            if key in self._parameters:
                print(f'Adding value {key} would overwrite a stored parameters. Skipping')
            else:
                self._parameters[key] = additional_information[key]

        with open(os.path.join(path, f'plasmid sim run {date_stamp}.yaml'), 'w') as fp:
            fp.write(yaml.dump(self._parameters))

        return date_stamp

    # ================================================================================== #
    # BIOLOGICAL FUNCTIONS
    # Biological simulation
    def s_m_phase(self):
        """Combined S-phase and M-phase duplicates plasmids, signals, and cells and divvies 
        up the resulting signals to daughter cells.
        S and M phases are combined here because we have to keep track of clusters to be
        partitioned separately
        Note: Ensures that all 0 values are at the end of each row

        We want to determine which cells are duplicated by asking whether a polynomial
        function, dependent on the number of plasmids in the cell, is less than a
        random number.

        Equation is c - a*x^2 - b*x

        Returns:
            bool: the success of the iteration of s- and m- phases

        """
        x = self.n_plasmids_per_cell[:self.n]
        c = self.average_cell_replication_prob

        # Determine which cells have duplicated
        # Account for positive selection
        if self.positive_selection_coefficient == 0:
            a = self.signal_selective_disadvantage_on_cell_replication_squared
            b = self.signal_selective_disadvantage_on_cell_replication
            cell_dup = -a*x*x + -b*x + c > np.random.random(self.n) # Checked the cell_dup is not 100%
        elif self.positive_selection_coefficient > 0:
            a = np.abs(self.positive_selection_coefficient)
            cell_dup = 2*c*(np.exp(a*x)/(1 + np.exp(a*x))) - 1 > np.random.random(self.n)
        else:
            raise ValueError('Positive selection coefficient must be positive')

        duplicated_cells = np.where(cell_dup > 0)
        logging.info('Percent of cells unduplicated: %.2f', 100.0 - 100*float(np.sum(cell_dup))/self.n)
        if len(duplicated_cells[0]) == 0: 
            return self.fail()

        # For each duplicated cell, duplicate plasmids
        for i in np.nditer(duplicated_cells):
            # Check if there are any plasmids in the cell, just in case
            if self.n_signals_per_cell[i] > 0:
                # Duplicate the plasmids and join aggregates in s phase
                if self.negative_cluster_selection_coefficient == 0:
                    clusters = self.sphase(i)
                else:
                    clusters = self.sphase_negative_selection(i)

                # Split the clusters if necessary
                clusters = self.s_phase_cluster_breakup(clusters)

                # Split the clusters in G2-phase and distribute
                done = self.mphase(i, clusters)

                if not done:
                    logging.info('\tWARNING: Hitting maximum cluster boundary in S-phase')
                    return False
            self.n += 1

        return True

    def sphase(self, i: int):
        """Simulate S-phase for a single cell and return duplicated aggregates and split cells.
        Duplicate plasmids and add to clusters
        Generate a list of random numbers to check for plasmid duplication of length
        n_plasmids_per_cell
        Gives a list of 0s and 1s which can be summed

        Args:
            i (int): the cell to run s-phase on
        
        Returns:
            np.ndarray: an array of the sizes of the clusters after s-phase

        """
        plas_dup = self.s_phase_duplication_prob > np.random.random(self.n_plasmids_per_cell[i])

        # Generate an array for each cluster duplicate and duplicate clusters
        clusters = np.zeros((2, self.n_signals_per_cell[i]), dtype=np.uint16)
        tally = 0
        for j in range(self.n_signals_per_cell[i]):
            clusters[0, j] = self.cells[i, j]
            clusters[1, j] = np.sum(plas_dup[tally:tally + self.cells[i, j]])
            tally += self.cells[i, j]

        # Note: clusters[0] has all of the original clusters, clusters[1] has the
        # duplicated versions which, by definition, are <= clusters[0] in size
        plas_aggregation = self.plasmid_repulsion_attraction > np.random.random(self.n_signals_per_cell[i])
        clusters[0, plas_aggregation] += clusters[1][plas_aggregation]
        clusters[1, plas_aggregation] = 0

        return clusters

    def sphase_negative_selection(self, i: int):
        """Simulate S-phase for a single cell and return duplicated aggregates and split cells.
        Duplicate plasmids and add to clusters
        Generate a list of random numbers to check for plasmid duplication of length
        n_plasmids_per_cell
        Gives a list of 0s and 1s which can be summed

        Args:
            i (int): the cell to run s-phase on
        
        Returns:
            np.ndarray: an array of the sizes of the clusters after s-phase

        """
        plas_dup = np.random.random(self.n_plasmids_per_cell[i])

        # Generate an array for each cluster duplicate and duplicate clusters
        clusters = np.zeros((2, self.n_signals_per_cell[i]), dtype=np.uint16)
        tally = 0
        for j in range(self.n_signals_per_cell[i]):
            clusters[0, j] = self.cells[i, j]
            dup_probs = self.s_phase_duplication_prob*np.exp(-self.negative_cluster_selection_coefficient*self.cells[i,j])*np.ones(self.cells[i,j])
            clusters[1, j] = np.sum(dup_probs > plas_dup[tally:tally + self.cells[i, j]])
            tally += self.cells[i, j]

        # Note: clusters[0] has all of the original clusters, clusters[1] has the
        # duplicated versions which, by definition, are <= clusters[0] in size
        plas_aggregation = self.plasmid_repulsion_attraction > np.random.random(self.n_signals_per_cell[i])
        clusters[0, plas_aggregation] += clusters[1][plas_aggregation]
        clusters[1, plas_aggregation] = 0

        return clusters

    def sphase_negative_selection_chain(self, i: int):
        """Simulate S-phase for a single cell and return duplicated aggregates and split cells.
        Assumes that the probability of replication is proportional PER PLASMID to its permission in a chain.
        Duplicate plasmids and add to clusters
        Generate a list of random numbers to check for plasmid duplication of length
        n_plasmids_per_cell
        Gives a list of 0s and 1s which can be summed

        Args:
            i (int): the cell to run s-phase on
        
        Returns:
            np.ndarray: an array of the sizes of the clusters after s-phase

        """
        plas_dup = np.random.random(self.n_plasmids_per_cell[i])

        # Generate an array for each cluster duplicate and duplicate clusters
        clusters = np.zeros((2, self.n_signals_per_cell[i]), dtype=np.uint16)
        tally = 0
        for j in range(self.n_signals_per_cell[i]):
            clusters[0, j] = self.cells[i, j]
            dup_probs = self.s_phase_duplication_prob*np.exp(-self.negative_cluster_selection_coefficient*np.arange(self.cells[i,j]))
            clusters[1, j] = np.sum(dup_probs > plas_dup[tally:tally + self.cells[i, j]])
            tally += self.cells[i, j]

        # Note: clusters[0] has all of the original clusters, clusters[1] has the
        # duplicated versions which, by definition, are <= clusters[0] in size
        plas_aggregation = self.plasmid_repulsion_attraction > np.random.random(self.n_signals_per_cell[i])
        clusters[0, plas_aggregation] += clusters[1][plas_aggregation]
        clusters[1, plas_aggregation] = 0

        return clusters

    def mphase(self, i: int, clusters: np.ndarray):
        """Jostle clusters before mitosis, then assign signals to daughter cells.
        
        Args:
            i (int): the cell in which to run m-phase
            cluster (np.ndarray): the current state of the clusters in the cell
            
        Returns:
            bool: the success of the m-phase process
            
        """
        # Reset cell rows
        self.cells[i, :] = 0

        # Figure out which cells each go to
        assign = np.random.randint(2, size=clusters[0].size)

        # Set counters for position in daughter cells 1 and 2
        c1, c2 = 0, 0

        # Iterate over clusters and assign
        for j in range(clusters[0].size):
            if clusters[assign[j], j] != 0:
                self.cells[i, c1] = clusters[assign[j], j]
                c1 += 1
            if clusters[1 - assign[j], j] != 0:
                self.cells[self.n, c2] = clusters[1 - assign[j], j]
                c2 += 1

            if c1 >= self.maximum_cluster_n - 1 or c2 >= self.maximum_cluster_n - 1:
                return False

        # Fixed self.counts which was non-specific. Do we need both measurements
        # calculated here?
        self.n_signals_per_cell[i] = c1
        self.n_signals_per_cell[self.n] = c2

        # Make sure to maintain matches for selection
        if self.selection_passed_partner[i] > -1:
            self.selection_passed_partner[i] = self.n
            self.selection_passed_partner[self.n] = -1

        return True

    def g1phase(self):
        """Form and break up clusters of plasmids during G1 phase.
        
        Returns:
            bool: the success of G1 phase
            
        """
        for i in range(self.n):
            if self.n_signals_per_cell[i] > 0:
                out = []
                for j in range(self.n_signals_per_cell[i]):
                    out += self.dissociate_cluster_and_jostle(self.cells[i, j], sphase=False)
                if len(out) >= self.maximum_cluster_n:
                    logging.info('\tWARNING: Hitting maximum cluster boundary in G1 phase')
                    return False

                self.cells[i, :len(out)] = out
                self.n_signals_per_cell[i] = len(out)
        
        return True

    # ---------------------------------------------------------------------------------- #
    # BIOLOGICAL CLUSTERING FUNCTIONS0
    @staticmethod
    def _pad(a: list, b: list):
        """Padding function to ensure that clusters stay associated
        The first cluster of a and b are paired, all others are unpaired

        Args:
            a (list): a list of integer cluster sizes
            b (list): a list of integer cluster sizes after replication

        Returns:
            tuple of lists of ints: the clusters of cell b

        """

        if a == [] and b == []: 
            return [], []
        if a == []: 
            return [0]*len(b), b
        if b == []: 
            return a, [0]*len(a)
        
        padto = len(a) + len(b) - 1
        a = a + [0]*(padto - len(a))
        b = [b[0]] + [0]*(padto - len(b)) + b[1:]
        return a, b

    def s_phase_cluster_breakup(self, clusters: np.ndarray):
        """In S-phase, takes a pair of lists of clusters and breaks them up.
        Clusters have already been combined, if necessary. They also have been lined up
        so that they can have shared fates

        Note: clusters[0] has all of the original clusters, clusters[1] has the
        duplicated versions which, by definition, are <= clusters[0] in size
        
        Args:
            clusters (np.ndarray): a list of clusters to duplicate and possibly break up

        Returns:
            np.ndarray: updated array of clusters

        """

        p1 = []
        p2 = []

        for j in range(clusters[0].size):
            a = self.dissociate_cluster_and_jostle(clusters[0, j], sphase=True)
            b = self.dissociate_cluster_and_jostle(clusters[1, j], sphase=True)
            a, b = self._pad(a, b)
            p1 += a
            p2 += b

        return np.array([p1, p2], dtype=np.uint16)

    def dissociate_cluster_and_jostle(self, n: int, sphase: bool = True):
        """Cluster dissociation is a two-part process. First, check for cluster 
        "jostling" from one parameter. If jostled, break up cluster using the 
        Chinese Restaurant Process.

        Args:
            n (int): an input cluster
            sphase (bool, optional): Whether the analysis is in s-0phase. Defaults to True.

        Returns:
            list of ints: a list of output cluster sizes
        
        """
        # Check if 1 or 0- those values have guaranteed outputs
        if n == 0: 
            return []
        if n == 1: 
            return [1]
            
        comparison = (self.cluster_jostling*self.cluster_jostling_s_vs_g1 if sphase else self.cluster_jostling)
        if np.random.random() > comparison: 
            return [n]

        return crp(n, self.cluster_crp_alpha)

    # ================================================================================== #
    # SET UP FUNCTIONS
    # Ancillary functions to set up the simulation

    def sample(self):
        """Draw a random sample of the cells so that the simulation takes a reasonable amount of time.
        Find the nonzero elements and only sample from them, shuffle, and set the new n."""
        nonzeros = np.nonzero(self.n_signals_per_cell[:self.n])
        np.random.shuffle(nonzeros[0])
        new_n = np.minimum(self.sample_to_n, nonzeros[0].size)
        
        self.zeros += (1.0 - self.zeros)*(1.0 - float(nonzeros[0].size)/self.n)
        self.cells[:new_n] = self.cells[np.sort(nonzeros[0][:new_n])]
        self.cells[new_n:self.n] = 0

        self.n = new_n
        self._reset_counts_per_cell()
        self._sample_just_run = True

        return self

    def _reset_counts_per_cell(self):
        """Reset the number of signals per cell and the number of plasmids per cell.

        Notes:
            Calculate the number of clusters per cell using a sum
            Could use the shuffled values from above, but it's harder to error check
            The slight hit in time per generation is worthwhile
        
        """
        n_clusters = np.copy(self.cells[:self.n])
        n_clusters[np.nonzero(n_clusters)] = 1
        self.n_signals_per_cell[:self.n] = np.sum(n_clusters, axis=1)
        self.n_signals_per_cell[self.n:] = 0

        # Reset the counts to plasmids per cell
        self.n_plasmids_per_cell[:self.n] = np.sum(self.cells[:self.n], axis=1)
        self.n_plasmids_per_cell[self.n:] = 0

    # ================================================================================== #
    # OUTPUT FUNCTIONS
    # Ancillary functions to return the results of the simulation

    def _current_signals_per_cell(self):
        """
        Calculate the number of signals per cell. Used to generate the initial population.
        
        """
        # Sample to make sure that we have the correct population size and zeros are where they should be
        self.sample()

        # Calculate the histogram of the n_signals_per_cell
        cluster_bins = np.bincount(self.n_signals_per_cell).astype(np.float64)/self.n

        # Set the zeros
        if cluster_bins.size > 0: 
            cluster_bins[0] = self.zeros
            cluster_bins[1:] *= 1.0 - self.zeros
        else: 
            cluster_bins = np.array([self.zeros])

        return cluster_bins

    def get_plasmids_per_cell(self):
        """
        Calculate the number of plasmids per cell.
        
        """
        # Sample to make sure that we have the correct population size and zeros are where they should be
        self.sample()

        # Calculate the histogram of the n_signals_per_cell
        cluster_bins = np.bincount(self.n_plasmids_per_cell).astype(np.float64)/self.n

        # Set the zeros
        if cluster_bins.size > 0: 
            cluster_bins[0] = self.zeros
            cluster_bins[1:] *= 1.0 - self.zeros
        else: 
            cluster_bins = np.array([self.zeros])

        return cluster_bins

    def get_signals_per_cell(self, init_or_final='final'):
        """Return the percentages of signals/cell in the initial or final population by signal size."""
        # Determine if it's the initial population or the final population and return the
        # results as percentages
        if init_or_final[0].lower() == 'i': 
            return np.copy(self.initpop)
        return np.copy(self._current_signals_per_cell())

    def get_plasmids_per_signal(self):
        """Return the percentages of plasmids/signal by number of plasmids."""
        # Calculate the histogram of the number of plasmids per signal and return as
        # percentage
        cluster_size_bins = np.bincount(self.cells[:self.n].flatten().astype(np.int64)).astype(np.float64)
        cluster_size_bins[1:] = cluster_size_bins[1:]/np.sum(cluster_size_bins[1:])
        return cluster_size_bins*100.0
        

if __name__ == '__main__':
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

    virus = 'kshv'
    gens = 50
    lpp = LatentPlasmidPopulation()
    lpp.set_virus(virus)
    lpp.population(5000, 1, 0)
    lpp.simulate(gens)
    print(lpp.get_plasmids_per_cell())
