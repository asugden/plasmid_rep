import numpy as np
from typing import Tuple

from plasmid_rep import lpp


class LPPParent(lpp.LatentPlasmidPopulation):
    """Subclass LatentPlasmidPopulation to be able to track average replication.

    """
    def __init__(self):
        self._parent = None
        self._total_parents = 0
        self._orig_mphase = super().mphase
        super().__init__()

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

        # Add in parent tracking
        self._parent = np.zeros(self.n*2, dtype=np.uint16)
        self._parent[:self.n] = np.arange(self.n)

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

    def mphase(self, i: int, clusters: np.ndarray):
        """Account for which parent is duplicated
        
        Args:
            i (int): the cell in which to run m-phase
            cluster (np.ndarray): the current state of the clusters in the cell
            
        Returns:
            bool: the success of the m-phase process

        """
        self._parent[self.n] = self._parent[i]
        return self._orig_mphase(i, clusters)

    def sample(self):
        """Draw a random sample of the cells so that the simulation takes a reasonable amount of time.
        Find the nonzero elements and only sample from them, shuffle, and set the new n."""
        nonzeros = np.nonzero(self.n_signals_per_cell[:self.n])
        np.random.shuffle(nonzeros[0])
        new_n = np.minimum(self.sample_to_n, nonzeros[0].size)
        
        sorted_nonzeros = np.sort(nonzeros[0][:new_n])
        self.zeros += (1.0 - self.zeros)*(1.0 - float(nonzeros[0].size)/self.n)
        self.cells[:new_n] = self.cells[sorted_nonzeros]
        self.cells[new_n:self.n] = 0
        self._parent[:new_n] = self._parent[sorted_nonzeros]
        self._parent[new_n:self.n] = -1

        self.n = new_n
        self._reset_counts_per_cell()
        self._sample_just_run = True

        return self

    def cell_parent_fractions(self) -> Tuple[np.ndarray, float, int]:
        """Measure the fractions of cells with each unique parent as a fraction of the initial population

        Returns:
            np.ndarray: sorted parent fraction from highest to lowest
            float: the fraction of cells  with zero plasmids
            int: the number of total parents
        
        """
        self.sample()
        _, counts = np.unique(self._parent, return_counts=True)
        return np.sort(counts).astype(float)[::-1], self.zeros, self._total_parents
