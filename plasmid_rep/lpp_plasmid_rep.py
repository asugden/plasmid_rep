import numpy as np

from plasmid_rep import lpp


class LPPPlasmidRep(lpp.LatentPlasmidPopulation):
    """Subclass LatentPlasmidPopulation to be able to track average replication.

    """
    def __init__(self):
        self._average_rep = []
        super().__init__()

    def sphase_negative_selection(self, i: int) -> np.ndarray:
        """Simulate S-phase for a single cell and return duplicated aggregates and split cells.
        Duplicate plasmids and add to clusters
        Generate a list of random numbers to check for plasmid duplication of length
        n_plasmids_per_cell
        Gives a list of 0s and 1s which can be summed
        This version adds saving of duplication rates

        Args:
            i (int): the cell to run s-phase on
        
        Returns:
            np.ndarray: an array of the sizes of the clusters after s-phase

        """
        plas_dup = np.random.random(self.n_plasmids_per_cell[i])
        n_plasmids_duplicated = [0, 0]

        # Generate an array for each cluster duplicate and duplicate clusters
        clusters = np.zeros((2, self.n_signals_per_cell[i]), dtype=np.uint16)
        tally = 0
        for j in range(self.n_signals_per_cell[i]):
            clusters[0, j] = self.cells[i, j]
            dup_probs = self.s_phase_duplication_prob*np.exp(-self.negative_cluster_selection_coefficient*np.arange(self.cells[i,j]))
            clusters[1, j] = np.sum(dup_probs > plas_dup[tally:tally + self.cells[i, j]])
            n_plasmids_duplicated[0] += clusters[1, j]
            n_plasmids_duplicated[1] += clusters[0, j]
            tally += self.cells[i, j]

        self._average_rep.append(n_plasmids_duplicated)

        # Note: clusters[0] has all of the original clusters, clusters[1] has the
        # duplicated versions which, by definition, are <= clusters[0] in size
        plas_aggregation = self.plasmid_repulsion_attraction > np.random.random(self.n_signals_per_cell[i])
        clusters[0, plas_aggregation] += clusters[1][plas_aggregation]
        clusters[1, plas_aggregation] = 0

        return clusters

    def mean_duplication_rate(self) -> float:
        """Measure the mean plasmid duplication rate with negative selection

        Returns:
            float: mean plasmid duplication rate
        
        """
        n_dup = np.sum([v[0] for v in self._average_rep])
        n_tot = np.sum([v[1] for v in self._average_rep])
        return float(n_dup)/n_tot
