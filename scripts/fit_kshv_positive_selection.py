import numpy as np

import plasmid_rep.LPP as LPP
import plasmid_rep.LPP_parent_tracking as LPPpt


def parent_fraction_for_top_n_percent(sim_counts: np.ndarray, 
                                      zero_fraction: float, 
                                      total_parents: int, 
                                      top_fraction=0.5) -> float:
    """Return the fraction of parental lines that produce e.g. 90% of all daughter cells

    Args:
        sim_fractions (np.ndarray): output of counts from each parent
        zero_fraction (float): the fraction of cells with zero plasmids
        total_parents (int): the total number of parents for the simulation
        top_fraction (float, optional): top fraction of parents to consider. Defaults to 0.9.

    Returns:
        float: the fraction of total parents that make up e.g. 90% of daugher cells
    
    """
    pass


def fraction_of_ebvpos_cells(pos_coefficient: float = 0.10, generations: int = 100) -> float:

    ebv_pop = LPP.LatentPlasmidPopulation()
    ebv_pop.set_virus('ebv')
    ebv_pop.positive_selection_coefficient = pos_coefficient

    # Set up a large population and burn in
    ebv_pop.population(5000, lambda_=2)
    ebv_pop.simulate(generations)
    ppc = ebv_pop.get_plasmids_per_cell()
    av_nonzero_ppc = np.average(np.arange(len(ppc) - 1) + 1, weights=ppc[1:])
    print(av_nonzero_ppc)

    return ppc[0]


def fit_positive_selection(base_path: str, pos_coefficient: float = 0.1, generations: int = 100):
    """Fit negative selection and confirm stability

    Args:
        base_path (str): base path into which parameters and output should be saved
        pos_coefficient (float): the positive selection coefficient
    
    """
    # ebv_zeros = fraction_of_ebvpos_cells(generations=generations)  # kshv 0.006154319463987479

    # Set the current variables
    kshv_pop = LPPpt.LPP_ParentTracking()
    kshv_pop.set_virus('kshv')
    kshv_pop.positive_selection_coefficient = pos_coefficient

    # Set up a large population and burn in
    kshv_pop.population(5000, lambda_=0.02)
    print(kshv_pop._total_parents)
    print(kshv_pop.zeros)
    kshv_pop.simulate(generations)
    frac, zeros, parents = kshv_pop.cell_parent_fractions()
    ppc = kshv_pop.get_plasmids_per_cell()[1:]
    av_nonzero_ppc = np.average(np.arange(len(ppc)) + 1, weights=ppc)

    topfrac = parent_fraction_for_top_n_percent(frac, zeros, parents)
    print(topfrac)

    kshv_pop.save(base_path, {'cell_parent_fractions': (0.9, topfrac)})
    np.set_printoptions(suppress=True)
    logging.info(frac)
    print(frac)


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
    
    fit_positive_selection('/Users/arthur/Desktop')
