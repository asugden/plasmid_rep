import numpy as np
import os.path

from plasmid_rep import lpp
import plasmid_rep.plotting.standard_plots


def kshv_pos_selection(base_path: str, pos_coefficient: float = 0.20, generations: int = 150) -> float:
    """Simulate KSHV with and without positive selection and plot the results.

    Args:
        base_path (str): path into which the data should be saved
        pos_coefficient (float, optional): the nonzero positive selection coefficient to use. Defaults to 0.10.
        generations (int, optional): the number of generations to simulate. Defaults to 50.

    Returns:
        float: the fraction of cells with zero plasmids per cell after generations
    """
    for sel in [0, pos_coefficient]:
        kshv_pop = lpp.LatentPlasmidPopulation()
        kshv_pop.set_virus('kshv')
        kshv_pop.positive_selection_coefficient = sel
        kshv_pop.negative_selection_coefficient = 0
        kshv_pop.s_phase_duplication_prob = 0.92

        # Set up a large population and burn in
        kshv_pop.population(5000, lambda_=0.02)
        kshv_pop.simulate(generations)
        ppc = kshv_pop.get_plasmids_per_cell()
        av_nonzero_ppc = np.average(np.arange(len(ppc) - 1) + 1, weights=ppc[1:])
        print(av_nonzero_ppc)

        ts = kshv_pop.save(base_path, {})
        plasmid_rep.plotting.standard_plots.plasmids_per_cell(kshv_pop,
                                                            os.path.join(base_path, f'{ts}-kshv_sel{sel}_plasmids_oer_cell.pdf'), 
                                                            xmax=50, 
                                                            add_text=f'KSHV 50 gens positive sel: {sel}')

    return ppc[0]


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
    
    kshv_pos_selection('/Users/arthur/Desktop')
