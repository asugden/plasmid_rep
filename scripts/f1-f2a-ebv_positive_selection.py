import numpy as np
import os.path

from plasmid_rep import lpp
import plasmid_rep.plotting.standard_plots


def ebv_pos_selection(base_path: str, generations: int = 50) -> float:
    """Simulate EBV with and without positive selection and plot the results.

    Args:
        base_path (str): path into which the data should be saved
        generations (int, optional): the number of generations to simulate. Defaults to 50.

    Returns:
        float: the fraction of cells with zero plasmids per cell after generations
    """
    for sel in [0, None]:
        ebv_pop = lpp.LatentPlasmidPopulation()
        ebv_pop.set_virus('ebv')
        ebv_pop.update_parameters({'positive_selection_coefficient': sel})
        sel = ebv_pop.positive_selection_coefficient  # In case of None, get default

        # Set up a large population and burn in
        ebv_pop.population(5000, lambda_=2)
        ebv_pop.simulate(generations)
        ppc = ebv_pop.get_plasmids_per_cell()
        av_nonzero_ppc = np.average(np.arange(len(ppc) - 1) + 1, weights=ppc[1:])
        print(av_nonzero_ppc)

        ts = ebv_pop.save(base_path, {})
        plasmid_rep.plotting.standard_plots.plasmids_per_cell(ebv_pop,
                                                            os.path.join(base_path, f'{ts}-ebv_sel{sel}_plasmids_oer_cell.pdf'))

    return ppc[0]


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
    
    pop = lpp.LatentPlasmidPopulation()
    base_path = pop.load_config()['base_path']
    ebv_pos_selection(base_path)
