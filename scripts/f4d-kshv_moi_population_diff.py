import os.path

from plasmid_rep import lpp
import plasmid_rep.plotting.standard_plots


def kshv_pos_selection(base_path: str, generations: int = 50, turn_off_neg: bool = False) -> float:
    """Simulate KSHV with and without positive selection and plot the results.

    Args:
        base_path (str): path into which the data should be saved
        pos_coefficient (float, optional): the nonzero positive selection coefficient to use. Defaults to 0.10.
        generations (int, optional): the number of generations to simulate. Defaults to 50.

    Returns:
        float: the fraction of cells with zero plasmids per cell after generations
    """
    # kshv_pop = lpp.LatentPlasmidPopulation()
    # kshv_pop.set_virus('kshv')
    # if turn_off_neg:
    #     kshv_pop.update_parameters({'negative_cluster_selection_coefficient': 0,
    #                                 's_phase_duplication_prob': 0.92})
    # # Set up a large population and burn in
    # kshv_pop.population(2000, mu=1, sd=0)
    # kshv_pop.simulate(generations)

    # ts = kshv_pop.save(base_path, {})

    big_pop = lpp.LatentPlasmidPopulation()
    big_pop.set_virus('kshv')
    if turn_off_neg:
        big_pop.update_parameters({'negative_cluster_selection_coefficient': 0,
                                   's_phase_duplication_prob': 0.92})
    # Set up a large population and burn in
    big_pop.population(2000, mu=2, sd=0)
    big_pop.simulate(generations)

    ts2 = big_pop.save(base_path, {})
    neg_text = 'noneg' if turn_off_neg else 'neg'
    plasmid_rep.plotting.standard_plots.multi_plasmids_per_cell([big_pop],
                                                                os.path.join(base_path, f'{ts2}-{neg_text}-kshv_small_big.pdf'), 
                                                                xmax=50, 
                                                                add_text=f'KSHV {generations} gens both selections')


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
    
    pop = lpp.LatentPlasmidPopulation()
    base_path = pop.load_config()['base_path']
    kshv_pos_selection(base_path)
