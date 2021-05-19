import os.path

import plasmid_rep.LPP_plasmid_replication as LPPpr
import plasmid_rep.plotting.standard_plots as plotting


def fit_negative_selection(base_path: str, neg_coefficient: float = 0.07):
    """Fit negative selection and confirm stability

    Args:
        base_path (str): base path into which parameters and output should be saved
        neg_coefficient (float): the negative selection coefficient
    
    """
    # Set the current variables
    kshv_pop = LPPpr.LPP_PlasmidReplication()
    kshv_pop.set_virus('kshv')
    kshv_pop.s_phase_duplication_prob = 0.98
    kshv_pop.negative_cluster_selection_coefficient = neg_coefficient

    # Set up a large population and burn in
    kshv_pop.population(5000, 10, 3)
    kshv_pop.simulate_burn_in(20, 0)
    kshv_pop._average_rep = []  # pylint: disable=protected-access
    kshv_pop.simulate(30)
    date_str = kshv_pop.save(base_path, {'mean_duplication_rate': kshv_pop.mean_duplication_rate()})
    logging.info(kshv_pop.mean_duplication_rate())

    # Rerun to confirm stability
    kshv_pop._average_rep = []  # pylint: disable=protected-access
    kshv_pop.simulate(20)
    date_str = kshv_pop.save(base_path, {'mean_duplication_rate': kshv_pop.mean_duplication_rate()})
    print(kshv_pop.mean_duplication_rate())

    # Plot the results
    plotting.plasmids_per_cell(kshv_pop, os.path.join(base_path, f'{date_str} plasmids per cell.pdf'))


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
    
    fit_negative_selection('/Users/arthur/Desktop')
