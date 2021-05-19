"""Dually infected cells are modeled as if KSHV and EBV are fully independent.
This accurately replicates the results we've seen.

"""
import matplotlib.pyplot as plt
import numpy as np
import os.path as opath
import pandas as pd
import seaborn as sns

import LatentPlasmidPopulation as lpp


def simulate(generations: int = 20):
    """Simulate KSHV and EBV distributions after infection with a single plasmid per cell

    Args:
        generations (int): the number of generations to simulate

    Returns:
        lpp.LatentPlasmidPopulation: the distribution of EBV amongst the cells
        lpp.LatentPlasmidPopulation: the distribution of KSHV amongst the cells

    """
    ebv_pop = lpp.LatentPlasmidPopulation()
    ebv_pop.set_virus('ebv')
    ebv_pop.population(2000, 1, 0)
    ebv_pop.simulate(generations)

    kshv_pop = lpp.LatentPlasmidPopulation()
    kshv_pop.set_virus('kshv')
    kshv_pop.population(2000, 1, 0)
    kshv_pop.simulate(generations)

    return ebv_pop, kshv_pop


def plot_dually_infected(ebv_sim: lpp.LatentPlasmidPopulation, 
                         kshv_sim: lpp.LatentPlasmidPopulation,
                         path: str,
                         plot_n: int = 1000):
    """Create a series of plots for dually infected cells

    Args:
        ebv_sim (lpp.LatentPlasmidPopulation): the simulated EBV population
        kshv_sim (lpp.LatentPlasmidPopulation): the simulated KSHV population
    
    """
    sns.set_theme(style='whitegrid')

    # Plasmids per signal
    pps_ebv = ebv_sim.get_plasmids_per_signal()
    pps_kshv = kshv_sim.get_plasmids_per_signal()
    pps_ebv = pps_ebv[1:]/np.sum(pps_ebv[1:])
    pps_kshv = pps_kshv[1:]/np.sum(pps_kshv[1:])

    # Plasmids per cell
    ppc_ebv = ebv_sim.get_plasmids_per_cell()
    ppc_kshv = kshv_sim.get_plasmids_per_cell()

    # Signals per cell
    spc_ebv = ebv_sim.get_signals_per_cell()
    spc_kshv = kshv_sim.get_signals_per_cell()

    # Plot plasmids per signal
    pps_ebv_sample = np.random.choice(np.arange(len(pps_ebv)) + 1, size=plot_n, p=pps_ebv)
    pps_kshv_sample = np.random.choice(np.arange(len(pps_kshv)) + 1, size=plot_n, p=pps_kshv)
    data = pd.DataFrame({'plasmids per signal': np.concatenate((pps_ebv_sample, pps_kshv_sample)), 
                         'category': ['ebv']*plot_n + ['kshv']*plot_n})
    ax = sns.violinplot(data=data, x='category', y='plasmids per signal')
    ax.set(ylim=(0, 10))
    plt.savefig(opath.join(path, 'plasmids_per_signal.pdf'))
    plt.clf()

    # Plot plasmids per cell zeros
    data = pd.DataFrame({'fraction of cells with nonzero plasmids': [1-ppc_ebv[0], 1-ppc_kshv[0]], 
                         'category': ['ebv', 'kshv']})
    sns.barplot(data=data, x='category', y='fraction of cells with nonzero plasmids')
    plt.savefig(opath.join(path, 'plasmids_per_cell_zeros.pdf'))
    plt.clf()

    # Plot mean plasmids per cell
    ppc_ebv_sample = np.random.choice(np.arange(len(ppc_ebv)) + 1, size=plot_n, p=ppc_ebv)
    ppc_kshv_sample = np.random.choice(np.arange(len(ppc_kshv)) + 1, size=plot_n, p=ppc_kshv)
    data = pd.DataFrame({'plasmids per cell': np.concatenate((ppc_ebv_sample, ppc_kshv_sample)), 
                         'category': ['ebv']*plot_n + ['kshv']*plot_n})
    ax = sns.violinplot(data=data, x='category', y='plasmids per cell')
    ax.set(ylim=(0, 10))
    plt.savefig(opath.join(path, 'mean_plasmids_per_cell.pdf'))
    plt.clf()

    # Plot plasmids per cell > 0
    points_to_plot = 20
    res_ppc_ebv = np.resize(ppc_ebv[1:], points_to_plot)
    res_ppc_ebv[len(ppc_ebv)-1:] = 0
    data = pd.DataFrame({'fraction of cells with n plasmids': np.concatenate([res_ppc_ebv, 
                                                                              np.resize(ppc_kshv[1:], points_to_plot)]), 
                         'n plasmids': [(i+1) for i in range(points_to_plot)] + [(i+1) for i in range(points_to_plot)],
                         'category': ['ebv']*points_to_plot + ['kshv']*points_to_plot})
    sns.barplot(data=data, x='n plasmids', y='fraction of cells with n plasmids', hue='category')
    plt.savefig(opath.join(path, 'plasmids_per_cell.pdf'))
    plt.clf()

    print(pps_ebv)
    print(pps_kshv)
    print(1)

    


if __name__ == '__main__':
    import logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

    ebv_sim, kshv_sim = simulate()
    plot_dually_infected(ebv_sim, kshv_sim, '/Users/arthur/Desktop/')

    # print(lpp.get_signals_per_cell())
    # print(lpp.get_plasmids_per_cell())
    # print(1)
    # lpp.save()

    # df = pd.DataFrame({'KSHV plasmids/cell': np.arange(20), 'Percentage of cells': np.resize(lpp.get_plasmids_per_cell(), 20)})
    # sns.lineplot(data=df, x='KSHV plasmids/cell', y='Percentage of cells').set_title('KSHV+EBV after 20 generations')
    # plt.savefig('/Users/arthur/Desktop/kshv-ppc-20gen.pdf')

    # df = pd.DataFrame({'KSHV signals/cell': np.arange(20), 'Percentage of cells': np.resize(lpp.get_plasmids_per_cell(), 20)})
    # sns.lineplot(data=df, x='KSHV signals/cell', y='Percentage of cells').set_title('KSHV+EBV after 20 generations')
    # plt.savefig('/Users/arthur/Desktop/kshv-spc-20gen.pdf')
