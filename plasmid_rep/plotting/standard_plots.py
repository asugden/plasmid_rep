import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os.path
import seaborn as sns

from plasmid_rep import lpp


def plasmids_per_cell(sim: lpp.LatentPlasmidPopulation,
                      path: str,
                      xmax: int = 20,
                      signals_per_cell: bool = False,
                      add_text=''):
    """Plot the number of plasmids per cell into path

    Args:
        sim (lpp.LatentPlasmidPopulation): a simulated population
        path (str): the output path for saving
        xmax (int, optional): the maximum X value
    
    """
    name, _ = os.path.splitext(path)

    if signals_per_cell:
        ppc = sim.get_signals_per_cell()
    else:
        ppc = sim.get_plasmids_per_cell()
    x = np.arange(len(ppc))
    txt = 'signals' if signals_per_cell else 'plasmids'

    # Plot zeros
    sns.set_style('whitegrid')
    (sns.barplot(y=[ppc[0]])
        .set(ylim=(0, 1),
             title='Fraction of cells with zero {txt}}',
             ylabel='cell fraction'))
    plt.savefig(name + '_zeros.pdf')
    plt.clf()

    # Get averages for plotting
    av_nonzero_ppc = np.average(np.arange(len(ppc) - 1) + 1, weights=ppc[1:])
    av_zero_ppc = np.average(np.arange(len(ppc)) + 1, weights=ppc)
        
    # Plot nonzeros
    divide_by = 1.0 - ppc[0]
    x, ppc = x[1:], ppc[1:]/divide_by

    extreme_values = [np.sum(ppc[101:]), np.sum(ppc[201:]), np.sum(ppc[501:]), np.sum(ppc[1001:])]
    x = x[:xmax]
    ppc = ppc[:xmax]
    x = np.concatenate([x, [100, 200, 500, 1000]])
    ppc = np.concatenate([ppc, extreme_values], axis=None)
    print(x[-4:], ppc[-4:])

    ymax = 1
    _ = (sns.barplot(x=x, y=ppc)
             .set(ylim=(0, ymax),
                  #xlim=(-0.5, xmax+0.5),
                  title=f'{txt.capitalize()} per cell{add_text}', 
                  xlabel=f'N {txt} per cell',
                  ylabel=f'fraction of nonzero cells with N {txt}'))
    plt.plot([av_nonzero_ppc, av_nonzero_ppc], [1, 0.9], color='gray')
    plt.plot([av_zero_ppc, av_zero_ppc], [1, 0.9], color='black')

    plt.savefig(name + '_nonzeros.pdf')
    plt.clf()


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
    plt.savefig(os.path.join(path, 'plasmids_per_signal.pdf'))
    plt.clf()

    # Plot plasmids per cell zeros
    data = pd.DataFrame({'fraction of cells with nonzero plasmids': [1-ppc_ebv[0], 1-ppc_kshv[0]], 
                         'category': ['ebv', 'kshv']})
    sns.barplot(data=data, x='category', y='fraction of cells with nonzero plasmids')
    plt.savefig(os.path.join(path, 'plasmids_per_cell_zeros.pdf'))
    plt.clf()

    # Plot mean plasmids per cell
    ppc_ebv_sample = np.random.choice(np.arange(len(ppc_ebv)) + 1, size=plot_n, p=ppc_ebv)
    ppc_kshv_sample = np.random.choice(np.arange(len(ppc_kshv)) + 1, size=plot_n, p=ppc_kshv)
    data = pd.DataFrame({'plasmids per cell': np.concatenate((ppc_ebv_sample, ppc_kshv_sample)), 
                         'category': ['ebv']*plot_n + ['kshv']*plot_n})
    ax = sns.violinplot(data=data, x='category', y='plasmids per cell')
    ax.set(ylim=(0, 10))
    plt.savefig(os.path.join(path, 'mean_plasmids_per_cell.pdf'))
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
    plt.savefig(os.path.join(path, 'plasmids_per_cell.pdf'))
    plt.clf()

    print(pps_ebv)
    print(pps_kshv)
    print(1)
