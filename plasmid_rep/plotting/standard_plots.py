import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set_style('whitegrid')

from plasmid_rep import LPP


def plasmids_per_cell(sim: LPP.LatentPlasmidPopulation,
                      path: str,
                      zeros: bool = True,
                      nonzeros: bool = True):
    assert(zeros or nonzeros)
    ppc = sim.get_plasmids_per_cell()
    x = np.arange(len(ppc))

    if not nonzeros:
        (sns.barplot(y=[ppc[0]])
            .set_ylim((0, 1))
            .set_title('Plasmids per cell')
            .set_ylabel('fraction of cells with zero plasmids'))
        plt.savefig(path)
        return
    elif not zeros:
        x, ppc = x[1:], ppc[1:]

    ymax = 1 if zeros else 0.2
    (sns.barplot(x=x, y=ppc)
        .set(ylim=(0, ymax),
             xlim=(x[0]-0.5, 20+0.5),
             title='Plasmids per cell', 
             xlabel='fraction of cells with N plasmids',
             ylabel='N plasmids per cell'))
    plt.savefig(path)


def plot_dually_infected(ebv_sim: LPP.LatentPlasmidPopulation, 
                         kshv_sim: LPP.LatentPlasmidPopulation,
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
