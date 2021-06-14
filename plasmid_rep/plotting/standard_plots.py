import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os.path
import seaborn as sns
from typing import List, Union

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
        signals_per_cell (bool, optional): extract signals per cell rather than plasmids per cell
        add_text (str, optional): additional text to add to plot title
    
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


def multi_plasmids_per_cell(sim: Union[lpp.LatentPlasmidPopulation, List[lpp.LatentPlasmidPopulation]],
                            path: str,
                            xmax: int = 20,
                            signals_per_cell: bool = False,
                            add_text=''):
    """Plot the number of plasmids per cell into path

    Args:
        sim (Union[lpp.LatentPlasmidPopulation, List[lpp.LatentPlasmidPopulation]]): a simulated population
        path (str): the output path for saving
        xmax (int, optional): the maximum X value
        signals_per_cell (bool, optional): extract signals per cell rather than plasmids per cell
        add_text (str, optional): additional text to add to plot title
    
    """
    path, _ = os.path.splitext(path)
    if not isinstance(sim, list):
        sim = [sim]

    ppc = ([v.get_signals_per_cell() for v in sim] if signals_per_cell 
            else [v.get_plasmids_per_cell() for v in sim])
    txt = 'signals' if signals_per_cell else 'plasmids'
        
    # Clean
    ppc_cleaned = None
    for i, pop in enumerate(ppc):
        x = np.arange(len(pop))
        divide_by = 1.0 - pop[0]
        x, pop = x[1:], pop[1:]/divide_by

        extreme_values = [np.sum(pop[101:]), np.sum(pop[201:]), np.sum(pop[501:]), np.sum(pop[1001:])]
        x = x[:xmax]
        pop = pop[:xmax]
        x = np.concatenate([x, [100, 200, 500, 1000]])
        pop = np.concatenate([pop, extreme_values], axis=None)
        print(x[-4:], pop[-4:])
        
        df = pd.DataFrame({'x': x, 'pop': pop, 'hue': [i]*len(x)})
        ppc_cleaned = pd.concat([ppc_cleaned, df], axis=0)

    (sns.catplot(data=ppc_cleaned, x='x', y='pop', hue='hue', kind='bar')
            .set(ylim=(0, 0.2),
                 title=f'{txt.capitalize()} per cell{add_text}', 
                 xlabel=f'N {txt} per cell',
                 ylabel=f'fraction of nonzero cells with N {txt}'))

    plt.savefig(path + '_comb.pdf')
    plt.clf()
