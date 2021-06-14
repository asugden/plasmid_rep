import matplotlib.pyplot as plt
import numpy as np

from plasmid_rep import lpp


def plot_crp(path: str, xmax: int = 50):
    """Plot the chinese restaurant process values

    Args:
        path (str): save location
        xmax (int, optional): the maximum for plotting. Defaults to 20.
    
    """
    for crp in [0.1, 0.5, 10]:
        xs = np.arange(xmax) + 1
        av_clusters = [1]
        for x in range(2, xmax+1):
            n_clust = []
            for _ in range(1000):
                out = lpp.crp(x, crp)
                n_clust.append(len(out))
            av_clusters.append(np.mean(n_clust))
        plt.plot(xs, av_clusters, label=crp)
    
    plt.legend()
    plt.xlabel('Cluster size (number of plasmids)')
    plt.ylabel('Number of resulting clusters')
    plt.title('Cluster breakup by CRP')
    plt.xlim((0, xmax))
    plt.savefig(f'{path}_resulting_clusters.pdf')
    plt.clf()

    for crp in [0.1, 0.5, 10]:
        xs = np.arange(xmax) + 1
        largest_cluster = [1]
        for x in range(2, xmax+1):
            large = []
            for _ in range(1000):
                out = lpp.crp(x, crp)
                large.append(max(out))
            largest_cluster.append(np.mean(large))
        plt.plot(xs, largest_cluster, label=crp)
    
    plt.legend()
    plt.xlabel('Cluster size (number of plasmids)')
    plt.ylabel('Largest resulting cluster')
    plt.title('First cluster size by CRP')
    plt.xlim((0, xmax))
    plt.savefig(f'{path}_resulting_cluster_size.pdf')


if __name__ == '__main__':
    import os.path
    pop = lpp.LatentPlasmidPopulation()
    base_path = pop.load_config()['base_path']
    plot_crp(os.path.join(base_path, 'crp'))
