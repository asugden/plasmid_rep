import matplotlib.pyplot as plt
import numpy as np


def plot_coeff(path: str, dup: float = 0.98, neg_coeff: float = 0.07, xmax: int = 20):
    """Plot the positive selection coefficient

    Args:
        path (str): save location
        coeff (float, optional): positive selection coefficient. Defaults to 0.2.
        xmax (int, optional): the maximum for plotting. Defaults to 20.
    """
    x = np.linspace(0, xmax, 200)
    y = dup*np.exp(-neg_coeff*x)
    plt.plot(x, y, label='coeff')

    neg_coeff = neg_coeff/2
    y = dup*np.exp(-neg_coeff*x)
    plt.plot(x, y, label='coeff/2')

    neg_coeff = neg_coeff*4
    y = dup*np.exp(-neg_coeff*x)
    plt.plot(x, y, label='coeff*2')

    plt.plot(x, [0.92]*len(x), label='no sel')

    plt.legend()

    plt.xlabel('Number of plasmids per cluster')
    plt.ylabel('Probability of duplication of each plasmid')
    plt.title('Negative selection')
    plt.xlim((0, xmax))
    plt.ylim((0, 1.05))
    plt.savefig(path)
    # plt.show()


if __name__ == '__main__':
    plot_coeff('/Users/arthur/Desktop/negative_selection_0.7.pdf')