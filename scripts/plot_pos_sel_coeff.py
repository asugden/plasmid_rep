import matplotlib.pyplot as plt
import numpy as np


def plot_coeff(path: str, coeff: float = 0.2, xmax: int = 50):
    """Plot the positive selection coefficient

    Args:
        path (str): save location
        coeff (float, optional): positive selection coefficient. Defaults to 0.2.
        xmax (int, optional): the maximum for plotting. Defaults to 20.
    """
    x = np.linspace(0, xmax, 200)
    y = 2*(np.exp(coeff*x)/(1.0 + np.exp(coeff*x))) - 1
    plt.plot(x, y, label='coeff')

    coeff = coeff/2
    y = 2*(np.exp(coeff*x)/(1.0 + np.exp(coeff*x))) - 1
    plt.plot(x, y, label='coeff/2')

    coeff = coeff*4
    y = 2*(np.exp(coeff*x)/(1.0 + np.exp(coeff*x))) - 1
    plt.plot(x, y, label='coeff*2')

    plt.plot(x, [1]*len(x), label='no sel')

    plt.legend()

    plt.xlabel('Number of plasmids per cell')
    plt.ylabel('Probability of cell replicating per "round"')
    plt.title('Positive selection')
    plt.xlim((0, xmax))
    plt.ylim((0, 1.05))
    plt.savefig(path)


if __name__ == '__main__':
    plot_coeff('/Users/arthur/Desktop/positive_selection_0.2.pdf')