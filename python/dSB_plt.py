import numpy as np
import matplotlib.pyplot as plt

def plot_ising_energy(energies, save_path=None, show=True):
    """
    Plot Ising energy as a function of stepcount.
    
    Parameters:
    energies (list or array): List of Ising energy values at each step
    save_path (str, optional): Path to save the figure
    show (bool): Whether to display the plot
    """
    steps = np.arange(len(energies))
    fig, ax = plt.subplots()
    ax.plot(steps, energies, color="C0")
    ax.set_xlabel("stepcount")
    ax.set_ylabel("Ising energy")
    ax.set_title("Ising energy vs stepcount")
    ax.grid(alpha=0.3)
    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)

def plot_oscillator_positions(positions, save_path=None, show=True):
    """
    Plot positions of oscillators as a function of stepcount.
    
    Parameters:
    positions (list of arrays or 2D array): Positions at each step, shape (steps, N)
    save_path (str, optional): Path to save the figure
    show (bool): Whether to display the plot
    """
    positions = np.array(positions)  # Ensure it's a numpy array
    steps = np.arange(positions.shape[0])
    fig, ax = plt.subplots()
    for i in range(positions.shape[1]):
        ax.plot(steps, positions[:, i], label=f"Oscillator {i+1}")
    ax.set_xlabel("stepcount")
    ax.set_ylabel("Position")
    ax.set_title("Positions of oscillators vs stepcount")
    ax.legend()
    ax.grid(alpha=0.3)
    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)

