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
    #ax.legend()
    ax.grid(alpha=0.3)
    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)

def plot_errors(energies, energies_q, positions, positions_q, save_path=None, show=True):
    """
    Plot the absolute errors between dSB and quantized dSB energies and positions.
    
    Parameters:
    energies (list or array): Ising energy values from dSB at each step
    energies_q (list or array): Ising energy values from quantized dSB at each step
    positions (list of arrays or 2D array): Positions from dSB at each step, shape (steps, N)
    positions_q (list of arrays or 2D array): Positions from quantized dSB at each step, shape (steps, N)
    save_path (str, optional): Path to save the figure
    show (bool): Whether to display the plot
    """
    energies = np.array(energies)
    energies_q = np.array(energies_q)
    positions = np.array(positions)
    positions_q = np.array(positions_q)
    
    steps = np.arange(len(energies))
    
    # Energy error
    energy_error = np.abs(energies - energies_q)
    
    # Position error: max absolute error per step
    position_error = np.abs(positions - positions_q)
    max_pos_error = np.max(position_error, axis=1)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    
    ax1.plot(steps, energy_error, color='C0')
    ax1.set_ylabel("Energy Error")
    ax1.set_title("Errors between dSB and quantized dSB")
    ax1.grid(alpha=0.3)
    
    ax2.plot(steps, max_pos_error, color='C1')
    ax2.set_xlabel("stepcount")
    ax2.set_ylabel("Max Position Error")
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)

def plot_energy_errors(energies, energies_qs, labels, save_path=None, show=True):
    """
    Plot the absolute energy errors for different configurations.
    
    Parameters:
    energies (list or array): Ising energy values from dSB at each step
    energies_qs (list of lists or arrays): List of Ising energy values from quantized dSB for each configuration
    labels (list): List of labels for each configuration
    save_path (str, optional): Path to save the figure
    show (bool): Whether to display the plot
    """
    energies = np.array(energies)
    steps = np.arange(len(energies))
    
    fig, ax = plt.subplots()
    for i, label in enumerate(labels):
        energies_q = np.array(energies_qs[i])
        error = np.abs(energies - energies_q)
        ax.plot(steps, error, label=label)
    ax.set_xlabel("stepcount")
    ax.set_ylabel("Energy Error")
    ax.set_title("Energy Errors for different bit configurations")
    ax.legend()
    ax.grid(alpha=0.3)
    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)

