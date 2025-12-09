import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
import re

def plot_energy(energies, A_val, M_val, save_path=None, show_energy=False):
    steps = np.arange(len(energies))
    fig, ax = plt.subplots()
    ax.plot(steps, energies, color="C0")
    ax.set_xlabel("steps m")
    ax.set_ylabel("Ising energy")
    ax.set_title(f"Energy vs steps (A={A_val}, M={M_val})")
    ax.grid(alpha=0.3)
    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    if show_energy:
        plt.show()
    else:
        plt.close(fig)

def plot_hamiltonian(hamiltonians, A_val, M_val, save_path=None, show_ham=False):
    steps = np.arange(len(hamiltonians))
    fig, ax = plt.subplots()
    ax.plot(steps, hamiltonians, color="C1", label="Hamiltonian")
    ax.set_xlabel("steps m")
    ax.set_ylabel("Hamiltonian")
    ax.set_title(f"Hamiltonian vs steps (A={A_val}, M={M_val})")
    ax.grid(alpha=0.3)
    ax.legend()
    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    if show_ham:
        plt.show()
    else:
        plt.close(fig)

def plot_heatmap(A_values, energy_matrix, heat_steps, heat_M, cmap="viridis_r", save_path=None, show_heat=False):
    plt.figure(figsize=(8, max(4, 0.2*len(A_values))))
    im = plt.imshow(energy_matrix, aspect="auto", origin="lower", cmap=cmap, extent=[0, heat_steps, A_values[0], A_values[-1]])
    plt.colorbar(im, label="ising energy")
    plt.xlabel("steps m")
    plt.ylabel("A")
    plt.title(f"Heatmap of energy over steps (M={heat_M})")
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
    if show_heat:
        plt.show()
    else:
        plt.close()

def make_save_prefix(base: str, args):
    base = re.sub(r'[<>:"/\\|?*]', '_', base)  # sanitize
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    pref = f"{base}_seed{args.seed}_{ts}"
    # If base included a folder, ensure it exists:
    out_dir = Path(pref).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    return pref