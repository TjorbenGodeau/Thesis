import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
from GbSB import GbSB
from Gset_extract import read_matrix_from_file, compute_cut_value
from plotting import plot_energy, plot_hamiltonian, plot_heatmap

def create_random_J(N, seed=0):
    rng = np.random.default_rng(seed)
    A = rng.normal(scale=1.0, size=(N,N))
    J = (A + A.T) / 2.0
    np.fill_diagonal(J, 0.0)
    for i in range(J.shape[0]):
        for j in range(J.shape[1]):
            if i > j:
                J[i,j] = 0
            else:
                continue
    return J

def parse_list(arg, cast):
    if arg is None:
        return []
    return [cast(x) for x in arg]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gset-j", action="store_true", help="Use interaction coefs matrix from G Set")
    parser.add_argument("--gset_path", nargs="?", default=r"C:\Users\tjorb\Documents\Thesis\benchmark\G_Set\G1.txt", help="path to G Set file")
    parser.add_argument("--dt", type=float, default=1.0, help="time step delta t")
    parser.add_argument("--M", type=int, default=200, help="total discrete time steps")
    parser.add_argument("--A", type=float, default=0.0, help="nonlinear control constant A")
    # parser.add_argument("A", nargs="*", type=float, default=list(np.linspace(0.0, 1.0, 11)), help="nonlinear control constants A (default: 0.0..1.0 inclusive)")
    parser.add_argument("-p", "--per-spin", action="store_true", help="use individual p_i for each oscillator")
    parser.add_argument("-e", "--energy", action="store_true", help="show energy plot")
    parser.add_argument("-hm", "--hamiltonian", action="store_true", help="show hamiltonian plot")
    parser.add_argument("--N", type=int, default=2, help="number of spins (used if not using G Set file)")
    parser.add_argument("--seed", type=int, default=0, help="random seed (used if not using G Set file)")
    args = parser.parse_args()

    if args.gset_j:
        J, N = read_matrix_from_file(args.gset_path)
    else:
        J = np.zeros((args.N,args.N))
        J[0,1] = 1.0
        # J = create_random_J(args.N, seed=args.seed)

    model = GbSB(J, dt=args.dt, M=args.M, A=args.A, per_spin=args.per_spin)

    initial_x = np.random.uniform(-1.0, 1.0, size=J.shape[0])
    initial_y = np.zeros(J.shape[0])
    if args.per_spin:
        initial_p = np.ones(J.shape[0])
    else:
        initial_p = 1.0

    print("Initial positions:", initial_x)
    print("Initial momenta:", initial_y)
    model.initialize(x0=initial_x, y0=initial_y, p0=initial_p)

    energies = []
    hamiltonians = []

    def record_callback(m, x, y, p):
        """Callback to record energy and hamiltonian at each step."""
        energies.append(model.energy())
        hamiltonians.append(model.hamiltonian())
        print("Positions:", model.x)
        print("Momenta:", model.y)
        
    model.run(steps=args.M, callback=record_callback)

    spins = model.spins()
    print("Final spins:", spins)
    energy_final = model.energy()
    hamiltonian_final = model.hamiltonian()
    cut_value = compute_cut_value(J, energy_final) if "compute_cut_value" in globals() else None

    # print(f"A={args.A}, M={args.M}, dt={args.dt}, seed={args.seed}, per_spin={args.per_spin} -> Energy={energy_final:.6g}, Hamiltonian={hamiltonian_final:.6g}, Cut value={cut_value}")

    plot_energy(energies, A_val=args.A, M_val=args.M, show_energy=args.energy)
    plot_hamiltonian(hamiltonians, A_val=args.A, M_val=args.M, show_ham=args.hamiltonian)




if __name__ == "__main__":
    main()