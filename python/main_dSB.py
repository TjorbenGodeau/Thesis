from dSB import dSB, dSB_quantized, quantize
import numpy as np
import dSB_plt as plt
from Gset_extract import read_matrix_from_file

def main():
    #coupling matrix J
    J, N = read_matrix_from_file(r"C:\Users\tjorb\Documents\Thesis\benchmark\G_Set\G7.txt")
    
    # Configurations to test different bit lengths separately
    configs = [
        {'bits_J': 2, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=2 bits'},
        {'bits_J': 4, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=4 bits'},
        {'bits_J': 8, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=8 bits'},
        {'bits_J': 16, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=16 bits'},
        {'bits_J': 64, 'bits_x': 2, 'bits_y': 64, 'bits_a': 64, 'label': 'x=2 bits'},
        {'bits_J': 64, 'bits_x': 4, 'bits_y': 64, 'bits_a': 64, 'label': 'x=4 bits'},
        {'bits_J': 64, 'bits_x': 8, 'bits_y': 64, 'bits_a': 64, 'label': 'x=8 bits'},
        {'bits_J': 64, 'bits_x': 16, 'bits_y': 64, 'bits_a': 64, 'label': 'x=16 bits'},
        {'bits_J': 64, 'bits_x': 64, 'bits_y': 2, 'bits_a': 64, 'label': 'y=2 bits'},
        {'bits_J': 64, 'bits_x': 64, 'bits_y': 4, 'bits_a': 64, 'label': 'y=4 bits'},
        {'bits_J': 64, 'bits_x': 64, 'bits_y': 8, 'bits_a': 64, 'label': 'y=8 bits'},
        {'bits_J': 64, 'bits_x': 64, 'bits_y': 16, 'bits_a': 64, 'label': 'y=16 bits'},
        {'bits_J': 64, 'bits_x': 64, 'bits_y': 64, 'bits_a': 2, 'label': 'a=2 bits'},
        {'bits_J': 64, 'bits_x': 64, 'bits_y': 64, 'bits_a': 4, 'label': 'a=4 bits'},
        {'bits_J': 64, 'bits_x': 64, 'bits_y': 64, 'bits_a': 8, 'label': 'a=8 bits'},
        {'bits_J': 64, 'bits_x': 64, 'bits_y': 64, 'bits_a': 16, 'label': 'a=16 bits'}
    ]
    
    # Initialize the dSB solver
    solver = dSB(J, dt=0.5, a0=1.0, c0=1.0, Nstep=50)
    solver.initialize()
    
    # Initialize quantized solvers for each configuration
    solvers_q = []
    for config in configs:
        solver_q = dSB_quantized(J, dt=0.5, a0=1.0, c0=1.0, Nstep=50, 
                                 bits_J=config['bits_J'], bits_x=config['bits_x'], 
                                 bits_y=config['bits_y'], bits_a=config['bits_a'], 
                                 J_range=1.0, y_range=10.0)
        solver_q.initialize()
        solvers_q.append(solver_q)
    
    # Track energy and positions during evolution
    energies = []
    energies_qs = [[] for _ in configs]
    positions = []
    positions_qs = [[] for _ in configs]
    
    # Run the schedule and collect data
    for _ in range(solver.Nstep):
        solver.step()
        energies.append(solver.energy())
        positions.append(solver.x.copy())

        for i, solver_q in enumerate(solvers_q):
            solver_q.step()
            energies_qs[i].append(solver_q.energy())
            positions_qs[i].append(solver_q.x.copy())

    
    # Plot results
    plt.plot_ising_energy(energies, show=False)
    for i, config in enumerate(configs):
        plt.plot_ising_energy(energies_qs[i], show=False)
    plt.plot_ising_energy(energies_qs[-1], show=False)  # Show the last one
    
    plt.plot_oscillator_positions(positions, show=False)
    for i, config in enumerate(configs):
        plt.plot_oscillator_positions(positions_qs[i], show=False)
    plt.plot_oscillator_positions(positions_qs[-1], show=False)  # Show the last one
    
    labels = [config['label'] for config in configs]
    plt.plot_energy_errors(energies, energies_qs, labels, show=True)


if __name__ == "__main__":
    main()