from dSB import dSB, dSB_quantized, quantize
import numpy as np
import dSB_plt as plt
from Gset_extract import read_matrix_from_file

def main():
    #coupling matrix J
    J, N = read_matrix_from_file(r"C:\Users\tjorb\Documents\Thesis\benchmark\G_Set\G7.txt")
    
    # Configurations to test different bit lengths separately
    configs = [
        #{'bits_J': 5, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=5 bits'},
        {'bits_J': 6, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=6 bits'},
        {'bits_J': 7, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=7 bits'},
        {'bits_J': 8, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=8 bits'},
        {'bits_J': 9, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=9 bits'},
        {'bits_J': 10, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=10 bits'},
        {'bits_J': 11, 'bits_x': 64, 'bits_y': 64, 'bits_a': 64, 'label': 'J=11 bits'},
        #{'bits_J': 64, 'bits_x': 5, 'bits_y': 64, 'bits_a': 64, 'label': 'x=5 bits'},
        #{'bits_J': 64, 'bits_x': 6, 'bits_y': 64, 'bits_a': 64, 'label': 'x=6 bits'},
        #{'bits_J': 64, 'bits_x': 7, 'bits_y': 64, 'bits_a': 64, 'label': 'x=7 bits'},
        #{'bits_J': 64, 'bits_x': 8, 'bits_y': 64, 'bits_a': 64, 'label': 'x=8 bits'},
        #{'bits_J': 64, 'bits_x': 64, 'bits_y': 5, 'bits_a': 64, 'label': 'y=5 bits'},
        #{'bits_J': 64, 'bits_x': 64, 'bits_y': 6, 'bits_a': 64, 'label': 'y=6 bits'},
        #{'bits_J': 64, 'bits_x': 64, 'bits_y': 7, 'bits_a': 64, 'label': 'y=7 bits'},
        #{'bits_J': 64, 'bits_x': 64, 'bits_y': 8, 'bits_a': 64, 'label': 'y=8 bits'},
        #{'bits_J': 64, 'bits_x': 64, 'bits_y': 64, 'bits_a': 5, 'label': 'a=5 bits'},
        #{'bits_J': 64, 'bits_x': 64, 'bits_y': 64, 'bits_a': 6, 'label': 'a=6 bits'},
        #{'bits_J': 64, 'bits_x': 64, 'bits_y': 64, 'bits_a': 7, 'label': 'a=7 bits'},
        #{'bits_J': 64, 'bits_x': 64, 'bits_y': 64, 'bits_a': 8, 'label': 'a=8 bits'},
        {'bits_J': 512, 'bits_x': 512, 'bits_y': 512, 'bits_a': 512, 'label': 'all=64 bits'}
    ]
    
    # Initialize the dSB solver
    solver = dSB(J, dt=0.25, a0=1.0, c0=1.0, Nstep=50)
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
    plt.plot_oscillator_positions(positions_qs[-1], show=True)  # Show the last one
    
    labels = [config['label'] for config in configs]
    plt.plot_energy_errors(energies, energies_qs, labels, show=True)


if __name__ == "__main__":
    main()