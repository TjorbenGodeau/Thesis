from dSB import dSB
import numpy as np
import dSB_plt as plt

def main():
    #coupling matrix J
    J = np.zeros((4, 4))
    J[0, 1] = 0.5
    J[0, 2] = 0.3
    J[0, 3] = 0.4
    J[1, 2] = 0.2
    J[1, 3] = 0.6
    J[2, 3] = 0.1
    
    # Initialize the dSB solver
    solver = dSB(J, dt=1.0, a0=1.0, c0=1.0, Nstep=50)
    solver.initialize()
    
    # Track energy and positions during evolution
    energies = []
    positions = []
    
    # Run the schedule and collect data
    for _ in range(solver.Nstep):
        solver.step()
        energies.append(solver.energy())
        positions.append(solver.x.copy())
    
    # Plot results
    plt.plot_ising_energy(energies, show=True)
    plt.plot_oscillator_positions(positions, show=True)

if __name__ == "__main__":
    main()