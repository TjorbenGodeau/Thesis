import argparse
import numpy as np
from GbSB import GbSB

def create_random_J(N, seed=0):
    rng = np.random.default_rng(seed)
    A = rng.normal(scale=1.0, size=(N,N))
    J = (A + A.T) / 2.0
    np.fill_diagonal(J, 0.0)
    return J

def main():
    parser = argparse.ArgumentParser(description="Run GbSB with configurable parameters.")
    parser.add_argument("--N", type=int, default=16, help="Number of spins/oscillators")
    parser.add_argument("--M", type=int, default=1000, help="Total discrete time steps in schedule")
    parser.add_argument("--dt", type=float, default=0.1, help="Time step Δt")
    parser.add_argument("--A", type=float, default=0.01, help="Nonlinear control constant A")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed for coupling matrix and init")
    parser.add_argument("--per-spin", action="store_true", help="Use per-spin p_i schedule")
    parser.add_argument("--steps", type=int, default=None, help="Number of integration steps to run (defaults to M)")
    args = parser.parse_args()

    J = create_random_J(args.N, seed=args.seed)
    model = GbSB(J, dt=args.dt, M=args.M, A=args.A, per_spin=args.per_spin)

    rng = np.random.default_rng(args.seed)
    x0 = 0.01 * rng.standard_normal(args.N)
    y0 = np.zeros(args.N)
    model.initialize(x0=x0, y0=y0)

    model.run(steps=args.steps)

    spins = model.spins()
    energy = model.energy()

    print(f"Parameters: N={args.N}, M={args.M}, dt={args.dt}, A={args.A}, seed={args.seed}, per_spin={args.per_spin}")
    print("Final spins (+1/-1):", spins)
    print("Ising energy:", energy)

if __name__ == "__main__":
    main()