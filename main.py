import argparse
import csv
import numpy as np
from GbSB import GbSB

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
    parser = argparse.ArgumentParser(description="Run GbSB with configurable parameters.")
    parser.add_argument("--N", type=int, default=16, help="Number of spins/oscillators")
    parser.add_argument("--M", type=int, nargs="+", default=[50_000], help="One or more M values (space separated)")
    parser.add_argument("--dt", type=float, default=0.1, help="Time step Δt")
    parser.add_argument("--A", type=float, nargs="+", default=[i/100.0 for i in range(0,101)], help="One or more A values (space separated)")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed for coupling matrix and init")
    parser.add_argument("--per-spin", action="store_true", default=True , help="Use per-spin p_i schedule")
    parser.add_argument("--steps", type=int, default=None, help="Number of integration steps to run (defaults to M for each run)")
    parser.add_argument("--regen-j", action="store_true", default=False , help="Regenerate J for each (A,M) run (otherwise same J reused)")
    parser.add_argument("--out", type=str, default=None, help="Optional CSV output file to save results")
    args = parser.parse_args()
    for name, val in vars(args).items():
        print(name, ":", val)

    J_base = create_random_J(args.N, seed=args.seed)
    print(J_base)

    results = []
    for A_val in args.A:
        for M_val in args.M:
            J = create_random_J(args.N, seed=args.seed) if args.regen_j else J_base
            model = GbSB(J, dt=args.dt, M=int(M_val), A=float(A_val), per_spin=args.per_spin)

            #re-seed initial state for fair comparisions across runs
            rng = np.random.default_rng(args.seed)
            x0 = 0.01 * rng.standard_normal(args.N)
            y0 = np.zeros(args.N)
            model.initialize(x0=x0, y0=y0)

            steps = args.steps if args.steps is not None else int(M_val)
            model.run(steps=steps)

            spins = model.spins()
            energy = model.energy()

            print(f"A={A_val}, M={M_val}, dt={args.dt}, seed={args.seed}, per_spin={args.per_spin} -> Energy={energy:.6g}")
            results.append({
                "A": A_val,
                "M": int(M_val),
                "dt": args.dt,
                "seed": args.seed,
                "per_spin": args.per_spin,
                "energy": float(energy),
                "spins": " ".join(str(int(s)) for s in spins)
            })
    
    if args.out:
        fieldnames = ["A", "M", "dt", "seed", "per_spin", "energy", "spins"]
        with open(args.out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow(r)

if __name__ == "__main__":
    main()