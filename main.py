import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
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

def plot_energy(energies, dt, A_val, M_val, save_path=None, show_energy=False):
    steps = np.arange(len(energies))
    plt.figure()
    plt.plot(steps, energies, color="C0")
    plt.xlabel("steps m")
    plt.ylabel("Ising energy")
    plt.title(f"Energy vs steps (A={A_val}, M={M_val})")
    plt.grid(alpha=0.3)
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
    if show_energy:
        plt.show()
    plt.close()

def plot_heatmap(A_values, energy_matrix, dit,  heat_steps, heat_M, cmap="viridis_r", save_path=None, show_heat=False):
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
        plt.close()

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
    parser.add_argument("--plot", choices=["none","energy","heatmap","both"], default="both", help="Which plots to produce")
    parser.add_argument("--heat-steps", type=int, default=50_000, help="Number of steps (m) to use along heatmap x-axis")
    parser.add_argument("--heat-M", type=int, default=None, help="Which M value to use for heatmap (defaults to first M)")
    parser.add_argument("--show-energy", action="store_true", default=False, help="Show plots interactively")
    parser.add_argument("--show-heat", action="store_true", default=True, help="Show plots interactively")
    parser.add_argument("--save-prefix", type=str, default=None, help="Prefix to save generated plots (files will be created)")
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
            
            energies = []
            def cb(m, x, y, p):
                energies.append(model.energy())
            
            model.run(steps=steps, callback=cb)

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
                "spins": " ".join(str(int(s)) for s in spins),
                "energies": list(energies)
            })

            if args.plot in ("energy","both"):
                save_path = None
                if args.save_prefix:
                    save_path = f"{args.save_prefix}_energy_A{A_val}_M{M_val}.png"
                plot_energy(energies, args.dt, A_val, M_val, save_path=save_path, show_energy=args.show_energy)
    
    # optional CSV output
    if args.out:
        fieldnames = ["A", "M", "dt", "seed", "per_spin", "energy", "spins"]
        with open(args.out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow({k: r[k] for k in fieldnames})

    # heatmap: use specified heat-M (or the first M) and A sweep
    if args.plot in ("heatmap","both"):
        heat_M = args.heat_M if args.heat_M is not None else args.M[0]
        # filter results for the chosen M and sort by A
        rows = [r for r in results if r["M"] == int(heat_M)]
        if len(rows) == 0:
            print(f"No runs found for M={heat_M}, cannot build heatmap.")
            return
        rows.sort(key=lambda r: r["A"])
        A_values = [r["A"] for r in rows]
        heat_steps = int(args.heat_steps)
        energy_matrix = np.zeros((len(rows), heat_steps), dtype=float)
        for i, r in enumerate(rows):
            e = r["energies"]
            if len(e) >= heat_steps:
                energy_matrix[i, :] = e[:heat_steps]
            else:
                # pad with last known energy (or nan if empty)
                if len(e) == 0:
                    energy_matrix[i, :] = np.nan
                else:
                    energy_matrix[i, :len(e)] = e
                    energy_matrix[i, len(e):] = e[-1]

        save_path = None
        if args.save_prefix:
            save_path = f"{args.save_prefix}_heatmap_M{heat_M}.png"
        plot_heatmap(A_values, energy_matrix, args.dt, heat_steps, heat_M, save_path=save_path, show_heat=args.show_heat)

if __name__ == "__main__":
    main()