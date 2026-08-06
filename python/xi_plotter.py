import numpy as np
import matplotlib.pyplot as plt
import argparse

# ----------------------------------------------------------------------
# Coupling matrix generators
def generate_coupling(structure, N, density=0.3, seed=42):
    """Return a symmetric zero-diagonal coupling matrix J (N x N)."""
    np.random.seed(seed)
    J = np.zeros((N, N))
    
    if structure == 'ferro2':
        # Exactly the 2-spin ferromagnetic problem from the paper
        if N != 2:
            raise ValueError("ferro2 structure requires N=2.")
        J[0,1] = 1.0
        J[1,0] = 1.0

    elif structure == 'ring':
        for i in range(N):
            j = (i + 1) % N
            J[i, j] = -1.0
            J[j, i] = -1.0

    elif structure == 'path':
        for i in range(N-1):
            J[i, i+1] = -1.0
            J[i+1, i] = -1.0

    elif structure == 'star':
        for i in range(1, N):
            J[0, i] = -1.0
            J[i, 0] = -1.0

    elif structure == 'complete_bipartite':
        if N % 2 != 0:
            raise ValueError(f"N={N} must be even for complete_bipartite.")
        half = N // 2
        for i in range(half):
            for j in range(half, N):
                J[i, j] = -1.0
                J[j, i] = -1.0

    elif structure == 'bipartite':
        if N % 2 != 0:
            raise ValueError(f"N={N} must be even for bipartite.")
        half = N // 2
        for i in range(half):
            for j in range(half, N):
                if np.random.rand() < density:
                    val = np.random.choice([-1.0, 1.0])
                    J[i, j] = val
                    J[j, i] = val

    elif structure == 'ladder':
        if N % 2 != 0:
            raise ValueError(f"N={N} must be even for ladder.")
        half = N // 2
        for i in range(half-1):
            J[i, i+1] = -1.0
            J[i+1, i] = -1.0
        for i in range(half, N-1):
            J[i, i+1] = -1.0
            J[i+1, i] = -1.0
        for i in range(half):
            J[i, i+half] = -1.0
            J[i+half, i] = -1.0

    elif structure == '2dgrid':
        side = int(round(np.sqrt(N)))
        if side*side != N:
            raise ValueError(f"N={N} is not a perfect square, required for 2dgrid.")
        for r in range(side):
            for c in range(side):
                idx = r*side + c
                if c+1 < side:
                    j = idx + 1
                    J[idx, j] = -1.0
                    J[j, idx] = -1.0
                if r+1 < side:
                    j = idx + side
                    J[idx, j] = -1.0
                    J[j, idx] = -1.0

    elif structure == 'chimera':
        if N % 2 != 0:
            raise ValueError(f"N={N} must be even for chimera.")
        half = N // 2
        for i in range(half):
            for j in range(half):
                if i != j:
                    J[i, j] = -1.0
                    J[j, i] = -1.0
        for i in range(half, N):
            for j in range(half, N):
                if i != j:
                    J[i, j] = -1.0
                    J[j, i] = -1.0
        for i in range(half):
            for j in range(half, N):
                J[i, j] = 1.0
                J[j, i] = 1.0

    elif structure == 'random':
        for i in range(N):
            for j in range(i+1, N):
                if np.random.rand() < density:
                    val = np.random.choice([-1.0, 1.0])
                    J[i, j] = val
                    J[j, i] = val

    elif structure == 'full':
        for i in range(N):
            for j in range(i+1, N):
                val = np.random.choice([-1.0, 1.0])
                J[i, j] = val
                J[j, i] = val

    else:
        raise ValueError(f"Unknown structure: {structure}")
    return J

# ----------------------------------------------------------------------
# SB simulation functions with correct inelastic wall condition
def run_abs(x0, y0, J, a0, c0, T, dt):
    x = x0.copy()
    y = y0.copy()
    history = [x.copy()]
    for step in range(1, T + 1):
        a = a0 * (step / T)          # a(t) linearly increases from 0 to a0
        Jx = J @ x
        dx = a0 * y
        dy = -x**3 - (a0 - a) * x + c0 * Jx
        y += dy * dt
        x += dx * dt
        history.append(x.copy())
    return np.array(history)

def run_bbs(x0, y0, J, a0, c0, T, dt):
    x = x0.copy()
    y = y0.copy()
    history = [x.copy()]
    for step in range(1, T + 1):
        a = a0 * (step / T)
        Jx = J @ x
        dx = a0 * y
        dy = -(a0 - a) * x + c0 * Jx
        y += dy * dt
        x += dx * dt
        # Inelastic walls: strictly outside -> clamp and zero momentum
        mask = (x > 1.0) | (x < -1.0)
        y[mask] = 0.0
        np.clip(x, -1.0, 1.0, out=x)
        history.append(x.copy())
    return np.array(history)

def run_dbs(x0, y0, J, a0, c0, T, dt):
    x = x0.copy()
    y = y0.copy()
    history = [x.copy()]
    for step in range(1, T + 1):
        a = a0 * (step / T)
        sgn_x = np.sign(x)
        Jsgn = J @ sgn_x
        dx = a0 * y
        dy = -(a0 - a) * x + c0 * Jsgn
        y += dy * dt
        x += dx * dt
        mask = (x > 1.0) | (x < -1.0)
        y[mask] = 0.0
        np.clip(x, -1.0, 1.0, out=x)
        history.append(x.copy())
    return np.array(history)

# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Simulate SB variants on different Ising graphs.")
    parser.add_argument('--structure', type=str, default='ferro2',
                        choices=['ferro2', 'ring', 'path', 'star', 'complete_bipartite',
                                 'bipartite', 'ladder', '2dgrid', 'chimera', 'random', 'full'],
                        help='Type of coupling matrix')
    parser.add_argument('--N', type=int, default=2, help='Number of oscillators (ignored for ferro2)')
    parser.add_argument('--T', type=int, default=1000, help='Number of time steps')
    parser.add_argument('--dt', type=float, default=0.1, help='Time step size')
    parser.add_argument('--density', type=float, default=0.3, help='Edge density for random/bipartite')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--out', type=str, default='sb_figure1.png', help='Output filename')
    args = parser.parse_args()

    # If using the paper's 2-spin ferromagnetic problem, force N=2 and set paper parameters
    if args.structure == 'ferro2':
        N = 2
        a0 = 1.0
        c0 = 0.2                  # as in the paper
        # a(t) = 0.01 t, so to reach a0=1 at t = 100, we set total time = 100
        # With dt=0.1, we need T=1000 steps
        # Override if user didn't change T/dt
        if args.T == 1000 and args.dt == 0.1:
            pass  # default matches
        else:
            print(f"Using custom T={args.T}, dt={args.dt} – total time = {args.T*args.dt}")
    else:
        N = args.N
        a0 = 1.0
        c0 = 0.5 / (np.sqrt(N) * 0.5)   # default scaling, can be overridden

    # Generate coupling matrix
    J = generate_coupling(args.structure, N, args.density, args.seed)
    nnz = np.count_nonzero(J) // 2
    print(f"Structure: {args.structure}, N={N}, edges={nnz}")

    # Initial conditions: small random values around zero (as in the paper)
    np.random.seed(args.seed + 1)
    x0 = np.random.uniform(-0.1, 0.1, N)
    y0 = np.random.uniform(-0.1, 0.1, N)

    # Run simulations
    print("Running aSB...")
    hist_a = run_abs(x0, y0, J, a0, c0, args.T, args.dt)
    print("Running bSB...")
    hist_b = run_bbs(x0, y0, J, a0, c0, args.T, args.dt)
    print("Running dSB...")
    hist_d = run_dbs(x0, y0, J, a0, c0, args.T, args.dt)

    time = np.arange(0, args.T + 1) * args.dt

    # ------------------------------------------------------------------
    # Plotting: three panels, each with x1 and x2 (if N=2, else generic)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    colors = ['#1f77b4', '#ff7f0e']  # blue, orange

    titles = ['(a) aSB – adiabatic following',
              '(b) bSB – ballistic with inelastic walls',
              '(c) dSB – tunneling-like escape']

    for ax, hist, title in zip(axes, [hist_a, hist_b, hist_d], titles):
        for i in range(N):
            label = f'$x_{i+1}$' if N <= 4 else f'$x_{i+1}$'
            ax.plot(time, hist[:, i], color=colors[i % len(colors)], 
                    linewidth=2, label=label)
        ax.set_xlabel('Time')
        ax.set_ylabel('Oscillator position $x_i$')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-1.2, 1.2)
        if N <= 4:
            ax.legend()

    plt.suptitle(f'Ising graph: {args.structure} (N={N}), $c_0$={c0:.2f}', fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(args.out, dpi=150)
    plt.show()

if __name__ == '__main__':
    main()