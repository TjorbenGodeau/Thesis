import numpy as np

class GbSB:
    def __init__(self, J, dt=1.0, M=1000, A=0.0, per_spin=False):
        """
        J: (N,N) coupling matrix (diagonal elements and below are zero)
        dt: time step Δt
        M: total discrete time steps of schedule (used for p schedule)
        A: nonlinear control constant for individuql p_i. A=0 --> conventional bSB
        per_spin: if True use individual p_i for each oscillator
        """

        self.J = np.asarray(J, dtype=float)
        if self.J.ndim != 2 or self.J.shape[0] != self.J.shape[1]:
            raise ValueError("J must be a square matrix")
        self.N = self.J.shape[0]
        self.dt = float(dt)
        self.M = int(M)
        self.A = float(A)
        self.per_spin = bool(per_spin)

        #position and momentum arrays
        self.x = np.zeros(self.N, dtype=float)
        self.y = np.zeros(self.N, dtype=float)

        #bifurcation parameter p
        if self.per_spin:
            self.p = np.ones(self.N, dtype=float)   #p_i(0)=1
        else:
            self.p = 1.0                            #p(0)=1

        #internal step counter
        self.m = 0

    def initialize(self, x0=None, y0=None, p0=None):
        if x0 is not None:
            x0 = np.asarray(x0, dtype=float)
            if x0.shape != (self.N,):
                raise ValueError("x0 must have shape (N,)")
            self.x = x0.copy()
        if y0 is not None:
            y0 = np.asarray(y0, dtype=float)
            if y0.shape != (self.N,):
                raise ValueError("y0 must have shape (N,)")
            self.y = y0.copy()
        if p0 is not None:
            if self.per_spin:
                p0 = np.asarray(p0, dtype=float)
                if p0.shape != (self.N,):
                    raise ValueError("p0 must have shape (N,) for per_spin mode")
                self.p = p0.copy()
            else:
                self.p = float(p0)
        self.m = 0

    def _compute_Jx_excluding_self(self, x):
        # sum_{j != i} J_ij x_j ; assumes J diagonal is zero
        return self.J.dot(x)
    
    def step(self):
        """
        Perform one discrete update step
        """
        if self.m >= self.M:
            return  # No more updates after M steps
        
        p_cur = self.p if self.per_spin else np.full(self.N, self.p, dtype=float)

        #Eq (3)
        Jx = self._compute_Jx_excluding_self(self.x)
        y_next = self.y - self.dt * (p_cur * self.x - Jx)

        #Eq (4)
        x_next = self.x + self.dt * y_next

        #Eq (5) - inelastic walls
        mask = np.abs(x_next) > 1.0
        if np.any(mask):
            s = np.sign(x_next[mask])
            x_next[mask] = s
            y_next[mask] = 0

        #commit updates
        self.x = x_next
        self.y = y_next

        #update p
        if self.per_spin:
            #Eq (7)
            denom = max(1, self.M - self.m) #avoid divide-by-zero
            xi_sq = self.x * self.x
            self.p = self.p - ((1.0 - self.A * xi_sq) * self.p) / denom
        else:
            self.m += 1
            self.p = max(0.0, 1.0 - float(self.m) / float(self.M))
            return #already increased m
        
        self.m += 1

    def run(self, steps=None, callback=None):
        """
        Run for `steps` iterations. If steps is None, run until schedule end (M - m).
        callback(m, x, y, p) called after each step if provided.
        """
        if steps is None:
            steps = max (0, self.M - self.m)
        for _ in range(steps):
            self.step()
            if callback is not None:
                callback(self.m, self.x.copy(), self.y.copy(), (self.p.copy() if self.per_spin else float(self.p)))
        return self.x
    
    def spins(self):
        """
        Return ±1 spins as sign(x). Zeros mapped to +1
        """
        s = np.sign(self.x)
        s[s == 0] = 1.0
        return s.astype(int)
    
    def energy(self):
        s = self.spins().astype(float)
        return -0.5 * s.dot(self.J.dot(s))
    
    def hamiltonian(self):
        c = 1 / (2 * np.sqrt(self.N))
        return 0.5 * np.sum(self.y ** 2) + 0.5 * np.sum(self.p * (self.x ** 2)) - (c / 2) * self.x.dot(self.J.dot(self.x))