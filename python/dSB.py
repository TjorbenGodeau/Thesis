import numpy as np

class dSB:
    def __init__(self, J, dt = 1.0, a0 = 1.0, c0 = 1.0, Nstep = 50):
        """
        J: (N,N) coupling matrix (diagonal elements and below are zero)
        dt: time step delta t
        a0: initial value of control parameter a
        c0: positive constant value
        Nstep: total number of steps in the schedule
        """
        self.J = np.asarray(J, dtype=float)
        if self.J.ndim != 2 or self.J.shape[0] != self.J.shape[1]:
            raise ValueError("J must be a square matrix")
        self.N = self.J.shape[0]
        self.dt = float(dt)
        self.a0 = float(a0)
        self.c0 = float(c0)
        self.Nstep = int(Nstep)

        #position and momentum arrays
        self.x = np.zeros(self.N, dtype=float)
        self.y = np.zeros(self.N, dtype=float)

        #internal step counter
        self.m = 0
    
    def initialize(self, x0=None, y0=None):
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
        self.m = 0
    
    def step(self):
        if self.m >= self.Nstep:
            return #complete schedule
        
        print(f"Step {self.m}: a={self.a0 * (1 - self.m / self.Nstep)},sign x={self.signs(self.x)}, x={self.x}, y={self.y}")
        #update control parameter a
        a = self.a0 * (1 - self.m / self.Nstep)

        #update momentum y
        Jx = np.dot(self.J, self.signs(self.x))
        self.y += self.dt * (-a * self.x + self.c0 * Jx)

        #update position x
        self.x += self.dt * self.y

        #inelastic walls
        mask = np.abs(self.x) > 1.0
        self.x[mask] = self.signs(self.x[mask]) * 1.0
        self.y[mask] = 0.0

        #increment step
        self.m += 1

    def run_schedule(self, callback=None):
        for _ in range(self.Nstep):
            self.step()

    def get_state(self):
        return self.x.copy(), self.y.copy(), self.m
    
    def energy(self):
        s = self.signs(self.x)
        return -np.dot(s, np.dot(self.J, s))
    
    def hamiltonian(self):
        potential = self.a0 * (1 - self.m / self.Nstep) * np.sum(self.x**2) / 2 - self.c0 * np.dot(self.signs(self.x), np.dot(self.J, self.x))
        return potential + np.sum(self.y**2) / 2
    
    def signs(self, x):
        return np.where(x == 0, 1, np.sign(x))