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
        self.x = np.random.uniform(-0.1, 0.1, self.N)
        self.y = np.random.uniform(-0.1, 0.1, self.N)

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
        
        #print(f"Step {self.m}: a={self.a0 * (1 - self.m / self.Nstep)},sign x={self.signs(self.x)}, x={self.x}, y={self.y}")
        #print(f"Step {self.m}: a={self.a0 * (1 - self.m / self.Nstep)},x={self.x}, y={self.y}")
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

    def run_schedule(self):
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

def quantize(val, bits, v_min=-1.0, v_max=1.0):
    """
    quantize 'val' to a fixed-point grid with 'bits' bits over [v_min, v_max].
    """
    if bits is None:    #None means full (float64) precision
        return val
    bits = int(bits)
    if bits < 1:
        raise ValueError("bits must be >= 1")
    levels = 2 ** bits - 1
    val = np.asarray(val, dtype=float)
    val_clipped = np.clip(val, v_min, v_max)
    val_q = np.round((val_clipped - v_min) / (v_max - v_min) * levels)
    return val_q / levels * (v_max - v_min) + v_min

class dSB_quantized:
    def __init__(self, J, dt=1.0, a0=1.0, c0=1.0, Nstep=50, bits_J=None, bits_x=None, bits_y=None, bits_a=None, J_range=1.0, y_range=10.0):
        self.J_raw = np.asarray(J, dtype=float)
        if self.J_raw.ndim != 2 or self.J_raw.shape[0] != self.J_raw.shape[1]:
            raise ValueError("J must be a square matrix")
        self.N = self.J_raw.shape[0]

        self.dt = float(dt)
        self.a0 = float(a0)
        self.c0 = float(c0)
        self.Nstep = int(Nstep)

        self.bits_J = bits_J
        self.bits_x = bits_x
        self.bits_y = bits_y
        self.bits_a = bits_a

        self.J_range = float(J_range)
        self.y_range = float(y_range)

        self.J = quantize(self.J_raw, bits_J, -self.J_range, self.J_range)

        self.x = quantize(np.random.uniform(-0.1, 0.1, self.N), bits_x, -1.0, 1.0)
        self.y = quantize(np.random.uniform(-0.1, 0.1, self.N), bits_y, -self.y_range, self.y_range)
        self.m = 0

    def _qx(self, v):
        return quantize(v, self.bits_x, -1.0, 1.0)
    
    def _qy(self, v):
        return quantize(v, self.bits_y, -self.y_range, self.y_range)
    
    def _qa(self, v):
        return quantize(v, self.bits_a, 0.0, self.a0)
    
    @staticmethod
    def signs(x):
        return np.where(x==0, 1, np.sign(x))
    
    def initialize(self, x0=None, y0=None):
        if x0 is not None:
            x0 = np.asarray(x0, dtype=float)
            if x0.shape != (self.N,):
                raise ValueError("x0 must have shape (N,)")
            self.x = self._qx(x0.copy())
        if y0 is not None:
            y0 = np.asarray(y0, dtype=float)
            if y0.shape != (self.N,):
                raise ValueError("y0 must have shape (N,)")
            self.y = self._qy(y0.copy())

        self.m = 0

    def step(self):
        if self.m >= self.Nstep:
            return
        
        a = self._qa(self.a0 * (1 - self.m / self.Nstep))

        Jx = np.dot(self.J, self.signs(self.x))
        self.y = self._qy(self.y + self.dt * (-a * self.x + self.c0 * Jx))

        self.x = self._qx(self.x + self.dt * self.y)

        mask = np.abs(self.x) >= 1.0
        self.x[mask] = self.signs(self.x[mask]) * 1.0
        self.y[mask] = 0.0

        self.m += 1
    
    def run_schedule(self, callback=None):
        for _ in range(self.Nstep):
            self.step()
            if callback is not None:
                callback(self)

    def get_state(self):
        """Return (x, y, m) — positions, momenta, step index."""
        return self.x.copy(), self.y.copy(), self.m

    def energy(self):
        """Ising energy E = -s^T J s  (uses quantized J)."""
        s = self.signs(self.x)
        return -np.dot(s, np.dot(self.J, s))

    def hamiltonian(self):
        """
        SB Hamiltonian:
            H = Σ y²/2  +  a(t)/2 Σ x²  −  c0 sign(x)^T J x
        """
        a = self._qa(self.a0 * (1.0 - self.m / self.Nstep))
        potential = (
            a * np.sum(self.x ** 2) / 2.0
            - self.c0 * np.dot(self.signs(self.x), np.dot(self.J, self.x))
        )
        return potential + np.sum(self.y ** 2) / 2.0