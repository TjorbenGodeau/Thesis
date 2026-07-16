#!/usr/bin/env python3
import os
import math
import random

# --- 1. Target Problems & Best Known Cuts ---
TARGET_PROBLEMS = [
    (1, "G1.txt", 800, 11624),
    (63, "G63.txt", 7000, 27047),
    (72, "G72.txt", 10000, 7008) 
    # Add G77 (14000, 9940) and G81 (20000, 14060) when ready
]

def to_signed(val, bits):
    """Simulates hardware register overflow/wraparound for N-bit signed integers."""
    val = val & ((1 << bits) - 1)
    if val >= (1 << (bits - 1)):
        val -= (1 << bits)
    return val

class GsetGraph:
    def __init__(self, filepath):
        self.edges = []
        self.adj = {}
        
        with open(filepath, 'r') as f:
            lines = f.readlines()
            header = lines[0].strip().split()
            self.num_nodes = int(header[0])
            self.num_edges = int(header[1])
            
            for i in range(self.num_nodes):
                self.adj[i] = []
                
            for line in lines[1:]:
                parts = line.strip().split()
                if len(parts) >= 3:
                    u, v, w = int(parts[0])-1, int(parts[1])-1, int(parts[2])
                    self.edges.append((u, v, w))
                    self.adj[u].append((v, w))
                    self.adj[v].append((u, w))

class DSB_GoldenModel:
    def __init__(self, graph, nstep, xy_frac, xy_int=2, a_bits=16):
        self.graph = graph
        self.nstep = nstep
        self.xy_frac = xy_frac
        self.xy_w = xy_int + xy_frac  # Total bit width (usually 16)
        self.a_bits = a_bits
        
        # Dynamically scale constants to match the fractional precision!
        self.ONE_FP = 1 << xy_frac
        # A0_FP is locked to Q14 scale because A_BITS dictates a fixed >> 15 shift in SV
        self.A0_FP = int(1.0 * (1 << 14)) 
        self.DT_FP = int(0.01 * self.ONE_FP)
        self.C0_FP = int(0.025 * self.ONE_FP)
        
        # Initialize x with small random perturbations (~0.006 magnitude)
        random.seed(12345)
        seed_limit = int(0.006 * self.ONE_FP)
        if seed_limit < 1: seed_limit = 1
        self.x = [random.randint(-seed_limit, seed_limit) for _ in range(graph.num_nodes)]
        self.y = [0] * graph.num_nodes
        self.signs = [1 if val >= 0 else -1 for val in self.x]

    def schedule_a(self, m):
        if m >= self.nstep:
            return 0
        numerator = self.A0_FP * (self.nstep - m)
        return int(numerator // self.nstep)

    def run(self):
        for m in range(self.nstep):
            a_m = self.schedule_a(m)
            self.signs = [1 if val >= 0 else -1 for val in self.x]
            
            new_x = [0] * self.graph.num_nodes
            new_y = [0] * self.graph.num_nodes
            
            for i in range(self.graph.num_nodes):
                Jx_i = sum(w * self.signs[v] for v, w in self.graph.adj[i])
                
                # Fixed-point math mimicking update_unit.sv
                ax_full = self.x[i] * a_m
                ax_shifted = ax_full >> (self.a_bits - 1)
                
                c0jx_full = Jx_i * self.C0_FP
                force_wide = -ax_shifted + c0jx_full
                
                dy_full = force_wide * self.DT_FP
                delta_y = dy_full >> self.xy_frac
                
                # Apply hardware register wrap-around
                y_pre = to_signed(self.y[i] + delta_y, self.xy_w)
                
                dx_full = y_pre * self.DT_FP
                x_pre = to_signed(self.x[i] + (dx_full >> self.xy_frac), self.xy_w)
                
                # Wall hits
                if x_pre > self.ONE_FP:
                    new_x[i] = self.ONE_FP
                    new_y[i] = 0
                elif x_pre < -self.ONE_FP:
                    new_x[i] = -self.ONE_FP
                    new_y[i] = 0
                else:
                    new_x[i] = x_pre
                    new_y[i] = y_pre
                    
            self.x = new_x
            self.y = new_y

    def compute_cut(self):
        final_signs = [1 if val >= 0 else -1 for val in self.x]
        cut = 0
        for u, v, w in self.graph.edges:
            if final_signs[u] != final_signs[v]:
                cut += w
        return cut

def calculate_hardware_cost(graph, k_max):
    total_chunks = 0
    for i in range(graph.num_nodes):
        degree = len(graph.adj[i])
        total_chunks += math.ceil(degree / k_max) if degree > 0 else 0
    # ~15 states in SLC FSM per chunk
    return total_chunks * 15 

def main():
    print("=== Python Golden Model: Hardware & Algorithmic Sweep ===")
    
    SWEEP_NSTEP = [30, 50]
    SWEEP_XY_FRAC = [14, 12]  
    SWEEP_K_MAX = [16, 32, 64] 
    
    header = f"{'Graph':<5} | {'Nstep':<5} | {'FracBit':<7} | {'K_MAX':<5} | {'Cut':<6} | {'Best':<6} | {'Accuracy':<8} | {'Est. Cycles/Step':<16}"
    print(header)
    print("-" * len(header))

    for gid, gfile, nodes, best_known in TARGET_PROBLEMS:
        if not os.path.exists(gfile):
            print(f"G{gid:<4} | [FILE NOT FOUND: {gfile}]")
            continue
            
        graph = GsetGraph(gfile)
        
        for nstep in SWEEP_NSTEP:
            for frac in SWEEP_XY_FRAC:
                # 1. Algorithmic Simulation
                model = DSB_GoldenModel(graph, nstep=nstep, xy_frac=frac)
                model.run()
                cut = model.compute_cut()
                accuracy = (cut / best_known) * 100
                
                # 2. Hardware Estimate
                for k_max in SWEEP_K_MAX:
                    est_cycles = calculate_hardware_cost(graph, k_max)
                    
                    row = (f"G{gid:<4} | {nstep:<5} | {frac:<7} | {k_max:<5} | "
                           f"{cut:<6} | {best_known:<6} | {accuracy:>6.2f}% | {est_cycles:<16,}")
                    print(row)
        print("-" * len(header))

if __name__ == "__main__":
    main()