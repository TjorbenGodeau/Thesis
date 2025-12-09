from __future__ import annotations
import sys
from typing import Optional, Tuple, List
import numpy as np

"""
Reads a text document describing a (square) matrix in coordinate form and returns
(or writes) the matrix.

Expected input format (example):
4 3
1 1 2.5
1 3 3.0
4 2 -1

First (non-empty, non-comment) line: two numbers: N M (N is dimension, M is number of non-zero entries)
Following lines: row_index column_index value

By default, indices are 1-based (use --index-base 0 if your file uses 0-based indices).
"""

def parse_first_line(tokens: List[str]) -> Tuple[int, Optional[int]]:
    """
    Parse first line tokens to extract dimension and optional nonzero count.
    
    Return (N, M or None)
    """
    if not tokens:
        raise ValueError("First line is empty or invalid.")
    if len(tokens) == 1:
        n = int(tokens[0])
        return n, None
    n = int(tokens[0])
    m = int(tokens[1])
    return n, m

def read_matrix_from_file(path: str, index_base: int = 1, dtype = float):
    """
    Read coordinate format matrix file and return a NxN matrix.

    - path: input file path
    - index_base: 1 for 1-based indexing, 0 for 0-based indexing
    - dtype: type for the matrix entries (float by default)
    """
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f]

    # Filter out empty lines and lines starting with '#'
    content_lines = [ln for ln in lines if ln and not ln.startswith("#")]
    if not content_lines: 
        raise ValueError("Input file has no non-empty, non-comment lines.")
    
    # First non-empty line is header: N and optional M
    header_tokens = [t for t in content_lines[0].replace(","," ").split() if t]
    n, m_expected = parse_first_line(header_tokens)

    # Initialize matrix of zeros
    matrix = np.zeros((n, n), dtype=dtype)

    # Remaining lines (coordinates)
    entries = []
    for line in content_lines[1:]:
        tokens = [t for t in line.replace(",", " ").split() if t]
        if len(tokens) < 3:
            raise ValueError(f"Line does not have 3 tokens: '{line}'")
        r = int(tokens[0]) - index_base
        c = int(tokens[1]) - index_base
        try: 
            v = dtype(tokens[2])
        except Exception:
            # in case of int
            v = dtype(float(tokens[2]))
        if r < 0 or r >= n or c < 0 or c >= n:
            raise IndexError(f"Index out of bounds in line: '{line}' (converted to 0-based: {r},{c})")
        entries.append((r, c, v))
        matrix[r, c] = v
    
    if m_expected is not None and m_expected != len(entries):
        print(f"Warning: header expected {m_expected} non-zero entries, " f"but {len(entries)} were read.", file=sys.stderr)
    
    return matrix, n
        
def write_matrix_to_csv(matrix, output_path: str):
    """
    Write matrix (numpy array) to CSV file.
    """
    fmt = "%s"
    np.savetxt(output_path, matrix, delimiter=",", fmt=fmt)


def write_matrix_to_txt(matrix, output_path: str):
    """
    Write matrix as whitespace-separated values
    """
    with open(output_path, "w", encoding="utf-8") as f:
        for row in matrix:
            f.write(" ".join(str(x) for x in row) + "\n")

def compute_cut_value(matrix, energy):
    """
    Computes the cut value of the maxCut problem with current IC-matrix (J) and energy.
    
    :param matrix: Square matrix with interaction coeficients
    :param energy: Ising energy
    """
    total_weight = np.sum(matrix)
    return (total_weight - energy) / 2

def main():
    matrix = read_matrix_from_file(r"C:\Users\tjorb\Documents\Thesis\benchmark\G_Set\G1.txt")
    print("Matrix has been extracted")
    write_matrix_to_txt(matrix, r"C:\Users\tjorb\Documents\Thesis\output_files\G1.txt")
    print("Matrix has been written to text file")

if __name__ == "__main__":
    main()