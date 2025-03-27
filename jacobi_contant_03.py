import numpy as np
import matplotlib.pyplot as plt
import sys

# Import equations from parts 1 and 2
sys.path.append('/Users/blakejohnson/Documents/Thesis/Three Body Problem')
from r3bp_calculations_02 import *

"""
Solving for the Jacobi Constant
"""

def jacobi_constant(state, mu):
    x, y, z, vx, vy, vz = state
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    U = 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2

    # Jacobi constant formula
    C = 2 * U - (vx**2 + vy**2 + vz**2)
    return C

# Calculate Jacobi constant for initial conditions
C = jacobi_constant(state0, mu)

if __name__ == '__main__':
    print(f"Jacobi constant: {C}")

