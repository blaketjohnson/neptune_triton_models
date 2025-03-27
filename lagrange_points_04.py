import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import root

# Import equations from parts 1 and 2
sys.path.append('/Users/blakejohnson/Documents/rtbp_neptune_triton/neptune_triton_models')
from jacobi_contant_03 import *

def find_Lagrange_points(mu):
    def find_L1(mu):
        initial_guess = np.array([0.5 - mu, 0.0, 0.0])
        sol = root(lambda coords: gradient_U(*coords, mu), initial_guess)
        if sol.success:
            return sol.x
        else:
            raise ValueError(f"Failed to find L1: {sol.message}")

    def find_L2(mu):
        initial_guess = np.array([1.5 - mu, 0.0, 0.0])
        sol = root(lambda coords: gradient_U(*coords, mu), initial_guess)
        if sol.success:
            return sol.x
        else:
            raise ValueError(f"Failed to find L2: {sol.message}")

    def find_L3(mu):
        initial_guess = np.array([-1.0 - mu, 0.0, 0.0])
        sol = root(lambda coords: gradient_U(*coords, mu), initial_guess)
        if sol.success:
            return sol.x
        else:
            raise ValueError(f"Failed to find L3: {sol.message}")

    def find_L4(mu):
        initial_guess = np.array([0.5 - mu, np.sqrt(3)/2, 0.0])
        sol = root(lambda coords: gradient_U(*coords, mu), initial_guess)
        if sol.success:
            return sol.x
        else:
            raise ValueError(f"Failed to find L4: {sol.message}")

    def find_L5(mu):
        initial_guess = np.array([0.5 - mu, -np.sqrt(3)/2, 0.0])
        sol = root(lambda coords: gradient_U(*coords, mu), initial_guess)
        if sol.success:
            return sol.x
        else:
            raise ValueError(f"Failed to find L5: {sol.message}")

    L1_position = find_L1(mu)
    L2_position = find_L2(mu)
    L3_position = find_L3(mu)
    L4_position = find_L4(mu)
    L5_position = find_L5(mu)

    return L1_position, L2_position, L3_position, L4_position, L5_position

# Calculate Lagrange points
L1, L2, L3, L4, L5 = find_Lagrange_points(mu)

if __name__ == '__main__':
    # Plot Neptune, Triton, and all Lagrange points in 2D
    neptune_position = np.array([0.0, 0.0])
    triton_position = np.array([1.0, 0.0])

    plt.figure(figsize=(10, 8))
    plt.scatter(neptune_position[0], neptune_position[1], color='blue', label='Neptune')
    plt.scatter(triton_position[0], triton_position[1], color='saddlebrown', label='Triton')
    plt.scatter(L1[0], L1[1], color='black', label='L1')
    plt.scatter(L2[0], L2[1], color='black', label='L2')
    plt.scatter(L3[0], L3[1], color='black', label='L3')
    plt.scatter(L4[0], L4[1], color='black', label='L4')
    plt.scatter(L5[0], L5[1], color='red', label='L5')

    plt.text(neptune_position[0], neptune_position[1], ' Neptune', fontsize=12, ha='center', va='bottom')
    plt.text(triton_position[0], triton_position[1], ' Triton', fontsize=12, ha='center', va='bottom')
    plt.text(L1[0], L1[1], 'L1', fontsize=12, ha='center', va='top')
    plt.text(L2[0], L2[1], 'L2', fontsize=12, ha='center', va='top')
    plt.text(L3[0], L3[1], 'L3', fontsize=12, ha='center', va='top')
    plt.text(L4[0], L4[1], 'L4', fontsize=12, ha='center', va='top')
    plt.text(L5[0], L5[1], 'L5', fontsize=12, ha='center', va='top')

    plt.xlabel('Distance in unitless dimensions')
    plt.ylabel('Distance in unitless dimensions')
    plt.title('Positions of Lagrange Points with Neptune and Triton')
    plt.legend()
    plt.grid(True)

    # Uncomment these lines if you need to zoom out
    # plt.xlim(-1.5, 2.5)
    # plt.ylim(-1.5, 1.5)

    plt.show()
    



