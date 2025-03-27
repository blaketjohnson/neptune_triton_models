import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp
import sys
from scipy.integrate import solve_ivp

# Import constants from your specific path
#sys.path.append('/Users/blakejohnson/Documents/Thesis/Three Body Problem')
from r2bp_calculations_01 import *


"""
Define the constants and calculate mu for the three body problem
"""

# Define constants
mu = m2 / (m1 + m2)  # Mass ratio
mu1 = 1 - mu
mu2 = mu  

"""
2.2 Finding the Equation for the Gradient

Define the symbols and potential energy function U. Compute the gradient.
"""

# Define symbols
x, y, z, u = sp.symbols('x y z u')

# Define distances r1 and r2
r1 = sp.sqrt((x + u)**2 + y**2 + z**2)
r2 = sp.sqrt((x + u - 1)**2 + y**2 + z**2)

# Define the effective potential function U
U = (1/2) * (x**2 + y**2 + z**2) + ((1 - u) / r1) + (u / r2)

# Calculate the partial derivatives
U_x = sp.diff(U, x)
U_y = sp.diff(U, y)
U_z = sp.diff(U, z)

if __name__ == '__main__':

    print("Effective Potential U:")
    sp.pretty_print(U)
    print("\nPartial derivatives:")
    print("U_x:")
    sp.pretty_print(U_x)
    print("\nU_y:")
    sp.pretty_print(U_y)
    print("\nU_z:")
    sp.pretty_print(U_z)

# Convert symbolic expressions to Python functions for numerical evaluation
U_x_func = sp.lambdify((x, y, z, u), U_x, 'numpy')
U_y_func = sp.lambdify((x, y, z, u), U_y, 'numpy')
U_z_func = sp.lambdify((x, y, z, u), U_z, 'numpy')

"""
2.3 Calculating the Equations of Motion

Define the equations of motion using the calculated gradients.
"""

# Potential Energy
def U(x, y, z, mu):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    return 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2

# Define the gradients of the effective potential numerically
def gradient_U(x, y, z, mu):
    Ux = U_x_func(x, y, z, mu)
    Uy = U_y_func(x, y, z, mu)
    Uz = U_z_func(x, y, z, mu)
    return np.array([Ux, Uy, Uz])

# Define the equations of motion
def equations(t, state, mu):
    x, y, z, vx, vy, vz = state
    
    # Calculate partial derivatives
    Ux, Uy, Uz = gradient_U(x, y, z, mu)
    
    # Acceleration components
    ax = Ux + 2 * vy
    ay = Uy - 2 * vx
    az = Uz
    
    return [vx, vy, vz, ax, ay, az]

"""
Set initial conditions and solve the system of equations.
"""

# Initial conditions
x0 = 0.5
y0 = 0
z0 = 0
vx0 = 0
vy0 = 1
vz0 = 0
state0 = [x0, y0, z0, vx0, vy0, vz0]

# Time span for the solution
t_span = (0, 100)
t_eval = np.linspace(*t_span, 1000)

# Solve the system of equations
solution = solve_ivp(equations, t_span, state0, t_eval=t_eval, method='RK45', args=(mu,))

# Extract the results
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]

"""
2.4 Plot The RTBP
Plot the satellite trajectory and mark Neptune and Triton.
"""

# Plot the trajectory
def plot_trajectory(x, y, z):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, color='green', label='Satellite Trajectory')
    ax.set_xlabel('x (non-dimensional units)')
    ax.set_ylabel('y (non-dimensional units)')
    ax.set_zlabel('z (non-dimensional units)')
    ax.set_title('Trajectory in the Restricted Three Body Problem')

    # Add markers for Neptune and Triton
    ax.scatter([-mu], [0], [0], color='blue', s=100, label='Neptune')
    ax.scatter([1 - mu], [0], [0], color='saddlebrown', s=50, label='Triton')

    # Add legend
    ax.legend()
    plt.show()

# Main function execution
if __name__ == '__main__':
    plot_trajectory(x, y, z)

