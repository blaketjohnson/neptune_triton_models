import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define constants for Neptune and Triton
m1 = 1.024e26  # Mass of Neptune in kg
m2 = 2.14e22   # Mass of Triton in kg
mu = m2 / (m1 + m2)  # Mass ratio

# Define symbolic variables
x, y, z = sp.symbols('x y z')

# Define the distances r1 and r2
r1 = sp.sqrt((x + mu)**2 + y**2 + z**2)
r2 = sp.sqrt((x + mu - 1)**2 + y**2 + z**2)

# Define the effective potential function Omega
Omega = 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2

# Calculate the partial derivatives
Omega_x = sp.diff(Omega, x)
Omega_y = sp.diff(Omega, y)
Omega_z = sp.diff(Omega, z)

# Solve for the Lagrange points (where the gradients are zero)
solutions = sp.solve([Omega_x, Omega_y, Omega_z], (x, y, z))

# Evaluate solutions numerically
lagrange_points = [(float(sol[0]), float(sol[1]), float(sol[2])) for sol in solutions]

# Add L4 and L5 points manually in 3D (z=0)
L4 = (0.5 - mu, np.sqrt(3) / 2, 0)
L5 = (0.5 - mu, -np.sqrt(3) / 2, 0)
lagrange_points.extend([L4, L5])

# Print Lagrange points
for i, point in enumerate(lagrange_points):
    print(f"L{i+1}: {point}")

# Define the equations of motion
def equations(t, state):
    x, y, z, vx, vy, vz = state
    
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    
    Omega_x = x - (1 - mu) * (x + mu) / r1**3 - mu * (x + mu - 1) / r2**3
    Omega_y = y - (1 - mu) * y / r1**3 - mu * y / r2**3
    Omega_z = z - (1 - mu) * z / r1**3 - mu * z / r2**3
    
    ax = Omega_x + 2 * vy
    ay = Omega_y - 2 * vx
    az = Omega_z
    
    return [vx, vy, vz, ax, ay, az]

# Initial conditions
x0 = 0.5
y0 = 0
z0 = 0
vx0 = 0
vy0 = 1
vz0 = 0
state0 = [x0, y0, z0, vx0, vy0, vz0]

# Time span for the solution
t_span = (0, 50)
t_eval = np.linspace(*t_span, 1000)

# Solve the system of equations
solution = solve_ivp(equations, t_span, state0, t_eval=t_eval, method='RK45')

# Extract the results
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]

# Plot the trajectory
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.plot(x, y, z, label='Trajectory')
ax1.scatter([0], [0], [0], color='b', s=100, label='Neptune')
ax1.scatter([1], [0], [0], color='g', s=50, label='Triton')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')
ax1.legend()
ax1.set_title('Trajectory in the Restricted Three Body Problem')
plt.show()

# Plot the 3D surface defined by C = 2U and Lagrange points
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')

# Generate mesh grid for 3D plot
X = np.linspace(-1.5, 1.5, 100)
Y = np.linspace(-1.5, 1.5, 100)
X, Y = np.meshgrid(X, Y)
Z = np.linspace(-1.5, 1.5, 100)

# Compute U for each point in the grid
R1 = np.sqrt((X + mu)**2 + Y**2 + Z[:, np.newaxis, np.newaxis]**2)
R2 = np.sqrt((X + mu - 1)**2 + Y**2 + Z[:, np.newaxis, np.newaxis]**2)
U = 0.5 * (X**2 + Y**2 + Z[:, np.newaxis, np.newaxis]**2) + (1 - mu) / R1 + mu / R2

# Plot the surface where C = 2U
C = 2 * U
ax2.contour3D(X, Y, C[50], 50, cmap='viridis')

# Plot Lagrange points
for point in lagrange_points:
    ax2.scatter(point[0], point[1], 2 * (0.5 * (point[0]**2 + point[1]**2 + point[2]**2) + (1 - mu) / np.sqrt((point[0] + mu)**2 + point[1]**2 + point[2]**2) + mu / np.sqrt((point[0] + mu - 1)**2 + point[1]**2 + point[2]**2)), color='r', s=100)

ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('z')
ax2.set_title('3D Surface of C = 2U with Lagrange Points')
plt.show()
