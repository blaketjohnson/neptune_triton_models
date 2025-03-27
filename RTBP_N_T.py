import numpy as np
from scipy.optimize import fsolve
from mpl_toolkits.mplot3d import Axes3D


""" Part 1 - Equations of Motion"""

# Constants for Neptune and Triton
mass_neptune = 1.024e26  # kg
mass_triton = 2.14e22    # kg
G = 6.67430e-11          # m^3 kg^-1 s^-2

# Compute mu, mu1, and mu2
mu = mass_triton / (mass_neptune + mass_triton)
mu1 = 1 - mu
mu2 = mu

print(f"mu: {mu}, mu1: {mu1}, mu2: {mu2}")

# Calculate r values
def calculate_distances(x, y, mu):
    r1 = np.sqrt((x + mu)**2 + y**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2)
    return r1, r2

# Calculate Effective potential
def effective_potential(x, y, mu):
    r1, r2 = calculate_distances(x, y, mu)
    U = -mu1 / r1 - mu2 / r2 - 0.5 * (x**2 + y**2)
    return U

# Calculate Acceleration
def acceleration(x, y, vx, vy, mu):
    r1, r2 = calculate_distances(x, y, mu)
    
    Ux = mu1 * (x + mu) / r1**3 + mu2 * (x - 1 + mu) / r2**3 + x
    Uy = mu1 * y / r1**3 + mu2 * y / r2**3 + y
    
    ax = 2 * vy + Ux
    ay = -2 * vx + Uy
    
    return ax, ay

import matplotlib.pyplot as plt

# Function to simulate the motion
def simulate_motion(mu, t_max, dt):
    # Initial conditions (starting at Triton's perigee)
    x, y = 0.3548e6, 0.0  # km
    vx, vy = 0.0, 4.39    # km/s (orbital speed at perigee)
    
    # Time array
    t = np.arange(0, t_max, dt)
    
    # Arrays to store positions
    x_vals = []
    y_vals = []
    
    for _ in t:
        ax, ay = acceleration(x, y, vx, vy, mu)
        
        # Update velocities
        vx += ax * dt
        vy += ay * dt
        
        # Update positions
        x += vx * dt
        y += vy * dt
        
        # Store positions
        x_vals.append(x)
        y_vals.append(y)
    
    return x_vals, y_vals


# Simulate the motion
t_max = 50000  # seconds
dt = 10        # time step in seconds

x_vals, y_vals = simulate_motion(mu, t_max, dt)

# Plot the orbits
plt.figure(figsize=(10, 10))
plt.plot(x_vals, y_vals, label='Satellite Orbit')
plt.scatter([0], [0], color='blue', label='Neptune')
plt.scatter([354800], [0], color='red', label='Triton')
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.title('Orbits in the RTBP')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()


"""Part 2 - Jacobi Constant"""

def jacobi_constant(x, y, vx, vy, mu):
    U = effective_potential(x, y, mu)
    C_j = 2 * U - (vx**2 + vy**2)
    return C_j

# Example calculation for the initial conditions
x, y = 0.3548e6, 0.0  # km
vx, vy = 0.0, 4.39    # km/s

C_j = jacobi_constant(x, y, vx, vy, mu)
print(f"Jacobi Constant: {C_j}")

# Find the Crit Values and Zero Velocity Curves

def lagrange_points(mu):
    # Equations for Lagrange points L1, L2, and L3
    def equations_L1L2L3(x):
        r1 = np.abs(x + mu)
        r2 = np.abs(x - 1 + mu)
        return x - mu1 * (x + mu) / r1**3 - mu2 * (x - 1 + mu) / r2**3

    L1 = fsolve(equations_L1L2L3, 0.8)[0]
    L2 = fsolve(equations_L1L2L3, 1.2)[0]
    L3 = fsolve(equations_L1L2L3, -1.0)[0]
    
    # Equations for L4 and L5
    L4_x = 0.5 - mu
    L4_y = np.sqrt(3) / 2
    L5_x = 0.5 - mu
    L5_y = -np.sqrt(3) / 2
    
    return (L1, 0), (L2, 0), (L3, 0), (L4_x, L4_y), (L5_x, L5_y)

L_points = lagrange_points(mu)
print(f"Lagrange Points: {L_points}")

# Calculate the critical Jacobi constants
Cj_values = [jacobi_constant(x, y, 0, 0, mu) for x, y in L_points]
print(f"Critical Jacobi Constants: {Cj_values}")

def zero_velocity_curve(Cj, mu, xlim, ylim, res=100):
    x = np.linspace(xlim[0], xlim[1], res)
    y = np.linspace(ylim[0], ylim[1], res)
    X, Y = np.meshgrid(x, y)
    
    Z = 2 * effective_potential(X, Y, mu) - Cj
    
    plt.contour(X, Y, Z, levels=[0], colors='r')
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    plt.title('Zero-Velocity Curve')
    plt.grid(True)
    plt.axis('equal')

# Plot zero-velocity curves for each critical Jacobi constant

plt.figure(figsize=(10, 10))
for Cj in Cj_values:
    zero_velocity_curve(Cj, mu, (-2, 2), (-2, 2))

plt.scatter([0], [0], color='blue', label='Neptune')
plt.scatter([0.3548e6], [0], color='red', label='Triton')
plt.legend()
plt.show()

""" Part 3 - Lagrange Points """

plt.figure(figsize=(10, 10))
for Cj in Cj_values:
    zero_velocity_curve(Cj, mu, (-2, 2), (-2, 2))

# Plot Lagrange points
L_points = lagrange_points(mu)
L_x, L_y = zip(*L_points)
plt.scatter(L_x, L_y, color='green', label='Lagrange Points')

# Plot Neptune and Triton
plt.scatter([0], [0], color='blue', label='Neptune')
plt.scatter([0.3548e6], [0], color='red', label='Triton')

# Orbit of Triton
theta = np.linspace(0, 2*np.pi, 100)
x_orbit = 0.3548e6 * np.cos(theta)
y_orbit = 0.3548e6 * np.sin(theta)
plt.plot(x_orbit, y_orbit, 'b--', label='Triton Orbit')

plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.title('Lagrange Points and Zero-Velocity Curves')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()



# Plot 3D
def plot_3d_surface(mu, xlim, ylim, res=100):
    x = np.linspace(xlim[0], xlim[1], res)
    y = np.linspace(ylim[0], ylim[1], res)
    X, Y = np.meshgrid(x, y)
    Z = effective_potential(X, Y, mu)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
    
    # Plot Lagrange points
    L_points = lagrange_points(mu)
    L_x, L_y = zip(*L_points)
    L_z = [effective_potential(x, y, mu) for x, y in L_points]
    ax.scatter(L_x, L_y, L_z, color='red', s=50, label='Lagrange Points')
    
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_zlabel('U')
    ax.set_title('3D Surface of Effective Potential')
    plt.legend()
    plt.show()

plot_3d_surface(mu, (-2, 2), (-2, 2))

""" Part 4 - Stability of the Equilibrium """

def stability_analysis(L_points, mu):
    # Calculate the Jacobian matrix at each Lagrange point
    for x, y in L_points:
        r1 = np.sqrt((x + mu)**2 + y**2)
        r2 = np.sqrt((x - 1 + mu)**2 + y**2)

        dUdx = mu1 * (x + mu) / r1**3 + mu2 * (x - 1 + mu) / r2**3 + x
        dUdy = mu1 * y / r1**3 + mu2 * y / r2**3 + y
        d2Udx2 = mu1 * (3 * (x + mu)**2 - r1**2) / r1**5 + mu2 * (3 * (x - 1 + mu)**2 - r2**2) / r2**5 + 1
        d2Udy2 = mu1 * (3 * y**2 - r1**2) / r1**5 + mu2 * (3 * y**2 - r2**2) / r2**5 + 1
        d2Udxdy = 3 * y * (mu1 * (x + mu) / r1**5 + mu2 * (x - 1 + mu) / r2**5)
        
        jacobian = np.array([[0, 2, 0, 0],
                             [-2, 0, 0, 0],
                             [d2Udx2, d2Udxdy, 0, 2],
                             [d2Udxdy, d2Udy2, -2, 0]])
        
        eigenvalues = np.linalg.eigvals(jacobian)
        print(f"Eigenvalues at ({x:.4f}, {y:.4f}): {eigenvalues}")

L_points = lagrange_points(mu)
stability_analysis(L_points, mu)

# plot the stability

def plot_stability_regions(L_points, mu, xlim, ylim, res=50):
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for i, (x, y) in enumerate(L_points):
        X, Y = np.meshgrid(np.linspace(x - 0.05, x + 0.05, res), np.linspace(y - 0.05, y + 0.05, res))
        Ux = mu1 * (X + mu) / ((X + mu)**2 + Y**2)**(3/2) + mu2 * (X - 1 + mu) / ((X - 1 + mu)**2 + Y**2)**(3/2) + X
        Uy = mu1 * Y / ((X + mu)**2 + Y**2)**(3/2) + mu2 * Y / ((X - 1 + mu)**2 + Y**2)**(3/2) + Y

        ax = axes[i]
        ax.quiver(X, Y, Ux, Uy, color='purple', alpha=0.6)
        ax.scatter([x], [y], color='red', label=f'L{i+1}')
        ax.set_xlim([x - 0.05, x + 0.05])
        ax.set_ylim([y - 0.05, y + 0.05])
        ax.set_title(f'Stability Region around L{i+1}')
        ax.grid(True)
        ax.legend()
    
    plt.tight_layout()
    plt.show()

plot_stability_regions(L_points, mu, (-2, 2), (-2, 2))
