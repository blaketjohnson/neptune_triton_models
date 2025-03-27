import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Co6nstants and parameters
m1 = 1.024e26 # Mass of Neptune in kg
m2 = 2.14e22   # Mass of Triton in kg
distance_km = 354800  # Distance between Neptune and Triton in km
au_conversion = 6.6846e-9  # Conversion factor from km to AU
distance = distance_km * au_conversion  # Distance in AU
mu = m2 / (m1 + m2)  # Mass ratio

# Define the effective potential and its gradients
def U(x, y, z, mu):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    return 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2

# Define the gradients of the effective potential
def gradient_U(x, y, z, mu):
    r1 = (x + mu)**2 + y**2 + z**2
    r2 = (x + mu - 1)**2 + y**2 + z**2
    
    Ux = x - (1 - mu) * (x + mu) / r1**1.5 - mu * (x + mu - 1) / r2**1.5
    Uy = y - (1 - mu) * y / r1**1.5 - mu * y / r2**1.5
    Uz = -(1 - mu) * z / r1**1.5 - mu * z / r2**1.5
    
    return np.array([Ux, Uy, Uz])

# Function to find Lagrange points
def find_Lagrange_points(mu):
    def find_L1(mu):
        initial_guess = np.array([0.5 - mu, 0.0, 0.0])
        sol = root(lambda coords: gradient_U(*coords, mu), initial_guess)
        if sol.success:
            return sol.x
        else:
            raise ValueError(f"Failed to find L1: {sol.message}")

    def find_L2(mu):
        initial_guess = np.array([0.5 - mu, 0.0, 0.0])  # Adjusted initial guess
        sol = root(lambda coords: gradient_U(*coords, mu), initial_guess)
        if sol.success:
            return sol.x
        else:
            raise ValueError(f"Failed to find L2: {sol.message}")

    def find_L3(mu):
        initial_guess = np.array([-1 - mu, 0.0, 0.0])  # Adjusted initial guess
        sol = root(lambda coords: gradient_U(*coords, mu), initial_guess)
        if sol.success:
            return sol.x
        else:
            raise ValueError(f"Failed to find L3: {sol.message}")

    def find_L4(mu):
        initial_guess = np.array([1 - mu, np.sqrt(3)/2, 0.0])
        sol = root(lambda coords: gradient_U(*coords, mu), initial_guess)
        if sol.success:
            return sol.x
        else:
            raise ValueError(f"Failed to find L4: {sol.message}")

    def find_L5(mu):
        initial_guess = np.array([1 - mu, -np.sqrt(3)/2, 0.0])
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

# Function to calculate the Jacobi constant
def jacobi_constant(x, y, z, mu):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    return x**2 + y**2 + 2 * (1 - mu) / r1 + 2 * mu / r2

# Calculate Jacobi constants for L1, L2, L3, L4, and L5
C_L1 = jacobi_constant(L1[0], L1[1], L1[2], mu)
C_L2 = jacobi_constant(L2[0], L2[1], L2[2], mu)
C_L3 = jacobi_constant(L3[0], L3[1], L3[2], mu)
C_L4 = jacobi_constant(L4[0], L4[1], L4[2], mu)
C_L5 = jacobi_constant(L5[0], L5[1], L5[2], mu)

# Function to compute zero-velocity curves
def zero_velocity_curve(x, y, C, mu):
    z = 0  # Zero-velocity curves in the x-y plane
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    return x**2 + y**2 + 2 * (1 - mu) / r1 + 2 * mu / r2 - C

# Grid for plotting zero-velocity curves
x_vals = np.linspace(-2, 2, 400)
y_vals = np.linspace(-2, 2, 400)
X, Y = np.meshgrid(x_vals, y_vals)

Z_L1 = zero_velocity_curve(X, Y, C_L1, mu)
Z_L2 = zero_velocity_curve(X, Y, C_L2, mu)
Z_L3 = zero_velocity_curve(X, Y, C_L3, mu)

# Plot Neptune, Triton, Lagrange points, and zero-velocity curves in 2D
neptune_position = np.array([0.0, 0.0])
triton_position = np.array([distance, 0.0])

plt.figure(figsize=(10, 8))
plt.scatter(neptune_position[0], neptune_position[1], color='b', label='M1')
plt.scatter(triton_position[0], triton_position[1], color='g', label='M2')
plt.scatter(L1[0], L1[1], color='pink', label='L1')
plt.scatter(L2[0], L2[1], color='orange', label='L2')
plt.scatter(L3[0], L3[1], color='purple', label='L3')
plt.scatter(L4[0], L4[1], color='black', label='L4')
plt.scatter(L5[0], L5[1], color='red', label='L5')

plt.text(neptune_position[0], neptune_position[1], ' M1', fontsize=12, ha='center', va='bottom')
plt.text(triton_position[0], triton_position[1], ' M2', fontsize=12, ha='center', va='bottom')
plt.text(L1[0], L1[1], ' L1', fontsize=12, ha='center', va='top')
plt.text(L2[0], L2[1], ' L2', fontsize=12, ha='center', va='top')
plt.text(L3[0], L3[1], ' L3', fontsize=12, ha='center', va='top')
plt.text(L4[0], L4[1], ' L4', fontsize=12, ha='center', va='top')
plt.text(L5[0], L5[1], ' L5', fontsize=12, ha='center', va='top')
plt.text(0, 0, 'O', fontsize=12, ha='right', va='top')  # Origin (center of mass)

plt.contour(X, Y, Z_L1, levels=[0], colors='r', linestyles='--', linewidths=0.5)
plt.contour(X, Y, Z_L2, levels=[0], colors='orange', linestyles='--', linewidths=0.5)
plt.contour(X, Y, Z_L3, levels=[0], colors='purple', linestyles='--', linewidths=0.5)

plt.xlabel('Distance in AU')
plt.ylabel('Distance in AU')
plt.title('Positions of Lagrange Points and Zero-Velocity Curves with Neptune and Triton')
plt.legend()
plt.grid(True)

plt.show()

# 3D Plot of Jacobi constant surface
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Define the range and resolution for the 3D plot
x_vals = np.linspace(-2, 2, 100)
y_vals = np.linspace(-2, 2, 100)
X, Y = np.meshgrid(x_vals, y_vals)

# Compute the Jacobi constant for the grid points
Z = jacobi_constant(X, Y, 0, mu)

# Plot the surface
ax.plot_surface(X, Y, -Z, cmap='gray', alpha=0.7)

# Plot Neptune and Triton
ax.scatter(0, 0, 0, color='b', s=100, label='M1')
ax.scatter(distance, 0, 0, color='g', s=100, label='M2')

# Plot the Lagrange points as scatter points
ax.scatter(L1[0], L1[1], -C_L1, color='pink', s=100, label='L1')
ax.scatter(L2[0], L2[1], -C_L2, color='orange', s=100, label='L2')
ax.scatter(L3[0], L3[1], -C_L3, color='purple', s=100, label='L3')
ax.scatter(L4[0], L4[1], -C_L4, color='black', s=100, label='L4')
ax.scatter(L5[0], L5[1], -C_L5, color='red', s=100, label='L5')

# Annotate the Lagrange points
ax.text(L1[0], L1[1], -C_L1, ' L1', color='pink', fontsize=12, ha='center', va='bottom')
ax.text(L2[0], L2[1], -C_L2, ' L2', color='orange', fontsize=12, ha='center', va='bottom')
ax.text(L3[0], L3[1], -C_L3, ' L3', color='purple', fontsize=12, ha='center', va='bottom')
ax.text(L4[0], L4[1], -C_L4, ' L4', color='black', fontsize=12, ha='center', va='bottom')
ax.text(L5[0], L5[1], -C_L5, ' L5', color='red', fontsize=12, ha='center', va='bottom')
ax.text(0, 0, 0, ' M1', color='b', fontsize=12, ha='center', va='bottom')
ax.text(distance, 0, 0, ' M2', color='g', fontsize=12, ha='center', va='bottom')

# Set labels and title
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('-Jacobi Constant')
ax.set_title('3D Surface of Jacobi Constant with Lagrange Points')

# Adjust the view angle for better visibility
ax.view_init(elev=30, azim=30)

plt.legend()
plt.show()





