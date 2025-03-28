import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

m1 = 1.024e26  # mass of Neptune
m2 = 2.14e22   # mass of Triton

# Define constants
mu = m2 / (m1 + m2)  # Mass ratio
mu1 = 1 - mu
mu2 = mu  

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
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z)
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_zlabel('z (km)')
ax.set_title('Trajectory in the Restricted Three Body Problem')

# Add markers for Neptune and Triton with specified colors
ax.scatter([-mu, 1-mu], [0, 0], [0, 0], color=['b', 'g'], label=['Neptune', 'Triton'])
ax.legend()

plt.show()


