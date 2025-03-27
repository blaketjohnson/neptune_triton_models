import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp

import sys

# Import constants from your specific path
sys.path.append('/Users/blakejohnson/Documents/Thesis/Three Body Problem')
from constants import *

"""
Part 1
"""

# (Curtis p. 63 2.21)
mu_2d = G * (m1 + m2)  # Gravitational parameter in m^3/s^2

# curtis p. 82 (2.73)
rp2_km = a2_km * (1 - e2)  # Radius at perigee in km
rp2_m = rp2_km * 1000  # Radius at perigee in meters

# curtis p. 82 (2.76)
b2_km = a2_km * math.sqrt(1 - e2**2)  # Semi-minor axis in km


#h2_p = math.sqrt(mu_2d * rp2_m * (1 + e2))  # Angular momentum of Triton in m^2/s

# curtis p. 81 (2.71)
h2_p = math.sqrt(mu_2d * a2_m * (1 - e2**2))  # Angular momentum of Triton in m^2/s


# curtis p 69 (2.31)
v2_p = h2_p / rp2_m  # Velocity of Triton at perigee in m/s

# curtis p.83 (2.80)
E2_sp = -mu_2d / (2 * a2_m)  # Specific orbital energy in J/kg

#curtis p.84 (2.83)
T2_orbit = (2 * math.pi / math.sqrt(mu_2d)) * (a2_m ** (3 / 2))  # Orbital period in seconds

def dfdt(f, t):
    x, y, z, vx, vy, vz = f
    r = np.sqrt(x**2 + y**2 + z**2)
    dx = vx
    dy = vy
    dz = vz
    ddx = (-mu_2d / r**3) * x
    ddy = (-mu_2d / r**3) * y
    ddz = (-mu_2d / r**3) * z
    return [dx, dy, dz, ddx, ddy, ddz]

# Initial conditions
t0 = 0  # Set initial time to 0
r0 = np.array([rp2_km * 1000, 0, 0])  # Initial position in meters
v0 = np.array([0, v2_p, 0])  # Initial velocity in m/s

# Initial conditions for odeint: [x, y, z, vx, vy, vz]
f0 = [r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]]

# Time span for the integration
tf = T2_s * 10  # Integrate for multiple orbits
dt = 100  # Time step for integration
tspan = np.arange(t0, tf, dt)

# Integrate the equations of motion
sol = odeint(dfdt, f0, tspan)

# Convert position from meters to kilometers for plotting
sol_km = sol[:, :3] / 1000

"""
Part 2
"""

# Define the effective potential and its gradients
def U(x, y, z, mu):
    # murray p. 66 (3.8)(3.9)
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)

    # return the potential energy U
    # murray p.67 (3.22)
    return 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2 # NOTE: The (x**2 + y**2 + z**2) is the centrifugal potential, 1/r1 and 1/r2 are the gravitational potential.


# Define the gradients of the effective potential
def gradient_U(x, y, z, mu):
    r1 = (x + mu)**2 + y**2 + z**2
    r2 = (x + mu - 1)**2 + y**2 + z**2
    
    Ux = x - (1 - mu) * (x + mu) / r1**1.5 - mu * (x + mu - 1) / r2**1.5
    Uy = y - (1 - mu) * y / r1**1.5 - mu * y / r2**1.5
    Uz = -(1 - mu) * z / r1**1.5 - mu * z / r2**1.5
    
    return np.array([Ux, Uy, Uz])


# Define the equations of motion
def equations(t, state, mu):
    x, y, z, vx, vy, vz = state
    
    # Calculate partial derivatives
    Ux, Uy, Uz = gradient_U(x, y, z, mu)
    
    # murray p. 67 (3.23)(3.24)(3.25)
    ax = Ux + 2 * vy
    ay = Uy - 2 * vx
    az = Uz
    
    return [vx, vy, vz, ax, ay, az]

# Initial conditions
x0 = 0.5
y0 = 0
z0 = 0
vx0 = 0
vy0 = 1
vz0 = 0
state0 = [x0, y0, z0, vx0, vy0, vz0]

# Time span for the solution I set 100 with 1000 intervals
t_span = (0, 100)
t_eval = np.linspace(*t_span, 1000)

# Solve the system of equations
solution = solve_ivp(equations, t_span, state0, t_eval=t_eval, method='RK45', args=(mu,))

# Extract the results
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]


"""
Part 3
"""

def jacobi_constant(state, mu):
    x, y, z, vx, vy, vz = state
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    U = 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2
    C = 2*U - (vx**2 + vy**2 + vz**2)
    return C

# Calculate Jacobi constant for initial conditions
C = jacobi_constant(state0, mu)