import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math

"""
Part 1: Tritons Orbit
1. Import values from appendix A Solar System Dynamics
2. Using these values to determine Tritons Orbit
3. Use the Two Body Problem to find the Equations of Motion and plot Tritons orbit around Neptune

"""


""" 1.1 - Input Values from Appendix A and convert as necessary"""


# Constants
k = 0.01720209895  # Gaussian Gravitational Constant
c = 2.99792458e8  # Speed of light in m/s
G = 6.672e-11  # Gravitational Constant in m^3 kg^-1 s^-2
AU = 1.495978707e11  # Astronomical unit in meters
day = 23.9345  # Hours (Appendix of Orbital Mechanics)

# Neptune
m1 = 1.0243e26  # Mass in kg
R1 = 2.5225e4  # Radius in km
R1_m = R1 * 10**3  # Radius in meters
a1_au = 30.06896348  # Semi-major axis in AU
e1 = 0.00858587  # Eccentricity
I0_1 = 1.76917  # Inclination in degrees

# Triton
m2 = 2.15e22  # Mass in kg
R2 = 1353  # Radius in km
R2_m = R2 * 10**3  # Radius in meters
a2_km = 354760  # Semi-major axis in km
a2_m = a2_km * 1000  # Semi-major axis in meters
a2_au = a2_m / AU  # Semi-major axis in AU
e2 = 0.0004  # Eccentricity
I0_2 = 156.834  # Inclination in degrees
T2 = 5.876854  # Orbital period in days (positive value)
T2_s = T2 * 24 * 3600  # Orbital period in seconds
print(f"T2_s: {T2_s:.2f} seconds")


""" 1.2 Using Orbital Mechanics for Engineering Students by Curtis - Evaluate Triton's Orbit around Neptune and Perigee Values"""


# Triton's elliptical orbit
mu_2d = G * (m1 + m2)  # Gravitational parameter in m^3/s^2

rp2_km = a2_km * (1 - e2)  # Radius at perigee in km
rp2_m = rp2_km * 1000  # Radius at perigee in meters

b2_km = a2_km * math.sqrt(1 - e2**2)  # Semi-minor axis in km

h2_p = math.sqrt(mu_2d * rp2_m * (1 + e2))  # Angular momentum of Triton in m^2/s

v2_p = h2_p / rp2_m  # Velocity of Triton at perigee in m/s

E2_sp = -mu_2d / (2 * a2_m)  # Specific orbital energy in J/kg

T2_orbit = (2 * math.pi / math.sqrt(mu_2d)) * (a2_m ** (3 / 2))  # Orbital period in seconds

# Output the results
print(f"Radius at perigee: {rp2_km:.2f} km")
print(f"Radius at perigee: {rp2_m:.2f} m")
print(f"Semi-minor axis: {b2_km:.2f} km")
print(f"Angular momentum: {h2_p:.2e} m^2/s")
print(f"Velocity at perigee: {v2_p:.2f} m/s")
print(f"Specific orbital energy: {E2_sp:.2e} J/kg")
print(f"Orbital period: {T2_orbit:.2f} s")


"""  1.3 Calculating the Two Body Problem """


# Define the two-body problem differential equations
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

# Plotting the trajectory
plt.figure()
plt.plot(sol_km[:, 0], sol_km[:, 1], color='saddlebrown', label="Triton's Orbit")
plt.scatter(0, 0, color='blue', s=50, label='Neptune')  # Neptune at the origin
plt.scatter(sol_km[0, 0], sol_km[0, 1], color='saddlebrown', s=30, label='Triton at Perigee')  # Triton at initial position

# Plotting settings
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.legend(loc='upper right')
plt.title("Triton's Orbit around Neptune")
plt.grid()
plt.axis()
plt.show()

'''
Part 2: Placing A Satellite Around Neptune
1.) Determine the Constants needed for the equations of motion
2.) Find the equations of motion and create a definition to solve for them
3.) Plot the orbit of a satellite around Neptune

'''


""" 2.1 - Calculate the adjusted gravitational parameter """

mu = m2/(m1+m2) # RTBP focuses on a different axis as the x,y,z axis are fixed and a new reference frame is used wrt rotation
mu1 = 1 - mu # gravitational parameter for Neptune ?
mu2 = mu # graviational parameter for Triton ?



