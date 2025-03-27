import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from constants import *

# Compute Triton's orbital parameters
mu_2d = G * (m1 + m2)  # Gravitational parameter (m^3/s^2)

rp2_km = a2_km * (1 - e2)  # Perigee radius (km)
rp2_m = rp2_km * 1000  # Convert to meters
b2_km = a2_km * math.sqrt(1 - e2**2)  # Semi-minor axis (km)

h2_p = math.sqrt(mu_2d * a2_m * (1 - e2**2))  # Angular momentum (m^2/s)
v2_p = h2_p / rp2_m  # Velocity at perigee (m/s)
E2_sp = -mu_2d / (2 * a2_m)  # Specific orbital energy (J/kg)
T2_orbit = (2 * math.pi / math.sqrt(mu_2d)) * (a2_m ** (3 / 2))  # Orbital period (s)

# Print computed values
if __name__ == '__main__':
    print(f"Radius at perigee: {rp2_km:.2f} km")
    print(f"Angular momentum: {h2_p:.2e} m^2/s")
    print(f"Velocity at perigee: {v2_p:.2f} m/s")
    print(f"Orbital period: {T2_orbit:.2f} s")

# Define the two-body equations of motion
def dfdt(f, t):
    x, y, z, vx, vy, vz = f
    r = np.sqrt(x**2 + y**2 + z**2)
    return [vx, vy, vz, (-mu_2d / r**3) * x, (-mu_2d / r**3) * y, (-mu_2d / r**3) * z]

# Initial conditions (Triton at perigee)
t0 = 0
r0 = np.array([rp2_km * 1000, 0, 0])  # Position (m)
v0 = np.array([0, v2_p, 0])  # Velocity (m/s)
f0 = [r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]]

# Time span for integration (10 orbits)
tf = T2_orbit * 10
dt = 100  # Time step (s)
tspan = np.arange(t0, tf, dt)

# Integrate equations of motion
sol = odeint(dfdt, f0, tspan)
sol_km = sol[:, :3] / 1000  # Convert to km

# Plot results
if __name__ == '__main__':
    plt.figure()
    plt.plot(sol_km[:, 0], sol_km[:, 1], color='saddlebrown', label="Triton's Orbit")
    plt.scatter(0, 0, color='blue', s=50, label='Neptune')
    plt.scatter(sol_km[0, 0], sol_km[0, 1], color='saddlebrown', s=30, label='Triton at Perigee')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.legend(loc='upper right')
    plt.title("Triton's Orbit around Neptune")
    plt.grid()
    plt.show()





