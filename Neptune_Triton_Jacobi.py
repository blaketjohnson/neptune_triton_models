import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
G = 6.672e-11  # Gravitational constant (m^3 kg^-1 s^-2)

# Neptune
m1 = 1.0243e26  # Mass of Neptune (kg)
R1 = 25225e3  # Radius of Neptune (m)

# Triton
m2 = 2.15e22  # Mass of Triton (kg)
a2_m = 354760e3  # Semi-major axis of Triton (m)
e2 = 0.0004  # Eccentricity of Triton

# Satellite
m3 = 721.9  # Mass of Voyager 1 (kg)
a3 = R1 + 20000e3  # Initial distance from Neptune (m) - increased altitude

# Initial conditions for Neptune (stationary at the origin)
x1, y1, z1 = 0.0, 0.0, 0.0
vx1, vy1, vz1 = 0.0, 0.0, 0.0

# Initial conditions for Triton at periapsis
r2_periapsis = a2_m * (1 - e2)  # Distance at periapsis (m)
v2_periapsis = np.sqrt(G * m1 * (1 + e2) / (a2_m * (1 - e2)))  # Velocity at periapsis (m/s)

x2, y2, z2 = r2_periapsis, 0.0, 0.0
vx2, vy2, vz2 = 0.0, v2_periapsis, 0.0

# Initial conditions for the satellite
x3, y3, z3 = a3, 0.0, 0.0
vx3, vy3, vz3 = 0.0, np.sqrt(G * m1 / a3), 0.0

# State vector: [x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2, x3, vx3, y3, vy3, z3, vz3]
state0 = [x1, vx1, y1, vy1, z1, vz1, x2, vx2, y2, vy2, z2, vz2, x3, vx3, y3, vy3, z3, vz3]

def der_state(t, state):
    """Compute the derivative of the given state"""
    # Extract positions and velocities
    r1 = np.array([state[0], state[2], state[4]])   # x1, y1, z1
    r2 = np.array([state[6], state[8], state[10]])  # x2, y2, z2
    r3 = np.array([state[12], state[14], state[16]]) # x3, y3, z3

    v1 = np.array([state[1], state[3], state[5]])   # vx1, vy1, vz1
    v2 = np.array([state[7], state[9], state[11]])  # vx2, vy2, vz2
    v3 = np.array([state[13], state[15], state[17]]) # vx3, vy3, vz3

    # Relative position vectors
    r21 = r2 - r1
    r31 = r3 - r1
    r12 = r1 - r2
    r32 = r3 - r2
    r13 = r1 - r3
    r23 = r2 - r3

    # Magnitudes of relative position vectors
    R21 = np.linalg.norm(r21)
    R31 = np.linalg.norm(r31)
    R12 = np.linalg.norm(r12)
    R32 = np.linalg.norm(r32)
    R13 = np.linalg.norm(r13)
    R23 = np.linalg.norm(r23)

    # Accelerations
    a1 = G * ((m2 * r21 / R21**3) + (m3 * r31 / R31**3))
    a2 = G * ((m1 * r12 / R12**3) + (m3 * r32 / R32**3))
    a3 = G * ((m1 * r13 / R13**3) + (m2 * r23 / R23**3))

    # Derivative of the state vector
    deriv = [
        state[1], a1[0], state[3], a1[1], state[5], a1[2],  # Neptune
        state[7], a2[0], state[9], a2[1], state[11], a2[2], # Triton
        state[13], a3[0], state[15], a3[1], state[17], a3[2] # Satellite
    ]

    return deriv

# Solve ODE for the satellite (10 days)
t_span_satellite = (0, 3600 * 24 * 10)  # Ten days
t_eval_satellite = np.linspace(*t_span_satellite, 500)

sol_satellite = solve_ivp(der_state, t_span_satellite, state0, t_eval=t_eval_satellite, method='RK45', rtol=1e-6, atol=1e-9)

# Extract solution for the satellite
x1_s, vx1_s, y1_s, vy1_s, z1_s, vz1_s, x2_s, vx2_s, y2_s, vy2_s, z2_s, vz2_s, x3_s, vx3_s, y3_s, vy3_s, z3_s, vz3_s = sol_satellite.y

# Convert back to km for plotting
x1_s_km, y1_s_km, z1_s_km = x1_s / 1e3, y1_s / 1e3, z1_s / 1e3
x2_s_km, y2_s_km, z2_s_km = x2_s / 1e3, y2_s / 1e3, z2_s / 1e3
x3_s_km, y3_s_km, z3_s_km = x3_s / 1e3, y3_s / 1e3, z3_s / 1e3

# Compute distances from the satellite to Neptune and Triton
r13_s = np.sqrt(x1_s**2 + y1_s**2 + z1_s**2)
r23_s = np.sqrt(x2_s**2 + y2_s**2 + z2_s**2)

# Compute velocities of the satellite
v3_s = np.sqrt(vx3_s**2 + vy3_s**2 + vz3_s**2)

# Compute the effective potential U
Omega = np.sqrt(G * m1 / a2_m**3)  # Angular velocity of the rotating frame
U = G * m1 / r13_s + G * m2 / r23_s + 0.5 * (Omega * np.sqrt(x3_s**2 + y3_s**2))**2

# Compute the Jacobi Integral
C = 2 * U - v3_s**2

# Plot the Jacobi Integral over time
fig = plt.figure()
plt.plot(t_eval_satellite / 3600, C)
plt.xlabel("Time (hours)")
plt.ylabel("Jacobi Integral (C)")
plt.title("Jacobi Integral over Time for the Satellite")
plt.grid(True)
plt.show()
