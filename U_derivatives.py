import sympy as sp

# Define the variables
x, y, z, u = sp.symbols('x y z u')

# Define the distances r1 and r2
r1 = sp.sqrt((x + u)**2 + y**2 + z**2)
r2 = sp.sqrt((x + u - 1)**2 + y**2 + z**2)

# Define the effective potential function Omega
Omega = (1/2) * (x**2 + y**2 + z**2) + ((1 - u) / r1) + (u / r2)

# Calculate the partial derivatives
Omega_x = sp.diff(Omega, x)
Omega_y = sp.diff(Omega, y)
Omega_z = sp.diff(Omega, z)

# Display the results
print("Omega_x:")
sp.pretty_print(Omega_x)
print("\nOmega_y:")
sp.pretty_print(Omega_y)
print("\nOmega_z:")
sp.pretty_print(Omega_z)
