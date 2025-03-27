

""" 
Murray has several constants listed in Appendix A. He also has several known values for Neptune and Triton in Appendix A.
I entered these values below as given. Additionally I went ahead and converted the units to SI units and common units such as km.

NOTE: I might not use all of these values, but I am listing all of them now in case I end up needing them later in the project.
"""

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

