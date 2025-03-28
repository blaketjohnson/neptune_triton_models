import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import root
from r2bp_calculations_01 import *
from jacobi_contant_03 import *
from constants import *

# Import equations from parts 1 and 2
sys.path.append('/Users/blakejohnson/Documents/rtbp_neptune_triton/neptune_triton_models')


'''
This section is used to determine the input values to run a pointcare map. 
The input values are read from a file and used to determine the range of values 
for the map.

Most of the necessary inputs are in the constants.py file. However I will need to
calulate the mass parameter for the Neptune-Triton system and the distance between the
two bodies.
'''

print(f"r_max (SOI): {SOI_neptune} km")
print(f"r1_min (Neptune): {R1} km")
print(f"r2_min (Triton): {R2} km")
mu_Neptune_triton = m2 / (m1 + m2)
print(f"mu_Neptune_triton: {mu_Neptune_triton}")
print(f"The Jacobi constant is: {C}")
distance_neptune_triton = R1+R2+rp2_km
print(f"Distance between Neptune and Triton: {distance_neptune_triton:.2f} km")



