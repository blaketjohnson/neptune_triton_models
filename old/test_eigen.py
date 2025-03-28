import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from scipy.optimize import root
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp
import sys


from part_5 import *

# Function to calculate and print the second partial derivatives Uxx, Uxy, Uyy
def dU():
    x, y, mu = sp.symbols('x y mu')
    r1 = sp.sqrt((x + mu)**2 + y**2)
    r2 = sp.sqrt((x + mu - 1)**2 + y**2)
    U_sym = 0.5 * (x**2 + y**2) + (1 - mu) / r1 + mu / r2
    
    Uxx = sp.diff(U_sym, x, x)
    Uxy = sp.diff(U_sym, x, y)
    Uyy = sp.diff(U_sym, y, y)
    
    return Uxx, Uxy, Uyy

# Function to calculate the second partial derivatives Uxx, Uxy, Uyy with given values
def ddU(x_val, y_val, mu_val):
    Uxx, Uxy, Uyy = dU()
    Uxx_val = Uxx.subs({sp.symbols('x'): x_val, sp.symbols('y'): y_val, sp.symbols('mu'): mu_val})
    Uxy_val = Uxy.subs({sp.symbols('x'): x_val, sp.symbols('y'): y_val, sp.symbols('mu'): mu_val})
    Uyy_val = Uyy.subs({sp.symbols('x'): x_val, sp.symbols('y'): y_val, sp.symbols('mu'): mu_val})
    
    return float(Uxx_val), float(Uxy_val), float(Uyy_val)

# Function to round small numbers to zero based on a threshold
def round_small_numbers_to_zero(value, threshold=1e-10):
    if abs(value) < threshold:
        return 0
    else:
        return value

# Function to round small components of eigenvalues
def round_eigenvalues(eigenvalues, threshold=1e-10):
    rounded_eigenvalues = []
    for eigenvalue in eigenvalues:
        real_part = round_small_numbers_to_zero(sp.re(eigenvalue), threshold)
        imag_part = round_small_numbers_to_zero(sp.im(eigenvalue), threshold)
        rounded_eigenvalues.append(real_part + imag_part * sp.I)
    return rounded_eigenvalues

# Function to calculate the eigenvalues of the system matrix
def eigenvalues(x_val, y_val, mu_val):
    Uxx_val, Uxy_val, Uyy_val = ddU(x_val, y_val, mu_val)
    
    # Print Uxx, Uxy, Uyy values for inspection
    print(f"Uxx: {Uxx_val}, Uxy: {Uxy_val}, Uyy: {Uyy_val}")
    
    # Define the matrix A
    A = sp.Matrix([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [Uxx_val, Uxy_val, 0, 2],
        [Uxy_val, Uyy_val, -2, 0]
    ])
    
    # Print matrix A for inspection
    print("Matrix A:")
    sp.pprint(A, use_unicode=True)
    
    # Calculate the eigenvalues of matrix A
    eigenvalues = A.eigenvals()
    
    return eigenvalues

# Example usage with specific values
x_val = 0.848
y_val = 0
mu_val = 0.01

eigenvalues_result = eigenvalues(x_val, y_val, mu_val)

# Round small components of eigenvalues
rounded_eigenvalues = round_eigenvalues(eigenvalues_result)

# Display the rounded eigenvalues with more precision
print("Rounded Eigenvalues:")
for eigenvalue in rounded_eigenvalues:
    sp.pprint(eigenvalue, use_unicode=True)


