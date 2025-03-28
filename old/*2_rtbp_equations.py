import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp
import sys

# Import constants from your specific path
sys.path.append('/Users/blakejohnson/Documents/Thesis/Three Body Problem')
from constants import *
from equations import *
from r2bp_calculations_01 import *




"""
Part 2: The Restricted Three Body Problem (RTBP)

Once the two body equations of motion have been solved, we can move to the restricted three body problem

"""


"""
2.1 Constants and Known Values:

I have imported the constants from part 1 above. However, in the three body equation we have a different mu.
Lets calculate mu for the three body problem
"""

# Define constants

# Murray p.64 (3.1)(3.2)
mu = m2 / (m1 + m2)  # Mass ratio
mu1 = 1 - mu
mu2 = mu  


"""
2.2 Finding the Equation for the Gradient

The RTBP utilizes potential energy and when we look at the Lagrange points we will need to determine the gradient of potential energy.
The lagrange points are at a position of equilibrium, meaning the gradient will be zero.
So before I get into the functions, I want to find an equation for the gradient of potential energy for x, y, and z.
"""

# Define symbols
# I want to find an equation and not a numerical value so I am stating that I want x,y,z to be char not variables.
x, y, z, u = sp.symbols('x y z u')

# Define distances r1 and r2
# murray p. 66 (3.8)(3.9)
r1 = sp.sqrt((x + u)**2 + y**2 + z**2)
r2 = sp.sqrt((x + u - 1)**2 + y**2 + z**2)

# Define the effective potential function U
# murray p.67 (3.22)
U = (1/2) * (x**2 + y**2 + z**2) + ((1 - u) / r1) + (u / r2)

# Calculate the partial derivatives
#This will return an equation for the gradient, not a numerical value
U_x = sp.diff(U, x)
U_y = sp.diff(U, y)
U_z = sp.diff(U, z)

# Display the results
print("Effective Potential U:")
sp.pretty_print(U)
print("\nPartial derivatives:")
print("U_x:")
sp.pretty_print(U_x)
print("\nU_y:")
sp.pretty_print(U_y)
print("\nU_z:")
sp.pretty_print(U_z)

# Convert symbolic expressions to Python functions for numerical evaluation
U_x_func = sp.lambdify((x, y, z, u), U_x, 'numpy')
U_y_func = sp.lambdify((x, y, z, u), U_y, 'numpy')
U_z_func = sp.lambdify((x, y, z, u), U_z, 'numpy')

"""
2.3 Calculating the Equations of Motion

Now that we have the necessary equations from murray and dervied above we can put them together to calculate the equations of motion.
First I calculated Potential Energy, then the gradient of potential energy, then the acceleration (equations of motion)
"""

# Potential Energy

def U(x, y, z, mu):
    # murray p. 66 (3.8)(3.9)
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)

    # return the potential energy U
    # murray p.67 (3.22)
    return 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2 # NOTE: The (x**2 + y**2 + z**2) is the centrifugal potential, 1/r1 and 1/r2 are the gravitational potential.

# Define the gradients of the effective potential numerically
def gradient_U(x, y, z, mu):
    Ux = U_x_func(x, y, z, mu)
    Uy = U_y_func(x, y, z, mu)
    Uz = U_z_func(x, y, z, mu)
    return np.array([Ux, Uy, Uz])

"""
Now that the potential energy has been calculated for each component, I can solve for the acceleration
I set up the function equations. This time I used 'state' instead of 'f' to differentiate between functions (I dont think this was 
necessary, but I did that to quickly identify any errors in my code that could come up). 

I then set the function to take the coordinate values, apply it to the potential energy equation calculated above to provide a potential 
energy at a given position.

I then set up the acceleration equations. The potential energy is used with the velocity at that position to give acceleration in the x,y,z
coordinates

The function will then return the acceleration and velocity in each coordinate
"""

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

"""
Set inital conditions when the satellite is at perigee. Here I placed the satellite 1/2 between Triton and Neptune. I do not have a 
velocity so I used a value 1 in the vy to represent 100% of the velocity will be in the y direction.
"""

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

"""
Next I used solve_ivp to solve the initial value problem using the Runge-kutta Method (RK45) to get the motion of the satellite
"""
# Solve the system of equations
solution = solve_ivp(equations, t_span, state0, t_eval=t_eval, method='RK45', args=(mu,))

# Extract the results
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]

"""
2.3 Plot The RTBP

Once the eqations of motion have been calculated it is time to plot.
I decided to make this plot in 3D just because.
ax.plot shows that we plot the satellite motion along the x, y, z axis.
the add_subplot makes the plot 3D
The scatter sets the positions of Neptune and Triton

"""
# Plot the trajectory
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, color='green', label='Satellite Trajectory')
ax.set_xlabel('x (non-dimensional units)')
ax.set_ylabel('y (non-dimensional units)')
ax.set_zlabel('z (non-dimensional units)')
ax.set_title('Trajectory in the Restricted Three Body Problem')

# Add markers for Neptune and Triton with specified colors and labels
ax.scatter([-mu], [0], [0], color='blue', s=100, label='Neptune')
ax.scatter([1 - mu], [0], [0], color='saddlebrown', s=50, label='Triton')

# Add legend with custom labels
ax.legend()

plt.show()