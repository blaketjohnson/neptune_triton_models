# r3bp_calculations_02.py

## Overview
This script models the Restricted Three-Body Problem (RTBP) for the Neptune-Triton system. It calculates the effective potential, derives equations of motion, and simulates a satellite trajectory in the RTBP framework.

## Functionality
- Computes the mass ratio \( \mu \) and related constants.
- Defines and calculates the effective potential \( U \).
- Computes the gradients of \( U \) and derives equations of motion.
- Numerically integrates the equations of motion using `solve_ivp`.
- Plots the satellite's trajectory in the RTBP, marking Neptune and Triton.

## Key Equations and References
- Mass ratio: \( \mu = \frac{m_2}{m_1 + m_2} \)
- Effective potential:
  \[
  U = \frac{1}{2} (x^2 + y^2 + z^2) + \frac{1 - \mu}{r_1} + \frac{\mu}{r_2}
  \]
  where \( r_1 \) and \( r_2 \) are the distances from Neptune and Triton, respectively.
- Equations of motion:
  \[
  \begin{aligned}
  \dot{x} &= v_x, \quad \dot{y} = v_y, \quad \dot{z} = v_z, \\
  \dot{v_x} &= U_x + 2v_y, \quad \dot{v_y} = U_y - 2v_x, \quad \dot{v_z} = U_z
  \end{aligned}
  \]

Equations referenced from *Orbital Mechanics for Engineering Students* by Curtis.

## How to Run the Script
Ensure `r2bp_calculations_01.py` is available for import.
```sh
python r3bp_calculations_02.py
```

## Dependencies
- Python 3.x
- NumPy
- SciPy
- Matplotlib

## Output
- Console output of the effective potential and its gradients.
- A 3D plot of the satellite's trajectory in the RTBP.

## Notes
- The initial conditions assume a satellite at (0.5, 0, 0) with an initial velocity of (0,1,0).
- The simulation runs for a non-dimensional time span of 100.
- Neptune and Triton are fixed at their respective locations in the rotating frame.

## Next Steps
- Refine initial conditions to simulate specific orbits around LP4.
- Extend analysis to include stability conditions near the Lagrange points.

