# lagrange_points_04.py

## Overview
This script calculates the **Lagrange Points (L1 - L5)** in the **Restricted Three-Body Problem (RTBP)** for the Neptune-Triton system. It determines the equilibrium points where the gravitational and centrifugal forces balance.

## Functionality
- Uses numerical root-finding (`scipy.optimize.root`) to locate the Lagrange points.
- Imports the potential gradient function from `jacobi_contant_03.py`.
- Computes the five Lagrange points using initial guesses near their expected locations.
- Plots Neptune, Triton, and the Lagrange points in a 2D figure.

## Key Equations and References
- **Gradient of the Effective Potential:**
  \[
  \nabla U = 0
  \]
  The Lagrange points are solutions to this equation, where the gradient of the effective potential is zero.
- **Equilibrium Conditions:**
  - L1, L2, L3 are collinear points found along the Neptune-Triton axis.
  - L4 and L5 form equilateral triangles with Neptune and Triton.

Equations referenced from *Orbital Mechanics for Engineering Students* by Curtis.

## How to Run the Script
Ensure `jacobi_contant_03.py` is available for import.
```sh
python lagrange_points_04.py
```

## Dependencies
- Python 3.x
- NumPy
- SciPy
- Matplotlib

## Output
- Console output of computed Lagrange point positions.
- A 2D plot of Neptune, Triton, and the five Lagrange points.

## Notes
- The script assumes a non-dimensional rotating frame.
- L4 and L5 are expected to be stable equilibrium points.
- The script currently does not analyze stability but can be extended to do so.

## Next Steps
- Extend analysis to verify stability using eigenvalues of the Jacobian matrix.
- Compute and visualize zero-velocity curves near the Lagrange points.

