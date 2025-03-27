# zero_velocity_curves_05.py

## Overview
This script calculates and visualizes **Zero-Velocity Curves** in the **Restricted Three-Body Problem (RTBP)** for the Neptune-Triton system. These curves define regions where a satellite cannot travel due to energy constraints, helping identify possible stable orbits and forbidden regions.

## Functionality
- Computes **Jacobi Constants** at the collinear Lagrange points (L1, L2, L3).
- Uses contour plotting to visualize the zero-velocity curves in the rotating frame.
- Marks Neptune, Triton, and the five Lagrange points in a 2D representation.

## What is a Zero-Velocity Curve?
- A **Zero-Velocity Curve** represents the boundary where the velocity of a satellite would be **zero** in the rotating reference frame.
- These curves arise from the **Jacobi Constant**, which is conserved in the RTBP.
- They determine **forbidden regions** where a spacecraft **cannot travel** for a given energy level.

## How is it Calculated?
- The Zero-Velocity Curve is derived from the **Jacobi Integral**:
  \[
  x^2 + y^2 + 2 \left( \frac{1 - \mu}{r_1} + \frac{\mu}{r_2} \right) = C
  \]
  where \( C \) is the **Jacobi Constant**, and \( r_1, r_2 \) are distances from Neptune and Triton.
- Given a specific \( C \), solving this equation for \( x, y \) provides the curve boundary.
- The script calculates contours of this equation and plots them.

## How to Run the Script
Ensure `lagrange_points_04.py` is available for import.
```sh
python zero_velocity_curves_05.py
```

## Dependencies
- Python 3.x
- NumPy
- SciPy
- Matplotlib

## Output
- Console output of **Jacobi Constants** at L1, L2, and L3.
- A **2D plot** showing:
  - Neptune, Triton, and the five Lagrange points.
  - **Zero-velocity curves** corresponding to the computed Jacobi Constants.

## Notes
- The script assumes motion in the **x-y plane** (z = 0).
- Zero-velocity curves enclose **forbidden regions** where no motion is possible.
- The plotted contours correspond to **L1, L2, and L3 Jacobi Constants**.

## Next Steps
- Extend analysis to include stability and escape trajectories.
- Compute **zero-velocity surfaces** in **3D** for better visualization.
- Compare results with actual satellite orbital paths around LP4.

