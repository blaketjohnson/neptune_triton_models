# jacobi_contant_03.py

## Overview
This script calculates the **Jacobi Constant** in the **Restricted Three-Body Problem (RTBP)** for the Neptune-Triton system. The Jacobi Constant is an important integral of motion that describes the total energy in the rotating frame.

## Functionality
- Computes the **Jacobi Constant** using the effective potential and velocity components.
- Imports initial conditions from the previous RTBP calculations (`r3bp_calculations_02.py`).
- Outputs the computed **Jacobi Constant** for given initial conditions.

## Key Equations and References
- **Effective Potential:**
  \[
  U = \frac{1}{2} (x^2 + y^2 + z^2) + \frac{1 - \mu}{r_1} + \frac{\mu}{r_2}
  \]
  where \( r_1 \) and \( r_2 \) are the distances to Neptune and Triton, respectively.
- **Jacobi Constant Definition:**
  \[
  C = 2U - (v_x^2 + v_y^2 + v_z^2)
  \]

Equations referenced from *Orbital Mechanics for Engineering Students* by Curtis.

## How to Run the Script
Ensure `r3bp_calculations_02.py` is available for import.
```sh
python jacobi_contant_03.py
```

## Dependencies
- Python 3.x
- NumPy
- Matplotlib

## Output
- Console output displaying the **Jacobi Constant** for the initial conditions.

## Notes
- The script assumes initial conditions imported from `r3bp_calculations_02.py`.
- The Jacobi Constant helps define zero-velocity curves in the RTBP.

## Next Steps
- Extend the analysis to compute Jacobi constants at Lagrange points.
- Visualize zero-velocity curves based on computed values.

