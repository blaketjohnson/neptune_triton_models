# stability_lagrange_06.py

## Overview
This script analyzes the **stability of the Lagrange points (L1 - L5)** in the **Restricted Three-Body Problem (RTBP)** for the Neptune-Triton system. It calculates the eigenvalues of the system matrix and determines whether each Lagrange point is stable or unstable.

## Functionality
- Computes the **second partial derivatives** of the effective potential.
- Constructs the **linearized equations of motion** in matrix form.
- Computes **eigenvalues** of the Jacobian matrix to assess stability.
- Determines whether each Lagrange point is **stable or unstable**.

## What is Stability Analysis of Lagrange Points?
- Stability is determined by analyzing **small perturbations** around the equilibrium points.
- If all eigenvalues have **non-positive real parts**, the point is **stable**.
- If any eigenvalue has a **positive real part**, the point is **unstable**.

## How is Stability Calculated?
1. Compute the **second derivatives** of the effective potential:
   \[
   U_{xx}, U_{xy}, U_{yy}
   \]
2. Construct the **system matrix (A):**
   \[
   A = \begin{bmatrix} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ U_{xx} & U_{xy} & 0 & 2 \\ U_{xy} & U_{yy} & -2 & 0 \end{bmatrix}
   \]
3. Compute **eigenvalues** of matrix A.
4. Analyze **real parts** of eigenvalues:
   - If all real parts \( \leq 0 \) → **Stable**
   - If any real part \( > 0 \) → **Unstable**

## How to Run the Script
Ensure `zero_velocity_curves_05.py` is available for import.
```sh
python stability_lagrange_06.py
```

## Dependencies
- Python 3.x
- NumPy
- SciPy
- SymPy
- Matplotlib

## Output
- Console output of **eigenvalues** for each Lagrange point.
- Stability classification for **L1, L2, L3, L4, and L5**.

## Notes
- The script assumes a **rotating reference frame**.
- Stability is determined by analyzing **linearized perturbations**.
- L4 and L5 are expected to be stable for **sufficiently small mass ratios**.

## Next Steps
- Extend analysis to visualize **eigenvectors** to understand stable/unstable directions.
- Compute **stability regions** for varying mass ratios.
- Compare numerical results with theoretical predictions.

