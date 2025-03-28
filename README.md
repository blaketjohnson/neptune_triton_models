# Neptune-Triton Thesis: Dynamics in the Restricted Three-Body Problem

## Overview
This repository contains my graduate thesis work, which analyzes the orbital dynamics of the Neptune-Triton system using advanced methods in celestial mechanics. The project models both the Restricted Two-Body and Restricted Three-Body Problems (R2BP & RTBP), computes key quantities like the Jacobi constant, identifies the Lagrange points, visualizes zero-velocity curves, and assesses the stability of these equilibrium points.

## Motivation
The primary goal of this research is to understand the complex dynamics of the Neptune-Triton system by:
- Simulating Triton’s orbit around Neptune using the two-body problem.
- Extending the analysis to include the effects of a third body (RTBP), which introduces additional dynamical features.
- Computing the Jacobi constant to define energy boundaries in the rotating frame.
- Locating the Lagrange points where gravitational and centrifugal forces balance.
- Visualizing zero-velocity curves that outline forbidden regions for satellite motion.
- Analyzing the stability of the Lagrange points through eigenvalue analysis.

## Repository Structure

neptune_triton_thesis/
├── src/
│   ├── r2bp_calculations.py       # Two-body simulation: Triton orbiting Neptune
│   ├── r3bp_calculations.py       # Three-body simulation: RTBP in the rotating frame
│   ├── jacobi_constant.py         # Jacobi constant calculator
│   ├── lagrange_points.py         # L1–L5 Lagrange point solver & plotter
│   ├── zero_velocity_curves.py    # Contour plots for Jacobi-based forbidden regions
│   ├── stability_lagrange.py      # Eigenvalue analysis of Lagrange point stability
│   ├── constants.py               # Neptune/Triton system constants
│   ├── equations.py               # Motion equations and integration methods
│   └── U_derivatives.py           # Symbolic partial derivatives of effective potential
├── docs/
│   ├── r2bp_calculations.md
│   ├── r3bp_calculations.md
│   ├── jacobi_constant.md
│   ├── lagrange_points.md
│   ├── zero_velocity_curves.md
│   └── stability_lagrange.md
├── old/                           # Archived test/unused code
├── README.md                      # Master README (this file)
└── LICENSE                        # License file (if applicable)


## Getting Started

### Prerequisites
- **Python 3.8+** (using a virtual environment is recommended)
- The following Python libraries:
  - NumPy
  - SciPy
  - Matplotlib
  - SymPy

To install the dependencies, run:
```bash
pip install numpy scipy matplotlib sympy
