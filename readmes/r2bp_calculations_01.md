# 01_r2bp_calculations.py

## Overview
This script calculates the motion of Triton around Neptune using the Restricted Two-Body Problem (R2BP). It computes key orbital parameters and integrates the equations of motion to generate a plot of Triton's orbit.

## Functionality
- Calculates Triton's semi-major axis, semi-minor axis, angular momentum, velocity at perigee, specific orbital energy, and orbital period using standard orbital mechanics equations.
- Defines the two-body equations of motion and integrates them over time using `odeint`.
- Plots Triton's orbit around Neptune, with Neptune at the origin and Triton starting at perigee.

## Key Equations and References
- Gravitational parameter: \( \mu = G (m_1 + m_2) \)
- Radius at perigee: \( r_p = a (1 - e) \)
- Semi-minor axis: \( b = a \sqrt{1 - e^2} \)
- Angular momentum: \( h = \sqrt{\mu a (1 - e^2)} \)
- Velocity at perigee: \( v_p = \frac{h}{r_p} \)
- Specific orbital energy: \( E = -\frac{\mu}{2a} \)
- Orbital period: \( T = \frac{2\pi}{\sqrt{\mu}} a^{3/2} \)

Equations referenced from *Orbital Mechanics for Engineering Students* by Curtis.

## How to Run the Script
Ensure `constants.py` is available for import.
```sh
python 01_r2bp_calculations.py
```

## Dependencies
- Python 3.x
- NumPy
- SciPy
- Matplotlib

## Output
- Numerical calculations printed to the console
- A 2D plot of Triton's orbit

## Notes
- Initial conditions assume Triton starts at perigee.
- The integration runs for multiple orbits to produce a complete trajectory.
- Neptune is fixed at the origin, and the system is modeled in kilometers.

## Next Steps
- Clean up inline comments in the Python script for clarity.
- Extend to the Restricted Three-Body Problem (R3BP) in future scripts.

