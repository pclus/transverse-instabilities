# Transverse-instabilities

These codes require [GNU Scientific Library](https://www.gnu.org/software/gsl/doc/html/index.html) (tested in version 2.7).

The `src` folder contains the source code:

- `netdynamics`: contains the main functions for monitoring the integration of the large-scale brain model. Can work with either `system` or `system_lyap_exp`.
- `system`: contains the main functions for integration of the large-scale brain model, including the Runge-Kutta algorithm and the velocity fields.
- `system_lyap_exp`: same as `system`, but including computation of Lyapunov exponents.
- `floquet`: independent code to compute the Floquet exponents using the Master Stability Function formalism. 

