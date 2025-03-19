# Interplay between intraspecific suppression and environment in shaping biodiversity
This repository is the official implementation of reserach ``Interplay between intraspecific suppression and environment in shaping biodiversity.''


## Simulation.c
This C code simulates the Generalized Consumer-Resource model (GCRM) with intraspecific suppression.

To run the C code, packages below are needed: 
- GSL (for random number generator)
- openmp (to boost matrix multiplication)


## SelfConsistent.py
This Python code gives the stable solution of the GCRM with intraspecific suppression by minimizing the difference between the given variable values and the self-consistent equation solutions.

To run the Python code, SciPy library is needed.
