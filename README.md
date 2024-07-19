# Enhancing biodiversity through intraspecific suppression in large ecosystems
This repository is the official implementation of reserach ``[Enhancing biodiversity through intraspecific suppression in large ecosystems](https://arxiv.org/abs/2305.12341).''


## Simulation.c
This C code simulates the Generalized Consumer-Resource model (GCRM) with intraspecific suppression.

To run the C code, packages below are needed: 
- gsl (for random number generator)
- openmp (to boost matrix multiplication)


## SelfConsistent.py
This python code obtains the stable solution of the GCRM with intraspecific suppression by minimizing the difference between the given variable values and the self-consistent equation solutions.

To run the python code, SciPy library is needed.
