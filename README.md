![DynamicalSystems.jl logo: The Double Pendulum](https://i.imgur.com/nFQFdB0.gif)

| **Documentation** | Gitter | Travis | Citing |
|:--------:|:-----:|:-----:|:----:|
|[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaDynamics.github.io/DynamicalSystems.jl/latest) | [![Gitter](https://img.shields.io/gitter/room/nwjs/nw.js.svg)](https://gitter.im/JuliaDynamics/Lobby) | [![Build Status](https://travis-ci.org/JuliaDynamics/DynamicalSystems.jl.svg?branch=master)](https://travis-ci.org/JuliaDynamics/DynamicalSystems.jl) | [![DOI](http://joss.theoj.org/papers/10.21105/joss.00598/status.svg)](https://doi.org/10.21105/joss.00598)

**DynamicalSystems.jl** is a Julia software library for the exploration of chaos and nonlinear dynamics. The current repository holds the documentation and exports *all* relevant packages.

---

### **For installation instructions, full content overview and detailed documentation, [click here](https://juliadynamics.github.io/DynamicalSystems.jl/latest/).**

---

## Brief Content Overview
### [DynamicalSystemsBase.jl](https://juliadynamics.github.io/DynamicalSystems.jl/latest/definition/general/)   
1. Intuitive, consistent APIs for the definition of general dynamical systems, both maps and flows. In fact we have implementations for 8 possible dynamical systems:
    * Continuous or Discrete.
    * In-place or out-of-place (large versus small systems).
    * Auto-differentiated or not (for the Jacobian function).

4. Dedicated interface for numerical data
5. Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user.
4. Robust implementations of all kinds of integrators, that evolve the system,
   many states of the system, or even deviation vectors. See the advanced documentation for this.
6. Library of predefined well-known dynamical systems that have been used extensively in scientific research.


### [ChaosTools.jl](https://juliadynamics.github.io/DynamicalSystems.jl/latest/chaos/overview/)

* Lyapunov Exponents
* Poincare SOS, Orbit Diagrams
* Entropies and Dimensions
* Delay Coordinates Embedding
* Neighborhood estimation
* Lyapunov exponent of a timeseries
* Finding Fixed Points of Maps
* Detecting Chaos

## Double pendulum video:
Checkout this *amazing* video by Cormullion, featuring the double pendulum (logo of our library)

[![Double Pendulums Video](http://img.youtube.com/vi/vLDpLxU2fEg/0.jpg)](
https://www.youtube.com/watch?v=vLDpLxU2fEg)
