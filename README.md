# DynamicalSystems.jl

A Julia package for the exploration of continuous and discrete dynamical systems.

| **Documentation**   | [**Package Evaluator**](http://pkg.julialang.org/?pkg=DynamicalSystems#DynamicalSystems) | **Travis**     | **AppVeyor** |
|:--------:|:-------------------:|:-----------------------:|:-----:|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://datseris.github.io/DynamicalSystems.jl/latest) | Not yet! | [![Build Status](https://travis-ci.org/Datseris/DynamicalSystems.jl.svg?branch=master)](https://travis-ci.org/Datseris/DynamicalSystems.jl) | [![Build status](https://ci.appveyor.com/api/projects/status/oabd7hgibx63bo1l?svg=true)](https://ci.appveyor.com/project/Datseris/dynamicalsystems-jl)

*WARNING: Only trust what is written in the documentation! Anything else might be incomplete, broken, or give wrong results!*

`DynamicalSystems.jl` aims to be a useful companion for students and scientists treading
on the field of Chaos, nonlinear dynamics and dynamical systems in general.

This is the (non-final) list of what this package aims to offer:

1. Intuitive, consistent APIs for the definition of dynamical systems.
2. Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user.
3. Lyapunov exponent estimation.
4. Entropy estimation.
5. Attractor dimension estimation.
6. Entropy/Attractor dimension/Lyapunov exponents for *numerical data*.
6. Attractor reconstruction, embedding and all that jazz.
7. Numeric Computation of Kolmogorov-Sinai entropy.
8. Definition of chaos, by Ott.
7. Chaos control, TBA.
8. Other stuff I have not yet decided upon, since this is like a pre-alpha version.
8. Suggest or Contribute more stuff!
