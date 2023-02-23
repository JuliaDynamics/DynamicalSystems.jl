# [Overarching tutorial for DynamicalSystems.jl](@id tutorial)

This page serves as a short, but to-the-point, introduction in the **DynamicalSystems.jl** library. It outlines the core components, and how they establish an interface that is used by the rest of the library. It also provides a couple of usage examples to connect the various packages of the library together.

## Installation

To install **DynamicalSystems.jl**, simply do:
```julia
using Pkg; Pkg.add("DynamicalSystems")
```

As discussed in the [contents](@ref contents) page, this installs several packages for the Julia language, that are all exported under a common name. To use them, simply do:
```julia
using DynamicalSystems
```
in your Julia session.

## Core components

The individual packages that compose `DynamicalSystems` interact flawlessly with each other because of the following two components:

1. The [`DynamicalSystem`](@ref) represents a dynamical system with known dynamic rule ``f``. The system can be in discrete time (often called a map), ``\vec{u}_{n+1} = \vec{f}(\vec{u}_n, p, n)``, or in continuous time (often called an ordinary differential equation) ``\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}, p, t)``. In both cases ``u`` is the _state_ of the dynamical system and ``p`` a parameter container. You should have a look at the page [Dynamical System Definition](@ref) for how to create this object. A list of several pre-defined systems exists in the [Predefined Dynamical Systems](@ref) page.
2. Numerical data, that can represent measured experiments, sampled trajectories of dynamical systems, or just sets in the state space, are represented by [`Dataset`](@ref), which is a container of equally-sized data points. Timeseries in **DynamicalSystems.jl** are represented by the already existing `Vector` type of the Julia language.

These core structures `DynamicalSystem, Dataset` are used throughout the package to do useful calculations often used in the field of nonlinear dynamics and chaos.
For example, using [`lyapunovspectrum`](@ref) and [`DynamicalSystem`](@ref) gives you the Lyapunov exponents of a dynamical system with known equations of motion.
Alternatively, by using [`lyapunov_from_data`](@ref) and [`Dataset`](@ref) you can approximate the maximum Lyapunov exponent of a measured trajectory or a reconstructed set resulting from [`embed`](@ref).

All things possible in **DynamicalSystems.jl** are listed in the [Contents](@ref) page.