# DynamicalSystems.jl

A Julia package for the exploration of continuous and discrete dynamical systems.

| **Documentation**   | [**Package Evaluator**](http://pkg.julialang.org/?pkg=DynamicalSystems#DynamicalSystems) | **Travis**     | **AppVeyor** |
|:--------:|:-------------------:|:-----------------------:|:-----:|
| Not yet! (use this README.md) | Not yet! | [![Build Status](https://travis-ci.org/Datseris/DynamicalSystems.jl.svg?branch=master)](https://travis-ci.org/Datseris/DynamicalSystems.jl) | [![Build status](https://ci.appveyor.com/api/projects/status/oabd7hgibx63bo1l?svg=true)](https://ci.appveyor.com/project/Datseris/dynamicalsystems-jl)


`DynamicalSystems.jl` aims to be a useful companion for students and scientists treading
on the field of Chaos, nonlinear dynamics and dynamical systems in general. The package
treats discrete and continuous systems of the forms:
```math
\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}) \quad \text{or}\quad \vec{x}_{n+1} = \vec{f}(\vec{x}_n)
```
and it **does not** accept non-autonomous systems, since these can be made into autonomous ones with 1 more dynamical variable.

This is the (non-final) list of what this package aims to offer:

1. Intuitive, consistent UI for the definition of dynamical systems.
2. Automatic "completion" of the dynamics of the system with numerically computed
  Jacobians, in case they are not provided by the user.
3. Lyapunov exponent estimation.
4. Entropy estimation.
5. Attractor dimension estimation.
6. Entropy/Attractor dimension/Lyapunov exponents for *numerical data*.
7. Chaos control.
8. Other stuff I have not yet decided upon, since this is like a pre-alpha version.
8. Suggest or Contribute more stuff!

## Lyapunov exponents of discrete systems
The discrete systems assume the following UI:
```julia
DiscreteDS(state, eom, jacob)
```
The equations of motion (`eom`) function takes the current `state` and returns an `SVector` (from the module `StaticArrays`) with the next state. Similarly, the Jacobian function (`jacob`) returns an `SMatrix` given the current state. The `state` is always stored internally as an `SVector` as well.

*Why `SVector` instead of the standard Julia ones? Because they are absurdly
fast and offer just as absurdly fast computation of automatic Jacobians!*

For example, let's calculate the lyapunov exponents of the Hénon system:
```julia
using DynamicalSystems, StaticArrays
henon_eom(x) = @SVector [1.0 - 1.4x[1]^2 + x[2], 0.3*x[1]]
### let's say I am too bored to write down a jacobian
henon = DiscreteDS([0.0, 0.0], henon_eom)
```
These lines created a `DiscreteDS` that represents the Hénon map. Since we did not pass
a Jacobian function to the constructor as the third argument, one was created numerically using the package `ForwardDiff`. Alternatively, you could use the already existing definition in this package: `henon = DynamicalSystems.Systems.henon()`.

The spectrum of the lyapunov exponents is given by the
function:
```julia
λspectrum(ds::DiscreteDS, N::Int; Ntr::Int = 100)
```
and are computed using the QR-decomposition method for `N` steps.
```julia
henon_λ = λspectrum(henon, 1000000)
# should give you: [0.42..., -1.62...]
```
One is positive and the other is negative, which is necessary for a strange attractor, but you can also test:
```julia
sum(henon_ls) ≈ log(0.3)
```
which is also true (and should be for mappings where the determinant of the Jacobian is constant). If you only needed the maximum lyapunov exponent of the system,
the method function `λmax` is **much** more efficient:
```julia
λ = λmax(henon, 1000000)
# gives 0.42018717281774165 or so
```
