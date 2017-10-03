![DynamicalSystems.jl logo: The Double Pendulum](https://i.imgur.com/nFQFdB0.gif)

# Introduction
`DynamicalSystems.jl` is a Julia package for the exploration of continuous and discrete dynamical systems. It aims to be a useful and powerful companion for students and scientists treading
on the field of Chaos, nonlinear dynamics and dynamical systems in general. The package
treats the following types of dynamical systems:

* [Discrete Maps](system_definition/#discrete-systems)
* [Continuous Flows](system_definition/#continuous-systems)
* [Numerical Data](system_definition/#numerical-data)

One of a major goals of this package is to be completely transparent as to what is
going on "under the hood". In scientific research, you never want to use *black boxes*,
e.g. functions that give a result without telling you how it was calculated. `DynamicalSystems.jl` battles this in 2 ways: Firstly, it is written entirely in Julia,
making the source code clear and easy to understand for even novice users. Secondly,
almost every documentation string gives
**direct references to the original papers** where the algorithm is taken from, in case some users don't understand (or simply don't want to read) the source code. For example,
the documentation string of [`?lyapunovs`](https://datseris.github.io/DynamicalSystems.jl/latest/lyapunovs/#DynamicalSystems.lyapunovs) will cite relevant publications for the definition and computation of the lyapunov spectrum.

Be sure to visit the [Contributor Guide](contributors_guide) page, because you can
help make this package better without having to write a single line of code!
Also, if you find this package helpful please consider staring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl)! This would give us an
accurate lower bound of users that this package has already helped!

## Contents
This is the (non-final) list of what this package aims to offer:

1. [Intuitive, consistent APIs for the definition of general dynamical systems](system_definition).
2. [Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user](system_definition).
3. [Lyapunov exponent estimation](lyapunovs).
4. [Entropy and Attractor Dimension estimation](entropies).
6. [Attractor reconstruction, embedding and all that jazz](nlts).
6. [Entropy/Attractor dimension/Lyapunov exponents for *numerical data*](nlts/#numerical-lyapunov-estimation).
* [Finding unstable periodic orbits of any period of Discrete maps](periodicity).

The following methods are currently "wanted features", that will be implemented soon:

* Numeric Computation of Kolmogorov-Sinai entropy.
* Definition of chaos, by Ott.
* GALI method for estimating chaotic behavior, TBA.
