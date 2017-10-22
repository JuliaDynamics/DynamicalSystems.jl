![DynamicalSystems.jl logo: The Double Pendulum](https://i.imgur.com/nFQFdB0.gif)

# Introduction
`DynamicalSystems.jl` is a Julia package for the exploration of continuous and discrete dynamical systems. It aims to be a useful and powerful companion for students and scientists treading
on the fields of Chaos, nonlinear dynamics and dynamical systems in general.

One of a major goals of this package is to be completely transparent as to what is
going on "under the hood". In scientific research, you never want to use *black boxes*,
e.g. functions that give a result without telling you how it was calculated. `DynamicalSystems.jl` battles this in 3 ways: Firstly, it is written entirely in Julia,
making the source code clear and easy to understand for even novice users. Secondly,
almost every documentation string gives
**direct references to the original papers** where the algorithm is taken from, in case some users don't understand (or simply don't want to read) the source code. For example,
the documentation string of [`lyapunovs`](@ref) will cite relevant publications for the definition and computation of the lyapunov spectrum. Thirdly, all documentation strings
for all exported names have very detailed descriptions of the algorithms (whenever
it is possible).

*You
can [join our chatroom](https://gitter.im/JuliaDynamics/Lobby) for discussions related
to dynamical systems and Julia as well as for asking questions about the packages of the
JuliaDynamics organization!*

Be sure to visit the [Contributor Guide](contributors_guide) page, because you can
help make this package better without having to write a single line of code!
Also, if you find this package helpful please consider staring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl)! This gives us an
accurate lower bound of users that this package has already helped!

## Installation
This package is registered. Simply use `Pkg.add("DynamicalSystems")` to install it.

Bug-fixes and upgrades are constantly fed to the master branch, which accessed with `Pkg.checkout("DynamicalSystems")` after installing. On the other hand the master
branch may also have breaking changes and therefore caution is advised!

*The [stable](https://JuliaDynamics.github.io/DynamicalSystems.jl/stable) documentation refers to the version of the package installed with `Pkg.add()`. The [latest](https://JuliaDynamics.github.io/DynamicalSystems.jl/latest) documentation refers to the version under development, obtained with `Pkg.checkout("DynamicalSystems")`.*

To ensure that your installation works perfectly, you can use `Pkg.test("DynamicalSystems")`.


## Contents

### [System Definition](system_definition)

1. Intuitive, consistent APIs for the definition of general dynamical systems. The currently supported system types are:

    * [Discrete Maps](system_definition/#discrete-systems)
    * [Continuous Flows](system_definition/#continuous-systems)
    * [Numerical Data](system_definition/#numerical-data)


2. Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user.
3. Interface for [`DifferentialEquations.jl`](http://docs.juliadiffeq.org/latest/index.html) for flexible integration of continuous system.
4. Well-defined functions for (numerically) evolving dynamical systems.
5. Dedicated interface ([`Dataset`](@ref)) for handling sets of data, in a way that feels familiar to scientists but is also fast.
6. Library of predefined well-known dynamical systems that have been used extensively in scientific research.

### [Lyapunov Exponents](lyapunovs)

The following treat systems where the equations of motion are known:

1. Maximum Lyapunov exponent for both discrete and continuous systems: [`lyapunov`](@ref).
2. Lyapunov *spectrum* for both discrete and continuous systems: [`lyapunovs`](@ref).
    * It is also possible to obtain the convergence timeseries of the above algorithms,
      in order to e.g. estimate better parameters for faster convergence.

### [Entropies and Dimensions](entropies)

1. Generalized (Renyi) entropy and all related entropies: [`genentropy`](@ref).

    * Ultra-fast and cheap method for computing entropies of large datasets without ever having to worry about memory overflow.


2. Generalized dimensions (e.g. capacity dimension, information dimension, etc.): [`generalized_dim`](@ref).
3. Kaplan-Yorke dimension: [`kaplanyorke_dim`](@ref).
4. Automated detection of best algorithmic parameters for calculating attractor dimensions.

And, in order to automatically deduce dimensions, we also offer methods for:

* Partitioning a function $y(x)$ vs. $x$ into regions where it is approximated by a straight line, using a flexible algorithm with a lot of control over the outcome. See [`linear_regions`](@ref).
* Detection of largest linear region of a function $y(x)$ vs. $x$ and extraction of the slope of this region.

### [Nonlinear Timeseries Analysis](nlts)

1. Flexible and abstracted [`Reconstruction`](@ref) interface, that creates the delay-coordinates reconstruction of a timeseries efficiently.
2. Methods for estimating good `Reconstruction` parameters (delay and dimension).
3. *Four* different algorithms for numerically determining the maximum Lyapunov exponent of a (e.g. experimentally) measured timeseries: [`numericallyapunov`](@ref).

    * Fast computation of the above algorithms made possible by the interaction of [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl), multiple dispatch and smart indexing (through the `Reconstruction` abstraction).

### [Periodicity](periodicity)

1. Numerical method to find unstable and stable fixed points of *any order* $n$ of a discrete map (of any dimensionality): [`periodicorbits`](@ref).

    * Convenience functions for defining and realizing all possible combinations of $\mathbf{\Lambda}_k$ matrices required in the above method.

### Wanted Features
The [wanted features GitHub page](https://github.com/JuliaDynamics/DynamicalSystems.jl/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Awanted_feature) lists features that are wanted by the `DynamicalSystems`, and are open to contributors.
