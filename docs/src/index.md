![DynamicalSystems.jl logo: The Double Pendulum](https://i.imgur.com/nFQFdB0.gif)

# Introduction
**DynamicalSystems.jl** is a Julia suite for the exploration of chaos and nonlinear dynamics.

*You
can [join our chatroom](https://gitter.im/JuliaDynamics/Lobby) for discussions related
to dynamical systems and Julia as well as for asking questions about the packages of the
JuliaDynamics organization!*

Be sure to visit the [Contributor Guide](contributors_guide) page, because you can
help make this package better without having to write a single line of code!
Also, if you find this package helpful please consider staring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl)! This gives us an
accurate lower bound of users that this package has already helped!

## Our Goals
Our aim is for the **DynamicalSystems.jl** ecosystem to be a useful and powerful companion for students and scientists working on chaos and nonlinear dynamics.

One of a major goals of this ecosystem is to be completely transparent as to what is
going on "under the hood". In scientific research, you never want to use *black boxes*,
e.g. functions that give a result without telling you how it was calculated. **DynamicalSystems.jl** battles this in 3 ways:

1. It is written entirely in Julia,
   making the source code clear and easy to understand for even novice users.
2. Almost every documentation string gives
   **direct references to the original papers** where the algorithm is taken from, in case some users don't understand (or simply don't want to read) the source code. For example,
   the documentation string of [`lyapunovs`](@ref) will cite relevant publications for the definition and computation of the lyapunov spectrum.
3. Documentation strings
   for exported names have summarized descriptions of the algorithms (whenever
   it is possible).

Another major goal is to offer code that is concise, intuitive, performant and **general**.
All functions work just as well with *any* `DynamicalSystem`, whether it is a simple
continuous chaotic system, like the Lorenz attractor, or a high dimensional discrete
map like 20 coupled standard maps!

For example, provided you have first defined a [`DynamicalSystem`](definition/general)
(which simply reduces to writing a function for the equations of motion),
you should be able to e.g. calculate the Lyapunov spectrum for it
in a single line:
```julia
lyapunovs(system, times_to_do_QR; keywords...)
```
The same function call works with any system, no discriminations here!

## Installation
Simply use `Pkg.add("DynamicalSystems")` to install *everything*.

We highly suggest our users to read the The  [latest](https://JuliaDynamics.github.io/DynamicalSystems.jl/latest) documentation
and not the [stable](https://JuliaDynamics.github.io/DynamicalSystems.jl/stable) one.
The reasoning is simple: the repository of `DynamicalSystems` is a package coordinator
and documentation host for the **DynamicalSystems.jl** ecosystem. It will often be
that a new tag will exist for one of the packages of the ecosystem but not for the repository of `DynamicalSystems` itself. Thus you can only read the documentation of the latest features by visiting the latest documentation version.

Notice however, that this does not have anything to do with the official release versions
of packages that do contain the actual code. Because our documentation is based on utilizing [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) and the documentation strings of exported function names, you can be assured that the documentation you read on the [latest](https://JuliaDynamics.github.io/DynamicalSystems.jl/latest) page reflects accurately the latest release versions.

### Low Dependency usage
All packages of the **DynamicalSystems.jl** ecosystem have a dependency on DynamicalSystemsBase.jl. By running `Pkg.add("DynamicalSystems")` you install all the packages of the ecosystem.
That is not necessary however, since **DynamicalSystems.jl** is a bridging package
that exports everything and hosts the documentation.

For example, if you only need the features of [ChaosTools.jl](chaos/overview) then
you can get away by doing only `Pkg.add("ChaosTools")` and all other dependencies
will be resolved accordingly.

## Contents

### DynamicalSystemsBase.jl

1. Intuitive, consistent APIs for the definition of general dynamical systems.
2. [Discrete Maps](definition/discrete)
3. [Continuous Flows](definition/continuous)
4. Dedicated interface for [Numerical Data](definition/dataset)
5. Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user.
4. Well-defined functions for (numerically) evolving dynamical systems.
6. Library of [predefined well-known dynamical systems](definition/predefined) that have been used extensively in scientific research.

### ChaosTools.jl
Please see the [overview section](chaos/overview) for a full list of features.

Quick summary:

* Lyapunov Exponents
* Entropies and Dimensions
* Delay Coordinates Embedding
* Neighborhood estimation
* Lyapunov exponent of a timeseries
* Finding Fixed Points of any Map of any order
* Detecting Chaos

### Wanted Features
The [wanted features GitHub page](https://github.com/JuliaDynamics/DynamicalSystems.jl/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Awanted_feature) lists features that are wanted by the **DynamicalSystems.jl** ecosystem and are open to contributors.
