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

!!! tip "Use latest documentation"
      We highly suggest our users to read the  [latest](https://JuliaDynamics.github.io/DynamicalSystems.jl/latest) documentation
      and not the [stable](https://JuliaDynamics.github.io/DynamicalSystems.jl/stable) one.

The current documentation was built with the following versions
```@example docs
Pkg.status("DynamicalSystemsBase") # hide
Pkg.status("ChaosTools") # hide
```


## Our Goals
Our aim is for the **DynamicalSystems.jl** ecosystem to be a useful and powerful companion for students and scientists working on chaos and nonlinear dynamics.

Our goals with this ecosystem can be summarized in the following three:

1. Be concise, intuitive, and general. All functions we offer work just as well with any system, whether it is a simple continuous chaotic system, like the [Lorenz attractor](definition/predefined/#DynamicalSystemsBase.Systems.lorenz), or a high dimensional discrete map like [coupled standard maps](definition/predefined/#DynamicalSystemsBase.Systems.coupledstandardmaps).
2. Be accurate, reliable and performant.
3. Be transparent with respect to what is happening "under the hood", i.e. be clear about exactly what each function call does. We take care of this aspect in many ways; by being well-documented, giving references to scientific papers and having clear source code.

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

### Low Dependency usage
By running `Pkg.add("DynamicalSystems")` you install all packages of the ecosystem.
That is not necessary however, since **DynamicalSystems.jl** is a bridging package
that exports everything and hosts the documentation.

For example, if you only need the features of [ChaosTools.jl](chaos/overview) then
you can get away by doing only `Pkg.add("ChaosTools")` and all other dependencies
will be resolved accordingly.

## Contents

### [DynamicalSystemsBase.jl](definition/general)

1. Intuitive, consistent APIs for the definition of general dynamical systems.
2. [Discrete Maps](definition/discrete)
3. [Continuous Flows](definition/continuous)
4. Dedicated interface for [Numerical Data](definition/dataset)
5. Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user.
4. Well-defined functions for (numerically) evolving dynamical systems.
6. Library of [predefined well-known dynamical systems](definition/predefined) that have been used extensively in scientific research.

### [ChaosTools.jl](chaos/overview)
Please see the [overview section](chaos/overview) for a full list of features.

Quick summary:

* Poincare S.O.S. and orbit diagrams
* Lyapunov Exponents
* Entropies and Dimensions
* Delay Coordinates Embedding
* Neighborhood estimation
* Lyapunov exponent of a timeseries
* Finding Fixed Points of any Map of any order
* Detecting Chaos

## Wanted Features
The following lists state features that are wanted by the **DynamicalSystems.jl** ecosystem and are open to contributors. These are structured in the form of GitHub Issues, with the label `wanted_feature`:

* [DynamicalSystemsBase.jl wanted features](https://github.com/JuliaDynamics/DynamicalSystemsBase.jl/issues?q=is%3Aissue+is%3Aopen+label%3Awanted_feature)
* [ChaosTools.jl wanted features](https://github.com/JuliaDynamics/ChaosTools.jl/issues?q=is%3Aissue+is%3Aopen+label%3Awanted_feature)
