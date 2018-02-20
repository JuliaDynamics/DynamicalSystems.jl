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
See the [News](news) page for recent updates!

## Our Goals
The ultimate goal for **DynamicalSystems.jl** is
to be a useful *library* for scientists working on chaos, nonlinear dynamics and
in general dynamical systems. We don't want to have "just code", but also detailed descriptions and references for as many methods as possible.

With **DynamicalSystems.jl** we try to

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

## Contents

### [DynamicalSystemsBase.jl](definition/general)

1. Intuitive, consistent APIs for the definition of general [dynamical systems](definition/general), both maps and flows. In fact we have implementations for 8 possible dynamical systems:
    * Continuous or Discrete.
    * In-place or out-of-place (large versus small systems).
    * Auto-differentiated or not (for the Jacobian function).

4. Dedicated interface for [numerical data](definition/dataset).
5. Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user.
4. Robust implementations of all kinds of integrators, that evolve the system,
   many states of the system, or even deviation vectors. See the [advanced documentation](advanced) for this.
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
* Finding Fixed Points of Maps
* Detecting Chaos

## Wanted Features
The following lists state features that are wanted by the **DynamicalSystems.jl** ecosystem and are open to contributors. These are structured in the form of GitHub Issues, with the label `wanted_feature`:

* [DynamicalSystemsBase.jl wanted features](https://github.com/JuliaDynamics/DynamicalSystemsBase.jl/issues?q=is%3Aissue+is%3Aopen+label%3Awanted_feature)
* [ChaosTools.jl wanted features](https://github.com/JuliaDynamics/ChaosTools.jl/issues?q=is%3Aissue+is%3Aopen+label%3Awanted_feature)
