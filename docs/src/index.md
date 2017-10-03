![DynamicalSystems.jl logo: The Double Pendulum](https://i.imgur.com/nFQFdB0.gif)

# Introduction
`DynamicalSystems.jl` is a Julia package for the exploration of continuous and discrete dynamical systems. It aims to be a useful and powerful companion for students and scientists treading
on the field of Chaos, nonlinear dynamics and dynamical systems in general. You
can [join our chatroom](https://gitter.im/JuliaDynamics/Lobby) for discussions related
to dynamical systems and Julia as well as for asking questions about the packages of the
JuliaDynamics organization!

Currently the package
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
Also, if you find this package helpful please consider staring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl)! This gives us an
accurate lower bound of users that this package has already helped!

## Contents
This is the list of what this package currently offer (updated very frequently):

1. [System Definition](system_definition)
      * Intuitive, consistent APIs for the definition of general dynamical systems
      * Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user.
      * Interface for [`DifferentialEquations.jl`](http://docs.juliadiffeq.org/latest/index.html) for
      flexible integration of continuous system.
      * Well-defined functions for evolving dynamical systems.
      * Dedicated interface (`Dataset`) for handling sets of data (numbers), in a way that
      feels familiar to scientists.
      * Library of predefined well-known dynamical systems that have been used
      extensively in scientific research.
3. [Lyapunov Exponents](lyapunovs)
      * Maximum Lyapunov exponent for both discrete and continuous systems.
      * Lyapunov *spectrum* for both discrete and continuous systems.
4. [Entropies and Dimensions](entropies)
      * Generalized (Renyi) entropy and all related entropies.
      * Ultra-fast and cheap method for computing entropies of large datasets
      without ever having to worry about memory overflow. It uses typically 3-6 orders
      of magnitude less time and memory than traditional histogram-based methods.
      * Generalized Dimensions (e.g. capacity dimension, information dimension, etc.).
      * Kaplan-Yorke dimension.
      * Partitioning of a function $y(x)$ vs. $x$ into regions where it is approximated
      by a straight line, using a flexible function with a lot of control over the outcome.
      * Detection of largest linear region of a function $y(x)$ vs. $x$ and extraction
      of the slope of this region (used e.g. in estimating dimensions of chaotic tractors).
      * Methods for detecting best algorithmic parameters for calculating attractor
      dimensions, including a fast implementation of minimum pairwise distance of a
      `Dataset`.
6. [Nonlinear Timeseries Analysis](nlts)
      * Flexible and abstracted `Reconstruction` interface, that creates
      the delay-coordinates reconstruction of a 1D timeseries efficiently.
      * Methods for estimating good `Reconstruction` parameters (delay and dimension).
      * *Four* different algorithms for numerically determining the maximum Lyapunov
      exponent of a (e.g. experimentally) measured timeseries.
      * Fast computation of the above algorithms made possible by the interaction of [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl), multiple dispatch and smart indexing (through the `Reconstruction` abstraction).
7. [Periodicity](periodicity)
      * Numerical method to find unstable and stable fixed points of *any order* of
      a discrete map (of any dimensionality). Fixed points of order $n>1$ are simply
      periodic orbits of order $n$.
      * Convenience functions for defining and realizing all possible combinations of
      $\mathbf{\Lambda}_k$ matrices required in the above method.


The [wanted features GitHub page](https://github.com/JuliaDynamics/DynamicalSystems.jl/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Awanted_feature) lists features that are wanted by the `DynamicalSystems`, and are open to contributors.
