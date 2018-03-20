![DynamicalSystems.jl logo: The Double Pendulum](https://i.imgur.com/nFQFdB0.gif)

# Introduction
**DynamicalSystems.jl** is a Julia software library for the exploration of chaos and nonlinear dynamics.

The current documentation was built with the following versions
```@example docs
Pkg.status("DynamicalSystemsBase") # hide
Pkg.status("ChaosTools") # hide
Pkg.status("TimeseriesPrediction") # hide
```
See the [News](news) page for recent updates!

!!! info "Introductory textbooks"
    Our library assumes basic some basic knowledge of nonlinear dynamics and complex systems.

    If you are new to the field but want to learn more, we can suggest the following textbooks as introductions:
    * Nonlinear Dynamics And Chaos - S. Strogatz
    * An Exploration of Dynamical Systems and Chaos - J. Argyris *et al.*
    * Chaos in Dynamical Systems - E. Ott

## Contents

### [Fundamentals](definition/general)

1. Intuitive, consistent APIs for the definition of general [dynamical systems](definition/general), both maps and flows. The following combinations are possible:
    * Continuous or Discrete systems. Continuous systems use [`DifferentialEquations.jl`](http://docs.juliadiffeq.org/latest/) for solving the ODE problem.
    * In-place or out-of-place (large versus small systems).
    * Auto-differentiated or not (for the Jacobian function).


2. Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user.
4. Robust implementations of all kinds of integrators, that evolve the system,
   many states of the system, or even deviation vectors. See the [advanced documentation](advanced) for this.
4. Dedicated interface for [numerical data](definition/dataset).
5. Efficient [`neighborhood`](@ref) estimation by interfacing [`NearestNeighbors`](https://github.com/KristofferC/NearestNeighbors.jl).
5. Delay Coordinates Embedding: flexible and abstracted [`Reconstruction`](@ref) interface, that creates the delay-coordinates reconstruction of a timeseries efficiently.
    * Supports multiple dimensions and multiple timescales.


6. Library of [predefined well-known dynamical systems](definition/predefined) that have been used extensively in scientific research.

### [ChaosTools](chaos/overview)
Please see the [overview section](chaos/overview) for a full list of features. Here
is a quick summary:

* Poincare S.O.S. and orbit diagrams
* Lyapunov Exponents
* Entropies and Dimensions
* Estimation of Reconstruction parameters
* Lyapunov exponent of a timeseries
* Finding Fixed Points of Maps
* Detecting Chaos


### [TimeseriesPrediction](tsprediction/localmodels)

* Predicting the future of one or multiple timeseries using average local models.


## Our Goals
The ultimate goal for **DynamicalSystems.jl** is
to be a useful *library* for students and scientists working on chaos, nonlinear dynamics and
in general dynamical systems. We don't want to have "just code", but also detailed descriptions and references for as many methods as possible.

With **DynamicalSystems.jl** we try to

1. Be concise, intuitive, and general. All functions we offer work just as well with any system, whether it is a simple continuous chaotic system, like the [Lorenz attractor](definition/predefined/#DynamicalSystemsBase.Systems.lorenz), or a high dimensional discrete map like [coupled standard maps](definition/predefined/#DynamicalSystemsBase.Systems.coupledstandardmaps).
2. Be accurate, reliable and performant.
3. Be transparent with respect to what is happening "under the hood", i.e. be clear about exactly what each function call does. We take care of this aspect in many ways; by being well-documented, giving references to scientific papers and having clear source code.

## Installation
Simply use `Pkg.add("DynamicalSystems")` to install *everything*.

## Citing
There is a (very small) paper associated with **DynamicalSystems.jl**. If we have helped
you in research that led to a publication, please be kind enough to cite it, using
the DOI `10.21105/joss.00598` or the following BiBTeX entry:
```
@article{Datseris2018,
  doi = {10.21105/joss.00598},
  url = {https://doi.org/10.21105/joss.00598},
  year  = {2018},
  month = {mar},
  volume = {3},
  number = {23},
  pages = {598},
  author = {George Datseris},
  title = {{DynamicalSystems}.jl: A Julia software library for chaos and nonlinear dynamics},
  journal = {Journal of Open Source Software}
}
```

## Contacting

You
can [join our chatroom](https://gitter.im/JuliaDynamics/Lobby) for discussions and/or questions about the packages of the JuliaDynamics organization!

Be sure to visit the [Contributor Guide](contributors_guide) page, because you can
help make this package better without having to write a single line of code!
Also, if you find this package helpful please consider staring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl)! This gives us an
accurate lower bound of users that this package has already helped!

!!! tip "Use latest documentation"
      We highly suggest our users to read the  [latest](https://JuliaDynamics.github.io/DynamicalSystems.jl/latest) documentation
      and not the [stable](https://JuliaDynamics.github.io/DynamicalSystems.jl/stable) one.
