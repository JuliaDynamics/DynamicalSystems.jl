For **DynamicalSystems.jl** a "system" is simple a structure that contains the system's state, the equations of motion and the Jacobian. The last two are *functions* that take as an input a state. It is highly advised to create a `DynamicalSystem`
using Functors (see below).

The above of course stand for systems where one already *knows* the equations of motion.
if instead, your "system" is in the form of [numerical data](definition/dataset), then see the appropriate section.

All core definitions in **DynamicalSystems.jl** are contained in the [DynamicalSystemsBase.jl](https://github.com/JuliaDynamics/DynamicalSystemsBase.jl) Julia package and are necessarily required in every other package of this ecosystem.
All system `sturct`s are also a subtype of the abstract type `DynamicalSystem`.

!!! info "Trajectory and Timeseries"
    The word "timeseries" can be very confusing, because it can mean a univariate (also called scalar or one-dimensional)
    timeseries or a multivariate (also called multi-dimensional) timeseries. To resolve this confusion, in
    **DynamicalSystems.jl** we have the following convention: **"timeseries"** always
    refers to a one-dimensional vector of numbers, which exists with respect to
    some other one-dimensional vector of numbers that corresponds to a time-vector.
    On the other hand,
    the word **"trajectory"** is used to refer to a *multi-dimensional* timeseries,
    which is of course simply a group/set of one-dimensional timeseries.

    Note that the data representation of a "trajectory" in Julia may vary: from
    a 2D Matrix to independent Vectors. In our package, a trajectory is always
    represented using a [`Dataset`](@ref), which is a `Vector` of `SVector`s, and
    each `SVector` represents a data-point (the values of the variables at a given
    time-point).

---

## General Functions
The following functions are defined for convenience for any dynamical system:
```@docs
dimension
state
jacobian
```
