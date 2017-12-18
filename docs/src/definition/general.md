For **DynamicalSystems.jl** a system is simple a structure that contains the system's state, the equations of motion and the Jacobian. The last two are *functions* that take as an input a state. It is highly advised to create a `DynamicalSystem`
using Functors (see below).

The above of course stand for systems where one already *knows* the equations of motion.
if instead, your "system" is in the form of [numerical data](definition/dataset), then see the appropriate section.

All core definitions in **DynamicalSystems.jl** are contained in the [DynamicalSystemsBase.jl](https://github.com/JuliaDynamics/DynamicalSystemsBase.jl) Julia package and are necessarily required in every other package of this ecosystem.
All system `sturct`s are also a subtype of the abstract type `DynamicalSystem`.


!!! warning "Non-autonomous systems"
    We do **not** accept non-autonomous systems. To use such systems with this ecosystem increase
    the dimensionality of your system by 1, by introducing an additional variable
    ``\tau`` such that ``d\tau/dt = 1`` (or ``\tau_{n+1} = \tau_n + 1``).
    This additional variable will serve as
    the "time" in your equations of motion.

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

## Coordination with other packages

*(this section assumes that you have read through the basic documentation on defining
a system)*

One advantage of using [Functors](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1), is that it offers a universal definition of your
equations of motion that fits both the expected structure of **DynamicalSystems.jl** as well as other packages.

For example, take a look at the following code:
```julia
mutable struct Rössler
    a::Float64
    b::Float64
    c::Float64
end
# E.o.m. as expected by DynamicalSystems.jl:
@inline @inbounds function (s::Rössler)(du::EomVector, u::EomVector)
    du[1] = -u[2]-u[3]
    du[2] = u[1] + s.a*u[2]
    du[3] = s.b + u[3]*(u[1] - s.c)
    return nothing
end
# E.o.m. as expected by e.g. LTISystems.jl or OrdinaryDiffEq.jl:
@inline @inbounds (s::Rössler)(t::Real, u::EomVector, du::EomVector) = s(du, u)
```
You could then give `s = Rössler(1, 2)` to `ContinuousDS` *as well as* the interfaces
of e.g. [LTISystems.jl](https://github.com/JuliaSystems/LTISystems.jl)
and everything will "just work"!
