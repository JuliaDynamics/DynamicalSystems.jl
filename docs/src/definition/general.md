All core definitions for **DynamicalSystems.jl** are contained in [DynamicalSystemsBase.jl](https://github.com/JuliaDynamics/DynamicalSystemsBase.jl).

For **DynamicalSystems.jl** a "dynamical system" is a simple structure with
three fundamental parts:

1. The state,
2. The equations of motion function and
3. The Jacobian function.

The last two are *functions* that take as an input a state as well as the parameters
of the model. Depending on the type, some
dynamical system types may also have some other fields that are of minor importance.

The above "definition" of course stands for systems where one already *knows* the equations of motion. if instead, your "system" is in the form of [numerical data](definition/dataset), then see the appropriate section.

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

## Definition Table

Here is a handy table that summarizes in what form should be the functions required for the equations of motion and the Jacobian, for each system type:

|   System Type   |   Equations of Motion  |         Jacobian         |
|:---------------:|:----------------------:|:------------------------:|
|  `ContinuousDS` |   `eom!(du, u, p, t)`  |   `jacob!(J, u, p, t)`   |
| `BigDiscreteDS` |   `eom!(xnew, x, p)`   |     `jacob!(J, x, p)`    |
|   `DiscreteDS`  | `eom(x, p) -> SVector` | `jacob(x, p) -> SMatrix` |
|  `DiscreteDS1D` |  `eom(x, p) -> Number` |  `deriv(x, p) -> Number` |


!!! tip "Use mutable containers for the parameters"
    It is highly suggested to use a subtype of `Array` or [`LMArray`](https://github.com/JuliaDiffEq/LabelledArrays.jl) for the container
    of the model's parameters. Some functions offered by **DynamicalSystems.jl**,
    like e.g. [`orbitdiagram`](@ref),
    assume that the parameters can be first accessed by `p[x]` with `x` some qualifier
    as well as that this value can be set by `p[x] = newvalue`.

    The [Labelled Arrays](https://github.com/JuliaDiffEq/LabelledArrays.jl) package
    offers `Array` implementations that can be accessed both by index as
    well as by some name.


## General Functions
The following functions are defined for convenience for any dynamical system:
```@docs
dimension
state
jacobian
set_state!
```
