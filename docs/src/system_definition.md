# System Definition
For `DynamicalSystems.jl` a system is simple a structure that contains the system's state, the equations of motion and the Jacobian. The last two are *functions* that take as an input a state.

This of course stands for systems where one already **knows the equations of motion**.
if instead, your "system" is in the form of [numerical data](#numerical-data), then see the appropriate section.


By taking advantage of the package [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) an automated Jacobian
function can always be supplemented by the package. More details are enclosed in the individual sections, however the documentation strings of all the constructors are
also self-contained.

!!! warning "Non-autonomous systems"
    This package does **not** accept non-autonomous systems. To use such systems with this package increase
    the dimensionality of your system by 1, by introducing an additional variable
    ``\tau`` such that ``d\tau/dt = 1`` (or ``\tau_{n+1} = \tau_n + 1``). This additional variable will serve as
    the "time" in your equations of motion.

!!! tip "User-defined Jacobian"
    Providing a user-defined Jacobian to the system constructors is faster than using
    the one generated from `ForwardDiff.jl`, albeit slightly. Currently, the only function
    that uses the Jacobian is `lyapunovs`, therefore it is advised to provided a
    user-defined one if you want to use it.

## Discrete Systems
Discrete systems are of the form:
```math
\vec{x}_{n+1} = \vec{f}(\vec{x}_n).
```
The Type representing such systems for dimensionality `D ≤ 10` is called `DiscreteDS`.

The constructor is:
```julia
DiscreteDS(state, eom [, jacob])
```
The documentation string of the constructor is perfectly self-contained, but for the sake of clarity we will go through all the steps in the following.

`state` is simply the state the system starts (a.k.a. initial conditions) and
`eom` is a *function* that takes a `state` as an input and returns the next state
as an output.

The `jacob` is also a *function* that takes a `state` as an input and returns the
Jacobian matrix of the system (at this state). This however is optional and if not provided by the user, will be calculated automatically using the package `ForwardDiff.jl`.

!!! tip "Return form of the `eom` function for small `D`"
    It is **heavily** advised that the equations of motion `eom` function returns an `SVector` from
    the julia package [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) and similarly the `jacob` function returns an `SMatrix`. [Numerous benchmarks](https://github.com/Datseris/DynamicalSystems.jl/tree/master/test/benchmarks) have been made in order to deduce the most efficient possible way to define
    a system, and this way was proved to be the best when the system's dimension is small.

For example, let's create one of the [Predefined Systems](#predefined-systems) offered by this package, the Hénon map:
```julia
using DynamicalSystems
using StaticArrays #only necessary when defining a system

eom_henon(x) = SVector{2}(1.0 - a*x[1]^2 + x[2], b*x[1])
jacob_henon(x) = @SMatrix [-2*a*x[1] 1.0; b 0.0]

ds = DiscreteDS(rand(2), eom_henon, jacob_henon)
```
If we did not want to write a Jacobian (due to e.g. unending laziness), we could
```julia
ds_nojac = DiscreteDS(rand(2), eom_henon)
```
and the Jacobian function would be created automatically.


### Large Discrete Systems
TBA

### 1-dimensional Discrete Systems
In the case of maps, there a special structure for one-dimensional systems, since
they are commonly used in scientific research. The syntax is `DiscreteDS1D(state, eom [, deriv])`. In this one-dimensional case, you don't need to worry about `StaticArrays.jl`
because everything is in plain numbers. For example:
```julia
using DynamicalSystems

@inline eom_logistic(r) = (x) -> r*x*(1-x)  # this is a closure
@inline deriv_logistic(r) = (x) -> r*(1-2x) # this is a closure
r = 3.7
logistic = DiscreteDS1D(rand(), eom_logistic(r), deriv_logistic(r))
```
Once again, if you skip the derivative functions it will be calculated automatically
using `ForwardDiff.jl`.

## Continuous Systems
Continuous systems of the form
```math
\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}),
```
are defined almost identically with the discrete systems. The documentation string
of the constructor `ContinuousDS` has all the necessary information:
```@docs
ContinuousDS
```
Once again, using `StaticArrays` is preferred.

For example, the continuous Rössler system can be defined as:
```julia
eom_roessler(u) =
SVector{3}(-u[2]-u[3], u[1] + a*u[2], b + u[3]*(u[1] - c))

function jacob_roessler(u)
  i = one(eltype(u))
  o = zero(eltype(u))
  @SMatrix [o     -i      -i;
            i      a       o;
            u[3]   o       u[1] - c]
end

ros = ContinuousDS(rand(3), eom_roessler, jacob_roessler)
```

## System evolution
`DynamicalSystems.jl` provides convenient interfaces for the evolution of systems. In general, these are the functions you want to use:
```@docs
evolve
timeseries
ODEProblem
evolve!
```
Especially in the continuous case, an interface is provided to the module `DifferentialEquations.jl`, with an approach that fits more the structuring of the present package (e.g. time is never passed to the equations of motion). Also, the function `timeseries` is the only one that stores the actual time-series of the system. All the other functions only keep the final state.


## Numerical Data
In the most general case, the numerical data representing the evolution of a system
are in the form of time-series. `DynamicalSystems.jl` accepts two forms of numerical data
in the most common function calls:
```julia
foo(dataset)
bar(vectors...)
```
where the `dataset` is an `N×D` matrix that contains `N` data points of a `D` dimensional
system. The `vectors... = v1, v2, ..., vD` are simply the individual columns of the `dataset` (each column corresponds to a dynamic variable), so that `dataset ≡ hcat(vectors...)`.


## Predefined Systems
Predefined systems exist in the `Systems` submodule exported by `DynamicalSystems.jl`, in the form of functions that return a `DynamicalSystem`.

All of these functions have very similar documentation strings:

1. Call signature (parameters of the system are always passed as keyword arguments).
1. Introductory text about what this system is and who introduced it first.
2. Couple of sentences that contain cool science info about the system.
3. Reference to the original papers.

For example, the documentation of the [Lorenz system](https://en.wikipedia.org/wiki/Lorenz_system) reads:
```@docs
Systems.lorenz
```
So far, the predefined systems that exist in the `Systems` sub-module are:
```@autodocs
Modules = [Systems]
Order   = [:function, :type]
```
