# System Definition
For `DynamicalSystems.jl` a system is simple a structure that contains the system's state, the equations of motion and the Jacobian. The last two are *functions* that take as an input a state.

This of course stands for systems where one already **knows the equations of motion**.
if instead, your "system" is in the form of [numerical data](#numerical-data), then see the appropriate section.


!!! warning "Non-autonomous systems"
    This package does **not** accept non-autonomous systems. To use such systems with this package increase
    the dimensionality of your system by 1, by introducing an additional variable
    `τ` such that `dτdt = 1` (or `τ_next = τ_prev + 1`).
    This additional variable will serve as
    the "time" in your equations of motion.


## Discrete Systems
Discrete systems are of the form:
```math
\vec{x}_{n+1} = \vec{f}(\vec{x}_n).
```
The Type representing such systems is called `DiscreteDS`:
```@docs
DiscreteDS
```

The documentation string of the constructor is perfectly self-contained, but for the sake of clarity we will go through all the steps in the following.

`state` is simply the state the system starts (a.k.a. initial conditions) and
`eom` is a *function* that takes a `state` as an input and returns the next state
as an output.
The `jacob` is also a *function* that takes a `state` as an input and returns the
Jacobian matrix of the system (at this state).
This however is optional and if not provided by the user, will be calculated automatically using the package `ForwardDiff.jl`.

!!! note "Return form of the `eom` function"
    It is **heavily** advised that the equations of motion `eom` function returns an `SVector` from
    the julia package [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) and similarly the `jacob` function returns an `SMatrix`. [Numerous benchmarks](https://github.com/Datseris/DynamicalSystems.jl/tree/master/test/benchmarks) have been made in order to deduce the most efficient possible way to define
    a system, and this way was proved to be the best when the system's dimension is small.

For example, let's create one of the [Predefined Systems](#predefined-systems) offered by this package, the Hénon map:
```julia
using DynamicalSystems
using StaticArrays #only necessary when defining a system

eom_henon(x) = SVector{2}(1.0 - a*x[1]^2 + x[2], b*x[1])
jacob_henon(x) = @SMatrix [-2*a*x[1] 1.0; b 0.0]

hen = DiscreteDS(rand(2), eom_henon, jacob_henon)
```
If we did not want to write a Jacobian, we could do
```julia
hen_nojac = DiscreteDS(rand(2), eom_henon)
```
and the Jacobian function would be created automatically.

### 1-dimensional Discrete Systems
In the case of maps, there a special structure for one-dimensional systems, since
they are commonly used in scientific research.
The syntax is `DiscreteDS1D(state, eom [, deriv])`.
In this one-dimensional case, you don't need to worry about `StaticArrays.jl`
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
$$
\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}),
$$
are defined in a similar manner with the discrete systems:
```@docs
ContinuousDS
```
There are two major differences compared to the discrete case:

1. The second field `eom!` ends with an `!` to remind users that it is an in-place
   function. This is necessary because the integration of continuous systems using
   [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl)
   is much better this way.
2. Automated Jacobian function evaluation is not yet supported due to the dissonance
   of the interfaces of [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) and [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)

Notice that providing a Jacobian is not necessary, since currently it is only used by
the function [`lyapunovs`](lyapunovs/#DynamicalSystems.lyapunovs).
If you do provide a Jacobian,
it is best if it returns an `SMatrix`, just like with the discrete systems case.

As an example, the continuous Rössler system can be defined as:
```julia
@inline @inbounds function eom_roessler!(du, u)
    a = 0.2; b = 0.2; c = 5.7
    du[1] = -u[2]-u[3]
    du[2] = u[1] + a*u[2]
    du[3] = b + u[3]*(u[1] - c)
end
@inline @inbounds function jacob_roessler(u)
    i = one(eltype(u))
    o = zero(eltype(u))
    @SMatrix [o     -i      -i;
              i      a       o;
              u[3]   o       u[1] - c]
end

ros = ContinuousDS(rand(3), eom_roessler!, jacob_roessler)
```


## System evolution
`DynamicalSystems.jl` provides convenient interfaces for the evolution of systems. Especially in the continuous case, an interface is provided to the module `DifferentialEquations.jl`, with an approach that fits more the structuring of the present package (e.g. time is never passed to the equations of motion).

Notice that if you want to do repeated evolutions of a system, use the
`evolve!(::ODEProblem)` interface, which does not create a new `ODEProblem` every time.

These are the functions related to system-evolution:
```@docs
evolve
evolve!
timeseries
```

In addition, interfaces are provided for usage directly with [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl), by giving additional constructors:
```@docs
ODEProblem
ODEIntegrator
```


## Numerical Data
Numerical data in `DynamicalSystems.jl` are represented by a structure called
`Dataset`:
```@docs
Dataset
```
In essence a `Dataset` is simply a container for a `Vector` of `Vector`s, but only for
cases where the all inner vectors are of equal size.
However, it
is visually represented as a matrix, similarly to how numerical data would be printed
on a spreadsheet (with time being the *column* direction). It also offers a lot more
functionality and pretty-printing. Besides the examples in the documentation string,
you can also do:
```julia
data = timeseries(hen, 10000)
for point in data
# do stuff with each datapoint (vector with as many elements as system dimension)
end
```

All functions that manipulate and use data are expecting a `Dataset` instance. If given
a matrix, they will first convert to `Dataset`. This means that you should first
convert your data to a `Dataset` if you want to call functions more than once, to avoid
constantly converting.

## Predefined Systems
Predefined systems exist in the `Systems` submodule exported by `DynamicalSystems`, in the form of functions that return a `DynamicalSystem`. They are accessed
like:
```julia
using DynamicalSystems
ds = Systems.lorenz(ρ = 32.0)
ts = timeseries(ds, 10.0)
```

All of these functions have very similar documentation strings:

1. Call signature (parameters of the system are always passed as keyword arguments).
1. Introductory text about what this system is and who introduced it first.
2. Couple of sentences that contain cool science info about the system.
3. Reference to the original papers.

So far, the predefined systems that exist in the `Systems` sub-module are:
```@autodocs
Modules = [Systems]
Order   = [:function, :type]
```
