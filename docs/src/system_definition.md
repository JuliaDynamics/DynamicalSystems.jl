# System Definition
For DynamicalSystems.jl a system is simple a structure that contains the system's state, the equations of motion and the Jacobian. The last two are *functions* that take as an input a state. It is highly advised to create a `DynamicalSystem`
using Functors (see below).

The above of course stand for systems where one already *knows* the equations of motion.
if instead, your "system" is in the form of [numerical data](#numerical-data), then see the appropriate section.


!!! warning "Non-autonomous systems"
    This package does **not** accept non-autonomous systems. To use such systems with this package increase
    the dimensionality of your system by 1, by introducing an additional variable
    `τ` such that `dτdt = 1` (or `τ_next = τ_prev + 1`).
    This additional variable will serve as
    the "time" in your equations of motion.

!!! info "Trajectory and Timeseries"
    The word "timeseries" can be very confusing, because it can mean a univariate (also called scalar or one-dimensional)
    timeseries or a multivariate (also called multi-dimensional) timeseries. To resolve this confusion, in
    DynamicalSystems.jl we have the following convention: **"timeseries"** always
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
## Continuous Systems
Continuous systems of the form
```math
\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}),
```
are defined almost identically with the [`BigDiscreteDS`](@ref) systems:
```@docs
ContinuousDS
```

Notice that the fields `eom!` and `jacob!` end with a `!`, to remind users
that these functions should operate in-place. Also notice that the type `ContinuousDS`
is actually immutable.

You can give any function that complies with the requirements stated by the documentation string.
However, it is highly advised to use [Functors](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1) for dynamical systems where the equations
of motion contain parameters.

### Defining a `DynamicalSystem` using Functors
A Functor is a shorthand for saying [Function-like objects](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1),
i.e. `struct`s that are also callable (see the linked documentation page). Using
such objects one can create both the equations of motion and a parameter container
under a single `struct` definition.

For example, let's take a look at the source code that generates continuous Rössler
system, from the [Predefined Systems](#predefined-systems):
```@example 1
using DynamicalSystems
mutable struct Rössler
    a::Float64
    b::Float64
    c::Float64
end
function (s::Rössler)(du::EomVector, u::EomVector)
    du[1] = -u[2]-u[3]
    du[2] = u[1] + s.a*u[2]
    du[3] = s.b + u[3]*(u[1] - s.c)
    return nothing
end
function (s::Rössler)(J::EomMatrix, u::EomVector)
    J[2,2] = s.a
    J[3,1] = u[3]; J[3,3] = u[1] - s.c
    return nothing
end
nothing # hide
```
The first code-block defines a `struct` that is simply a container for the
parameters of the Rössler system. The second code-block defines the equations
of motion of the system, by taking advantage of the fact that you can *call*
this `struct` as if it was a function:
```@example 1
s = Rössler(1,2,3)
u = rand(3); du = copy(u)
s(du, u)
du == u
```
The third code-block then defines the Jacobian function using multiple dispatch.
This allows us to use the *same* instance of `Rössler` for *both* the equations
of motion *and* the Jacobian function!

!!! important "Use `EomVector` and `EomMatrix`"
    The Types `EomVector` and `EomMatrix` are simple aliases exported by our package. They represent
    a `Union{}` over all possible vectors or matrices that are involved in the
    equations of motion as used in DynamicalSystems.jl (namely `Vector`, `SubArray`
    and `SVector`).

    **You must define your equations of motion / Jacobian functions using these
    Types instead of a simple `Vector` or `Matrix`, otherwise functions like
    `lyapunovs` won't work properly.**

The possibility of providing an initialized
Jacobian to the `ContinuousDS` constructor allows us to "cheat".
Notice that the Jacobian function only accesses
fields that depend on the parameters and/or the state variables, because the other
fields are constants and will be initialized properly later.

Next, we define a "set-up" function, that returns a `ContinuousDS`:
```@example 1
# this is in fact the function Systems.roessler()
function roessler(u0=rand(3); a = 0.2, b = 0.2, c = 5.7)
    # Initialize Jacobian matrix:
    i = one(eltype(u0))
    o = zero(eltype(u0))
    J = zeros(eltype(u0), 3, 3)
    J[1,:] .= [o, -i,      -i]
    J[2,:] .= [i,  a,       o]
    J[3,:] .= [u0[3], o, u0[1] - c]
    s = Rössler(a, b, c)
    # Pass the same system to both fields!
    return ContinuousDS(u0, s, s, J)
end

ds = roessler()
```
Then, it is trivial to change a parameter of the system by e.g. doing `ds.eom!.c = 2.2`.
Notice that this parameter change will affect both the equations of motion as well
as the Jacobian function, making everything concise and easy-to-use!


!!! info "`Vectors` vs. `SVectors` for continuous systems."
    There is no distinction based on the size of the system for continuous because using `SVectors` or in-place operations with normal `Vectors` yield almost no speed differences in conjunction with [DifferentialEquations.jl](http://docs.juliadiffeq.org/stable/index.html) for small
    dimensions.

An example of a system definition without Functors is [also shown here](#defining-a-dynamicalsystem-without-functors).

---

## Discrete Systems
Discrete systems are of the form:
```math
\vec{x}_{n+1} = \vec{f}(\vec{x}_n).
```
DynamicalSystems.jl categorizes discrete systems in three cases, due to the
extremely performant handling that [`StaticArrays`](https://github.com/JuliaArrays/StaticArrays.jl) offers for small dimensionalities.


### High-Dimensional Discrete Systems
At around `D=10` dimensions, Static Arrays start to become less efficient than Julia's
base Arrays, provided that the latter use in-place operations. For cases of
discrete systems with much high dimensionality, we offer a
type called `BigDiscreteDS`:
```@docs
BigDiscreteDS
```

This system is identical to [`ContinuousDS`](@ref) as far as definition is concerned.
All operations are done in place, and the type is immutable. The same suggestions
about using Functors and initialized Jacobians also apply here.

See the source code of the pre-defined [coupled standard maps](https://juliadynamics.github.io/DynamicalSystems.jl/latest/system_definition/#DynamicalSystems.Systems.coupledstandardmaps) for an example of a `BigDiscreteDS` definition.

---

### Low-dimensional Discrete Systems
The definition of low-dimensional discrete systems differs fundamentally from
high dimensional ones, because everything is *much* more efficiently done with
statically sized vectors. The Type representing such systems is called `DiscreteDS`:
```@docs
DiscreteDS
```

The documentation string of the constructor is perfectly self-contained, but for the sake of clarity we will go through all the steps in the following.

`state` is simply the state the system starts (a.k.a. initial conditions) and
is always of type `SVector` from [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl).
`eom` is a *function* that takes a `state` as an input and returns the next state
as an output.
The `jacob` is also a *function* that takes a `state` as an input and returns the
Jacobian matrix of the system (at this state). So far this is actually
different than `BigDiscreteDS` where the functions where in-place.

!!! note "Return form of the `eom` function"
    It is **heavily** advised that the equations of motion `eom` function returns an `SVector` from
    the julia package [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) and similarly the `jacob` function returns an `SMatrix`. [Numerous benchmarks](https://github.com/Datseris/DynamicalSystems.jl/tree/master/test/benchmarks) have been made in order to deduce the most efficient way to define
    a system, and this way was proved to be the best when the system's dimension is small.

Something important to note when defining a `DiscreteDS` using Functors: since
the function calls take only one argument (always a state), it is impossible to
use multiple dispatch to differentiate between a call to the e.o.m. or the Jacobian
functions.

However, it is very easy to still define both function calls
using a single `struct`, by using
a 2 argument function given to the constructor. For example:
```@example henon
using DynamicalSystems
using StaticArrays # only necessary when defining a system

function henon(u0=zeros(2); a = 1.4, b = 0.3)
    he = HénonMap(a,b)
    # The jacobian function: (still uses the HénonMap)
    @inline jacob_henon(x) = he(x, nothing)
    return DiscreteDS(u0, he, jacob_henon)
end
mutable struct HénonMap
    a::Float64
    b::Float64
end
(f::HénonMap)(x::EomVector) = SVector{2}(1.0 - f.a*x[1]^2 + x[2], f.b*x[1])
(f::HénonMap)(x::EomVector, no::Void) = @SMatrix [-2*f.a*x[1] 1.0; f.b 0.0]

ds = henon()
```
In this case, doing `ds.eom.a = 2.5` would still affect *both* the equations
of motion as well as the Jacobian, making everything work perfectly!


#### Defining a `DynamicalSystem` without Functors
As an example of defining a system without
using Functors, let's create another one of the [Predefined Systems](#predefined-systems) offered by this package, the Folded Towel map.

Because this map doesn't have any parameters, it is unnecessary to
associate a Functor with it.
```@example 2
using DynamicalSystems
using StaticArrays # only necessary when defining a system

@inline @inbounds function eom_towel(x)
x1, x2, x3 = x[1], x[2], x[3]
SVector( 3.8*x1*(1-x1) - 0.05*(x2+0.35)*(1-2*x3),
0.1*( (x2+0.35)*(1-2*x3) - 1 )*(1 - 1.9*x1),
3.78*x3*(1-x3)+0.2*x2 )
end

@inline @inbounds function jacob_towel(x)
    @SMatrix [3.8*(1 - 2x[1]) -0.05*(1-2x[3]) 0.1*(x[2] + 0.35);
    -0.19((x[2] + 0.35)*(1-2x[3]) - 1)  0.1*(1-2x[3])*(1-1.9x[1])  -0.2*(x[2] + 0.35)*(1-1.9x[1]);
    0.0  0.2  3.78(1-2x[3]) ]
end

u0=[0.085, -0.121, 0.075]
towel =  DiscreteDS(u0, eom_towel, jacob_towel)
```
If we did not want to write a Jacobian for it, we could do
```julia
towl_nojac = DiscreteDS(rand(3), eom_towel)
```
and the Jacobian function is created automatically.

### 1-dimensional Discrete Systems
In the case of maps, there a special structure for one-dimensional systems.
The syntax is `DiscreteDS1D(state, eom [, deriv])`.
In this one-dimensional case, you don't need to worry about StaticArrays.jl
because everything is in plain numbers. For example:
```@example 3
using DynamicalSystems

@inline eom_logistic(r) = (x) -> r*x*(1-x)  # this is a closure
@inline deriv_logistic(r) = (x) -> r*(1-2x) # this is a closure
r = 3.7
logistic = DiscreteDS1D(rand(), eom_logistic(r), deriv_logistic(r))
```
Once again, if you skip the derivative functions it will be calculated automatically
using ForwardDiff.jl.


---

## System evolution
DynamicalSystems.jl provides convenient interfaces for the evolution of systems.
```@docs
evolve
trajectory
```

Especially in the continuous case, an API is provided for usage directly with [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl), by giving additional constructors:
```@docs
ODEProblem
ODEIntegrator
variational_integrator
```
Notice that if you want to do repeated evolutions of different states of a
continuous system,
you should use the
`ODEIntegrator(ds::DynamicalSystem)` in conjunction with `DifferentialEquations.reinit!(integrator, newstate)` to avoid the intermediate initializations of the integrator each time.
---

## Coordination with other packages
One more advantage of using [Functors](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1), is that it offers a universal definition of your
equations of motion that fits both the expected structure of DynamicalSystems.jl as well as other packages.

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

---

## Numerical Data
Numerical data in DynamicalSystems.jl is represented by a structure called
`Dataset`:
```@docs
Dataset
```
---
In essence a `Dataset` is simply a container for a `Vector` of `SVector`s.
However, it
is visually represented as a matrix, similarly to how numerical data would be printed
on a spreadsheet (with time being the *column* direction). It also offers a lot more
functionality than just pretty-printing.
Besides the examples in the documentation string,
you can also do:
```julia
using DynamicalSystems
hen == Systems.henon()
data = trajectory(hen, 10000)
for point in data
# do stuff with each datapoint (vector with as many elements as system dimension)
end
```

All functions of our package that manipulate and use data are expecting an `AbstractDataset` subtype. This allows us to define efficient methods that coordinate
well with other packages, like e.g. [`neighborhood`](@ref).

If given a matrix, we first convert to `Dataset`. This means that you should first
convert your data to a `Dataset` if you want to call functions more than once, to avoid
constantly converting.

---

## Predefined Systems
Predefined systems exist in the `Systems` submodule exported by `DynamicalSystems`, in the form of functions that return a `DynamicalSystem`. They are accessed
like:
```julia
using DynamicalSystems
ds = Systems.lorenz(ρ = 32.0)
typeof(ds) # ContinuousDS
ts = trajectory(ds, 10.0)
```

So far, the predefined systems that exist in the `Systems` sub-module are:
```@autodocs
Modules = [Systems]
Order   = [:function]
```
