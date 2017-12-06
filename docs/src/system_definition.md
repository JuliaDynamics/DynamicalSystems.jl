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




## Continuous Systems
Continuous systems of the form
```math
\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}),
```
are defined almost identically with the [`BigDiscreteDS`](@ref) systems:
```@docs
ContinuousDS
```
---
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
nothing #hide
```
The first code-block defines a `struct` that is simply a container for the
parameters of the Rössler system. The second code-block defines the equations
of motion of the system, by taking advantage of the fact that you can *call*
this `struct` as if it was a function:
```@example 1
s = Rössler(1,2,3)
u = rand(3); du = copy(du)
s(du, u)
du == u
```
The third code-block then defines the Jacobian function using multiple dispatch.
This allows us to use the *same* instance of `Rössler` for *both* the equations
of motion *and* the Jacobian function.

!!! important "Use `EomVector` and `EomMatrix`"
    The Types `EomVector` and `EomMatrix` are simple aliases exported by our package. They represent
    a `Union{}` over all possible vectors or matrices that are involved in the
    equations of motion as used in DynamicalSystems.jl (namely `Vector`, `SubArray`
    and `SVector`).

    **You must define your equations of motion / Jacobian functions using these
    Types instead of a simple `Vector` or `Matrix`, otherwise functions like
    `lyapunovs` won't work properly.**

Also notice that the Jacobian function only accesses
fields that depend on the parameters and/or the state variables, because the other
fields are constants and will be initialized properly later.

Next, we define a "set-up" function, that returns a `ContinuousDS`
```@example 1
# this is in fact the function accessed by Systems.roessler()
function roessler(u0=rand(3); a = 0.2, b = 0.2, c = 5.7)
    # Initialize Jacobian matrix:
    i = one(eltype(u0))
    o = zero(eltype(u0))
    J = zeros(eltype(u0), 3, 3)
    J[1,:] .= [o, -i,      -i]
    J[2,:] .= [i,  a,       o]
    J[3,:] .= [u0[3], o, u0[1] - c]
    s = Rössler(a, b, c)
    name = "Rössler system"
    # Pass the same system to both fields!
    return ContinuousDS(u0, s, s, J; name = name)
end

ds = roessler()
```
Then, it is trivial to change a parameter of the system by e.g. doing `ds.eom!.c = 2.2`.

!!! info "`Vectors` vs. `SVectors` for equations of motion."
    There is no distinction based on the size of the system for the continuous case because using `SVectors` or in-place operations with normal `Vectors` yield almost no speed differences in conjunction with [DifferentialEquations.jl](http://docs.juliadiffeq.org/stable/index.html) for small
    dimensions.

An example of a system definition without Functors is [also shown here](#asd).

## Discrete Systems
Discrete systems are of the form:
```math
\vec{x}_{n+1} = \vec{f}(\vec{x}_n).
```
DynamicalSystems.jl categorizes discrete systems in three cases, due to the
extremely performant handling that [`StaticArrays`](https://github.com/JuliaArrays/StaticArrays.jl) offers for small dimensionalities.


### High-Dimensional Discrete system.
At around `D=10` dimensions, Static Arrays start to become less efficient than Julia's
base Arrays, provided that the latter use in-place operations. For cases of
discrete systems with much high dimensionality, we offer a
type called `BigDiscreteDS`:
```@docs
BigDiscreteDS
```
---
This system is identical to [`ContinuousDS`](@ref) as far as definition is concerned.
All operations are done in place, and the type is immutable. The same suggestion
about using Functors also applies here.

For example, this is the source code that defines a `BigDiscreteDS` representing
coupled standard maps:

In addition, the possibility of providing an initialized
Jacobian allows one to "cheat". For example, let's look at the definition of the
function for the Jacobian for the [coupled standard maps](system_definition/#DynamicalSystems.Systems.henonhelies):
```julia
### The following are inside a local scope!
J = zeros(eltype(u0), 2M, 2M) #u0 is the state of the system
# Set ∂/∂p entries (they are eye(M,M))
# And they don't change, they are constants
for i in idxs
    J[i, i+M] = 1
    J[i+M, i+M] = 1
end

@inbounds function jacob_coupledsm!(J, x)
    # x[i] ≡ θᵢ
    # x[[idxsp1[i]]] ≡ θᵢ+₁
    # x[[idxsm1[i]]] ≡ θᵢ-₁
    for i in idxs
        cosθ = cos(x[i])
        cosθp= cos(x[idxsp1[i]] - x[i])
        cosθm= cos(x[idxsm1[i]] - x[i])
        J[i+M, i] = ks[i]*cosθ + Γ*(cosθp + cosθm)
        J[i+M, idxsm1[i]] = - Γ*cosθm
        J[i+M, idxsp1[i]] = - Γ*cosθp
        J[i, i] = 1 + J[i+M, i]
        J[i, idxsm1[i]] = J[i+M, idxsm1[i]]
        J[i, idxsp1[i]] = J[i+M, idxsp1[i]]
    end
end
```
The function that evaluates the Jacobian (in-place) only accesses half
of the matrix elements, since the other half is constant and correctly initialized.
Afterwards this function as well as `J` are passed into the constructor
with `BigDiscreteDS(u0, eom_coupledsm!, jacob_coupledsm!, J; name = "something")`.



### Low-dimensional Discrete System
The Type representing such systems is called `DiscreteDS`:
```@docs
DiscreteDS
```
---

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
In the case of maps, there a special structure for one-dimensional systems.
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
using ForwardDiff.jl.




## Dimension of a System
The dimension of any sub-type of `DynamicalSystem` is obtained by `D = dimension(ds)`.

## System evolution
DynamicalSystems.jl provides convenient interfaces for the evolution of systems. Because
all system types are immutable, the new states (initial conditions) can be set with
```@docs
set_state
```
The following functions are related to system evolution:
```@docs
evolve
trajectory
```
---
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
struct Lorenz
```
Notice that dispatch will work even without type annotations in this case,
because one call takes 2 arguments and the other 3.
This way you can simply pass the object `Lorenz96` to the constructor of `ContinuousDS`:
```julia
using DynamicalSystems
lor = Lorenz96(0.01, 5) # create struct
u0 = rand(5)
# pass the struct as "equations of motion":
ds = ContinuousDS(u0, lor; name="Lorenz 96 (chain of 5)")
traj = trajectory(ds, 100.0) # works!
```

Notice that this *is not necessary* if you want to incorporate only DifferentialEquations.jl and DynamicalSystems.jl, since we provide interfaces for
`ODEProblem` and `ODEIntegrator`.

!!! info "Closures and Functors"
    Julia handles [Closures](https://docs.julialang.org/en/stable/devdocs/functions/#Closures-1) like Functor objects. For our predefined systems we use closures instead of functors,
    but you can see for yourself that these 2 approaches are almost the same.

## Numerical Data
Numerical data in DynamicalSystems.jl is represented by a structure called
`Dataset`:
```@docs
Dataset
```
---
In essence a `Dataset` is simply a container for a `Vector` of `SVector`s, but only for
cases where the all inner vectors are of equal size.
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

All functions of our package that manipulate and use data are expecting a `Dataset` instance. This allows us to define efficient methods that coordinate
well with other packages, like e.g. [`neighborhood`](@ref).

If given
a matrix, we first convert to `Dataset`. This means that you should first
convert your data to a `Dataset` if you want to call functions more than once, to avoid
constantly converting.

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
