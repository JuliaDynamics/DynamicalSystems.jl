# Continuous Systems
Continuous systems of the form
```math
\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}),
```
are defined using the `ContinuousDS` structure:
```@docs
ContinuousDS
```
---
Notice that the fields `eom!` and `jacob!` end with a `!`, to remind users
that these functions should operate in-place. Also notice that the type `ContinuousDS`
is actually immutable.

You can use any function that complies with the requirements stated by the documentation string.
However, it is highly advised to use [Functors](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1) for dynamical systems where the equations
of motion contain parameters.

## Defining a `DynamicalSystem` using Functors
A Functor is a shorthand for saying [Function-like objects](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1),
i.e. `struct`s that are also callable (see the linked documentation page). Using
such objects one can create both the equations of motion and a parameter container
under a single `struct` definition.

For example, let's take a look at the source code that generates continuous Rössler
system, from the [Predefined Systems](definition/predefined):
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
    equations of motion as used in **DynamicalSystems.jl** (namely `Vector`, `SubArray`
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

An example of a system definition without Functors is [also shown here](definition/discrete#defining-a-dynamicalsystem-without-functors).
