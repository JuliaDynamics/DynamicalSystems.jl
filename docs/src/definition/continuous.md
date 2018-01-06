# Continuous Systems
Continuous systems of the form
```math
\frac{d\vec{u}}{dt} = \vec{f}(t, \vec{u}),
```
are defined using the `ContinuousDS` structure:
```@docs
ContinuousDS
```
---
You can use any function that complies with the requirements stated by the documentation string. However, it is highly advised to use [Functors](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1) for dynamical systems where the equations
of motion contain parameters.

In the following examples we will demonstrate how one can use both constructors.

## Defining a `DynamicalSystem` using Functors
A Functor is a shorthand for saying [Function-like objects](https://docs.julialang.org/en/stable/manual/methods/#Function-like-objects-1),
i.e. `struct`s that are also callable (see the linked documentation page). Using
such objects one can create both the equations of motion and a parameter container
under a single `struct` definition.

Here we will use the constructor
```julia
ContinuousDS(state, eom! [, jacob! [, J]]; tspan = (0.0, 100.0))
```
and create the continuous Rössler
system, from our [Predefined Systems](predefined):
```julia
using DynamicalSystems
mutable struct Rössler
    a::Float64
    b::Float64
    c::Float64
end
@inline @inbounds function (s::Rössler)(t, u::AbstractVector, du::AbstractVector)
    du[1] = -u[2]-u[3]
    du[2] = u[1] + s.a*u[2]
    du[3] = s.b + u[3]*(u[1] - s.c)
    return nothing
end
@inline @inbounds function (s::Rössler)(t, u::AbstractVector, J::AbstractMatrix)
    J[2,2] = s.a
    J[3,1] = u[3]; J[3,3] = u[1] - s.c
    return nothing
end
```
The first code-block defines a `struct` that is simply a container for the
parameters of the Rössler system. The second code-block defines the equations
of motion of the system, by taking advantage of the fact that you can *call*
this `struct` as if it was a function:
```julia
s = Rössler(1,2,3)
u = rand(3); du = copy(u)
s(0, u, du)
```

The third code-block then defines the Jacobian function using multiple dispatch.
This allows us to use the *same* instance of `Rössler` for *both* the equations
of motion *and* the Jacobian function!

!!! important "Use `AbstractVector` and `AbstractMatrix`"
    You must define your equations of motion / Jacobian functions using `Abstract`
    Types, and **not** types like `Vector` or `Matrix`, otherwise functions like
    `lyapunovs` won't work properly.

The possibility of providing an initialized
Jacobian to the `ContinuousDS` constructor allows us to "cheat".
Notice that the Jacobian function only accesses
fields that depend on the parameters and/or the state variables, because the other
fields are constants and will be initialized properly later.

Next, we define a "set-up" function, that returns a `ContinuousDS`:
```julia
function roessler(u0=rand(3); a = 0.2, b = 0.2, c = 5.7)
    i = one(eltype(u0))
    o = zero(eltype(u0))
    J = zeros(eltype(u0), 3, 3)
    J[1,:] .= [o, -i,      -i]
    J[2,:] .= [i,  a,       o]
    J[3,:] .= [u0[3], o, u0[1] - c]

    s = Rössler(a, b, c)
    return ContinuousDS(u0, s, s, J)
end

ds = roessler()
# Equivalent with our predefined system:
ds = Systems.roessler()
```
```
3-dimensional continuous dynamical system:
state: [0.021655, 0.530449, 0.0227049]
e.o.m.: DynamicalSystemsBase.Systems.Rössler(0.2, 0.2, 5.7)
```
Then, it is trivial to change a parameter of the system by e.g. doing
`ds.prob.f.c = 2.2`. The equations of motion for a continuous system are stored in the `ODEProblem` struct, the field `f`.

Notice that this parameter change will affect both the equations of motion as well
as the Jacobian function, making everything concise and easy-to-use!

## Using `ODEProblem` to define a `ContinuousDS`
Here we will show how one can take advantage of the callback capabilities of [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/) to define
a system.

!!! danger "Callbacks do not propagate in variation vector methods!"
    Methods that evolve variation vectors in time (currenlty [`gali`](@ref) and
    [`lyapunovs`](@ref)) do not inherit callbacks present in the definition
    of a `ContinuousDS`.

We will make a Hénon–Heiles that also satisfies energy conservation.
We first write the equations of motion and the Jacobian functions in the instructed form:
```julia
function hheom!(t, u::AbstractVector, du::AbstractVector)
    du[1] = u[3]
    du[2] = u[4]
    du[3] = -u[1] - 2u[1]*u[2]
    du[4] = -u[2] - (u[1]^2 - u[2]^2)
    return nothing
end
function hhjacob!(t, u::AbstractVector, J::AbstractMatrix)
    J[3,1] = -1 - 2u[2]; J[3,2] = -2u[1]
    J[4,1] = -2u[1]; J[4,2] =  -1 + 2u[2]
    return nothing
end
```
The Jacobian matrix will be initialized properly later. Now, we are going to use a
`Callback` to conserve energy. First, define the energy functions
```julia
@inline V(q1, q2) = 1//2 * (q1^2 + q2^2 + 2q1^2 * q2 - 2//3 * q2^3)
@inline T(p1, p2) = 1//2 * (p1^2 + p2^2)
@inline H(q1, q2, p1, p2) = T(p1, p2) + V(q1, q2)
@inline H(u::AbstractVector) = H(u...)
```

Then, create a "residual" function, used in the [`ManifoldProjection`](http://docs.juliadiffeq.org/latest/features/callback_library.html#Manifold-Conservation-and-Projection-1) callback:

```julia
u0 = [0.1, 0, 0, 0.5]
const E = H(u0[1],u0[2],u0[3],u0[4])

function g(u, resid)
    resid[1] = H(u[1],u[2],u[3],u[4]) - E
    resid[2:4] .= 0
end
```

Next we create the [`Callback`](http://docs.juliadiffeq.org/latest/features/callback_functions.html),
the [`ODEProblem`](http://docs.juliadiffeq.org/latest/types/ode_types.html) and then
dynamical system structure, `ContinuousDS`:

```julia
# Pkg.add("DiffEqCallbacks")
using DiffEqCallbacks, OrdinaryDiffEq

cb = ManifoldProjection(g, nlopts=Dict(:ftol=>1e-13), save = false)
prob = ODEProblem(hheom!, u0, (0., 100.0),  callback=cb)

# Initialize Jacobian
o = 0.0; i = 1.0; J = zeros(4,4)
J[1,:] = [o,    o,     i,    o]
J[2,:] = [o,    o,     o,    i]
J[3,:] = [ -i - 2u0[2],   -2u0[1],   o,   o]
J[4,:] = [-2u0[1],  -1 + 2u0[2],  o,   o]

ds = ContinuousDS(prob, hhjacob!, J)
```

Notice that using the argument `save = false` in the `ManifoldProjection` is crucial, because otherwise any data taken from the system,
using e.g. [`trajectory`](@ref) will necessarily have saved points at every
callback realization (which you *do not* want if you want timeseries of equi-sampled
points).

Let's see now if our system does indeed conserve energy!
```julia
a1 = trajectory(ds, 1000.0)
a2 = trajectory(ds, 1000.0, diff_eq_kwargs = Dict(:solver => Vern9(),
:abstol => 1e-9, :reltol => 1e-9))

energies1 = [H(p) for p in a1]
energies2 = [H(p) for p in a2]

println("Default solver: ΔE = ", std(energies1))
println("High-Accuracy solver: ΔΕ = ", std(energies2))
```
```
Default solver: ΔE = 2.7868013699842223e-5
High-Accuracy solver: ΔΕ = 1.2345950190344736e-12
```
By combining a high-accuracy solver with a callback one can get an incredibly low
energy error.
