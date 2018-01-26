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

In the following examples we will demonstrate how one can use both constructors.

## Defining a `ContinuousDS` using `eom!`
Here we will use the constructor
```julia
ContinuousDS(state, eom! [, jacob! [, J]]; tspan = (0.0, 100.0))
```
and create the continuous Rössler
system, from our [Predefined Systems](predefined):
```julia
using DynamicalSystems

@inline @inbounds function roessler_eom(du, u, p, t)
    a, b, c = p
    du[1] = -u[2]-u[3]
    du[2] = u[1] + a*u[2]
    du[3] = b + u[3]*(u[1] - c)
    return nothing
end

@inline @inbounds function roessler_jacob(J, u, p, t)
    J[2,2] = p[1]
    J[3,1] = u[3]; J[3,3] = u[1] - p[3]
    return nothing
end
```
The first code-block defines defines the equations
of motion of the system. The second code-block then defines
the Jacobian function of the system.

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

    return ContinuousDS(u0, roessler_eom, roessler_jacob, J; parameters = [a, b, c])
end

ds = roessler()
# Equivalent with our predefined system:
ds = Systems.roessler()
```
```
3-dimensional continuous dynamical system:
state: [0.021655, 0.530449, 0.0227049]
e.o.m.: DynamicalSystemsBase.Systems.roessler_eom
```
Then, it is trivial to change a parameter of the system by e.g. doing
`ds.prob.p[3] = 2.2`.
Notice that this parameter change will affect both the equations of motion as well
as the Jacobian function, making everything concise and easy-to-use!



## Defining a `ContinuousDS` using `ODEProblem`
Here we will show how one can take advantage of the callback capabilities of [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/) to define
a system.

!!! danger "Callbacks do not propagate in variation vector methods!"
    Methods that evolve variation vectors in time (currenlty [`gali`](@ref) and
    [`lyapunovs`](@ref)) do not inherit callbacks present in the definition
    of a `ContinuousDS`. The issue that keeps track of this is [here](https://github.com/JuliaDynamics/DynamicalSystemsBase.jl/issues/8).

We will make a Hénon–Heiles that also satisfies energy conservation. This is also available in the [predefined systems](predefined), but we will use it here as an example.

We first write the equations of motion and the Jacobian functions in the instructed form:
```julia
function hheom!(du, u, p, t)
    du[1] = u[3]
    du[2] = u[4]
    du[3] = -u[1] - 2u[1]*u[2]
    du[4] = -u[2] - (u[1]^2 - u[2]^2)
    return nothing
end
function hhjacob!(J, u, p, t)
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

function g!(resid, u)
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

cb = ManifoldProjection(g!, nlopts=Dict(:ftol=>1e-13), save = false)
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
callback realization (which you *do not* want if you want timeseries sampled at
regular intervals, which is also the whole purpose of [`trajectory`](@ref)).

Let's see now if our system does indeed conserve energy!
```julia
a1 = trajectory(ds, 1000.0)
energies1 = [H(p) for p in a1]
maxer = maximum(@. energies1 - E)
println("Default accuracy: max(ΔE) = $maxer")
```
```
Default accuracy: max(ΔE) = 9.926393040871062e-12
```
Reminder: by default **DynamicalSystems.jl** uses the solver `Vern9()` and
error tolerances of `1e-9`.
