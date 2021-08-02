# Advanced documentation
This section overviews the various integrators available from `DynamicalSystemsBase`, as well as gives some insight into the internals, so that other developers that want to use this library can build upon it.

## Integrators
```@docs
integrator
parallel_integrator
tangent_integrator
```
---
Notice that the state type `integrator.u` of each integrator is quite different and *does change* between the possible versions of a [`DynamicalSystem`](@ref)!

## Integrator state functions
There are four functions associated with the integrators that we export:
```@docs
get_state
set_state!
get_deviations
set_deviations!
```
!!! note
    These functions work with *any* possible integrator and it is best to use the to change states robustly!

## Re-initializing an integrator
It is more efficient to re-initialize an integrator using `reinit!`
than to create a new one.
This can be very helpful when looping over initial conditions and/or parameter values.

All high-level functions from `ChaosTools` have a set-up part that creates an integrator, and a low-level part that does the computation. The low level part is your friend! Use it! See the [Using GALI](@ref) page for an example as well as the section below.

The `reinit!` call signature is the same for continuous and discrete systems.
In the following, `state` is supposed to be a `D` dimensional vector (state of the dynamical system).

1. `reinit!(integ, state)` : to be used with standard [`integrator`](@ref).
3. `reinit!(integ, Vector_of_states)` : to be used with the [`parallel_integrator`](@ref).
2. `reinit!(integ, state, Q0::AbstractMatrix)` : to be used with the [`tangent_integrator`](@ref). This three argument version of `reinit!` is exported from `DynamicalSystemsBase`.

### Re-init of continuous tangent integrator
Here we compute the [`lyapunovspectrum`](@ref) for many different initial conditions.
```julia
ds = Systems.lorenz()
tinteg = tangent_integrator(ds, 2)
ics = [rand(3) for i in 1:100]
for ic in ics
  reinit!(tinteg, ic, orthonormal(3, 2))
  λ = lyapunovspectrum(tinteg, 1000, 0.1, 10.0)
  # reminder: lyapunovspectrum(tinteg, N, Δt::Real, Ttr::Real = 0.0)
end
```


### Re-init of discrete parallel integrator
Here we compute the [`lyapunov`](@ref) for many different parameters.
```julia
ds = Systems.henon()
u0 = rand(SVector{2})
ps = 1.2:0.01:1.4
pinteg = parallel_integrator(ds, [u0, u0 + 1e-9rand(SVector{2})])
for p in ps
  set_parameter!(ds, 1, p)
  reinit!(pinteg, [u0, u0 + 1e-9rand(SVector{2})])
  λ = lyapunov(pinteg, 1000, 10, 1, 1e-9, 1e-6, 1e-12)
  # reminder: lyapunov(pinteg, T, Ttr, Δt, d0, ut, lt)
end
```

## Using callbacks with integrators
For the case of continuous systems you can add callbacks from the event handling of **DifferentialEquations.jl**. This is done simply as a keyword argument to the initializers.

In this example we use a simple `SavingCallback` to save the distance between the two states of a [`parallel_integrator`](@ref).

```@example MAIN
using DynamicalSystems, DiffEqCallbacks
using LinearAlgebra: norm

kwargs = (abstol=1e-14, reltol=1e-14, maxiters=1e9)
ds = Systems.lorenz()
d0 = 1e-9
T = 100.0

save_func(u, t, integrator) = norm(u[1] - u[2])
saved_values = SavedValues(eltype(ds.t0), eltype(get_state(ds)))
cb = SavingCallback(save_func, saved_values)

u0 = get_state(ds)
pinteg = parallel_integrator(ds, [u0, u0 + rand(SVector{3})*d0*√3];
kwargs..., callback = cb)
step!(pinteg, T)
t = saved_values.t
n = saved_values.saveval
```
As expected you can see that the recorded distance between two states is increasing.

## `DynamicalSystem` implementation
```julia
abstract type DynamicalSystem{
        IIP,     # is in place , for dispatch purposes and clarity
        S,       # state type
        D,       # dimension
        F,       # equations of motion
        P,       # parameters
        JAC,     # jacobian
        JM,      # jacobian matrix
        IAD}     # is auto-differentiated
    # one-liner: {IIP, S, D, F, P, JAC, JM, IAD}
    # Subtypes of DynamicalSystem have fields:
    # 1. f
    # 2. u0
    # 3. p
    # 4. t0
    # 5. jacobian (function)
    # 6. J (matrix)
end
```
The `DynamicalSystem` stores only the absolutely necessary information. Every other functionality of **DynamicalSystems.jl** initializes an integrator.

The final type-parameter `IAD` is useful when creating the `tangent_integrator`, so that the vector field is not computed twice!

## Choosing a solver

`ContinuousDynamicalSystem`s are evolved using solvers from [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/). In this page we discuss the importance of which solver to choose.

### Default Solver
The default solver is:
```@example MAIN
using DynamicalSystems
DynamicalSystemsBase.DEFAULT_SOLVER
```
which is a Runge-Kutta-like solver. The number in the solver's name is the "order" of the solver.

### Speed of a solver
Estimating a given solver's performance for a particular problem is not trivial. The following are general rules of thumb:

1. Higher order solvers call the equations of motion function more times per step.
2. Higher order solvers can cover larger timespans per step.
3. Higher order solvers do better at small tolerances.

This means that there is a delicate balance between how expensive is your function and how large of a step a solver can take while it is still efficient. In general you want to strike a point of taking large steps but also not calling the function exceedingly often.

### How do I pick?
The answer to this question is easy: **benchmarks!**

Here is a simple case: let's compute the Lyapunov spectrum of the Lorenz system using [`lyapunovspectrum`](@ref):
```@example MAIN
ds = Systems.lorenz()
tols = (abstol = 1e-6, reltol = 1e-6)
lyapunovspectrum(ds, 2000; Ttr = 100.0, tols...)
```

The above uses the default solver. Let's now benchmark using two different solvers, `SimpleATsit5` and `Vern9`. Since the `SimpleATsit5` case is of lower order, naively one might think it is faster because it makes less function calls. This argument is not necessarily true though.

It is important to understand that when calling `lyapunovspectrum(ds, 2000)` you want the system (and the tangent space) to be evolved so that it reaches a total time of `2000*Δt`, which by default is `2000.0` units of time. Even though `SimpleATsit5` requires less function calls per step, `Vern9` can cover larger timespans per step.

Here are the numbers:
```julia
using BenchmarkTools, OrdinaryDiffEq, SimpleDiffEq, Statistics
b1 = @benchmark lyapunovspectrum(ds, 2000; alg = SimpleATsit5(), Ttr = 100.0, tols...);
b2 = @benchmark lyapunovspectrum(ds, 2000; alg = Vern9(),        Ttr = 100.0, tols...);
println("Timing for SimpleATsit5:")
println(mean(b1))
println("Timing for Vern9:")
println(mean(b2))
```

```
Timing for SimpleATsit5:
TrialEstimate(53.517 ms)
Timing for Vern9:
TrialEstimate(27.511 ms)
```

As you can see `Vern9` is faster in doing the _entire_ computation! Of course this does not have to be universally true. It is true for the Lorenz system, but for your specific system you should do dedicated benchmarks!

### DifferentialEquations.jl

For more info about the possible solvers be sure to head over to the documentation of [DifferentialEquations.jl](https://diffeq.sciml.ai/latest/)!
