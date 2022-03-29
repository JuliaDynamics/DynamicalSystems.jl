# Integrators
When writing algorithms for nonlinear dynamics in code, it is almost always more convenient to work with integrators. DynamicalSystems.jl provides several different integrators of dynamical systems and a unified API to treat all of them.
The API is based on the same API that DifferentialEquations.jl uses.

## Available integrators

* [`integrator`](@ref)
* [`parallel_integrator`](@ref)
* [`tangent_integrator`](@ref)
* [`projected_integrator`](@ref)
* [`stroboscopicmap`](@ref)
* [`poincaremap`](@ref)

```@docs
integrator
parallel_integrator
tangent_integrator
projected_integrator
stroboscopicmap
poincaremap
GeneralizedDynamicalSystem
```

## Integrator API
After you have initialized any integrator, you can use the same functions to handle them. The most important function is [`step!`](@ref), that will simply progress the integrator. [`get_state`](@ref) will return its current state. [`reinit!`](@ref) can be used to re-start the integrator at a possibly different new state.

Especially for the [`tangent_integrator`](@ref), there are also two more functions: [`get_deviations`](@ref), [`set_deviations`](@ref).

```@docs
step!(ds::DynamicalSystem)
get_state
reinit!(ds::DynamicalSystem)
get_deviations
set_deviations!
```


### Re-initializing an integrator
It is more efficient to re-initialize an integrator using `reinit!`
than to create a new one.
This can be very helpful when looping over initial conditions and/or parameter values.

All high-level functions from `ChaosTools` have a set-up part that creates an integrator, and a low-level part that does the computation. The low level part is your friend! Use it! See the [Using GALI](@ref) page for an example as well as the section below.

#### Example: Re-init of continuous tangent integrator
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


## Using callbacks with integrators
For the case of continuous systems you can add callbacks from the event handling of **DifferentialEquations.jl**. This is done simply as a keyword argument to the initializers.

In this example we use a simple `SavingCallback` to save the distance between the two states of a [`parallel_integrator`](@ref).

```@example MAIN
using DynamicalSystems, DiffEqCallbacks
using LinearAlgebra: norm

ds = Systems.lorenz()
d0 = 1e-9
T = 100.0

save_func(u, t, integrator) = norm(u[1] - u[2])
saved_values = SavedValues(eltype(ds.t0), eltype(get_state(ds)))
cb = SavingCallback(save_func, saved_values)

diffeq = (abstol=1e-14, reltol=1e-14, maxiters=1e9, callback = cb)

u0 = get_state(ds)
pinteg = parallel_integrator(ds, [u0, u0 + rand(SVector{3})*d0*√3]; diffeq)
step!(pinteg, T)
t = saved_values.t
n = saved_values.saveval
```
As expected you can see that the recorded distance between two states is increasing.

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
lyapunovspectrum(ds, 2000; Ttr = 100.0, diffeq = tols)
```

The above uses the default solver. Let's now benchmark using two different solvers, `SimpleATsit5` and `Vern9`. Since the `SimpleATsit5` case is of lower order, naively one might think it is faster because it makes less function calls. This argument is not necessarily true though.

It is important to understand that when calling `lyapunovspectrum(ds, 2000)` you want the system (and the tangent space) to be evolved so that it reaches a total time of `2000*Δt`, which by default is `2000.0` units of time. Even though `SimpleATsit5` requires less function calls per step, `Vern9` can cover larger timespans per step.

Here are the numbers:
```julia
using BenchmarkTools, OrdinaryDiffEq, SimpleDiffEq, Statistics
b1 = @benchmark lyapunovspectrum(ds, 2000; diffeq = (alg = SimpleATsit5(), tols...), Ttr = 100.0);
b2 = @benchmark lyapunovspectrum(ds, 2000; diffeq = (alg = Vern9(), tols...), Ttr = 100.0);
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
