# Choosing a solver

`ContinuousDynamicalSystem`s are evolved using solvers from [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/). In this page we discuss the importance of which solver to choose.

## Default Solver
The default solver is:
```@example solver
using DynamicalSystems
DynamicalSystemsBase.DEFAULT_SOLVER
```
which is a Runge-Kutta-like solver. The number in the solver's name is the "order" of the solver.

## Speed of a solver
Estimating a given solver's performance for a particular problem is not trivial. The following are general rules of thumb:

1. Higher order solvers call the equations of motion function more times per step.
2. Higher order solvers can cover larger timespans per step.
3. Higher order solvers do better at small tolerances.

This means that there is a delicate balance between how expensive is your function and how large of a step a solver can take while it is still efficient. In general you want to strike a point of taking large steps but also not calling the function exceedingly often.

## How do I pick?
The answer to this question is easy: **benchmarks!**

Here is a simple case: let's compute the Lyapunov spectrum of the Lorenz system using [`lyapunovs`](@ref):
```@example solver
ds = Systems.lorenz()
tols = (abstol = 1e-6, reltol = 1e-6)
lyapunovs(ds, 2000; Ttr = 100.0, tols...)
```

The above uses the default solver. Let's now benchmark using two different solvers, `SimpleATsit5` and `Vern9`. Since the `SimpleATsit5` case is of lower order, naively one might think it is faster because it makes less function calls. This argument is not necessarily true though.

It is important to understand that when calling `lyapunovs(ds, 2000)` you want the system (and the tangent space) to be evolved so that it reaches a total time of `2000*dt`, which by default is `2000.0` units of time. Even though `SimpleATsit5` requires less function calls per step, `Vern9` can cover larger timespans per step.

Here are the numbers:
```@example solver
using BenchmarkTools, OrdinaryDiffEq, SimpleDiffEq
@btime lyapunovs(ds, 2000; alg = SimpleATsit5(), Ttr = 100.0, tols...);
@btime lyapunovs(ds, 2000; alg = Vern9(),        Ttr = 100.0, tols...);
```

As you can see `Vern9` is faster in doing the _entire_ computation! Of course this does not have to be universally true. It is true for the Lorenz system, but for your specific system you should do dedicated benchmarks!

## DifferentialEquations.jl

For more info about the possible solvers be sure to head over to the documentation of [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/)!
