# Time Evolution of Systems
**DynamicalSystems.jl** provides convenient interfaces for the evolution of systems.
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
---
Notice that if you want to do repeated evolutions of different states of a
continuous system,
you should use the
`ODEIntegrator(ds::DynamicalSystem)` in conjunction with `DifferentialEquations.reinit!(integrator, newstate)` to avoid the intermediate initializations of the integrator each time.
