# Time Evolution of Systems
**DynamicalSystems.jl** provides convenient interfaces for the evolution of systems.
```@docs
evolve
evolve!
trajectory
get_sol
```
---
as well as some basic interfacing for [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/index.html):
```@docs
ODEProblem
```
---
Notice that if you want to do repeated evolutions of different states of a
continuous system, you should use the
`DiffEqBase.ODEIntegrator` in conjunction with `DiffEqBase.reinit!(integrator, newstate)` to avoid the intermediate initializations of the integrator each time.

## Solution precision for continuous systems
A numerical solution of an ODE is not the "true" solution, uniquely defined by a (well-defined) ODE and an initial condition. Especially for chaotic systems, where deviations are amplified exponentially, one is left worried if the numerical solutions truly are part of the system and can truly give insight in understanding the system.

DifferentialEquations.jl offers a tool, called [Uncertainty Quantification](http://docs.juliadiffeq.org/latest/analysis/uncertainty_quantification.html),
which allows users to asses up to what time-scales the numerical solution is close
to the "true" solution. For example, using the default solving parameters of
**DynamicalSystems.jl**, the Lorenz system is accurate up to time `t = 50.0`.

However, fortunately for us, there is not too much worry about the numerical solution diverging from the true solution. That is because of the [shadowing theorem](http://mathworld.wolfram.com/ShadowingTheorem.html) (or
[shadowing lemma](http://www.scholarpedia.org/article/Shadowing_lemma_for_flows)):

!!! quote "Shadowing Theorem"
    Although a numerically computed chaotic trajectory diverges exponentially from the true trajectory with the same initial coordinates, there exists an errorless trajectory with a slightly different initial condition that stays near ("shadows") the numerically computed one.

This simply means that one can always numerically study chaos not only qualitatively but also quantitatively. For more information, see the book *Chaos in Dynamical Systems* by E. Ott, or the
[scholarpedia](http://www.scholarpedia.org/article/Shadowing_lemma_for_flows) entry.
