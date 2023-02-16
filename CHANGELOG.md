The package DynamicalSystems.jl does not have any actual code.
It re-exports other packages. The changelog here only lists major changes to the overarching DynamicalSystems.jl ecosystem.

# DynamicalSystems.jl v3.0

Super duper major update.

1. [Modularization](#modularization)
2. [New documentation](#new-documentation)
3. [Re-write of `DynamicalSystem`](#re-write-of-dynamicalsystem)
4. [Continuous time systems changes](#continuous-time-systems-changes)
5. [New package StateSpaceSets.jl](#new-package-statespacesetsjl)
6. [New package FractalDimensions.jl](#new-package-fractaldimensionsjl)
7. [New amazing package ComplexityMeasures.jl](#new-amazing-package-complexitymeasuresjl)

## Modularization

The packages that compose DynamicalSystems.jl have been split up into more packages, making a more modular library, and also allowing a much smaller dependency tree for users that want only a specific part of the functionality.

## New documentation

This modularization also lead to a documentation overhaul: now every package of the library builds and hosts its own documentation. The main documentation of DynamicalSystems.jl now follows an approach similar to [DiffEqDocs](https://docs.sciml.ai/DiffEqDocs/stable/) which "dispatches" to documentations of other packages.

## Re-write of `DynamicalSystem`

The inner library DynamicalSystemsBase.jl has had a complete re-write. This has the following major impacts:

1. Concept of what is a "dynamical system" is much more general. now `DynamicalSystem` defines a proper, extendable interface. The set of functions that the interface satisfies are listed in `DynamicalSystem` documentation string.

2. More variants of a dynamical system are offered. Any arbitrary steppable thing, like an agent based model, can become a dynamical system.

3. Now **_all concrete implementations of `DynamicalSystem` can be iteratively evolved in time via the `step!` function_**. For example, `ContinuousDynamicalSystem` is now a light wrapper around `ODEIntegrator`. After working with this library for years I realized that the "middle man" of representing a dynamical system or an ode problem was not really necessary in downstream functions, since all of them exclusively use an integrator.

4. Dynamical systems are completely detached from a Jacobian creation, which now is the exclusive task of `TangentDynamicalSystem`. As a result all predefined dynamical systems do not include a hand-coded Jacobian; rather the Jacobian needs to be given from `Systems` submodule to the `TangentDynamicalSystem`.

5. There is no more a special type for 1D discrete time systems.
It led to a really large amount of lines of code that were simply not worth
saving the user the effort of transforming a 1-dimensional static vector to a number. `reinterpret` does it for free as well.


## Continuous time systems changes

- **Majorly breaking: There are no longer `diffeq` keyword arguments to downstream functions**. All functions such as `lyapunovspectrum, basins_of_attraction`, etc. no longer accept the `diffeq` keyword, that used to specify arguments for DifferentialEquations.jl. Instead, **the DiffEq solver options are decided during the creation of `CoupledODEs`** (what was before `ContinuousDynamicalSystem`).

- The default integrator used in continuous time Ordinary Differential Equations (what was before called `ContinuousDynamicalSystem` and now called `CoupledODEs`) is now `Tsit5()`. This means that OrdinaryDiffEq.jl is a dependency of DynamicalSystemsBase.jl. It is a better option than the previous `SimpleATsit5()`. Additional, the insane compile and first-use time improvements in latest Julia versions had made the impact of loading OrdinaryDiffEq.jl much smaller.

## New package StateSpaceSets.jl

It has all functionality surrounding `StateSpaceSet` (what was previously known as `Dataset`). It therefore detaches the infastructure for handing numeric data from delay embeddings, and now DelayEmbeddings.jl is a package dedicated to creating and optimizing delay coordinates embeddings.

https://juliadynamics.github.io/StateSpaceSets.jl/dev/

## New package FractalDimensions.jl

All functionality related to computing fractal dimensions has been detached from ChaosTools.jl into a new package: https://juliadynamics.github.io/FractalDimensions.jl/dev/

## New amazing package ComplexityMeasures.jl

This is a complete overhaul and massive enhancement of the previous Entropies.jl and deserves its own announcement post: