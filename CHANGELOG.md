The package DynamicalSystems.jl does not have any actual code.
It re-exports other packages. The changelog here only lists major changes to the overarching DynamicalSystems.jl ecosystem.

# DynamicalSystems.jl v3.0

Super duper major update.

1. [Modularization](#modularization)
2. [Overhauling](#overhauling)
3. [New documentation](#new-documentation)
4. [Re-write of `DynamicalSystem`](#re-write-of-dynamicalsystem)
5. [Continuous time systems changes](#continuous-time-systems-changes)
6. [Package split: StateSpaceSets.jl](#package-split-statespacesetsjl)
7. [Package split: FractalDimensions.jl](#package-split-fractaldimensionsjl)
8. [Partial overhaul: RecurrenceAnalysis.jl](#partial-overhaul-recurrenceanalysisjl)
9. [New amazing package ComplexityMeasures.jl](#new-amazing-package-complexitymeasuresjl)
10. [New amazing package Attractors.jl](#new-amazing-package-attractorsjl)
11. [Co-maintainer(s) needed](#co-maintainers-needed)

## Modularization

The packages that compose DynamicalSystems.jl have been split up into more packages, making a more modular library, and also allowing a much smaller dependency tree for users that want only a specific part of the functionality.

## Overhauling

Practically all packages that are part of DynamicalSystems.jl have been overhauled in some way or another. Some, have been completely overhauled, like `DynamicalSystemsBase` and `ComplexityMeasures`. Some others have been partially overhauled like `RecurrenceAnalysis` and `ChaosTools`.

DynamicalSystems.jl was the first library I've written and after progressing as a software developer I have learned so much more about good design. I've put all this knowledge to good use and decided it was a good point to just do a mass-overhaul on the package.

## New documentation

This modularization also lead to a documentation overhaul: now every package of the library builds and hosts its own documentation. The main documentation of DynamicalSystems.jl now follows an approach similar to [DiffEqDocs](https://docs.sciml.ai/DiffEqDocs/stable/) which includes an overarching tutorial and overview, and then "dispatches" to documentations of other packages.

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

- The default integrator used in continuous time Ordinary Differential Equations (what was before called `ContinuousDynamicalSystem` and now called `CoupledODEs`) is now `Tsit5()`. This means that OrdinaryDiffEq.jl is a dependency of DynamicalSystemsBase.jl. It is a better option than the previous `SimpleATsit5()`. Additionally, the insane compile and first-use time improvements in latest Julia versions had made the impact of loading OrdinaryDiffEq.jl much smaller.

## Package split: StateSpaceSets.jl

It has all functionality surrounding `StateSpaceSet` (what was previously known as `Dataset`). It therefore detaches the infastructure for handing numeric data from delay embeddings, and now DelayEmbeddings.jl is a package dedicated to creating and optimizing delay coordinates embeddings.

https://juliadynamics.github.io/StateSpaceSets.jl/dev/

## Package split: FractalDimensions.jl

All functionality related to computing fractal dimensions has been detached from ChaosTools.jl into a new package: https://juliadynamics.github.io/FractalDimensions.jl/dev/

## Partial overhaul: RecurrenceAnalysis.jl

RecurrenceAnalysis.jl, now in v2.0, has had its core type (`RecurrenceMatrix`) overhauled for better design, more clarity, extendibility, and overall more Julian programming utilizing multiple dispatch and dedicated types to specify options.

Now a subtype of `AbstractRecurrenceType` is used to specify how to create a recurrence matrix. More details in the docpage of the package:

## New amazing package ComplexityMeasures.jl

This is a complete overhaul and massive enhancement of the previous Entropies.jl and deserves its own announcement post:

## New amazing package Attractors.jl

## Co-maintainer(s) needed

Some parts have co-maintainers (Attractors.jl have , complexityeasures have ..) but this whole ecosystem has become too large and I don't have the energy to maintain it alone, so I am kindly asking for people with expertise in both NLD and programming to consider coming on board and helping with answering questions, addressing issues, fixing bugs, and testing out the ecosystem.