# Changelog for DynamicalSystems

The package DynamicalSystems.jl has little actual code.
It mainly re-exports other packages. The source code contained here is only regarding graphical interactive applications for dynamical systems.

The changelog here therefore lists either major changes to the overarching DynamicalSystems.jl ecosystem, or major changes to plotting infrastructure.

The changelogs of individual sub-packages are self-contained for each package.


# v3.5

- Improved documentation regarding versions
- Added SignalDecomposition.jl into DynamicalSystems.jl

# v3.4

Some important fixes on the `interactive_trajectory` GUI:

- Changed the internal handling for continuous time systems. Now they are stepped for exactly `Δt` by giving `true` as third input to `step!`. This increases consistency with discrete systems by not altering the integration protocol.
- The documentation around the "run" button of the GUI has been clarified.

# v3.3

DynamicalSystems.jl now integrates with ModelingToolkit.jl and allows using symbolic variables to access/observe state and parameter.

At a low level, this happens via the functions `observe_state`, `set_state!`,
`current_parameter` and `set_parameter!`.

Additionally, `interactive_trajectory_timeseries` allows symbolic indexing
for state space plot, timeseries plots, or parameter sliders.
Everything is also automatically named and limits are also automatically deduced for everything! Super convenient!

# v3.2

- Significant improvement of `interactive_trajectory` by allowing the states or parameters of the system to be reset to their original values!

# v3.1

DynamicalSystems.jl moved to Julia 1.9 and also ported all code related to interactive visualizations here (utilizing the new package extensions system).
Now the `InteractiveDynamics` package is obsolete.

# DynamicalSystems.jl v3.0

Super duper major update.

1. [Overhauling](#overhauling)
2. [Modularization](#modularization)
3. [New documentation](#new-documentation)
4. [Notable breaking changes](#notable-breaking-changes)
   1. [Renames](#renames)
5. [Re-write of `DynamicalSystem`](#re-write-of-dynamicalsystem)
6. [ODE solver change](#ode-solver-change)
7. [Package split: StateSpaceSets.jl](#package-split-statespacesetsjl)
8. [Package split: FractalDimensions.jl](#package-split-fractaldimensionsjl)
9. [Partial overhaul: RecurrenceAnalysis.jl](#partial-overhaul-recurrenceanalysisjl)
10. [New amazing package ComplexityMeasures.jl](#new-amazing-package-complexitymeasuresjl)
11. [New amazing package Attractors.jl](#new-amazing-package-attractorsjl)
12. [Co-maintainer(s) needed](#co-maintainers-needed)

## Overhauling

Practically all packages that are part of DynamicalSystems.jl have been overhauled in some way or another. Some have been completely overhauled, like `DynamicalSystemsBase` and `ComplexityMeasures`. Some others have been partially overhauled like `RecurrenceAnalysis` and `ChaosTools`.

DynamicalSystems.jl was the first library I've written and after progressing as a software developer I have learned so much more about good design. I've put all this knowledge to good use and decided it was a good point to just do a mass-overhaul on the package.

## Modularization

The packages that compose DynamicalSystems.jl have been split up into more packages, making a more modular library, and also allowing a much smaller dependency tree for users that want only a specific part of the functionality.
The modularization itself did not bring any breaking changes.

## New documentation

This modularization also lead to a documentation overhaul: now every package of the library builds and hosts its own documentation. The DynamicalSystems.jl brings everything together with an overarching tutorial and points to the individual docs.

The main documentation of DynamicalSystems.jl now includes an overarching tutorial and overview, and then "dispatches" to documentations of other packages.
In the near future the documentation will be in the style similar to [DiffEqDocs](https://docs.sciml.ai/DiffEqDocs/stable/), however, I am still exploring alternative and stable ways to achieve similar results.

## Notable breaking changes

- The `trajectory` function now returns both the state space set and the time vector. In sort, all code that did `A = trajectory(...)` is now broken and should be replaced by `A, tvec = trajectory(...)`.
- There are no longer `diffeq` keyword arguments to downstream functions. All functions such as `lyapunovspectrum, basins_of_attraction`, etc. no longer accept the `diffeq` keyword, that used to specify arguments for DifferentialEquations.jl. Instead, the DiffEq solver options are decided during the creation of `CoupledODEs` (what was before `ContinuousDynamicalSystem`).
- The basic dynamical system constructors do not accept a Jacobian. Giving a fourth argument will error. The dedicated `TangentDynamicalSystem` should be used for this.
- The special type for 1D discrete dynamical systems has been removed. 1D systems are now represented by a 1D static vector that needs to be reshaped accordingly.

### Renames

Some things have been renamed to have a clearer and more specific name. These aren't breaking changes, only deprecations. From the top of my head, the most important ones are:

- Dataset -> StateSpaceSet
- DiscreteDynamicalSystem -> DeterministicIteratedMap
- ContinuousDynamicalSystem -> CoupledODEs

## Re-write of `DynamicalSystem`

The core library DynamicalSystemsBase.jl has had a complete re-write. This has the following major impacts:

1. Concept of what is a "dynamical system" is much more general. now `DynamicalSystem` defines a proper, extendable interface. The set of functions that the interface satisfies are listed in `DynamicalSystem` documentation string.
2. More variants of a dynamical system are offered. Any arbitrary steppable thing, like an agent based model, can become a dynamical system.
3. Overall much better and rigorous tests.
3. Now **_all concrete implementations of `DynamicalSystem` can be iteratively evolved in time via the `step!` function_**. For example, `ContinuousDynamicalSystem` is now a light wrapper around `ODEIntegrator`. After working with this library for years I realized that the "middle man" of representing a dynamical system (like an ode problem) was not really necessary in downstream functions, since all of them exclusively use an integrator.
4. Dynamical systems are completely detached from a Jacobian creation, which now is the exclusive task of `TangentDynamicalSystem`. As a result all predefined dynamical systems do not include a hand-coded Jacobian; rather the Jacobian needs to be given from `Systems` submodule to the `TangentDynamicalSystem`.

Also, the `Systems` submodule has been removed from the DynamicalSystemsBase.jl package and became a new package PredefinedDynamicalSystems.jl. It is untested and it is not recommended to use in any other setting other a demonstrational example. (this was always the case, it's just that so far we weren't transparent that nothing there was tested...)

Read more here: https://juliadynamics.github.io/DynamicalSystemsBase.jl/dev/

## ODE solver change

The default integrator used in continuous time Ordinary Differential Equations (what was before called `ContinuousDynamicalSystem` and now called `CoupledODEs`) is now `Tsit5()`. This means that OrdinaryDiffEq.jl is a dependency of DynamicalSystemsBase.jl. It is a better option than the previous `SimpleATsit5()`. Additionally, the insane compile and first-use time improvements in latest Julia versions had made the impact of loading OrdinaryDiffEq.jl much smaller.

## Package split: StateSpaceSets.jl

It has all functionality surrounding `StateSpaceSet` (what was previously known as `Dataset`). It therefore detaches the infrastructure for handing numeric data from delay embeddings, and now DelayEmbeddings.jl is a package dedicated to creating and optimizing delay coordinates embeddings.

https://juliadynamics.github.io/StateSpaceSets.jl/dev/

## Package split: FractalDimensions.jl

All functionality related to computing fractal dimensions has been detached from ChaosTools.jl into a new package: https://juliadynamics.github.io/FractalDimensions.jl/dev/

## Partial overhaul: RecurrenceAnalysis.jl

RecurrenceAnalysis.jl, now in v2.0, has had its core type (`RecurrenceMatrix`) overhauled for better design, more clarity, extendibility, and overall more Julian programming utilizing multiple dispatch and dedicated types to specify options.

Now a subtype of `AbstractRecurrenceType` is used to specify how to create a recurrence matrix. More details in the docpage of the package:

https://juliadynamics.github.io/RecurrenceAnalysis.jl/dev/

## New amazing package ComplexityMeasures.jl

This is a complete overhaul and massive enhancement of the previous Entropies.jl and deserves its own announcement post: https://discourse.julialang.org/t/complexitymeasures-jl/95198

## New amazing package Attractors.jl

Incredibly exciting stuff that we (@awage @kalelr @Datseris) have been working on for the last year. It features an interface for finding attractors and their basin fractions, _as well as a new approach to "continuation"_. We are currently writing a paper for this, so the package doesn't have an announcement yet. Once we have the preprint though, we will write an announcement post for the package and share the pre-print!

For now, feel free to see the docs and give us feedback!!!

https://juliadynamics.github.io/Attractors.jl/dev/

## Co-maintainer(s) needed

Some parts of the library have co-maintainers (Attractors.jl has @awage and @kalelr, ComplexityMeasures.jl and TimeseriesSurrogates.jl have @kahaaga) but this whole ecosystem around DynamicalSystems.jl has become too large and I think I just don't have the time to maintain it alone anymore. So, I am kindly asking for people with expertise in both NLD and programming to consider coming on board and helping with answering questions, addressing issues, fixing bugs, and testing out the ecosystem. Feel free to contact me via a PM on Slack or here on Discourse!