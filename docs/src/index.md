```@docs
DynamicalSystems
```

!!! tip "Star us on GitHub!"
    If you have found this library useful, please consider starring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl).
    This gives us an accurate lower bound of the (satisfied) user count.

## Introduction

Welcome to the documentation of **DynamicalSystems.jl**!

- If you have not used the library before, and would like to get started, then please read the [overarching tutorial](@ref tutorial) for the library.
- The [contents](@ref contents) page gives a summary of all packages that are part of the library.
- See the [learning resources](@ref learning) below to find out more resources about learning the library and using it in scientific research and/or education.
- Besides the formal algorithmic/scientific content of **DynamicalSystems.jl** (those in the [contents](@ref contents)) page, the library also provides basic functionality for interactive or offline animations and visualizations. These are found in the [visualizations](@ref visualization) page.
- The remaining of this introduction page discusses our goals with the library, how to participate as a user or developer, how to cite, and other relevant information (see the sections of the sidebar on the left).


## Latest news

All subpackages composing DynamicalSystems.jl have their own human-readable CHANGELOG files in their corresponding GitHub repositories. The CHANGELOGs log all meaningful changes in the software.

Notable news for the **DynamicalSystems.jl** library are also posted on the official Julia language Discourse, and users may subscribe to this particular Topic to get notified of updates:

https://discourse.julialang.org/t/dynamicalsystems-jl-news-updates-and-announcements/122079

## [Learning resources](@id learning)

### Textbook with DynamicalSystems.jl

We have written an undergraduate level textbook as an introduction to nonlinear dynamics. The text is written in an applied, hands-on manner, while still covering all fundamentals. The book pages are interlaced with real Julia code that uses DynamicalSystems.jl and is published in the Undergraduate Lecture Notes in Physics by Springer Nature:
* [Nonlinear Dynamics: A concise introduction interlaced with code](https://link.springer.com/book/10.1007/978-3-030-91032-7) by G. Datseris & U. Parlitz.


Additional textbooks on nonlinear dynamics with practical focus are:
* Chaos in Dynamical Systems - E. Ott
* Nonlinear Time series Analysis - H. Kantz & T. Schreiber
* Nonlinear Dynamics and Chaos - S. Strogatz

### Course on applied nonlinear dynamics and complex systems

We are developing a full course (targeting a graduate or undergraduate semester long course) on applied nonlinear dynamics, nonlinear timeseries analysis, and complex systems, using the packages of [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/). **DynamicalSystems.jl** is part of this course.

The materials of the course are on GitHub: <https://github.com/JuliaDynamics/NonlinearDynamicsComplexSystemsCourses>


## How to cite

If using **DynamicalSystems.jl** resulted in a publication with references, we kindly ask that you give appropriate credit in three ways:

1. Cite the whole **DyamicalSystems.jl** ecosystem by citing the associated DOI `10.21105/joss.00598` (and see below for a BiBTeX entry). Even if you did not use **DynamicalSystems.jl** directly, but only a subcomponent, we still strongly encourage to cite **DynamicalSystems.jl** as well. More citations enable us to obtain funding to maintain and further develop **DynamicalSystems.jl** and the central citation associated with this effort is that assigned to the top-level `DynamicalSystems` module.
2. Cite the submodule (or subpackage, or subcomponent) that you used directly, if it has an associated publication (see example below).
3. Cite the specific algorithm(s) you used to also give credit to the originators of the methods.

For example, a typical citation that gives proper credit to using **DynamicalSystems.jl** could be:

> For our work we used the Sample Entropy `\cite{@SampleEntropy}` that is implemented in the ComplexityMeasures.jl component `\cite{ComplexityMeasures}` of the DynamicalSystems.jl library `\cite{DynamicalSystems}`.

or, one more example:

> For the surrogates for the atmospheric temperature timeseries we used the random cascade surrogates `\cite{Palus2008} that is implemented in the TimeseriesSurrogates.jl component `\cite{TimeseriesSurrogates}` of the DynamicalSystems.jl library `\cite{DynamicalSystems}`.

To find the specific references associated with the subpackages visit their respective dedicated documentation pages.

Besides the library, you may find useful the Nonlinear Dynamics textbook that utilizes **DynamicalSystems.jl**, which you can cite as:

```
@book{DatserisParlitz2022,
  doi = {10.1007/978-3-030-91032-7},
  url = {https://doi.org/10.1007/978-3-030-91032-7},
  year = {2022},
  publisher = {Springer Nature},
  author = {George Datseris and Ulrich Parlitz},
  title     = "Nonlinear dynamics: A concise introduction interlaced with code",
  address   = "Cham, Switzerland",
  language  = "en",
}
```

## [Asking questions](@id ask_questions)

There are three options for asking questions:

1. As a new post in the official [Julia discourse](https://discourse.julialang.org/) and ask a question under the category Specific Domains > Modelling & Simulations, also using `dynamical-systems` as a tag. This option is preferred for any meaningfully involved question, as the answer there will be future-searchable.
2. As a message in our channel `#dynamics-bridged` in the [Julia Slack](https://julialang.org/slack/) workplace. This option is preferred for a brief question with (expected) simple answer, or to get an opinion about something, or to chat about something.
3. By opening an issue directly on the [GitHub page of DynamicalSystems.jl](https://github.com/JuliaDynamics/DynamicalSystems.jl) while providing a Minimal Working Example. This option is preferred when you encounter unexpected behavior.

## Contributing

Be sure to visit the [Contributor Guide](@ref) page, because you can help make this package better without having to write a single line of code! Also, if you find this package helpful please consider staring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl)! This gives us an accurate lower bound of users that this package has already helped!

## Maintainers and Contributors

The **DynamicalSystems.jl** library is maintained by [George Datseris](https://github.com/Datseris), who is also curating and writing this documentation.
The software code however is built from the contributions of several individuals.
The list is too long to write and constantly update, so the best way to find out these contributions is to visit the GitHub page of each of the subpackages and checkout their respective "contributors" pages.

## Version numbers and SemVer

The version of `DynamicalSystems` by itself is a bit meaningless, because the module does not have any source code, besides re-exporting other modules and offering some visualization functionality.
For transparency, the packages and versions used to build the documentation you are reading now are:

```@setup MAIN
using CairoMakie, DynamicalSystems
```

```@example MAIN
using Pkg
Pkg.status([
    "DynamicalSystems",
    "StateSpaceSets", "DynamicalSystemsBase", "RecurrenceAnalysis", "FractalDimensions", "DelayEmbeddings", "ComplexityMeasures", "TimeseriesSurrogates", "PredefinedDynamicalSystems", "Attractors", "TransitionsInTimeseries", "SignalDecomposition", "ChaosTools", "CairoMakie",
    ];
    mode = PKGMODE_MANIFEST
)
```

!!! warn "Version numbers do not strictly follow SemVer2.0"
    Because of the nature of the **DynamicalSystems.jl** library, the exported API contains hundreds of algorithm implementations, most of which are independent of each other. Our development approach is that mildly breaking changes to these individual algorithms (due to e.g., better API design or better performance implementations or better default keyword arguments) can be done **without incrementing any major version numbers**. We increment major version numbers only for breaking changes that have wide impact over most of the **DynamicalSystems.jl** library.

    Every single subpackage of DynamicalSystems.jl has a human-written CHANGELOG.md file that details all changes done in each version. You should consult this package if you want to know what changed from version to version.


## Other NLD-relevant packages

Besides DynamicalSystems.jl, the Julia programming language has a thriving ecosystem with plenty of functionality that is relevant for nonlinear dynamics. We list some useful references below:

* [DifferentialEquations.jl](https://diffeq.sciml.ai/dev/index.html) - Besides providing solvers for standard ODE systems (infastructure already used in DynamicalSystems.jl), it also has much more features like SDE solvers or uncertainty quantification.
* [SciMLSensitivity.jl](https://github.com/SciML/SciMLSensitivity.jl) - Discrete and continuous local sensitivity analysis, i.e., derivatives of the solutions of ODEs, or functions of the solutions, versus parameters, hosting [various forward and adjoint methods as well as methods tailored to chaotic systems](https://docs.sciml.ai/SciMLSensitivity/stable/tutorials/chaotic_ode/).
* [GlobalSensitivity.jl](https://github.com/SciML/GlobalSensitivity.jl) - Global sensitivity analysis assessing the effect of any input variables over a larger domain on the output.
* [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl) - Featureful toolkit for automated bifurcation analysis.
* [NetworkDynamics.jl](https://github.com/PIK-ICoNe/NetworkDynamics.jl) - Simulating dynamics on networks and transforming network systems into `ODEProblem` (that can be made directly into a `ContinuousDynamicalSystem`).
* [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) - Agent based modelling.
* [EasyModelAnalysis.jl](https://github.com/SciML/EasyModelAnalysis.jl) - Analysis tools for conveniently analysing solutions of DiffEq systems.
* [ARFIMA.jl](https://github.com/JuliaDynamics/ARFIMA.jl) - generate ARFIMA process timeseries.
* [ConcurrentSim.jl](https://github.com/JuliaDynamics/ConcurrentSim.jl) - discrete event process oriented simulation framework.
* [Associations.jl](https://github.com/JuliaDynamics/Associations.jl) - hundreds of algorithms for relational/causal timeseries analysis and causal graphs.
