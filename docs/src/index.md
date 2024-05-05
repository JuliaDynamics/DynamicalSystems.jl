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

DynamicalSystems.jl now integrates with ModelingToolkit.jl and allows using symbolic variables to access/observe state and parameter.

At a low level, this happens via the functions `observe_state`, `set_state!`,
`current_parameter` and `set_parameter!`.

Additionally, `interactive_trajectory_timeseries` allows symbolic indexing
for state space plot, timeseries plots, or parameter sliders.
Everything is also automatically named and limits are also automatically deduced for everything! Super convenient!

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


## Citing

There is a (small) paper associated with **DynamicalSystems.jl**. If we have helped
you in research that led to a publication, please cite it using
the DOI `10.21105/joss.00598` or the following BiBTeX entry:
```
@article{Datseris2018,
  doi = {10.21105/joss.00598},
  url = {https://doi.org/10.21105/joss.00598},
  year  = {2018},
  month = {mar},
  volume = {3},
  number = {23},
  pages = {598},
  author = {George Datseris},
  title = {DynamicalSystems.jl: A Julia software library for chaos and nonlinear dynamics},
  journal = {Journal of Open Source Software}
}
```

Irrespectively of **DynamicalSystems.jl**, _please also cite the specific algorithm that you used from the library_. The documentation of the function used will point you to the correct reference.

Besides the library, we would also appreciate it if you cited the textbook we wrote that **DynamicalSystems.jl** accompanies:

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

## Contributing & Donating

Be sure to visit the [Contributor Guide](@ref) page, because you can help make this package better without having to write a single line of code! Also, if you find this package helpful please consider staring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl)! This gives us an accurate lower bound of users that this package has already helped!

Finally, you can donate for the development of **DynamicalSystems.jl**. You can do that by adding bounties to existing issues on the GitHub repositories (you can open new issues as well). Every issue has an automatic way to create a bounty using [Bountysource](https://www.bountysource.com/), see the first comment of each issue.

### Issues with Bounties

Money that **DynamicalSystems.jl** obtains from awards, sponsors, or donators are converted into bounties for GitHub issues. The full list of issues that have a bounty is [available here](https://github.com/issues?utf8=%E2%9C%93&q=is%3Aopen+is%3Aissue+org%3AJuliaDynamics+label%3Abounty).

By solving these issues you not only contribute to open source, but you also get some pocket money to boot :)


## Maintainers and Contributors

The **DynamicalSystems.jl** library is maintained by [George Datseris](https://github.com/Datseris), who is also curating and writing this documentation.
The software code however is built from the contributions of several individuals.
The list is too long to write and constantly update, so the best way to find out these contributions is to visit the GitHub page of each of the subpackages and checkout the "contributors" pages there.

## Version numbers

The version of `DynamicalSystems` by itself is a bit meaningless, because the module does not have any source code, besides re-exporting other modules.
For transparency, the packages and versions used to build the documentation you are reading now are:

```@setup MAIN
using CairoMakie, DynamicalSystems
```

```@example MAIN
using Pkg
Pkg.status([
    "DynamicalSystems",
    "StateSpaceSets", "DynamicalSystemsBase", "RecurrenceAnalysis", "FractalDimensions", "DelayEmbeddings", "ComplexityMeasures", "TimeseriesSurrogates", "PredefinedDynamicalSystems", "Attractors", "ChaosTools", "CairoMakie",
    ];
    mode = PKGMODE_MANIFEST
)
```

!!! warn "Version numbers do not strictly follow SemVer2.0"
    Because of the nature of the **DynamicalSystems.jl** library, the exported API contains hundreds of algorithm implementations, most of which are independent of each other. Our development approach is that breaking changes to these individual algorithms (due to e.g., better API design or better performance implementations or better default keyword arguments) can be done **without incrementing any major version numbers**. We increment major version numbers only for breaking changes that have wide impact over most of the **DynamicalSystems.jl** library.


## Other NLD-relevant packages

Besides DynamicalSystems.jl, the Julia programming language has a thriving ecosystem with plenty of functionality that is relevant for nonlinear dynamics. We list some useful references below:

* [DifferentialEquations.jl](https://diffeq.sciml.ai/dev/index.html) - Besides providing solvers for standard ODE systems (infastructure already used in DynamicalSystems.jl), it also has much more features like SDE solvers or uncertainty quantification.
* [DiffEqSensitivity.jl](https://github.com/SciML/DiffEqSensitivity.jl) - Discrete and continuous local sensitivity analysis, i.e., derivatives of the solutions of ODEs, or functions of the solutions, versus parameters, hosting [various forward and adjoint methods as well as methods tailored to chaotic systems](https://diffeq.sciml.ai/stable/analysis/sensitivity/).
* [GlobalSensitivity.jl](https://github.com/SciML/GlobalSensitivity.jl) - Global sensitivity analysis assessing the effect of any input variables over a larger domain on the output.
* [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl) - Featureful toolkit for automated bifurcation analysis.
* [NetworkDynamics.jl](https://github.com/PIK-ICoNe/NetworkDynamics.jl) - Simulating dynamics on networks and transforming network systems into `ODEProblem` (that can be made directly into a `ContinuousDynamicalSystem`).
* [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) - Agent based modelling.
* [EasyModelAnalysis.jl](https://github.com/SciML/EasyModelAnalysis.jl) - Analysis tools for conveniently analysing solutions of DiffEq systems.
* [SignalDecomposition.jl](https://github.com/JuliaDynamics/SignalDecomposition.jl) - Decompose a signal/timeseries into structure and noise or seasonal and residual components.
* [ARFIMA.jl](https://github.com/JuliaDynamics/ARFIMA.jl) - generate ARFIMA process timeseries.
* [ConcurrentSim.jl](https://github.com/JuliaDynamics/ConcurrentSim.jl) - discrete event process oriented simulation framework.
* [CausalityTools.jl](https://github.com/JuliaDynamics/CausalityTools.jl) - hundreds of algorithms for relational/causal timeseries analysis and causal graphs.