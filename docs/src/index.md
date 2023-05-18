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
- Besides the formal algorithmic/scientific content of **DynamicalSystems.jl** (those in the [contents](@ref contents)) page, the library also provides an extensive suite for interactive or offline animations and visualizations dynamical systems. These are found in the [visualizations](@ref visualization) page.
- The remaining of this introduction page discusses our goals with the library, how to participate as a user or developer, how to cite, and other relevant information (see the sections of the sidebar on the left).


## Latest news

The new major version, v3, of **DynamicalSystems.jl** has been released! This is a massive release, with a huge amount of new content, and a re-write of many of the fundamentals of the library!

The full details of this release are too long to describe here, so please read more about it in our [changelog](https://github.com/JuliaDynamics/DynamicalSystems.jl/blob/main/CHANGELOG.md)!


## [Learning resources](@id learning)

### Textbook with DynamicalSystems.jl

We have written an undergraduate level textbook as an introduction to nonlinear dynamics. The text is written in an applied, hands-on manner, while still covering all fundamentals. The book pages are interlaced with real Julia code that uses DynamicalSystems.jl and is published in the Undergraduate Lecture Notes in Physics by Springer Nature:
* [Nonlinear Dynamics: A concise introduction interlaced with code](https://link.springer.com/book/10.1007/978-3-030-91032-7) by G. Datseris & U. Parlitz.


Additional textbooks on nonlinear dynamics worth having a look are:
* Chaos in Dynamical Systems - E. Ott
* Nonlinear Time series Analysis - H. Kantz & T. Schreiber

### Course on applied nonlinear dynamics and complex systems

We are developing a full course (targeting a graduate or undergraduate semester long course) on applied nonlinear dynamics, nonlinear timeseries analysis, and complex systems, using the packages of [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/). **DynamicalSystems.jl** is part of this course.

The materials of the course are on GitHub: https://github.com/JuliaDynamics/NonlinearDynamicsComplexSystemsCourses


## [Our Goals](@id goals)

**DynamicalSystems.jl** was created with three goals in mind.
The first was to fill the missing gap of a high quality and general purpose software for nonlinear dynamics, which can make the field of nonlinear dynamics accessible and reproducible.
The second goal was to create a useful _library_ where students and scientists from different fields may come and learn about methods of nonlinear dynamics.

The third goal was to fundamentally change the perception of the role of code in both scientific education as well as research.
It is rarely the case that real, _runnable_ code is shown in the classroom, because it is often long and messy.
This is especially hurtful for nonlinear dynamics, a field where computer-assisted exploration is critical.
And published work in this field fares even worse, with the overwhelming majority of published research not sharing the code used to create the paper.
This makes reproducing these papers difficult, while some times straight-out impossible.

To achieve these goals we made **DynamicalSystems.jl** so that it is:

1. **Transparent**: extra care is taken so that the source code of all functions is clear and easy to follow, while remaining as small and concise as possible.
2. **Intuitive**: a software simple to use and understand makes experimentation easier.
3. **Easy to extend**: this makes contributions more likely, and can motivate researchers to implement their method here, instead of leaving it in a cryptic script stored in some data server, never-to-be-published with the paper.
4. **Reliable**: the algorithm implementations are tested extensively.
5. **Well-documented**: all implemented algorithms provide a high-level scientific description of their functionality in their documentation string as well as references to scientific papers.
6. **General**: all algorithms work just as well with any system, whether it is a simple continuous chaotic system, like the Lorenz model, or a high dimensional discrete system like coupled standard maps.
7. **Performant**: written entirely in Julia, and taking advantage of some of the best packages within the language, **DynamicalSystems.jl** is _really fast_.

## Citing

There is a (small) paper associated with **DynamicalSystems.jl**. If we have helped
you in research that led to a publication, please be kind enough to cite it, using
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

however, we would really appreciate it if you also cited the textbook we wrote that **DynamicalSystems.jl** accompanies:

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


## Asking questions

There are three options for asking questions:

1. Join the official [Julia discourse](https://discourse.julialang.org/) and ask a question under the category Specific Domains > Modelling & Simulations.
2. Join our channel `#dynamics-bridged` in the [Julia Slack](https://julialang.org/slack/) workplace.
3. Open an issue directly on the GitHub page of DynamicalSystems.jl while providing a Minimal Working Example.

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
using CairoMakie, InteractiveDynamics, DynamicalSystems
include("../style.jl")
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
* [GlobalSensitivity.jl](https://github.com/SciML/GlobalSensitivity.jl) Global sensitivity analysis assessing the effect of any input variables over a larger domain on the output.
* [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl) - Featureful toolkit for automated bifurcation analysis.
* [NetworkDynamics.jl](https://github.com/PIK-ICoNe/NetworkDynamics.jl) - Package for easily simulating dynamics on networks and transforming network systems into `ODEProblem` (that can be made directly into a `ContinuousDynamicalSystem`).
* [Agents.jl](https://github.com/JuliaDynamics/Agents.jl) for agent based modelling.