![DynamicalSystems.jl logo: The Double Pendulum](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/dynamicalsystems_logo.gif?raw=true)

**DynamicalSystems.jl** is an [award-winning](https://dsweb.siam.org/The-Magazine/Article/winners-of-the-dsweb-2018-software-contest) Julia software library for the exploration of chaos and nonlinear dynamics.
It is part of [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/), an organization dedicated to creating high quality scientific software.

To learn how to use this library please see [Getting started](@ref) below, and subsequently, the [Contents](@ref) page to get an overview of all offered functionality of **DynamicalSystems.jl**.

!!! tip "Latest news"
    Rework and improvement of optimal embedding dimension: [`optimal_traditional_de`](@ref)!

!!! info "Star us on GitHub!"
    If you have found this library useful, please consider starring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl).
    This gives us an accurate lower bound of the (satisfied) user count.

## Getting started
**DynamicalSystems.jl** is a collection of Julia packages bundled together under a single package `DynamicalSystems`. To install this bundle you can do:
```julia
using Pkg; Pkg.add("DynamicalSystems")
```

The individual packages that compose `DynamicalSystems` interact flawlessly with each other because of the following two structures:

1. The [`DynamicalSystem`](@ref) represents a dynamical system with known dynamic rule ``f``. The system can be in discrete time (often called a map), ``\vec{u}_{n+1} = \vec{f}(\vec{u}_n, p, n)``, or in continuous time (often called an ordinary differential equation) ``\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}, p, t)``. In both cases ``u`` is the _state_ of the dynamical system and ``p`` a parameter container. You should have a look at the page [Dynamical System Definition](@ref) for how to create this object. A list of several pre-defined systems exists in the [Predefined Dynamical Systems](@ref) page.
2. Numerical data, that can represent measured experiments, or sampled trajectories of dynamical systems, are represented by [`Dataset`](@ref), which is a container of equally-sized data points. Timeseries in **DynamicalSystems.jl** are represented by the already existing `Vector` type of the Julia language.

These core structures `DynamicalSystem, Dataset` are used throughout the package to do useful calculations often used in the field of nonlinear dynamics and chaos.
For example, using [`lyapunovspectrum`](@ref) and [`DynamicalSystem`](@ref) gives you the Lyapunov exponents of a dynamical system with known equations of motion.
Alternatively, by using [`numericallyapunov`](@ref) and [`Dataset`](@ref) you can approximate the maximum Lyapunov exponent of a measured trajectory.

All things possible in **DynamicalSystems.jl** are listed in the [Contents](@ref) page.

### Tutorials
Tutorials for **DynamicalSystems.jl** exist in the form of [Jupyter notebooks](https://github.com/JuliaDynamics/JuliaDynamics/tree/master/tutorials).

In addition, a full 2-hours YouTube tutorial is available below:

```@raw html
<iframe width="560" height="400" src="https://www.youtube.com/embed/A8g9rdEfdNg" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

!!! info "Introductory textbooks"
    Our library assumes some basic knowledge of nonlinear dynamics and complex systems.

    If you are new to the field but want to learn more, we can suggest the following textbooks as introductions:
    * Chaos in Dynamical Systems - E. Ott
    * Nonlinear Time series Analysis - H. Kantz & T. Schreiber


### Advanced installation
For more advanced users, you can choose which packages to install and use at a high level.
All packages depend on `DelayEmbeddings` which defines core numeric data structures and methods. For example `RecurrenceAnalysis` and `TimeseriesPrediction` depend only on `DelayEmbeddings`. Packages that require equations of motion also depend on `DynamicalSystemsBase`, like for example `ChaosTools`.

If you only need functionality of a specific package you can install only that one, e.g. `]add RecurrenceAnalysis` and only the minimum amount of requirements will be installed.

The documentation you are reading now was built with the following stable versions:
```@example versions
using Pkg
Pkg.status([
    "DelayEmbeddings", "RecurrenceAnalysis",
    "DynamicalSystemsBase", "ChaosTools",
    "Entropies",
])
```


## Our Goals
**DynamicalSystems.jl** was created with three goals in mind.
The first was to fill the missing gap of a software for nonlinear dynamics and chaos of the highest quality (none exist in any programming language).
The second was to create a useful _library_ where students and scientists from different fields may come and learn about methods of nonlinear dynamics and chaos.

The third was to fundamentally change the perception of the role of code in both scientific education as well as research.
It is rarely the case that real, _runnable_ code is shown in the classroom, because it is often long and messy.
This is especially hurtful for nonlinear dynamics, a field where computer-assisted exploration is critical.
But published work in this field fares even worse, with the overwhelming majority of published research not sharing the code used to create the paper.
This makes reproducing these papers difficult, while some times straight-out impossible.

To achieve these goals we made **DynamicalSystems.jl** so that it is:

1. Transparent: extra care is taken so that the source code of all functions is clear and easy to follow, while remaining as small and concise as possible.
1. Intuitive: a software simple to use and understand makes experimentation easier.
1. Easy to install, easy to extend: This makes contributions more likely, and can motivate researchers to implement their method here, instead of leaving it in a cryptic script stored in some data server, never-to-be-published with the paper.
1. Reliable: the algorithm implementations are tested extensively.
1. Well-documented: all implemented algorithms provide a high-level scientific description of their functionality in their documentation string as well as references to scientific papers.
1. General: all algorithms work just as well with any system, whether it is a simple continuous chaotic system, like the Lorenz model, or a high dimensional discrete system like coupled standard maps.
1. Performant: written entirely in Julia, and taking advantage of some of the best packages within the language, **DynamicalSystems.jl** is _really fast_.

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

## Issues with Bounties
Money that **DynamicalSystems.jl** obtains from awards, sponsors or donators are converted into bounties for GitHub issues. The full list of issues that have a bounty is [available here](https://github.com/issues?utf8=%E2%9C%93&q=is%3Aopen+is%3Aissue+org%3AJuliaDynamics+label%3Abounty).

By solving these issues you not only contribute to open source, but you also get some pocket money to boot :)

## Contacting

Feel free to open issues on GitHub if you have questions and/or suggestions.
You can also [join our chatroom](https://gitter.im/JuliaDynamics/Lobby) for discussions and/or questions about the packages of the JuliaDynamics organization! If you are using the Julia Slack workplace, please join the channel `#dynamics-bridged`.

## Contributing & Donating

*TL;DR: See ["good first issues"](https://github.com/issues?q=is%3Aopen+is%3Aissue+repo%3AJuliaDynamics%2FChaosTools.jl+repo%3AJuliaDynamics%2FDynamicalSystemsBase.jl+repo%3AJuliaDynamics%2FDelayEmbeddings.jl+repo%3AJuliaDynamics%2FRecurrenceAnalysis.jl+repo%3AJuliaDynamics%2FDynamicalSystems.jl+label%3A%22good+first+issue%22+) or ["wanted features"](https://github.com/issues?q=is%3Aopen+is%3Aissue+repo%3AJuliaDynamics%2FChaosTools.jl+repo%3AJuliaDynamics%2FDynamicalSystemsBase.jl+repo%3AJuliaDynamics%2FDelayEmbeddings.jl+repo%3AJuliaDynamics%2FRecurrenceAnalysis.jl+repo%3AJuliaDynamics%2FDynamicalSystems.jl+label%3A%22wanted+feature%22+).*

---

Be sure to visit the [Contributor Guide](@ref) page, because you can help make this package better without having to write a single line of code! Also, if you find this package helpful please consider staring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl)! This gives us an accurate lower bound of users that this package has already helped!

Finally, you can donate for the development of **DynamicalSystems.jl**. You can do that by adding bounties to existing issues on the GitHub repositories (you can open new issues as well). Every issue has an automatic way to create a bounty using [Bountysource](https://www.bountysource.com/), see the first comment of each issue.

## Maintainers and Contributors
The DynamicalSystems.jl software is maintained by [George Datseris](https://github.com/Datseris), who is also curating and writing this documentation page.

The software code however is built from the contributions of several individuals. For an accurate list of the names as well as contributions of each one, please visit the GitHub's contributor list for the sub-packages of DynamicalSystems.jl, e.g. [ChaosTools.jl](https://github.com/JuliaDynamics/ChaosTools.jl/graphs/contributors).
