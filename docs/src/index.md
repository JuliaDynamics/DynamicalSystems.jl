![DynamicalSystems.jl logo: The Double Pendulum](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/dynamicalsystems_logo.gif?raw=true)

## Introduction
**DynamicalSystems.jl** is an award-winning Julia software library for the exploration of chaos and nonlinear dynamics.

!!! info "DSWeb Competition 1st place winner"
    **DynamicalSystems.jl** has won the first place in the 2018 Software Contest held by the dynamical systems division of SIAM. For more details see the [official statement](https://dsweb.siam.org/The-Magazine/Article/winners-of-the-dsweb-2018-software-contest)

!!! info "JuliaDynamics"
    **DynamicalSystems.jl** is part of [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/), check out our [website](https://juliadynamics.github.io/JuliaDynamics/) for more cool stuff!



The documentation you are reading now was built with the following stable versions:
```@setup versions
using Pkg.API: installed
ins = installed()
function f()
for pkg in ["DelayEmbeddings", "RecurrenceAnalysis", "DynamicalSystemsBase", "ChaosTools", "InteractiveChaos"]
  println(rpad(" * $(pkg) ", 30, "."), " $(ins[pkg])")
end
end
```
```@example versions
f() # hide
```
See the [News](news) page for recent updates!

!!! info "Introductory textbooks"
    Our library assumes some basic knowledge of nonlinear dynamics and complex systems.

    If you are new to the field but want to learn more, we can suggest the following textbooks as introductions:
    * Nonlinear Dynamics And Chaos - S. Strogatz
    * An Exploration of Dynamical Systems and Chaos - J. Argyris *et al.*
    * Chaos in Dynamical Systems - E. Ott

!!! example "Jupyter Notebooks / Tutorials"
    [In this repository](https://github.com/JuliaDynamics/JuliaDynamics/tree/master/tutorials) you can find various Jupyter notebooks that have been used as introductory tutorials for **DynamicalSystems.jl**!

## YouTube tutorial

You can find a tutorial on **DynamicalSystems.jl** hosted on the official YouTube channel of the Julia language:

<iframe width="560" height="315" src="https://www.youtube.com/embed/A8g9rdEfdNg" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Contents
The _module_ `DynamicalSystems` re-exports all following functionality, grouped into different packages.

### [`DynamicalSystemsBase`](ds/general)
* Intuitive, consistent APIs for the definition of general [dynamical systems](definition/general), under a unified struct [`DynamicalSystem`](@ref). The following combinations are possible:
  * Continuous or Discrete systems. Continuous systems use [`DifferentialEquations.jl`](http://docs.juliadiffeq.org/latest/) for solving the ODE problem.
  * In-place or out-of-place (large versus small systems).
  * Auto-differentiated or not (for the Jacobian function).


* Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user.
* Robust implementations of all kinds of integrators, that evolve the system, many states of the system, or even deviation vectors. See the [advanced documentation](advanced) for this.
* Library of [predefined well-known dynamical systems](ds/predefined) that have been used extensively in scientific research.

### [`DelayEmbeddings`](embedding/reconstruction)
Is a package for performing delay coordinate embeddings and finding optimal parameters for doing so.
* Unified & dedicated interface for [numerical data](embedding/dataset): [`Dataset`](@ref).
* Simple and extendable [`neighborhood`](@ref) estimation by interfacing [`NearestNeighbors`](https://github.com/KristofferC/NearestNeighbors.jl).
* Flexible and abstracted [`reconstruct`](@ref) interface, that creates the delay-coordinates reconstruction of a timeseries efficiently.
    * Supports multiple dimensions and multiple timescales.

* Methods that estimate optimal embedding parameters: the delay time ([`estimate_delay`](@ref)) and the number of temporal neighbors  ([`estimate_dimension`](@ref)).

* Fast calculation of mutual information: [`mutualinformation`](@ref).

### [`ChaosTools`](chaos/overview)
Is a package that has many algorithms for chaotic dynamical systems. All algorithms are independent of each other but they are also not expansive enough to be a standalone package.

Please see the [overview section](chaos/overview) for a full list of features. Here is a quick summary:

* Poincare S.O.S. and orbit diagrams
* Lyapunov Exponents
* Entropies and Dimensions
* Lyapunov exponent of a timeseries (numerical data)
* Finding Fixed Points of Maps
* Categorizing and distinguishing chaotic and regular behavior
* Nonlinear timeseries analysis

### [`RecurrenceAnalysis`](rqa/rplots)
`RecurrenceAnalysis` offers tools to compute and analyze [Recurrence Plots](https://en.wikipedia.org/wiki/Recurrence_plot), a field called [Recurrence Quantification Analysis](https://en.wikipedia.org/wiki/Recurrence_quantification_analysis).

* [Recurrence, cross-recurrence and joint-recurrence "plots" (they are matrices)](rqa/rplots).
* [Recurrence quantification analysis](rqa/quantification) (RQA):
  * Recurrence rate, determinism, average/maximum diagonal length, divergence, laminarity, trend, entropy, trapping time, average/maximum vertical length.
  * Fine-tuning of the algorithms that compute the above (e.g. Theiler window and many more)
  * [Windowed versions](rqa/windowed) of the above

### `InteractiveChaos`
`InteractiveChaos` is a package that provides applications for interactively exploring dynamical systems. It is an _extension_ of **DynamicalSystems.jl** and builds upon existing code by hooking up to the `DynamicalSystem` structure.

`InteractiveChaos` **is not installed with DynamicalSystems.jl**. To install it do `]add InteractiveChaos Makie`. `Makie` is necessary for providing a plotting backend, since `InteractiveChaos` does not install one by default.

The functionality of `InteractiveChaos` is contained within individual functions, all of which launch a dedicated interactive application. Here is their list:

* [`interactive_orbitdiagram`](@ref)
* [`interactive_poincaresos`](@ref)
* [`trajectory_highlighter`](@ref)

!!! info "Videos & Animations"
    Besides the documentation strings, each interactive function is accompanied with an animation (`.gif` or `.mp4` file) displayed after the docstring, as well as a video tutorial demonstrating its use. See the individual pages for the video links (by clicking the documentation string links)!

## Our Goals
The ultimate goal for **DynamicalSystems.jl** is to be a useful **software library** for students and scientists working on chaos, nonlinear dynamics and in general dynamical systems. The word "library" is intended in the literal sense: a place where people go to learn things.

With **DynamicalSystems.jl** we try to

1. Be concise, intuitive, and general. All functions we offer work just as well with any system, whether it is a simple continuous chaotic system, like the [Lorenz attractor](definition/predefined/#DynamicalSystemsBase.Systems.lorenz), or a high dimensional discrete map like [coupled standard maps](definition/predefined/#DynamicalSystemsBase.Systems.coupledstandardmaps).
2. Be accurate, reliable and performant.
3. Be transparent with respect to what is happening "under the hood", i.e. be clear about exactly what each function call does. We take care of this aspect in many ways; by being well-documented, giving references to scientific papers and having clear source code.

## Installation
Simply use `]add DynamicalSystems` to install everything. Alternatively you can also do `using Pkg; Pkg.add("DynamicalSystems")`.

---

For more advanced users, you can choose which packages to install and use at a high level. The *package* `DynamicalSystems` serves two purposes: it re-exports everything under a single module `DynamicalSystems` and it also builds the documentation.

All packages depend on `DelayEmbeddings` which defines core numeric data structures and methods. For example `RecurrenceAnalysis` and `TimeseriesPrediction` depend only on `DelayEmbeddings`. Packages that require equations of motion also depend on `DynamicalSystemsBase`, like for example `ChaosTools`.

If you only need functionality of a specific package you can install only that one, e.g. `]add RecurrenceAnalysis` and only the minimum amount of requirements will be installed.

## Citing
There is a (very small) paper associated with **DynamicalSystems.jl**. If we have helped
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

You can [join our chatroom](https://gitter.im/JuliaDynamics/Lobby) for discussions and/or questions about the packages of the JuliaDynamics organization! If you are using the Julia Slack workplace, please join the channel `#dynamics-bridged`.

## Contributing & Donating

Be sure to visit the [Contributor Guide](contributors_guide) page, because you can help make this package better without having to write a single line of code! Also, if you find this package helpful please consider staring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl)! This gives us an accurate lower bound of users that this package has already helped!

Finally, you can donate for the development of **DynamicalSystems.jl**. You can do that by adding bounties to existing issues on the GitHub repositories (you can open new issues as well). Every issue has an automatic way to create a bounty using [Bountysource](https://www.bountysource.com/), see the first comment of each issue.
