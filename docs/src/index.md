![DynamicalSystems.jl logo: The Double Pendulum](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/chaos/dynamicalsystems_logo.gif?raw=true)

**DynamicalSystems.jl** is an [award-winning](https://dsweb.siam.org/The-Magazine/Article/winners-of-the-dsweb-2018-software-contest) Julia software library for the exploration of chaos and nonlinear dynamics.
It is part of [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/), an organization dedicated to creating high quality scientific software.

!!! tip "Latest news"
    Expansion entropy from Ott and Hunt implemented as [`expansionentropy`](@ref)!


The documentation you are reading now was built with the following stable versions:
```@setup versions
using Pkg.API: installed
ins = installed()
function f()
for pkg in ["DelayEmbeddings", "RecurrenceAnalysis", "DynamicalSystemsBase", "ChaosTools"]
  println(rpad(" * $(pkg) ", 30, "."), " $(ins[pkg])")
end
end
```
```@example versions
f() # hide
```

!!! info "Introductory textbooks"
    Our library assumes some basic knowledge of nonlinear dynamics and complex systems.

    If you are new to the field but want to learn more, we can suggest the following textbooks as introductions:
    * Nonlinear Dynamics And Chaos - S. Strogatz
    * An Exploration of Dynamical Systems and Chaos - J. Argyris *et al.*
    * Chaos in Dynamical Systems - E. Ott
    * Nonlinear Time series Analysis - H. Kantz & T. Schreiber
    * Chaos and Integrability in Nonlinear Dynamics - M. Tabor

## Installation
Simply use `]add DynamicalSystems` to install everything. Alternatively you can also do `using Pkg; Pkg.add("DynamicalSystems")`.

For more advanced users, you can choose which packages to install and use at a high level. The *package* `DynamicalSystems` serves two purposes: it re-exports everything under a single module `DynamicalSystems` and it also builds the documentation.

All packages depend on `DelayEmbeddings` which defines core numeric data structures and methods. For example `RecurrenceAnalysis` and `TimeseriesPrediction` depend only on `DelayEmbeddings`. Packages that require equations of motion also depend on `DynamicalSystemsBase`, like for example `ChaosTools`.

If you only need functionality of a specific package you can install only that one, e.g. `]add RecurrenceAnalysis` and only the minimum amount of requirements will be installed.

## Tutorials
Tutorials for **DynamicalSystems.jl** exist in the form of [Jupyter notebooks](https://github.com/JuliaDynamics/JuliaDynamics/tree/master/tutorials).

In addition, a full 2-hours YouTube tutorial is available below:

```@raw html
<iframe width="560" height="315" src="https://www.youtube.com/embed/A8g9rdEfdNg" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```


## Our Goals
The ultimate goal for **DynamicalSystems.jl** is to be a useful **software library** for students and scientists working on chaos, nonlinear dynamics and in general dynamical systems. The word "library" is intended in the literal sense: a place where people go to learn things.

With **DynamicalSystems.jl** we try to

1. Be concise, intuitive, and general. All functions we offer work just as well with any system, whether it is a simple continuous chaotic system, like the Lorenz attractor, or a high dimensional discrete map like coupled standard maps.
2. Be accurate, reliable and performant.
3. Be transparent with respect to what is happening "under the hood", i.e. be clear about exactly what each function call does. We take care of this aspect in many ways; by being well-documented, giving references to scientific papers and having clear source code.

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

Be sure to visit the [Contributor Guide](@ref) page, because you can help make this package better without having to write a single line of code! Also, if you find this package helpful please consider staring it on [GitHub](https://github.com/JuliaDynamics/DynamicalSystems.jl)! This gives us an accurate lower bound of users that this package has already helped!

Finally, you can donate for the development of **DynamicalSystems.jl**. You can do that by adding bounties to existing issues on the GitHub repositories (you can open new issues as well). Every issue has an automatic way to create a bounty using [Bountysource](https://www.bountysource.com/), see the first comment of each issue.
