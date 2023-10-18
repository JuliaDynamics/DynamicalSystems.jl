![DynamicalSystems.jl logo: The Double Pendulum](https://i.imgur.com/nFQFdB0.gif)

[![](https://img.shields.io/badge/docs-online-blue.svg)](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/)
[![DocBuild](https://github.com/juliadynamics/DynamicalSystems.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/DynamicalSystems.jl/actions)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00598/status.svg)](https://doi.org/10.21105/joss.00598)
[![Textbook](https://img.shields.io/badge/Textbook-10.1007%2F978--3--030--91032--7-purple)](https://link.springer.com/book/10.1007/978-3-030-91032-7)
[![Package Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/DynamicalSystems)](https://pkgs.genieframework.com?packages=DynamicalSystems)

**DynamicalSystems.jl** is an [award-winning](https://dsweb.siam.org/The-Magazine/Article/winners-of-the-dsweb-2018-software-contest) Julia software library for nonlinear dynamics and nonlinear timeseries analysis.

To install **DynamicalSystems.jl**, run `import Pkg; Pkg.add("DynamicalSystems")`.
To learn how to use it and see its contents visit the documentation, which you can either find [online](https://juliadynamics.github.io/DynamicalSystems.jl/dev/) or build locally by running the `docs/make.jl` file.

**DynamicalSystems.jl** is part of [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/), an organization dedicated to creating high quality scientific software.

## Highlights

Aspects of **DynamicalSystems.jl** that make it stand out among other codebases for nonlinear dynamics or nonlinear timeseries analysis are:

- **Exceptional documentation**. All implemented algorithms provide a high-level scientific description of their functionality in their documentation string as well as references to scientific papers. The documentation features hundreds of tutorials and examples ranging from introductory to expert usage.
- **Accessible source code**. One of the main priorities of the library is that the source code of (almost) all implementations is small, simple, easy to understand and modify. This increases confidence, reduces bugs, and allows users to become developers without unnecessary effort.
- **Open source community project**. Built from the ground up entirely on GitHub, **DynamicalSystems.jl** is 100% open source and built by community contributions. Anyone can be a developer of the library. Everyone is welcomed.
- **Extensive content**. It aims to cover the entire field of nonlinear dynamics. It has functionality for complexity measures, delay embeddings, stability and bifurcation analysis, chaos, surrogate testing, recurrence quantification analysis, and much more. Furthermore, all algorithms are "general" and work for any dynamical system applicable. Missing functionality that falls under nonlinear dynamics is welcomed to be part of the library!
- **Well tested**. All implemented functionality is extensively tested. Each time any change in the code base is done, the extensive test suite is run and checked before merging the change in.
- **Extendable**. **DynamicalSystems.jl** is a living, evolving project. New contributions can become part of the library and be accessed by all users in the next release. Most importantly, all parts of the library follow professional standards in software design and implement extendable interfaces so that it is easy to contribute new functionality.
- **Active development**. Since the start of the project (May 2017) there is activity each month: new features, bugfixes, and the developer team answers users questions on Discourse/Slack.
- **Performant**. Written entirely in Julia, heavily optimized and parallelized, and taking advantage of some of the best packages within the language, **DynamicalSystems.jl** is _really fast_.

## Goals

The primary goal of **DynamicalSystems.jl** is to be a library in the literal sense: where people go to learn something (here in particular for nonlinear dynamics). That is why the main priority is that the documentation is detailed and references articles and why the source code is written as clearly as possible, so that it is examinable by any user.

The second goal is to fill the missing gap of high quality general purpose software for nonlinear dynamics which can be easily extended with new functionality. The purpose of this is to make the field of nonlinear dynamics accessible and reproducible.

The third goal is to fundamentally change the perception of the role of code in both scientific education as well as research.
It is rarely the case that real, _runnable_ code is shown in the classroom, because it is often long and messy.
This is especially hurtful for nonlinear dynamics, a field where computer-assisted exploration is critical.
And published scientific work in this field fares even worse, with the overwhelming majority of published research not sharing the code used to create the paper.
This makes reproducing these papers difficult, while some times straight-out impossible.
**DynamicalSystems.jl** can change this situation, because it is high level (requires writing little code to get lots of results) while offering extensive and well-tested functionality.
