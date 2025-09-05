![DynamicalSystems.jl logo: The Double Pendulum](https://github.com/JuliaDynamics/JuliaDynamics/blob/master/videos/dynamicalsystems/juliadynamics_logo_anim_dark.gif?raw=true)

[![](https://img.shields.io/badge/docs-online-blue.svg)](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/dev/)
[![DocBuild](https://github.com/juliadynamics/DynamicalSystems.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/DynamicalSystems.jl/actions)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00598/status.svg)](https://doi.org/10.21105/joss.00598)
[![Textbook](https://img.shields.io/badge/Textbook-10.1007%2F978--3--030--91032--7-purple)](https://link.springer.com/book/10.1007/978-3-030-91032-7)
[![Package Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FDynamicalSystems&query=total_requests&label=Downloads)](http://juliapkgstats.com/pkg/DynamicalSystems)

**DynamicalSystems.jl** is an [award-winning](https://dsweb.siam.org/The-Magazine/Article/winners-of-the-dsweb-2018-software-contest) Julia-based general-purpose software library for the whole of nonlinear dynamics and nonlinear timeseries analysis.

To install **DynamicalSystems.jl**, run `import Pkg; Pkg.add("DynamicalSystems")` as a Julia language command.
To learn how to use it and see its contents visit the documentation, which you can either find [online](https://juliadynamics.github.io/DynamicalSystems.jl/dev/) or build locally by running the `docs/make.jl` file.

**DynamicalSystems.jl** is part of [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/), an organization dedicated to creating high quality scientific software.

## Highlights

Aspects of **DynamicalSystems.jl** that make it stand out among other codebases for nonlinear dynamics or nonlinear timeseries analysis are:

- **Exceptional documentation**. All implemented algorithms provide a high-level scientific description of their functionality in their documentation string as well as references to scientific papers. The documentation features hundreds of tutorials and examples ranging from introductory to expert usage.
- **Accessible source code**. One of the main priorities of the library is that the source code of (almost) all implementations is small, simple, easy to understand and modify. This increases confidence, reduces bugs, and allows users to become developers without unnecessary effort.
- **Open source community project**. Built from the ground up entirely on GitHub, **DynamicalSystems.jl** is 100% open source and based on community contributions. Anyone can be a developer of the library. Everyone is welcomed.
- **Extensive content**. It aims to cover the entire field of nonlinear dynamics and nonlinear timeseries analysis. It has functionality for complexity measures, delay embeddings, periodic orbits, nonlocal stability analysis, continuation, chaos, fractal dimensions, surrogate testing, recurrence quantification analysis, and much more. Furthermore, all algorithms are "general" and work for any dynamical system applicable. Missing functionality that falls under this wide category of content is welcomed to be part of the library!
- **Well tested**. All implemented functionality is extensively tested. Each time any change in the code base is done, the extensive test suite is run and checked before merging the change in.
- **Extendable**. New contributions can become part of the library and be accessed by all users in the next release. Most importantly, all parts of the library follow professional standards in software design and implement extendable interfaces so that it is easy to contribute new functionality.
- **Active development**. It is a living, evolving project. Since its beginning in May 2017, **DynamicalSystems.jl** has had some activity every single month: new features, new packages, bugfixes. The developer team routinely answers users questions on official Julia language forums.
- **Performant**. Written entirely in Julia, heavily optimized and parallelized, and taking advantage of some of the best packages within the language, **DynamicalSystems.jl** is _really fast_.

## Goals

The **DynamicalSystems.jl** library started as a vision with three main goals;
These same goals now are the core pillars guiding development, and are largely the source of where the aforementioned unique highlights stem from.

### Goal 1: Accessible and reproducible nonlinear dynamics

The first goal of the library is to make this beautiful field **accessible and reproducible**.

**Accessible** means that if you read on some sorts of fancy algorithm online in a scientific article, you should be able to use it instantly. You shouldn't have to put in the work to code it yourself. The authors of the paper already did that.
_So why should you do it again?!_ To resolve this problem we developed, and continue to develop, a library that has an incredibly low threshold of entry: contributing to **DynamicalSystems.jl** and making your code available to all is truly _easier_ than coding your own algorithms from scratch, due to the well thought out and generic interfaces it provides for dynamical systems.

**Reproducible** means that given some sorts of dynamical systems analysis in a scientific article, you should be able to do _exactly the same analysis_ and get _exactly the same results_ (within some numeric precision) as the article.
After all, computers are deterministic constructs.
**DynamicalSystems.jl** allows this by (1) being written in a modern programming language with incredible environment and reproducibility support, (2) being well tested, and (3) by providing thousands of algorithms out of the box, allowing most dynamical systems analysis to be done instantly while implementing only as little new stuff as necessary.

### Goal 2: Library in the literal sense

**DynamicalSystems.jl** is not just a software library. It is also a library in the literal sense: _where people go to learn something new_ (here in particular for nonlinear dynamics).
That is why the documentation is of exceptionally high quality: detailed descriptions and explanations of algorithms, with references to the scientific articles articles. It is also partly a reason for the source code to be written as clearly as possible, so that it is examinable by any user.

### Goal 3: A general purpose software

The third goal is to fill the missing gap of a high quality _general purpose software_ for nonlinear dynamics which can be easily extended with new functionality. This can be particularly impactful in teaching.
You see, it is unfortunately rarely the case that real, _runnable_ code is shown in the classroom, because it is often long and messy. This is especially hurtful for nonlinear dynamics, a field where computer-assisted exploration is critical.

**DynamicalSystems.jl** provides teachers with a framework capable of demonstrating actual, real-world nonlinear dynamics code and its output, without having to invest the weeks to code the internal infrastructure themselves.
Its high level syntax requires writing little code to get lots of meaningful analysis done, while its extensive functionality covers most typical classroom applications.

### Goal 4: Stopping the endless re-invention of wheel

Because Nonlinear Dynamics as a field lacks a general purpose and "widely accepted" software, almost every software implementation starts from scratch.
While doing so much of the code written actually implements functionality that already exists in some other codebase for nonlinear dynamics.
This is astonishingly, and shamefully, prevalent in nonlinear dynamics, where up to 90% of the functionality of a codebase may already exist somewhere else.
Needless to say this is just an absolute waste of time!

**DynamicalSystems.jl** hopes to establish itself as the central software for nonlinear dynamics, from which new algorithms can be implementing.
Re-using all the well-thought out interfaces and functionality means that one does not have to waste time writing code for functionality that already exists.
Rather, they can focus on coding only the _new_ stuff!