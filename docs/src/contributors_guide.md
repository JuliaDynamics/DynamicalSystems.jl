# Contributor Guide

The ultimate goal for **DynamicalSystems.jl** is
to be a useful *library* for scientists working on nonlinear dynamics and to make nonlinear dynamics accessible and reproducible.

Of course, for such an ambitious goal to be achieved, many of us should try to
work together to improve the library! If you want to help the cause, there are many ways to contribute to the **DynamicalSystems.jl** library:

1. Just *use it*! Share it with your colleagues if it was useful for you, and report unexpected behaviour if you find any (see [here](@ref ask_questions) for how).
2. Suggest methods that you think should be included in our library. This should be
   done by opening a new issue that describes the method, gives references to papers
   using the method and also justifies why the method should be included. Please open an issue to the GitHub page of the submodule of **DynamicalSystems.jl** that you feel is the most related to the method.
3. Contribute code by solving existing issues. The easiest issues to tackle are the ones with label "good first issue". Here is a list of all such issues from all submodules of **DynamicalSystems.jl**: [link](https://github.com/issues?page=1&q=is%3Aopen+is%3Aissue+repo%3AJuliaDynamics%2FChaosTools.jl+repo%3AJuliaDynamics%2FDynamicalSystemsBase.jl+repo%3AJuliaDynamics%2FDelayEmbeddings.jl+repo%3AJuliaDynamics%2FRecurrenceAnalysis.jl+repo%3AJuliaDynamics%2FDynamicalSystems.jl+repo%3AJuliaDynamics%2FAttractors.jl+repo%3AJuliaDynamics%2FComplexityMeasures.jl+repo%3AJuliaDynamics%2FFractalDimensions.jl+repo%3AJuliaDynamics%2FStateSpaceSets.jl+label%3A%22good+first+issue%22).
4. Contribute code by implementing new methods! That is by far the most impactful way to contribute to the library. The individual packages that compose **DynamicalSystems.jl** have plenty of issues that outline new methods wanted by the library, that are likely not tagged as "good first issues" because they will likely require familiarity that goes beyond a complete beginner. You can tackle one of these if you want to contribute! Additionally, we strongly welcome contributions of brand new algorithms that have just been developed during research in nonlinear dynamics. In fact, **DynamicalSystems.jl** started with the vision that researchers would add their newly developed methods directly to the library when they publish the methods.

## Contributing Code

When contributing code, principles for writing good scientific code apply. We recommend the [Good Scientific Code workshop] material for teaching you this.

You should keep these things in mind:

* In general, the
  speed of the implementation is important, but not as important as the
  **clarity of the implementation**. One of cornerstones of all of
  **DynamicalSystems.jl** is to have clear and readable source code. Fortunately,
  Julia allows you to have perfectly readable code but also super fast ;)
  If necessary add comments to the code, so that somebody that knows the method, can also understand the code immediately.
* Try to design general, extendable functions instead of unnecessarily specialized to the case at hand.
* The documentation strings of the new API functions you contribute are the most important to make as good as possible. Please follow the convention of the documentation strings of DynamicalSystems.jl outlined below.

## Documentation string style

Documentation strings are the most important thing in your pull request/code. The number 1 priority of DynamicalSystems.jl is highest possible quality of documentation and utmost transparency, and the best way to achieve this is with good documentation strings. In DynamicalSystems.jl we recommend that documentation strings are structured in the following way (and this is also the recommendation we give in the [Good Scientific Code Workshop](https://youtu.be/x3swaMSCcYk?t=11087)).

1. Clear call signature in code syntax, including expected input types if necessary. The call signature should ONLY include only the most important information, not list out in detail every keyword!
1. Brief summary of the function
1. [Optional] Return value and type if not obvious (almost always it is not obvious!)
1. [Optional] References to related functions if sensible with the `@ref` command.
1. [Optional] Keyword arguments list if the function has some with a `## Keyword arguments` subsection.
1. [Optional] Detailed discussion of functionality if function behavior is scientifically involved with a `## Description` subsection.
1. [Optional] Citations to relevant scientific papers with the `@cite` command.

The syntax of the documentation strings follows Documenter.jl protocol. Please see the documentation string of the [`lyapunov`](https://github.com/JuliaDynamics/ChaosTools.jl/blob/a7ba7f559e24bd6e32d270b9a6281d4b919a20a1/src/chaosdetection/lyapunovs/lyapunov.jl#L4-L65) function and use the same structure.
