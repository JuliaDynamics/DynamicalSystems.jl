# Contributor Guide

*TL;DR: To contribute via Pull Requests see ["good first issues"](https://github.com/issues?q=is%3Aopen+is%3Aissue+repo%3AJuliaDynamics%2FChaosTools.jl+repo%3AJuliaDynamics%2FDynamicalSystemsBase.jl+repo%3AJuliaDynamics%2FDelayEmbeddings.jl+repo%3AJuliaDynamics%2FRecurrenceAnalysis.jl+repo%3AJuliaDynamics%2FEntropies.jl+repo%3AJuliaDynamics%2FDynamicalSystems.jl+label%3A%22good+first+issue%22+) or ["wanted features"](https://github.com/issues?q=is%3Aopen+is%3Aissue+repo%3AJuliaDynamics%2FChaosTools.jl+repo%3AJuliaDynamics%2FDynamicalSystemsBase.jl+repo%3AJuliaDynamics%2FDelayEmbeddings.jl+repo%3AJuliaDynamics%2FRecurrenceAnalysis.jl+repo%3AJuliaDynamics%2FEntropies.jl+repo%3AJuliaDynamics%2FDynamicalSystems.jl+label%3A%22wanted+feature%22+). *

---

The ultimate goal for **DynamicalSystems.jl** is
to be a useful *library* for scientists working on nonlinear dynamics and to make nonlinear dynamics accessible and reproducible.

Of course, for such an ambitious goal to be achieved, many of us should try to
work together to improve the library! If you want to help the cause, there are many ways to contribute to the **DynamicalSystems.jl** library:

1. Just *use it*. If you encountered unexpected behavior simply report it either on
   our [gitter chatroom](https://gitter.im/JuliaDynamics/Lobby) or using the
   [DynamicalSystems.jl Issues](https://github.com/JuliaDynamics/DynamicalSystems.jl/issues) page.
2. Suggest methods that you think should be included in our library. This should be
   done by opening a new issue that describes the method, gives references to papers
   using the method and also justifies why the method should be included.
3. Contribute code by solving issues. The easiest issues to tackle are the ones with label ["good first issue"](https://github.com/issues?q=is%3Aopen+is%3Aissue+repo%3AJuliaDynamics%2FChaosTools.jl+repo%3AJuliaDynamics%2FDynamicalSystemsBase.jl+repo%3AJuliaDynamics%2FDelayEmbeddings.jl+repo%3AJuliaDynamics%2FRecurrenceAnalysis.jl+repo%3AJuliaDynamics%2FDynamicalSystems.jl+label%3A%22good+first+issue%22+).
4. Contribute code by implementing new methods! That is the most **awesome** way to
   contribute! The individual packages that compose **DynamicalSystems.jl** have plenty of issues with the tag ["wanted feature"](https://github.com/issues?q=is%3Aopen+is%3Aissue+repo%3AJuliaDynamics%2FChaosTools.jl+repo%3AJuliaDynamics%2FDynamicalSystemsBase.jl+repo%3AJuliaDynamics%2FDelayEmbeddings.jl+repo%3AJuliaDynamics%2FRecurrenceAnalysis.jl+repo%3AJuliaDynamics%2FDynamicalSystems.jl+label%3A%22wanted+feature%22+), which can get you started on a big contribution!
5. Contribute code by defining a new pre-defined dynamical system that you found useful.

## Contributing Code
When contributing code, you should keep these things in mind:

* In general, the
  speed of the implementation is important, but not as important as the
  **clarity of the implementation**. One of cornerstones of all of
  **DynamicalSystems.jl** is to have clear and readable source code. Fortunately,
  Julia allows you to have perfectly readable code but also super fast ;)
  If necessary add comments to the code, so that somebody that knows the method, can also understand the code immediately.
* Try to design general, extendable functions instead of unnecessarily specialized to the case at hand.
* For the documentation strings of new methods and systems please follow the convention of the documentation strings of DynamicalSystems.jl. Specifically, the first section should describe the function in a couple of sentences, its positional arguments and its return value. The next section `## Keyword Arguments` describes the keywords. The next section `## Description` describes the algorithm in detail if need be. Lastly, papers that are relevant to the method must be cited. Have a look at the documentation strings of `lyapunov` and `lyapunovspectrum` to get an idea.

## Documentation string style
Documentation strings are the most important thing in your pull request/code. The number 1 priority of DynamicalSystems.jl is highest possible quality of documentation and utmost transparency, and the best way to achieve this is with good documentation strings. In DynamicalSystems.jl we recommend that documentation strings are structured in the following way (and this is also the recommendation we give in the [Good Scientific Code Workshop](https://youtu.be/x3swaMSCcYk?t=11087)). 

1. Clear call signature in code syntax, including expected input types if necessary. The call signature should ONLY include only the most important information, not list out in detail every keyword!
1. Brief summary of the function
1. [Optional] Return value and type if not obvious
1. [Optional] References to related functions if sensible
1. [Optional] Keyword arguments list if the function has some
1. [Optional] Detailed discussion of functionality if function behavior is scientifically involved
1. [Optional] Citations to relevant scientific papers!

The syntax of the documentation strings follows Documenter.jl protocol. Please see the documentation string of the [lyapunov](https://github.com/JuliaDynamics/ChaosTools.jl/blob/a7ba7f559e24bd6e32d270b9a6281d4b919a20a1/src/chaosdetection/lyapunovs/lyapunov.jl#L4-L65) function and use the same structure.
