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
* For the documentation strings of new methods and systems please follow the convention of the documentation strings of DynamicalSystems.jl. Specifically, the first section should describe the function in a couple of sentences, its positional arguments and its return value. The next section `## Keyword Arguments` describes the keywords. The next section `## Description` describes the algorithm in detail if need be.
* Always have a reference to the original work that introduces the method
  or the system that you are using. You should put this reference
  to the main function's documentation string.
  See the existing documentation strings and do
  it in a similar manner.
