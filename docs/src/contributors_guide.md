# Contributor Guide
*You can contribute to this package even if you are very good with
coding with Julia.*

## Reporting Bugs and other Issues
The easiest and most common way to improve this package is
simply by *using it* and reporting
any unexpected behavior! You can use the
[DynamicalSystems.jl Issues](https://github.com/JuliaDynamics/DynamicalSystems.jl/issues) page or, if you think it is something minor not worth opening an issue,
just come over to our [gitter chatroom]((https://gitter.im/JuliaDynamics/Lobby))
and let us know what the problem is.

## Contributing New Methods and Algorithms

The ultimate goal for DynamicalSystems.jl is
to be a useful tool for scientists working on chaos, nonlinear dynamics and
in general dynamical systems.

For such a feat to be accomplished, many different methods across this interdisciplinary
field have to be not only implemented but suggested in the first place!

For a something to be implemented
in this package, the following steps have to happen:

1. *A suggestion that a method should be included* has to be brought upon notice
   of the developers. Since the current amount of developers actively maintaining
   the package is small, so is the amount of knowledge of important methods.
2. An algorithm that describes how the method will be implemented has to be formulated.
   This
   algorithm most probably already exist in the papers that first introduce the method,
   however it may not be trivial to transform this algorithm from a mathematical
   abstraction to something realistic and applicable in a computational manner.
3. The source code for the above has to be implemented in Julia. In general, the
   speed of the implementation is important, but not as important as the
   **reliability of the implementation**.

It is clear that one can contribute to `DynamicalSystems.jl` by contributing in steps
(1) and (2). Neither of those require any hardcore coding knowledge with Julia.

For step (1), you can open a new issue at the [DynamicalSystems.jl Issues](https://github.com/JuliaDynamics/DynamicalSystems.jl/issues) page. All issues
that refer methods that we would want to have in our package are labeled as
"wanted_feature". You can view the current wanted features [here](https://github.com/JuliaDynamics/DynamicalSystems.jl/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Awanted_feature) and see for yourself if you can contribute
to some of them!

If you have any idea about how to improve
this package please do not hesitate to [join our chatroom](https://gitter.im/JuliaDynamics/Lobby) and share your ideas!



## Examples of new things you could contribute

* Any method that calculates a quantity that has been used in at least one scientific
  (and peer-reviewed) journal.
* Any kind of new **Type** of Dynamical system, provided it is also used in research.
  If you do want to make something like this, please make it a subtype
  of `DynamicalSystem`. I have created the discrete and continuous general types, but
  more specialized types would allow for specialized methods.
* Any kind of existing discrete or continuous system that have been used in published
  literature at least once and you find it useful (put this in the
  `famous_systems.jl` file).

Notice that the above are not conclusive, but only examples!

## How you should contribute **code**

* For new methods and systems please follow the convention
  of the documentation strings (do e.g. `?lyapunov` to see how they are
  structured).
* Have enough comments in your code so that somebody that knows the method,
  can also understand the code immediately.
* Always have a reference to the original work that first introduces the method
  or the system that you are using. You should put this reference
  to the main function's documentation string.
  See the existing documentation strings and do
  it in an identical manner.
* If you are adding a new system type accompanied by new methods and quantities, please
  introduce a new file in the `/src` folder, or even better create your own subfolder, instead of adding code to the existing files.

When enhancing already existing code, make sure to:

* Have enough comments at parts that are not easily understood, so that somebody
  else may continue your work in the future.
