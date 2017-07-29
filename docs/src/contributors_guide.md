# Contributor Guide

## What you can contribute

* Any method that calculates a quantity that has been used in at least one published
  (and peer-reviewed) journal.
* Any kind of new **Type** of Dynamical system, provided it is also used in research.
  If you do want to make something like this, please make it a subtype
  of `DynamicalSystem`. I have created the discrete and continuous general types, but
  more specialized types would allow for specialized methods.
* Any kind of existing discrete or continuous system that have been used in published
  literature at least once and you find it useful.
  Put this in the `famous_systems.jl` file.
* If you are adding a new system type accompanied by new methods and quantities, please
  introduce a new file in the `/src` folder, or even better create your own subfolder, instead of adding code to the existing files.

Also, you can contribute in the enhancement of the existing package by:
* Solving the [github issues](https://github.com/Datseris/DynamicalSystems.jl/issues)
* Improving the speed of the existing methods. I am not good at writing extremely
  fast code by inspecting assembly and all that.

## How you should contribute
This section will be updated as the package matures.

* For new methods and systems please always have very clear and self-contained
  documentation strings.
* Always have a reference to the original work that first introduces the method
  or the system that you are using. See the existing documentation strings and do
  it in an identical manner.

When enhancing already existing code, make sure to:
* Have enough comments at parts that are not easily understood, so that somebody
  else may continue your work in the future.
