**DynamicalSystems.jl** offers some basic methods for visualizing chaotic systems in
the form of the functions described in this documentation page.

All plotting is done through the [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl)
package. However, this is not a dependency of DynamicalSystems.jl. Instead, all
functions described here are brought into scope as soon as the user executes
`using PyPlot`, which works regardless if `PyPlot` module was loaded before or
after `DynamicalSystems`. This is possible through the [Requires.jl](https://github.com/MikeInnes/Requires.jl) package.

Use the help mode (press `?` and then the function name) to access the documentation
strings for e.g. using keyword arguments.

## Visualization Library
* `phasespace` : Plots the phasespace of a 2D `DiscreteDynamicalSystem`.
* `plot_linear_regions` : Plots the results of `lnear_regions`.
* `plot_dataset` : Plots each column of a dataset.
* `plot_trajectory` : Produces a trajectory and plots each column.
