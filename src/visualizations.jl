export interactive_trajectory, interactive_cobweb, interactive_orbitdiagram, scaleod,
  interactive_poincaresos_scan, interactive_poincaresos, interactive_trajectory_timeseries

"""
    interactive_trajectory(ds::DynamicalSystem [, u0s]; kwargs...) → fig, dsobs

Create a panel (a new Makie figure) that contains an axis showing the evolving trajectories
of `ds` while optionally allowing for interactively starting/stopping the time evolution
and/or changing the system parameters via sliders. Return the newly created figure `fig`
that may include interactivity-related buttons, and a `dsobs::DynamicalSystemObservable`
that facilities the creation of custom animations and/or interactive
applications, see the custom animations section below.

The given dynamical system is always cast into a `ParallelDynamicalSystem`.
The trajectories from the initial conditions in `u0s` (a vector of vectors)
are all evolved and visualized in parallel. By default
only the current state of the system is used.

## Interactivity and time stepping keywords

- `add_controls = true`: If `true`, below the main axis containing the trajectories
  some controls for animating the trajectories live are added: run, pause, and how many
  steps to evolve for per animation step. The plotted trajectories can always be evolved
  manually using the custom animations etup that we describe below; `add_controls` only
  concerns the buttons and interactivity added to the created panel.
- `parameter_sliders = nothing`: If given, it must be a dictionary, mapping
  parameter indices (of the parameter container of `ds`) into ranges. Each combination
  of index and range becomes a slider that can be interactively controlled to alter
  a system parameter on the fly during time evolution. Note that all parameter systems
  can always be altered on the fly using the custom animation setup that we describe below;
  `parameter_sliders` only conserns the buttons and interactivity added to the created panel.
- `pnames = Dict(keys(ps) .=> keys(ps))` : Dictionary mapping parameter keys to labels.
  Only valid if `parameter_sliders` is given.
- `Δt`: Default time step of time evolution. 1 for discrete time,
  0.01 for continuous time systems.
- `pause = nothing`: If given, it must be a real number. This number is given to the `sleep`
  function, which is called between each plot update. This functionality exists because
  in most cases the dynamical evolution is so fast that we need to actively pause before
  each update to visually track the trajectory. For interactive applications this should
  be `0.001` or some other small number.

## Visualization keywords

- `colors`: The color for each initial condition (and resulting trajectory).
- `idxs = 1:min(length(transform(u0s[1])), 3)`: Which variables to plot
  (up to three can be chosen).
* `tail = 1000`: Length of plotted trajectory (in step units of the integrator).
* `fade = true`: For continuous time system, the trajectories in state space are faded
  to full transparency if `true`.
- `markersize = 15`: Size of markers of trajectory endpoints. For discrete systems
  half of that is used for the trajectory tail.
* `plotkwargs = NamedTuple()`: A named tuple of keyword arguments propagated to
  the state space plot (`lines` for continuous, `scatter` for discrete systems).
  `plotkwargs` can also be a vector of named tuples, in which case each initial condition
  gets different arguments.

## Layouting keywords

- `idxs = 1:min(length(u0s[1]), 3)`: Which variables of `ds` should be plotted.
  If three indices are given, the plot is also 3D, otherwise 2D. Note that no
  transformation of the dynamical system is done, you can use a `ProjectedDynamicalSystem`
  if you want to visualize an arbitrary transformation of `ds`.
- `lims`: A tuple of tuples (min, max) for the axis limits. If not given, they are
  automatically deduced by evolving each of `u0s` 100 units and picking most extreme
  values (limits are _not_ adjusted by default during the live animations)
- `figure, axis`: both can be named tuples with arbitrary keywords
  propagated to the generation of the `Figure` and state space `Axis` instances.

## Custom animations

The second return argument `dsobs` is a `DynamicalSystemObservable`.
The trajectories plotted in the main panel are linked to observables that are fields
of the `dsobs`. Specifically, the field `dsobs.state_obserable` is an observable containing the
final state of each of the trajectories, i.e., a vector of vectors like `u0s`.
`dsobs.param_observable` is an observable of the system parameters. These observables
are triggered by the interactive GUI buttons (the first two when the system is stepped in
time, the last one when the parameters are updated). However, these observables,
and hence the corresponding plotted trajectories that are `map`ed from these observables,
can be updated via the formal API of `DynamicalSystem`:

```
step!(dsobs, n::Int = 1)
```
will step the system for `n` steps of `Δt` time,
and only update the plot on the last step. `set_parameter!(dsobs, index, value)` will
update the system parameter and then trigger the parameter observable.
Lastly, `set_state!(dsobs, new_u [, i])` will set the `i`-th system state and clear the
trajectory plot to the new initial condition.

This information can be used to create custom animations and/or interactive apps.
In principle, the only thing a user has to do is create new observables from
the existing ones using e.g. the `on` function and plot these new observables.
Various examples are provided in the online documentation.
"""
function interactive_trajectory end


"""
    interactive_trajectory_timeseries(ds::DynamicalSystem, fs, [, u0s]; kwargs...) → fig, dsobs

An extension to [ `interactive_trajectory`](@ref), which adds a new panel to the right
of the original trajectory panel. This new panel contains timeseries of various observed
quantities from the state of the dynamical system. These observed timeseries are given in
`fs`. It is a vector of integers and/or functions. Integers mean to simply observe the
specified variable. Functions are all functions of the full state of the dynamical system
that return a real number.

## Keyword arguments

- `linekwargs = NamedTuple()`: Extra keywords propagated to the timeseries plots.
  Can also be a vector of named tuples, each one for each unique initial condition.
- `timeseries_names`: A vector of strings with length equal to `fs` giving names to
  the y-labels of the timeseries plots.

All other keywords are propagated to [`interactive_trajectory`](@ref).
"""
function interactive_trajectory_timeseries end



"""
    interactive_cobweb(ds::DiscreteDynamicalSystem, prange, O::Int = 3; kwargs...)
Launch an interactive application for exploring cobweb diagrams of 1D discrete
dynamical systems. Two slides control the length of the plotted trajectory and the
current parameter value. The parameter values are obtained from the given `prange`.

In the cobweb plot, higher order iterates of the dynamic rule `f` are plotted as well,
starting from order 1 all the way to the given order `O`.
Both the trajectory in the cobweb, as well as any iterate `f` can be turned off by
using some of the buttons.

## Keywords
* `fkwargs = [(linewidth = 4.0, color = randomcolor()) for i in 1:O]`: plotting keywords
  for each of the plotted iterates of `f`
* `trajcolor = :black`: color of the trajectory
* `pname = "p"`: name of the parameter slider
* `pindex = 1`: parameter index
* `xmin = 0, xmax = 1`: limits the state of the dynamical system can take
* `Tmax = 1000`: maximum trajectory length
* `x0s = range(xmin, xmax; length = 101)`: Possible values for the x0 slider.
"""
function interactive_cobweb end



"""
    interactive_orbitdiagram(
        ds::DynamicalSystem, p_index, pmin, pmax, i::Int = 1;
        u0 = nothing, parname = "p", title = ""
    )

Open an interactive application for exploring orbit diagrams (ODs) of discrete
dynamical systems. Requires `DynamicalSystems`.

In essense, the function presents the output of `orbitdiagram`
of the `i`th variable of the `ds`, and allows interactively zooming into it.

Keywords control the name of the parameter, the initial state (used for _any_ parameter)
or whether to add a title above the orbit diagram.

## Interaction

The application is separated in the "OD plot" (left) and the "control panel" (right).
On the OD plot you can interactively click
and drag with the left mouse button to select a region in the OD. This region is then
**re-computed** at a higher resolution.

The options at the control panel are straight-forward, with
* `n` amount of steps recorded for the orbit diagram (not all are in the zoomed region!)
* `t` transient steps before starting to record steps
* `d` density of x-axis (the parameter axis)
* `α` alpha value for the plotted points.

Notice that at each update `n*t*d` steps are taken.
You have to press `update` after changing these parameters.
Press `reset` to bring the OD in the original
state (and variable). Pressing `back` will go back through the history of your exploration
History is stored when the "update" button is pressed or a region is zoomed in.

You can even decide which variable to get the OD for
by choosing one of the variables from the wheel!
Because the y-axis limits can't be known when changing variable, they reset to the size
of the selected variable.

## Accessing the data

What is plotted on the application window is a _true_ orbit diagram, not a plotting
shorthand. This means that all data are obtainable and usable directly.
Internally we always scale the orbit diagram to [0,1]² (to allow `Float64` precision
even though plotting is `Float32`-based). This however means that it is
necessary to transform the data in real scale. This is done through the function
[`scaleod`](@ref) which accepts the 5 arguments returned from the current function:
```julia
figure, oddata = interactive_orbitdiagram(...)
ps, us = scaleod(oddata)
```
"""
function interactive_orbitdiagram end

"""
    scaleod(oddata) -> ps, us
Given the return values of [`interactive_orbitdiagram`](@ref), produce
orbit diagram data scaled correctly in data units. Return the data as a vector of
parameter values and a vector of corresponding variable values.
"""
function scaleod end



"""
    interactive_poincaresos_scan(A::StateSpaceSet, j::Int; kwargs...)
    interactive_poincaresos_scan(As::Vector{StateSpaceSet}, j::Int; kwargs...)

Launch an interactive application for scanning a Poincare surface of section of `A`
like a "brain scan", where the plane that defines the section can be arbitrarily
moved around via a slider. Return `figure, ax3D, ax2D`.

The input dataset must be 3 dimensional, and here the crossing plane is always
chosen to be when the `j`-th variable of the dataset crosses a predefined value.
The slider automatically gets all possible values the `j`-th variable can obtain.

If given multiple datasets, the keyword `colors` attributes a color to each one, e.g.
`colors = [JULIADYNAMICS_COLORS[mod1(i, 6)] for i in 1:length(As)]`.

The keywords `linekw, scatterkw` are named tuples that are propagated as keyword arguments
to the line and scatter plot respectively, while the keyword `direction = -1` is propagated
to the function `DyamicalSystems.poincaresos`.
"""
function interactive_poincaresos_scan end



"""
    interactive_poincaresos(cds, plane, idxs, complete; kwargs...)
Launch an interactive application for exploring a Poincaré surface of section (PSOS)
of the continuous dynamical system `cds`. Requires `DynamicalSystems`.

The `plane` can only be the `Tuple` type accepted by `DynamicalSystems.poincaresos`,
i.e. `(i, r)` for the `i`th variable crossing the value `r`. `idxs` gives the two
indices of the variables to be displayed, since the PSOS plot is always a 2D scatterplot.
I.e. `idxs = (1, 2)` will plot the 1st versus 2nd variable of the PSOS. It follows
that `plane[1] ∉ idxs` must be true.

`complete` is a three-argument **function** that completes the new initial state
during interactive use, see below.

The function returns: `figure, laststate` with the latter being
an observable containing the latest initial `state`.

## Keyword Arguments
* `direction, rootkw` : Same use as in `DynamicalSystems.poincaresos`.
* `tfinal = (1000.0, 10.0^4)` : A 2-element tuple for the range of values
  for the total integration time (chosen interactively).
* `color` : A **function** of the system's initial condition, that returns a color to
  plot the new points with. The color must be `RGBf/RGBAf`.
   A random color is chosen by default.
* `labels = ("u₁" , "u₂")` : Scatter plot labels.
* `scatterkwargs = ()`: Named tuple of keywords passed to `scatter`.
* `diffeq = NamedTuple()` : Any extra keyword arguments are passed into `init` of DiffEq.

## Interaction
The application is a standard scatterplot, which shows the PSOS of the system,
initially using the system's `u0`. Two sliders control the total evolution time
and the size of the marker points (which is always in pixels).

Upon clicking within the bounds of the scatter plot your click is transformed into
a new initial condition, which is further evolved and its PSOS is computed and then
plotted into the scatter plot.

Your click is transformed into a full `D`-dimensional initial condition through
the function `complete`. The first two arguments of the function are the positions
of the click on the PSOS. The third argument is the value of the variable the PSOS
is defined on. To be more exact, this is how the function is called:
```julia
x, y = mouseclick; z = plane[2]
newstate = complete(x, y, z)
```
The `complete` function can throw an error for ill-conditioned `x, y, z`.
This will be properly handled instead of breaking the application.
This `newstate` is also given to the function `color` that
gets a new color for the new points.
"""
function interactive_poincaresos end