export interactive_trajectory_panel

"""
    interactive_trajectory_panel(ds::DynamicalSystem [, u0s]; kwargs...) → fig, dso

Create a panel (a new Makie figure) that contains an axis showing the evolving trajectories
of `ds` while optionally allowing for interactively starting/stopping the time evolution
and/or changing the system parameters via sliders. Return the newly created figure `fig`
that may include interactivity-related buttons, and a `dso::DynamicalSystemObservable`
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
- `force_non_adaptive = true`: Because trajectories are not pre-computed and
  interpolated, but rather calculated on the fly step by step, the integrator of continuous
  time systems is casted into its non-adaptive version if this keyword is `true`.
  This results to smoother plotted curves.

## Visualization keywords

- `colors`: The color for each initial condition (and resulting trajectory).
- `idxs = 1:min(length(transform(u0s[1])), 3)`: Which variables to plot
  (up to three can be chosen).
* `tail = 1000`: Length of plotted trajectory (in step units of the integrator).
* `fade = true`: For continuous time system, the trajectories in state space are faded
  to full transparency if `true`.
- `markersize = 15`: Size of markers of trajectory endpoints. For discrete systems
  half of that is used for the trajectory tail.
* `plotkwargs = NamedTuple()` : A named tuple of keyword arguments propagated to
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

The second return argument `dso` is a `DynamicalSystemObservable`.
The trajectories plotted in the main panel are linked to observables that are fields
of the `dso`. Specifically, the field `dso.state_obserable` is an observable containing the
final state of each of the trajectories, i.e., a vector of vectors like `u0s`.
`dso.param_observable` is an observable of the system parameters. These observables
are triggered by the interactive GUI buttons (the first two when the system is stepped in
time, the last one when the parameters are updated). However, these observables,
and hence the corresponding plotted trajectories that are `map`ed from these observables,
can be updated via the formal API of `DynamicalSystem`.
`step!(dso, n::Int, [, Δt])` will step the system for `n` steps of `Δt` time,
and only update the plot on the last step. `set_parameter!(dso, index, value)` will
update the system parameter and then trigger the parameter observable.
Lastly, `set_state!(dso, new_u [, i])` will set the `i`-th system state and clear the
trajectory plot to the new initial condition.

This information can be used to create custom animations and/or interactive apps.
In principle, the only thing a user has to do is create new observables from
the existing ones using e.g. the `on` function and plot these new observables.
Various examples are provided in the online documentation.
"""
function interactive_trajectory_panel end