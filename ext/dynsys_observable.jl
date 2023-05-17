using DynamicalSystems.DataStructures
# This is a very light struct that contains the trajectory end state
# and the trajectory tail. Note that for convenience it always contains the state
# as a vector of observables, each observable containing each of the
# parallel states of the dynamical system
"""
    DynamicalSystemObservable

Struct containing a reference to a `ParallelDynamicalSystem`, as well as various observables
that keep track of the system state(s) and the parameters. The observables are linked
with the panel generated by `interactive_trajectory_panel` and allow setting up custom
animations and/or interactive applications. See `interactive_trajectory_panel` for more.
"""
struct DynamicalSystemObservable
    pds::ParallelDynamicalSystem # reference to the dynamical system
    state_observable::Observable
    tail_observables::Vector{Observable}
    param_observable::Observable
    idxs::AbstractVector{Int}
    Δt::Real # a default value for `step!`
end

# TODO: Somehow extract this from an online repo...?
COLORS = [
    "#7143E0",
    "#191E44",
    "#0A9A84",
    "#AF9327",
    "#791457",
    "#6C768C",
]


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
function DynamicalSystems.interactive_trajectory_panel(
        ds::DynamicalSystems.DynamicalSystem, u0s = [DynamicalSystems.current_state(ds)];
        # Selection of what to plot
        idxs = 1:min(length(u0s[1]), 3),
        # Time evolution
        tail = 1000,
        Δt = DynamicalSystems.isdiscretetime(ds) ? 1 : 0.01,
        force_non_adaptive = true,
        pause = nothing,
        # Visualization
        colors = [COLORS[mod1(i, 6)] for i in 1:length(u0s)],
        plotkwargs = NamedTuple(), markersize = 15,
        linekwargs = DynamicalSystems.isdiscretetime(ds)  ? () : (linewidth = 4,),
        fade = true,
        # parameters
        ps = nothing,
        pnames = isnothing(ps) ? nothing : Dict(keys(ps) .=> keys(ps)),
        add_controls = true,
        # figure and axis
        figure = NamedTuple(),
        axis = NamedTuple(),
        lims = nothing,
    )

    if ds isa CoupledODEs && force_non_adaptive
        newdiffeq = (ds.diffeq..., adaptive = false, dt = Δt)
        ds = CoupledODEs(ds, newdiffeq)
    end

    pds = DynamicalSystems.ParallelDynamicalSystem(ds, u0s)
    fig = Figure(; figure...)
    # Set up trajectrory plot
    statespacelayout = fig[1,1] = GridLayout()
    lims = isnothing(lims) ? _traj_lim_estimator(ds, u0s, idxs) : lims
    statespaceax, tailobs, finalpoints = _init_statespace_plot!(statespacelayout, ds, idxs,
        lims, pds, colors, plotkwargs, markersize, tail, axis, fade,
    )
    # Set up layouting and add controls
    if add_controls # Notice that `run` and `step` are already observables
        run, step, stepslider = _trajectory_plot_controls!(statespacelayout)
        display(fig) # attemp to display by default in interactive scenarios
    else
        # So that we can leave the interactive UI code as is
        run = Observable(0); step = Observable(0); stepslider = Observable(1)
    end

    # Create the dynamical system observable now with these linked
    po = Observable(deepcopy(current_parameters(ds)))
    dso = DynamicalSystemObservable(pds, finalpoints, tailobs, po, SVector(idxs...), Δt)

    # Functionality of live evolution. This links all observables with triggers.
    # The run button just triggers the step button perpetually
    isrunning = Observable(false)
    on(run) do c; isrunning[] = !isrunning[]; end
    on(run) do c
        @async while isrunning[]
            step[] = step[] + 1
            isopen(fig.scene) || break # ensures computations stop if closed window
            isnothing(pause) ? yield() : sleep(pause)
        end
    end
    # while the actual observables triggering happens from the step button
    on(step) do clicks
        n = stepslider[]
        # which of course calls the stepping function on the observable
        step!(dso, n)
    end

    # Live parameter changing
    # note here `ps` are parameters to have a slider; all parameters
    # can be changed after creation of `dso` via `set_parameter!`
    if !isnothing(ps)
        paramlayout = fig[2, :] = GridLayout(tellheight = true, tellwidth = false)
        slidervals = _add_ds_param_controls!(paramlayout, ps, pnames)
        update = Button(fig, label = "update", tellwidth = false, tellheight = true)
        paramlayout[2, 1] = update
        on(update.clicks) do clicks
            for l in keys(slidervals)
                v = slidervals[l][]
                set_parameter!(dso, l, v)
            end
        end
    end

    return fig, dso
end


# Main panels of animation
"Create the state space axis and evolution controls. Return the axis."
function _init_statespace_plot!(
        layout, ds, idxs, lims, pds, colors, plotkwargs, markersize, tail, axis, fade,
    )
    tailobs, finalpoints = _init_trajectory_observables(pds, tail, idxs)
    is3D = length(idxs) == 3
    statespaceax = !is3D ? Axis(layout[1,1]; xlabel = "x1", ylabel = "x2", axis...) :
        Axis3(layout[1,1]; xlabel = "x1", ylabel = "x2", zlabel = "x3", axis...)

    # Initialize trajectories plotted element
    for (i, ob) in enumerate(tailobs)
        pk = plotkwargs isa Vector ? plotkwargs[i] : plotkwargs
        x = to_color(colors[i])
        if fade
            x = [RGBAf(x.r, x.g, x.b, i/tail) for i in 1:tail]
        end
        if !DynamicalSystems.isdiscretetime(ds)
            Makie.lines!(statespaceax, ob;
                color = x, linewidth = 3.0, transparency = true, pk...
            )
        else
            Makie.scatter!(statespaceax, ob; color = x,
                markersize = markersize/2, transparency = true, strokewidth = 0.0, pk...
            )
        end
    end
    finalargs = if !DynamicalSystems.isdiscretetime(ds)
        (marker = :circle, )
    else
        (marker = :diamond, )
    end
    Makie.scatter!(statespaceax, finalpoints;
        color = colors, markersize = markersize, finalargs...)
    !isnothing(lims) && (statespaceax.limits = lims)
    is3D && (statespaceax.protrusions = 50) # removes overlap of labels
    return statespaceax, tailobs, finalpoints
end
function _init_trajectory_observables(pds, tail, idxs)
    N = length(DynamicalSystems.current_states(pds))
    tailobs = Observable[]
    T = length(idxs) == 2 ? Point2f : Point3f
    for i in 1:N
        cb = CircularBuffer{T}(tail)
        fill!(cb, T(DynamicalSystems.current_state(pds, i)[idxs]))
        push!(tailobs, Observable(cb))
    end
    finalpoints = Observable([x[][end] for x in tailobs])
    return tailobs, finalpoints
end
function _trajectory_plot_controls!(layout)
    layout[2, 1] = controllayout = GridLayout(tellwidth = false)
    run = Button(controllayout[1, 1]; label = "run")
    step = Button(controllayout[1, 2]; label = "step")
    slider_vals = vcat(1:10, 100:100:1000)
    sg = SliderGrid(controllayout[1,3],
        (label = "steps =", range = slider_vals, startvalue = 1),
    )
    return run.clicks, step.clicks, sg.sliders[1].value
end
function _traj_lim_estimator(ds, u0s, idxs)
    ds = deepcopy(ds)
    Δt = DynamicalSystems.isdiscretetime(ds) ? 1 : 0.1
    tr = DynamicalSystems.trajectory(ds, 100, u0s[1]; Δt)[1]
    _mi, _ma = DynamicalSystems.minmaxima(tr)
    mi, ma = _mi[idxs], _ma[idxs]
    for i in 2:length(u0s)
        tr = DynamicalSystems.trajectory(ds, 100, u0s[i]; Δt)[1]
        _mii, _maa = DynamicalSystems.minmaxima(tr)
        mii, maa = _mii[idxs], _maa[idxs]
        mi = min.(mii, mi)
        ma = max.(maa, ma)
    end
    # Alright, now we just have to put them into limits and increase a bit
    mi = mi .- 0.1mi
    ma = ma .+ 0.1ma
    lims = [(mi[i], ma[i]) for i in 1:length(idxs)]
    lims = (lims...,)
end

# TODO: also add extension for set_state!

# stepping code
function DynamicalSystems.step!(dso::DynamicalSystemObservable, n::Int, Δt = dso.Δt)
    N = length(dso.tail_observables)
    # Always store values, but only update observables after loop
    for _ in 1:n
        step!(dso.pds, Δt, true)
        for i in 1:N
            ob = dso.tail_observables[i]
            last_state = current_state(dso.pds, i)[dso.idxs]
            push!(ob[], last_state)
        end
    end
    # Here the observables are updated with their current values
    notify.(dso.tail_observables)
    dso.state_observable[] = [x[][end] for x in dso.tail_observables]
    return nothing
end

# Parameter handling
function _add_ds_param_controls!(paramlayout, ps, pnames)
    slidervals = Dict{keytype(ps), Observable}() # directly has the slider observables
    tuples_for_slidergrid = []
    for (i, (l, vals)) in enumerate(ps)
        startvalue = p0[l]
        label = string(pnames[l])
        push!(tuples_for_slidergrid, (;label, range = vals, startvalue))
    end
    sg = SliderGrid(paramlayout[1,1], tuples_for_slidergrid...; tellheight = true)
    for (i, (l, vals)) in enumerate(ps)
        slidervals[l] = sg.sliders[i].value
    end
    return slidervals
end


function DynamicalSystems.set_parameter!(dso::DynamicalSystemObservable, index, value)
    dso.param_observable[][index] = value
    set_parameter!(dso.pds, index, value)
    notify(dso.param_observable)
    return
end
