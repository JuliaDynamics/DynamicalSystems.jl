using DataStructures
export interactive_evolution

########################################################################################
# Dynamical system stepping structure
########################################################################################
# TODO:
struct DynamicalSystemAnimator
    pds
    states::Observable{Vector} # vector of vectors
    slidervals
    Δt::Real
end

########################################################################################
# GUI launcing function
########################################################################################
# TODO: Allow plotted timeseries to be arbitrary functions of state

"""
    interactive_evolution(ds::DynamicalSystem [, u0s]; kwargs...) → fig, obs, step

Launch an interactive GUI application that can evolve the initial conditions `u0s`
(vector of vectors) of the given dynamical system.
All initial conditions are evolved in parallel and at exactly the same time
using a `ParallelDynamicalSystem` of `ds`.

Added controls allow you to step/run/pause the evolution and to control after
how many integrator steps the plots are updated.
The application can run forever (trajectories are computed on demand).

By default the GUI window displays statespace and timeseries plots.
It also allows changing the parameters of `ds` live during the
system evolution, see keyword `ps` below in "Parameter Keywords".

The function returns `fig, obs, step, paramvals`. `fig` is the overarching figure
(the entire GUI) and can be recorded with `Makie.record`.
`obs` is a vector of observables, each containing the current state of each trajectory.
`step` is an observable that triggers the stepping of each trajectory and the update
of the plots. Do `step[] = 0` (or any other integer), to trigger an update.
`paramvals` is an observable containing current parameter values, and is only valid
if `ps` is not nothing, see keyword `ps` below in "Parameter Keywords".

The figure layout is as follows:
1. `fig[1,1]` = state space plot (`fig[1,1][1,1]`) and time evolution controls
2. `fig[1,2]` = timeseries plots
3. `fig[2,:]` = parameter controls (if `ps` is given)

This means that you can make any kind of composite plots and videos starting from
the figure returned from `interactive_evolution`. See the documentation online for
more such examples.

## State Space Keywords

* `transform = identity`: Transformation applied to the state of the dynamical system
  before plotting. Can even return a vector that is of higher dimension than `ds`.
* `idxs = 1:min(length(transform(u0s[1])), 3)`: Which variables to plot
  (up to three can be chosen).
  Variables are selected after `transform` has been applied.
* `colors`: The color for each initial condition (and resulting trajectory).
* `lims`: A tuple of tuples (min, max) for the axis limits. If not given, they are
  automatically deduced by evolving each of `u0s` 1000 units and picking most extreme
  values (limits are very hard to adjust after application is launched).
* `m = 1.0`: The trajectory endpoints have a marker. A heuristic is done to choose
  appropriate marker size given the trajectory size. `m` is a multiplier that scales
  the marker size.
* `tail = 1000`: Length of plotted trajectory (in step units of the integrator).
* `fade = true`: For continuous time system, the trajectories in state space are faded
  to full transparency if `true`.
* `plotkwargs = NamedTuple()` : A named tuple of keyword arguments propagated to
  the state space plot (`lines` for continuous, `scatter` for discrete systems).
  `plotkwargs` can also be a vector of named tuples, in which case each initial condition
  gets different arguments.
* `Δt`: Time step of time evolution. 1 for discrete time, 0.01 for continuous time systems.
  Because trajectories are not pre-computed and interpolated,
  but rather calculated on the fly step by step, a constant step size equal to `Δt`
  is enforced internally for continuous time systems.
* `add_controls = true`: Whether to add buttons and sliders for interactively
  controlling the trajectory evolution. Should be `false` only if composite
  videos are intended to be produced using the returned `step`. If `false`, the keyword
  `steps_per_update = 1` decides how many steps to take before updating plots.

## Timeseries Keywords

* `tsidxs = idxs`: Indices selecting variables to be plotted as timeseries. You can
  pass `nothing` instead and no timeseries will be plotted.
* `total_span`: How much the x-axis of the timeseries plots should span (in real time units)
* `linekwargs = NamedTuple()`: Extra keywords propagated to the timeseries plots.

## Parameter Keywords

* `ps = nothing`: If `ps` is not nothing, then it must be a dictionary, mapping keys
  of the system parameter container (`ds.p`) to possible ranges of values. The app then will
  add some additional controls on the bottom of the GUI which allow one to interactively change
  system parameters and then click the "update" button to translate the new parameters to
  system evolution. This can be done without stopping the live system evolution.
  Notice that in this scenario it is recommended to provide the `lims` keyword manually.
  An extra argument is returned in this case: a dictionary mapping parameter keys
  to _observables_ containing their current values. You can use this to generate additional
  plot elements that may depend on system parameters and thus need to be changed
  if the sliders are changed.
* `pnames = Dict(keys(ps) .=> keys(ps))` : Dictionary mapping parameter keys to labels.
  Only valid if `params` is a dictionary and not `nothing`.

In addition the keywords `figure, axis` can be named tuples with arbitrary keywords
propagated to the generation of the `Figure` and state space `Axis` instances.
"""
function interactive_evolution(
        ds::DynamicalSystems.DynamicalSystem, u0s = [DynamicalSystems.current_state(ds)];
        transform = identity, idxs = 1:min(length(transform(u0s[1])), 3), tsidxs = idxs,
        colors = [CYCLIC_COLORS[i] for i in 1:length(u0s)], tail = 1000,
        lims = nothing, Δt = DynamicalSystems.isdiscretetime(ds) ? 1 : 0.01,
        plotkwargs = NamedTuple(), m = 1.0,
        total_span = DynamicalSystems.isdiscretetime(ds) ? 50 : 10,
        linekwargs = DynamicalSystems.isdiscretetime(ds)  ? () : (linewidth = 4,),
        ps = nothing,
        pnames = isnothing(ps) ? nothing : Dict(keys(ps) .=> keys(ps)),
        add_controls = true, steps_per_update = 1,
        figure = (resolution = (isnothing(tsidxs) ? 800 : 1600, 800), ),
        axis = NamedTuple(),
        fade = true,
    )

    N = length(u0s)
    @assert length(colors) ≥ length(u0s) "You need to provide enough colors!"
    fig = Figure(; figure...)

    # Setup plots and integrator stuff
    pds = DynamicalSystems.ParallelDynamicalSystem(ds, u0s)
    statespacelayout = fig[1,1] = GridLayout()
    @assert length(idxs) ≤ 3 "Only up to three variables can be plotted!"
    idxs = DynamicalSystems.SVector(idxs...)
    lims = isnothing(lims) ? traj_lim_estimator(ds, u0s, idxs, transform) : lims
    statespaceax, obs, finalpoints = _init_statespace_plot!(statespacelayout, ds, idxs,
        lims, pds, colors, plotkwargs, m, tail, transform, axis, fade,
    )
    if !isnothing(tsidxs)
        timeserieslayout = fig[1,2] = GridLayout()
        allts, ts_axes = _init_timeseries_plots!(
            timeserieslayout, pds, tsidxs, colors, linekwargs, transform, tail, lims,
        )
        update_ts = true
    else
        update_ts = false
    end
    if !isnothing(ps)
        paramlayout = fig[2, :] = GridLayout(tellheight = true, tellwidth = false)
    end
    if add_controls # Notice that `run` and `step` are already observables
        run, step, stepslider = _trajectory_plot_controls!(statespacelayout)
        display(fig) # attemp to display by default in interactive scenarios
    else
        run = Observable(0); step = Observable(0); stepslider = Observable(steps_per_update)
    end


    # Functionality of live evolution. This links all observables with triggers.
    # The run button just triggers the step button perpetually
    isrunning = Observable(false)
    on(run) do c; isrunning[] = !isrunning[]; end
    on(run) do c
        @async while isrunning[]
            step[] = step[] + 1
            isopen(fig.scene) || break # ensures computations stop if closed window
            yield()
        end
    end
    # while the actual observables triggering happens from the step button
    on(step) do clicks
        n = stepslider[]
        # Always store values, but only update observables after loop
        for _ in 1:n
            DynamicalSystems.step!(pds, Δt, true)
            for i in 1:N
                ob = obs[i]
                last_state = transform(DynamicalSystems.current_state(pds, i))[idxs]
                push!(ob[], last_state)
                if update_ts
                    for k in 1:length(tsidxs)
                        push!(allts[k][i][], Point2f(DynamicalSystems.current_time(pds), last_state[tsidxs[k]]))
                    end
                end
            end
        end
        # Here the observables are updated with their current values
        notify.(obs)
        finalpoints[] = [x[][end] for x in obs]
        if update_ts
            for k in 1:length(tsidxs); notify.(allts[k]); end
            xlims!(ts_axes[end], max(0, DynamicalSystems.current_time(pds) - total_span), max(DynamicalSystems.current_time(pds), total_span))
        end
    end

    # Live parameter changing
    if !isnothing(ps)
        slidervals, paramvals = _add_ds_param_controls!(paramlayout, ps, DynamicalSystems.initial_parameters(ds), pnames)
        update = Button(fig, label = "update", tellwidth = false, tellheight = true)
        paramlayout[2, 1] = update
        on(update.clicks) do clicks
            _update_ds_parameters!(ds, slidervals, paramvals)
        end
    else
        paramvals = nothing
    end

    return fig, obs, step, paramvals
end




"Create the state space axis and evolution controls. Return the axis."
function _init_statespace_plot!(
        layout, ds, idxs, lims, pds, colors, plotkwargs, m, tail, transform, axis, fade,
    )
    obs, finalpoints = init_trajectory_observables(pds, tail, idxs, transform)
    is3D = length(idxs) == 3
    markersize = 15
    statespaceax = !is3D ? Axis(layout[1,1]; xlabel = "x1", ylabel = "x2", axis...) :
        Axis3(layout[1,1]; xlabel = "x1", ylabel = "x2", zlabel = "x3", axis...)

    # Initialize trajectories plotted element
    for (i, ob) in enumerate(obs)
        pk = plotkwargs isa Vector ? plotkwargs[i] : plotkwargs
        if !DynamicalSystems.isdiscretetime(ds)
            x = to_color(colors[i])
            if fade
                x = [RGBAf(x.r, x.g, x.b, i/tail) for i in 1:tail]
            end
            Makie.lines!(statespaceax, ob;
                color = x, linewidth = 2.0, transparency = true, pk...
            )
        else
            Makie.scatter!(statespaceax, ob; color = colors[i],
                markersize = markersize/2, strokewidth = 0.0, pk...
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

    return statespaceax, obs, finalpoints
end
function init_trajectory_observables(pds, tail, idxs, transform)
    N = length(DynamicalSystems.current_states(pds))
    obs = Observable[]
    T = length(idxs) == 2 ? Point2f : Point3f
    for i in 1:N
        cb = CircularBuffer{T}(tail)
        fill!(cb, T(transform(DynamicalSystems.current_state(pds, i))[idxs]))
        push!(obs, Observable(cb))
    end
    finalpoints = Observable([x[][end] for x in obs])
    return obs, finalpoints
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



function _init_timeseries_plots!(
        layout, pds, idxs, colors, linekwargs, transform, tail, lims
    )
    N = length(DynamicalSystems.current_states(pds))
    # Initialize timeseries data:
    allts = [] # each entry is a vector of observables; the contained observables
    # correspond to the timeseries of a given axis. So `length(ts)` == amount of axis.
    # However, `length(allts[i])` == amount of initial conditions.
    for i in 1:length(idxs)
        individual_ts = Observable[]
        for j in 1:N
            cb = CircularBuffer{Point2f}(tail)
            fill!(cb, Point2f(
                DynamicalSystems.current_time(pds), transform(DynamicalSystems.current_state(pds, j))[idxs][i])
            )
            push!(individual_ts, Observable(cb))
        end
        push!(allts, individual_ts)
    end
    # Initialize timeseries axis and plots:
    ts_axes = []
    for i in 1:length(idxs)
        ax = Axis(layout[i, 1]; xticks = LinearTicks(5))
        push!(ts_axes, ax)
        individual_ts = allts[i]
        for j in 1:N
            lines!(ax, individual_ts[j]; color = colors[j], linekwargs...)
            if DynamicalSystems.isdiscretetime(pds)
                scatter!(ax, individual_ts[j]; color = colors[j])
            end
        end
        tight_xticklabel_spacing!(ax)
        ax.ylabel = "x$i"
        ylims!(ax, lims[i])
    end
    linkxaxes!(ts_axes...)
    for i in 1:length(idxs)-1
        hidexdecorations!(ts_axes[i]; grid = false)
    end
    return allts, ts_axes
end

function traj_lim_estimator(ds, u0s, idxs, transform)
    Δt = DynamicalSystems.isdiscretetime(ds) ? 1 : 0.1
    _tr = DynamicalSystems.trajectory(ds, 1000, u0s[1]; Δt)[1]
    tr = DynamicalSystems.StateSpaceSet(transform.(_tr.data))
    _mi, _ma = DynamicalSystems.minmaxima(tr)
    mi, ma = _mi[idxs], _ma[idxs]
    for i in 2:length(u0s)
        _tr = DynamicalSystems.trajectory(ds, 1000, u0s[i]; Δt)[1]
        tr = DynamicalSystems.StateSpaceSet(transform.(_tr.data))
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




function _add_ds_param_controls!(paramlayout, ps, p0, pnames)
    # fig = paramlayout.parent.parent
    slidervals = Dict{keytype(ps), Observable}() # directly has the slider observables
    paramvals = Dict{keytype(ps), Observable}()  # will only get updated on button
    tuples_for_slidergrid = []
    for (i, (l, vals)) in enumerate(ps)
        startvalue = p0[l]
        label = string(pnames[l])
        # old code:
        # sll = labelslider!(fig, label, vals; sliderkw = Dict(:startvalue => startvalue))
        # slidervals[l] = sll.slider.value # directly add the observable
        # paramvals[l] = Observable(sll.slider.value[]) # will only get updated on button
        # paramlayout[i, :] = sll.layout
        # new code:
        push!(tuples_for_slidergrid, (;label, range = vals, startvalue))
    end
    sg = SliderGrid(paramlayout[1,1], tuples_for_slidergrid...; tellheight = true)
    for (i, (l, vals)) in enumerate(ps)
        slidervals[l] = sg.sliders[i].value
        paramvals[l] = Observable(sg.sliders[i].value[])
    end
    return slidervals, paramvals
end

function _update_ds_parameters!(ds, slidervals, paramvals)
    for l in keys(slidervals)
        v = slidervals[l][]
        paramvals[l][] = v
        DynamicalSystems.set_parameter!(ds, l, v)
    end
end
