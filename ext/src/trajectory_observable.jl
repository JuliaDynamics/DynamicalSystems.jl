using DynamicalSystems.DataStructures

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
    current_step::Observable{Int}
    Δt::Real # a default value for `step!`
end

function DynamicalSystems.interactive_trajectory(
        ds::DynamicalSystems.DynamicalSystem, u0s = [DynamicalSystems.current_state(ds)];
        # Selection of what to plot
        idxs = 1:min(length(u0s[1]), 3),
        # Time evolution
        tail = 1000,
        Δt = DynamicalSystems.isdiscretetime(ds) ? 1 : 0.01,
        pause = nothing,
        # Visualization
        colors = [COLORS[i] for i in 1:length(u0s)],
        plotkwargs = NamedTuple(), markersize = 15,
        fade = true,
        # parameters
        parameter_sliders = nothing,
        pnames = isnothing(parameter_sliders) ? nothing : Dict(keys(parameter_sliders) .=> keys(parameter_sliders)),
        add_controls = true,
        # figure and axis
        figure = (resolution = (800, 800),),
        axis = NamedTuple(),
        lims = nothing,
    )

    if ds isa CoupledODEs # force time evolution into non-adaptive
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
    else
        # So that we can leave the interactive UI code as is
        run = Observable(0); step = Observable(0); stepslider = Observable(1)
    end

    # Create the dynamical system observable now with these linked
    po = Observable(deepcopy(current_parameters(ds)))
    dso = DynamicalSystemObservable(pds, finalpoints, tailobs, po, Observable(0), Δt)

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
    # note here `parameter_sliders` are parameters to have a slider; all parameters
    # can be changed after creation of `dso` via `set_parameter!`
    if !isnothing(parameter_sliders)
        paramlayout = fig[2, :] = GridLayout(tellheight = true, tellwidth = false)
        slidervals = _add_ds_param_controls!(paramlayout, parameter_sliders, pnames, current_parameters(ds))
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
    tailobs, finalpoints = _init_trajectory_observables(pds, tail)
    is3D = length(idxs) == 3
    statespaceax = !is3D ? Axis(layout[1,1]; xlabel = "x1", ylabel = "x2", axis...) :
        Axis3(layout[1,1]; xlabel = "x1", ylabel = "x2", zlabel = "x3", axis...)

    # Here we make two more observables for the plotted tails and plotted final
    # states, so that the stored observables in `dsobs` are the full system state;
    # This simplifies drastically making custom animations
    T = length(idxs) == 2 ? Point2f : Point3f
    plotted_tailobs = [
        map(x -> T[y[idxs] for y in x], ob) for ob in tailobs
    ]
    plotted_finalpoints = map(x -> T[y[idxs] for y in x], finalpoints)

    # Initialize trajectories plotted element
    for (i, ob) in enumerate(plotted_tailobs)
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
    Makie.scatter!(statespaceax, plotted_finalpoints;
        color = colors, markersize = markersize, finalargs...)
    !isnothing(lims) && (statespaceax.limits = lims)
    is3D && (statespaceax.protrusions = 50) # removes overlap of labels
    return statespaceax, tailobs, finalpoints
end
function _init_trajectory_observables(pds, tail)
    N = length(current_states(pds))
    tailobs = Observable[]
    T = typeof(current_state(pds))
    for i in 1:N
        cb = CircularBuffer{T}(tail)
        fill!(cb, current_state(pds, i))
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

# Parameter handling
function _add_ds_param_controls!(paramlayout, parameter_sliders, pnames, p0)
    slidervals = Dict{keytype(parameter_sliders), Observable}() # directly has the slider observables
    tuples_for_slidergrid = []
    for (i, (l, vals)) in enumerate(parameter_sliders)
        startvalue = p0[l]
        label = string(pnames[l])
        push!(tuples_for_slidergrid, (;label, range = vals, startvalue))
    end
    sg = SliderGrid(paramlayout[1,1], tuples_for_slidergrid...; tellheight = true)
    for (i, (l, vals)) in enumerate(parameter_sliders)
        slidervals[l] = sg.sliders[i].value
    end
    return slidervals
end

###########################################################################################
# Extension of DynamicalSystems API
###########################################################################################
function DynamicalSystems.step!(dso::DynamicalSystemObservable, n::Int = 1)
    Δt = dso.Δt
    N = length(dso.tail_observables)
    # Always store values, but only update observables after loop
    for _ in 1:n
        step!(dso.pds, Δt)
        for i in 1:N
            ob = dso.tail_observables[i]
            last_state = current_state(dso.pds, i)
            push!(ob[], last_state)
        end
    end
    dso.current_step.val = dso.current_step[] + n
    # Here the observables are updated with their current values
    notify.(dso.tail_observables)
    dso.state_observable[] = [x[][end] for x in dso.tail_observables]
    return nothing
end

function DynamicalSystems.set_state!(dso::DynamicalSystemObservable, u, i::Int = 1)
    dso.current_step.val = 0
    set_state!(dso.pds, u, i)
    fill!(dso.tail_observables[i][], u)
    notify(dso.tail_obsrvables[i])
    dso.state_observable[].val[i] = u
    notify(dso.state_observable)
    return nothing
end

function DynamicalSystems.set_parameter!(dso::DynamicalSystemObservable, index, value)
    dso.param_observable[][index] = value
    set_parameter!(dso.pds, index, value)
    notify(dso.param_observable)
    return
end

###########################################################################################
# Timeseries extension
###########################################################################################
function DynamicalSystems.interactive_trajectory_timeseries(
    ds::DynamicalSystem, fs::Vector, u0s = [current_state(ds)];
    linekwargs = isdiscretetime(ds)  ? (linewidth = 1,) : (linewidth = 3,),
    timeseries_names = [_timeseries_name(f) for f in fs],
    colors = [COLORS[i] for i in 1:length(u0s)],
    timeseries_ylims = [(0, 1) for f in fs],
    kwargs...)

    fig, dsobs = interactive_trajectory(ds, u0s; colors, figure = (resolution = (1600, 800),), kwargs...)
    timeserieslayout = fig[1,2] = GridLayout()
    _init_timeseries_plots!(
        timeserieslayout, dsobs, fs, colors, linekwargs, timeseries_names, timeseries_ylims,
    )

    return fig, dsobs
end

function _init_timeseries_plots!(
        layout, dsobs, fs, colors, linekwargs, tsnames, tslims,
    )

    # First, create axis
    axs = [Axis(layout[i, 1]; ylabel = tsnames[i]) for i in 1:length(fs)]
    for i in 1:length(fs); ylims!(axs[i], tslims[i]); end
    linkxaxes!(axs...)
    for i in 1:length(fs)-1; hidexdecorations!(axs[i]; grid = false); end
    axs[end].xlabel = "time"

    # Create and plot the observables of the timeseries
    T = length(dsobs.tail_observables[1][])
    for (j, f) in enumerate(fs)
        for (i, tail) in enumerate(dsobs.tail_observables)
            observed_data = map(tail, dsobs.current_step) do x, n
                [Point2f(max(0, dsobs.Δt*(n - T + k)), _obtain_data(x[k], f)) for k in 1:length(x)]
            end
            # plot them
            lk = linekwargs isa AbstractVector ? linekwargs[i] : linekwargs
            plotf! = isdiscretetime(dsobs.pds) ? scatterlines! : lines!
            plotf!(axs[j], observed_data; color = colors[i], lk...)
        end

        # Add a last observable trigger that changes the axis xspan
        on(dsobs.tail_observables[end]) do x
            n = dsobs.current_step[]
            xlims!(axs[end], max(0, dsobs.Δt*(n - T)), max(T*dsobs.Δt, dsobs.Δt*n))
        end
    end
    return
end

_obtain_data(x::AbstractVector, f::Int) = x[f]
_obtain_data(x::AbstractVector, f::Function) = f(x)
_timeseries_name(f::Int) = "x"*subscript(f)
_timeseries_name(f) = string(f)