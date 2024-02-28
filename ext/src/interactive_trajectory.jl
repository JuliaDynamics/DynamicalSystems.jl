using DynamicalSystems.DataStructures

function DynamicalSystems.interactive_trajectory(
        ds::DynamicalSystems.DynamicalSystem, u0s = [DynamicalSystems.current_state(ds)];
        # Selection of what to plot
        idxs = 1:min(length(u0s[1]), 3),
        # Time evolution
        tail = 1000,
        Δt = DynamicalSystems.isdiscretetime(ds) ? 1 : 0.01,
        pause = nothing,
        # Visualization
        colors = collect(cgrad(COLORSCHEME, length(u0s); categorical = true)),
        plotkwargs = NamedTuple(), markersize = 15,
        fade = true,
        # parameters
        parameter_sliders = nothing,
        # `pnames` is deprecated
        pnames = isnothing(parameter_sliders) ? nothing : Dict(keys(parameter_sliders) .=> parameter_name.(keys(parameter_sliders))),
        parameter_names = pnames,
        add_controls = true,
        # figure and axis
        figure = (size = (800, 800),),
        axis = NamedTuple(),
        lims = nothing,
        statespace_axis = true,
    )

    if length(idxs) > dimension(ds)
        throw(ArgumentError("More indices given than the system has dimension! Change `idxs`."))
    elseif length(idxs) > 3
        throw(ArgumentError("State space plot can be up to 3 dimensional! Change `idxs`."))
    end

    if ds isa CoupledODEs # force time evolution into non-adaptive
        newdiffeq = (ds.diffeq..., adaptive = false, dt = Δt)
        ds = CoupledODEs(ds, newdiffeq)
    end

    u00s = deepcopy(u0s)
    p0 = initial_parameters(ds)
    pds = DynamicalSystems.ParallelDynamicalSystem(ds, u0s)
    fig = Figure(; figure...)
    # Set up trajectrory plot
    statespacelayout = fig[1,1] = GridLayout()
    lims = isnothing(lims) ? _traj_lim_estimator(ds, u0s, idxs, [1], Δt)[1] : lims
    tailobs, finalpoints = _init_statespace_plot!(statespacelayout, ds, idxs,
        lims, pds, colors, plotkwargs, markersize, tail, axis, fade, statespace_axis,
    )
    # Set up layouting and add controls
    if add_controls # Notice that `run` and `step` are already observables
        reset, run, step, stepslider = _trajectory_plot_controls!(statespacelayout, statespace_axis)
    else
        # So that we can leave the interactive UI code as is
        reset = Observable(0); run = Observable(0); step = Observable(0); stepslider = Observable(1)
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
    # Resetting system to initial states
    on(reset) do clicks
        for j in eachindex(u00s)
            set_state!(dso, copy(u00s[j]), j)
        end
    end

    # Live parameter changing
    # note here `parameter_sliders` are parameters to have a slider; all parameters
    # can be changed after creation of `dso` via `set_parameter!`
    if !isnothing(parameter_sliders)
        paramlayout = fig[2, :] = GridLayout(tellheight = true, tellwidth = false)
        slidervals, sliders = _add_ds_param_controls!(
            ds, paramlayout, parameter_sliders, parameter_names, current_parameters(ds)
        )
        update = Button(fig, label = "update", tellwidth = false, tellheight = true)
        resetp = Button(fig, label = "reset p", tellwidth = false, tellheight = true)
        gl = paramlayout[2, :] = GridLayout()
        gl[1,1] = update
        gl[1,2] = resetp
        on(update.clicks) do clicks
            for l in keys(slidervals)
                v = slidervals[l][]
                set_parameter!(dso, l, v)
            end
        end
        on(resetp.clicks) do clicks
            set_parameters!(pds, p0)
            # Also **visually** reset sliders to initial parameters
            for k in keys(p0)
                haskey(sliders, k) || continue
                set_close_to!(sliders[k], p0[k])
            end
        end
    end

    return fig, dso
end

vector_idx_observe(ds, u, idxs::AbstractVector{<:Int}) = u[idxs]
vector_idx_observe(ds, u, idxs) = [observe_state(ds, i, u) for i in idxs]

# Main panels of animation
"Create the state space axis and evolution controls. Return the axis."
function _init_statespace_plot!(
        layout, ds, idxs, lims, pds, colors, plotkwargs, markersize, tail, axis, fade,
        statespace_axis # whether to show the statespace axis
    )
    tailobs, finalpoints = _init_trajectory_observables(pds, tail)
    is3D = length(idxs) == 3
    axisposition = statespace_axis ? layout[1,1] : Figure()[1,1]
    xlabel = state_name(idxs[1])
    ylabel = state_name(idxs[2])
    statespaceax = if is3D
        zlabel = state_name(idxs[3])
        Axis3(axisposition; xlabel, ylabel, zlabel, axis...)
    else
        Axis(axisposition; xlabel, ylabel, axis...)
    end
    # Here we make two more observables for the plotted tails and plotted final
    # states, so that the stored observables in `dsobs` are the full system state;
    # This simplifies drastically making custom animations
    T = length(idxs) == 2 ? Point2f : Point3f
    plotted_tailobs = [
        # here x is a tail of the full dynamical system state (an observable!)
        # while y is a point on the tail, i.e., a state space point
        map(x -> T[vector_idx_observe(ds, y, idxs) for y in x], ob) for ob in tailobs
    ]
    plotted_finalpoints = map(x -> T[vector_idx_observe(ds, y, idxs) for y in x], finalpoints)

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
        color = colors, markersize, finalargs...
    )
    !isnothing(lims) && (statespaceax.limits = lims)
    is3D && (statespaceax.protrusions = 50) # removes overlap of labels
    return tailobs, finalpoints
end

function _init_trajectory_observables(pds, tail)
    N = length(current_states(pds))
    tailobs = Observable[]
    # Copy ensures that each state in the circular buffer is independent
    # of the state of the `pds`, if it is inplace
    # (copy should be free for out of place anyways)
    T = typeof(copy(current_state(pds)))
    for i in 1:N
        cb = CircularBuffer{T}(tail)
        fill!(cb, current_state(pds, i))
        push!(tailobs, Observable(cb))
    end
    finalpoints = Observable([x[][end] for x in tailobs])
    return tailobs, finalpoints
end
function _trajectory_plot_controls!(gl, statespace_axis)
    position = statespace_axis ? [2,1] : [1,1]
    controllayout = setindex!(gl, GridLayout(tellwidth = false, tellheight = statespace_axis), position...)
    reset = Button(controllayout[1, 0]; label = "reset")
    run = Button(controllayout[1, 1]; label = "run")
    step = Button(controllayout[1, 2]; label = "step")
    slider_vals = vcat(1:10, 100:100:1000)
    sg = SliderGrid(controllayout[1,3],
        (label = "steps =", range = slider_vals, startvalue = 1),
    )
    return reset.clicks, run.clicks, step.clicks, sg.sliders[1].value
end
function _traj_lim_estimator(ds, u0s, idxs, observables, dt)
    ds = deepcopy(ds)
    mi = fill(Inf, length(idxs))
    ma = fill(-Inf, length(idxs))
    omi = fill(Inf, length(observables))
    oma = fill(-Inf, length(observables))
    for i in 1:length(u0s)
        tr, = DynamicalSystems.trajectory(ds, 1000dt, u0s[i]; Δt = dt)
        _tr = StateSpaceSet(map(u -> [observe_state(ds, i, u) for i in idxs], tr))
        mii, maa = DynamicalSystems.minmaxima(_tr)
        mi = min.(mii, mi)
        ma = max.(maa, ma)
        # do same but now by transforming the `tr` set into the observed set
        otr = StateSpaceSet(map(u -> [observe_state(ds, f, u) for f in observables], tr))
        omii, omaa = DynamicalSystems.minmaxima(otr)
        omi = min.(omii, omi)
        oma = max.(omaa, oma)
    end
    # Alright, now we just have to put them into limits and increase a bit
    lims = [(mi[i]-0.1mi[i], ma[i]+0.1ma[i]) for i in 1:length(idxs)]
    lims = (lims...,)
    observable_lims = [(omi[i]-0.1omi[i], oma[i]+0.1oma[i]) for i in 1:length(observables)]
    return lims, observable_lims
end

# Parameter handling
function _add_ds_param_controls!(ds, paramlayout, parameter_sliders, pnames, p0)
    slidervals = Dict{keytype(parameter_sliders), Observable}() # directly has the slider observables
    sliders = Dict{keytype(parameter_sliders), Any}() # for updating via reset parameters
    tuples_for_slidergrid = []
    for (l, vals) in parameter_sliders # dictionary
        startvalue = current_parameter(ds, l, p0)
        label = string(pnames[l])
        push!(tuples_for_slidergrid, (;label, range = vals, startvalue))
    end
    sg = SliderGrid(paramlayout[1,1], tuples_for_slidergrid...; tellheight = true)
    for (i, (l, vals)) in enumerate(parameter_sliders)
        slidervals[l] = sg.sliders[i].value
        sliders[l] = sg.sliders[i]
    end
    return slidervals, sliders
end

###########################################################################################
# Timeseries extension
###########################################################################################
function DynamicalSystems.interactive_trajectory_timeseries(
    ds::DynamicalSystem, fs::Vector, u0s = [current_state(ds)];
    linekwargs = isdiscretetime(ds)  ? (linewidth = 1,) : (linewidth = 3,),
    timeseries_names = [state_name(f) for f in fs],
    colors = collect(cgrad(COLORSCHEME, length(u0s); categorical = true)),
    timeseries_ylims = nothing,
    timelabel = "time", timeunit = 1,
    Δt = DynamicalSystems.isdiscretetime(ds) ? 1 : 0.01,
    idxs = 1:min(length(u0s[1]), 3),
    lims = nothing,
    kwargs...)

    # automatic limits
    if isnothing(timeseries_ylims) || isnothing(lims)
        _lims, _timeseries_ylims = _traj_lim_estimator(ds, u0s, idxs, fs, Δt)
    end
    isnothing(lims) && (lims = _lims)
    isnothing(timeseries_ylims) && (timeseries_ylims = _timeseries_ylims)

    fig, dsobs = interactive_trajectory(ds, u0s; Δt, idxs, lims, colors, figure = (size = (1600, 800),), kwargs...)

    timeserieslayout = fig[:, 2] = GridLayout()
    _init_timeseries_plots!(
        timeserieslayout, dsobs, fs, colors, linekwargs, timeseries_names,
        timeseries_ylims, timelabel, timeunit
    )

    return fig, dsobs
end

function _init_timeseries_plots!(
        layout, dsobs, fs, colors, linekwargs, tsnames, tslims, timelabel, timeunit
    )

    # First, create axis
    axs = [Axis(layout[i, 1]; ylabel = tsnames[i]) for i in 1:length(fs)]
    for i in 1:length(fs); ylims!(axs[i], tslims[i]); end
    linkxaxes!(axs...)
    for i in 1:length(fs)-1; hidexdecorations!(axs[i]; grid = false); end
    axs[end].xlabel = timelabel

    # Create and plot the observables of the timeseries
    T = length(dsobs.tail_observables[1][])
    for (j, f) in enumerate(fs)
        for (i, tail) in enumerate(dsobs.tail_observables)
            observed_data = map(tail, dsobs.current_step) do x, n
                [Point2f(max(0, dsobs.Δt*(n - T + k)/timeunit), observe_state(dsobs.pds, f, x[k])) for k in 1:length(x)]
            end
            # plot them
            lk = linekwargs isa AbstractVector ? linekwargs[i] : linekwargs
            plotf! = isdiscretetime(dsobs.pds) ? scatterlines! : lines!
            plotf!(axs[j], observed_data; color = colors[i], lk...)
        end

        # Add a last observable trigger that changes the axis xspan
        on(dsobs.tail_observables[end]) do x
            n = dsobs.current_step[]
            xlims!(axs[end],
                max(0, dsobs.Δt*(n - T)/timeunit),
                max(T*dsobs.Δt/timeunit, dsobs.Δt*n/timeunit)
            )
        end
    end
    return
end

_obtain_data(x::AbstractVector, f::Int) = x[f]
_obtain_data(x::AbstractVector, f::Function) = f(x)
