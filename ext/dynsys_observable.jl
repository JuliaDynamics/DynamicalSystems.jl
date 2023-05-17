using DynamicalSystems.DataStructures
# This is a very light struct that contains the trajectory end state
# and the trajectory tail. Note that for convenience it always contains the state
# as a vector of observables, each observable containing each of the
# parallel states of the dynamical system
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

function DynamicalSystems.interactive_trajectory_panel(
        ds::DynamicalSystems.DynamicalSystem, u0s = [DynamicalSystems.current_state(ds)];
        # Selection of what to plot
        idxs = 1:min(length(u0s[1]), 3),
        # Time evolution
        tail = 1000,
        total_span = DynamicalSystems.isdiscretetime(ds) ? 50 : 10,
        Δt = DynamicalSystems.isdiscretetime(ds) ? 1 : 0.01,
        steps_per_update = 1,
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
        run = Observable(0); step = Observable(0); stepslider = Observable(steps_per_update)
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
            yield()
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
        slidervals = _add_ds_param_controls!(paramlayout, ps, po, pnames)
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
        if !DynamicalSystems.isdiscretetime(ds)
            x = to_color(colors[i])
            if fade
                x = [RGBAf(x.r, x.g, x.b, i/tail) for i in 1:tail]
            end
            Makie.lines!(statespaceax, ob;
                color = x, linewidth = 3.0, transparency = true, pk...
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
    tr = DynamicalSystems.trajectory(ds, 1000, u0s[1]; Δt)[1]
    _mi, _ma = DynamicalSystems.minmaxima(tr)
    mi, ma = _mi[idxs], _ma[idxs]
    for i in 2:length(u0s)
        tr = DynamicalSystems.trajectory(ds, 1000, u0s[i]; Δt)[1]
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
function _add_ds_param_controls!(paramlayout, ps, pobs, pnames)
    slidervals = Dict{keytype(ps), Observable}() # directly has the slider observables
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
    end
    return slidervals
end


function DynamicalSystems.set_parameter!(dso::DynamicalSystemObservable, index, value)
    dso.param_observable[][l] = v
    set_parameter!(dso.pds, l, v)
    notify(dso.param_observable)
    return
end
