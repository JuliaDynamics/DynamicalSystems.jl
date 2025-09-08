function DynamicalSystems.interactive_2d_clicker(ds;
        # DynamicalSystems kwargs:
        times = 100:10_000,
        Δt = 1,
        # Makie kwargs:
        color = randomcolor,
        plotkwargs = ()
    )

    figure = Figure(size = (1000, 800))

    T_slider, m_slider = _add_clicker_controls!(figure, times)
    ax = figure[1, :] = Axis(figure; tellheight = true)

    # Compute the initial plot
    u0 = DynamicalSystems.current_state(ds)
    data, = trajectory(ds, T_slider[]; Δt)
    positions_node = Observable(data)
    colors = (c = color(u0); [c for _ in 1:length(data)])
    colors_node = Observable(colors)

    if isdiscretetime(ds)
        scatter!(
            ax, positions_node, color = colors_node,
            markersize = lift(o -> o*px, m_slider), marker = :circle, plotkwargs...
        )
    else
        scatterlines!(
            ax, positions_node, color = colors_node,
            markersize = lift(o -> o*px, m_slider), marker = :circle, plotkwargs...
        )
    end

    # Interactive clicking on the phase space:
    laststate = Observable(u0)
    Makie.deactivate_interaction!(ax, :rectanglezoom)
    spoint = select_point(ax.scene)
    on(spoint) do newstate
        data, = trajectory(ds, T_slider[], newstate; Δt)
        pushfirst!(vec(data), fill(NaN, dimension(data))) # ensures break for scatterlines
        positions = positions_node[]; colors = colors_node[]
        append!(positions, data)
        c = color(newstate)
        append!(colors, fill(c, length(data)))
        # Update all the observables with Array as value:
        positions_node[], colors_node[], laststate[] = positions, colors, newstate
    end

    display(figure)
    return figure, laststate
end

function _add_clicker_controls!(figure, times)
    sg1 = SliderGrid(figure[2, :][1, 1],
        (label = "T", range = times,
        format = x -> string(round(x)),
        startvalue = times[1])
    )
    sg2 = SliderGrid(figure[2, :][1, 2],
        (label = "ms", range = 10.0 .^ range(0, 2, length = 100),
        format = x -> string(round(x)), startvalue = 10)
    )
    return sg1.sliders[1].value, sg2.sliders[1].value
end

# interactive psos is based in the 2D clicker
function DynamicalSystems.interactive_poincaresos(ds, plane, idxs, complete;
        # PSOS kwargs:
        direction = -1,
        rootkw = (xrtol = 1e-6, atol = 1e-6),
        Tmax = 1e3,
        kw...
    )

    # Basic sanity checks on the method arguments
    @assert typeof(plane) <: Tuple
    @assert length(idxs) == 2
    @assert eltype(idxs) == Int
    @assert plane[1] ∉ idxs

    i = DynamicalSystems.SVector{2, Int}(idxs)

    # Construct a new `PoincareMap` structure with the given parameters
    pmap = DynamicalSystems.DynamicalSystemsBase.PoincareMap(ds, plane;
        direction, rootkw, Tmax)

    # construct a 2d projected system compatible with the clicker
    z = plane[2] # third variable comes from plane
    complete_state = u -> complete(u..., z)
    project = i
    newds = ProjectedDynamicalSystem(pmap, project, complete_state)
    return interactive_2d_clicker(newds; kw...)
end
