ChaosTools = DynamicalSystems.ChaosTools

function DynamicalSystems.interactive_poincaresos(ds, plane, idxs, complete;
        # PSOS kwargs:
        direction = -1,
        tfinal = (1000.0, 10.0^4),
        rootkw = (xrtol = 1e-6, atol = 1e-6),
        # Makie kwargs:
        color = randomcolor,
        scatterkwargs = (),
        labels = ("u₁", "u₂")
    )

    # Basic sanity checks on the method arguments
    @assert typeof(plane) <: Tuple
    @assert length(idxs) == 2
    @assert eltype(idxs) == Int
    @assert plane[1] ∉ idxs
    u0 = DynamicalSystems.get_state(ds)

    i = DynamicalSystems.SVector{2, Int}(idxs)

    figure = Figure(size = (1000, 800), backgroundcolor = :white)

    T_slider, m_slider = _add_psos_controls!(figure, tfinal)
    ax = figure[0, :] = Axis(figure)

    # Construct a new `PoincareMap` structure with the given parameters
    pmap = DynamicalSystems.DynamicalSystemsBase.PoincareMap(ds, plane;
        direction, u0, rootkw, Tmax = tfinal[2])

    # Compute the initial section
    psos, = trajectory(pmap, T_slider[]; t0 = 0)
    data = psos[:, i]
    length(data) == 0 && error(ChaosTools.PSOS_ERROR)

    positions_node = Observable(data)
    colors = (c = color(u0); [c for i in 1:length(data)])
    colors_node = Observable(colors)
    scatter!(
        ax, positions_node, color = colors_node,
        markersize = lift(o -> o*px, m_slider), marker = :circle, scatterkwargs...
    )

    ax.xlabel, ax.ylabel = labels
    laststate = Observable(u0)

    # Interactive clicking on the psos:
    Makie.deactivate_interaction!(ax, :rectanglezoom)
    spoint = select_point(ax.scene)
    on(spoint) do pos
        x, y = pos; z = plane[2] # third variable comes from plane
        newstate = try
           complete(x, y, z)
        catch err
           @error "Could not get state, got error: " exception=err
           return
        end

        psos, = trajectory(pmap, T_slider[], newstate; t0 = 0)
        data = psos[:, i]
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

function _add_psos_controls!(figure, tfinal)
    sg1 = SliderGrid(figure[1, :][1,1],
        (label = "T", range = range(tfinal[1], tfinal[2], length = 1000),
        format = x -> string(round(x)), )
    )
    sg2 = SliderGrid(figure[1, :][1,2],
        (label = "ms", range = 10.0 .^ range(0, 2, length = 100),
        format = x -> string(round(x)), startvalue = 10)
    )
    return sg1.sliders[1].value, sg2.sliders[1].value
end
