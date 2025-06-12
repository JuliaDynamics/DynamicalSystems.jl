function DynamicalSystems.interactive_clicker(ds;
        # DynamicalSystems kwargs:
        tfinal = (1000.0, 10.0^4),
        complete = (x, y) -> [x, y],
        project = identity,
        # Makie kwargs:
        color = randomcolor,
        labels = ("x", "y"),
        scatterkwargs = ()
    )

    u0 = DynamicalSystems.get_state(ds)

    figure = Figure(size = (1000, 800), backgroundcolor = :white)

    T_slider, m_slider = _add_clicker_controls!(figure, tfinal)
    ax = figure[0, :] = Axis(figure)

    # Compute the initial section
    tr, = trajectory(ds, T_slider[]; t0 = 0)
    length(tr) == 0 && error("Initial computed trajectory is empty")

    data = project(tr)
    length(data[1]) != 2 && error("(Projected) trajectory is not 2D")

    positions_node = Observable(data)
    colors = (c = color(u0); [c for _ in 1:length(data)])
    colors_node = Observable(colors)
    scatter!(
        ax, positions_node, color = colors_node,
        markersize = lift(o -> o*px, m_slider), marker = :circle, scatterkwargs...
    )

    ax.xlabel, ax.ylabel = labels
    laststate = Observable(u0)

    # Interactive clicking on the phase space:
    Makie.deactivate_interaction!(ax, :rectanglezoom)
    spoint = select_point(ax.scene)
    on(spoint) do pos
        x, y = pos;
        newstate = try
           complete(x, y)
        catch err
           @error "Could not complete state, got error: " exception=err
           return
        end

        tr, = trajectory(ds, T_slider[], newstate; t0 = 0)
        data = project(tr)

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

function _add_clicker_controls!(figure, tfinal)
    sg1 = SliderGrid(figure[1, :][1, 1],
        (label = "T", range = range(tfinal[1], tfinal[2], length = 1000),
        format = x -> string(round(x)), )
    )
    sg2 = SliderGrid(figure[1, :][1, 2],
        (label = "ms", range = 10.0 .^ range(0, 2, length = 100),
        format = x -> string(round(x)), startvalue = 10)
    )
    return sg1.sliders[1].value, sg2.sliders[1].value
end
