function DynamicalSystems.interactive_clicker(dds;
        # DynamicalSystems kwargs:
        tfinal = (1000.0, 10.0^4),
        # TODO: project and complete, to be reusable from psos
        # Makie kwargs:
        color = randomcolor,
        scatterkwargs = (),
        labels = ("x", "y")
    )

    u0 = DynamicalSystems.get_state(dds)

    figure = Figure(size = (1000, 800), backgroundcolor = :white)

    T_slider, m_slider = _add_clicker_controls!(figure, tfinal)
    ax = figure[0, :] = Axis(figure)

    # Compute the initial section
    tr, = trajectory(dds, T_slider[]; t0 = 0)
    length(tr) == 0 && error("Initial computed trajectory is empty!")
    length(tr[1]) != 2 && error("Trajectory is not 2D")

    positions_node = Observable(tr)
    colors = (c = color(u0); [c for i in 1:length(tr)])
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
        newstate = [x, y]

        tr, = trajectory(dds, T_slider[], newstate; t0 = 0)
        positions = positions_node[]; colors = colors_node[]
        append!(positions, tr)
        c = color(newstate)
        append!(colors, fill(c, length(tr)))
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
