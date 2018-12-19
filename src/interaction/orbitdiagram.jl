using DynamicalSystems, Makie
using DynamicalSystems, BenchmarkTools
using DynamicalSystemsBase: DDS, MDI

# put this to makie
onlyleftclick(scplot) = (ispressed(scplot, Mouse.left) && !ispressed(scplot, Keyboard.space) &&
                   AbstractPlotting.is_mouseinside(scplot))

function interactive_orbitdiagram(ds::DDS, i::Int, p_index, p_min, p_max;
    density = 1000, u0 = get_state(ds), Ttr = 100, n = 100)

    # Initialization
    positions = Vector{Point2f0}(undef, n * density)
    integ = integrator(ds)

    # positions_node = Node(positions)

    pmin, pmax = p_min, p_max
    populate_orbitdiagram!(positions, integ, i, p_index, pmin, pmax, density, n, Ttr)

    scplot = scatter(positions, markersize = 0.01)

    on(scplot.events.mousebuttons) do buttons
        if onlyleftclick(scplot)
            pmin, xmin = mouseposition(scplot)
            # Get xmax, pmax with un-click
            pmax = 1.5pmin; xmax = 1.5xmin
            populate_orbitdiagram!(positions, integ, i,
                                   p_index, pmin, pmax, density, n, Ttr)
            # positions_node[] = positions
            limits = FRect(pmin,xmin,pmax-pmin,xmax-xmin)
            scplot = scatter(positions, markersize = 0.01, limits = limits)
        end
        display(scplot)
    end
    display(scplot)
    return scplot
end

function populate_orbitdiagram!(positions, integ, i, p_index, pmin, pmax, density, n, Ttr)
    pvalues = range(pmin, stop = pmax, length = density)

    k = 1
    for p in pvalues
        DynamicalSystemsBase.reinit!(integ)
        integ.p[p_index] = p
        DynamicalSystemsBase.step!(integ, Ttr)
        for z in 1:n
            DynamicalSystemsBase.step!(integ)
            @inbounds positions[k] = Point2f0(p, integ.u[i])
            k += 1
        end
    end
end



ds = Systems.henon()
i = 1
p_index = 1
p_min = 0.8
p_max = 1.4

interactive_orbitdiagram(ds, i, p_index, p_min, p_max)
