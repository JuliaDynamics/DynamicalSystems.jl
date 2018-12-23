using DynamicalSystemsBase, Makie

# add @inbounds when this is done.

tslider = AbstractPlotting.textslider
onlyleftclick(scplot) = (ispressed(scplot, Mouse.left) && !ispressed(scplot, Keyboard.space) &&
                   AbstractPlotting.is_mouseinside(scplot))

function interactive_orbitdiagram(ds::DiscreteDynamicalSystem,
    i::Int, p_index, p_min, p_max;
    density = 1000, u0 = get_state(ds), Ttr = 100, n = 100
    )

    # Initialization
    integ = integrator(ds)
    scene = Scene(resolution = (1000, 600))
    pmin, pmax = p_min, p_max
    od = minimal_od(integ, i, p_index, pmin, pmax, density, n, Ttr, u0)
    od_node = Node(od)

    ui_n, n = tslider(
        round.(Int, range(10, stop=10000, length=1000)),"n", start=n)
    ui_T, Ttr = tslider(
        round.(Int, range(10, stop=10000, length=1000)),"Ttr", start=Ttr)
    ui_d, density = tslider(
        round.(Int, range(10, stop=10000, length=1000)),"density", start=density)

    scplot = scatter(od_node, markersize = 0.01)

    on(scplot.events.mousebuttons) do buttons
        if onlyleftclick(scplot)
            pmin, xmin = mouseposition(scplot)
            # Get xmax, pmax with un-click
            pmax = pmin + 0.2; xmax = xmin + 0.4
            od_node[] = minimal_od(
                integ, i,  p_index, pmin, pmax,
                density[], n[], Ttr[], u0, xmin, xmax
            )

            limits = FRect(pmin,xmin,pmax-pmin,xmax-xmin)
            # AbstractPlotting.update_limits!(scplot, limits)
        end
    end
    hbox(vbox(ui_n, ui_T, ui_d), scplot, parent=scene)
    display(scene)
    return od_node
end

function minimal_od(integ, i, p_index, pmin, pmax,
                    density::Int, n::Int, Ttr::Int, u0, xmin = -Inf, xmax = Inf)

    pvalues = range(pmin, stop = pmax, length = density)
    od = Vector{Point2f0}()
    for p in pvalues
        DynamicalSystemsBase.reinit!(integ, u0)
        integ.p[p_index] = p
        DynamicalSystemsBase.step!(integ, Ttr)
        for z in 1:n
            DynamicalSystemsBase.step!(integ)
            x = integ.u[i]
            if xmin ≤ x ≤ xmax
                push!(od, Point2f0(p, integ.u[i]))
            end
        end
    end
    return od
end

ds = Systems.henon()
i = 1
p_index = 1
p_min = 0.8
p_max = 1.4

od_node = interactive_orbitdiagram(ds, i, p_index, p_min, p_max)
