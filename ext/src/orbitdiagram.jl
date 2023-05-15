export interactive_orbitdiagram, scaleod
# TODO: Allow initial state to be a function of parameter (define function `get_u(f, p)`)

"""
    interactive_orbitdiagram(
        ds::DynamicalSystem, p_index, pmin, pmax, i::Int = 1;
        u0 = nothing, parname = "p", title = ""
    )

Open an interactive application for exploring orbit diagrams (ODs) of discrete
dynamical systems. Requires `DynamicalSystems`.

In essense, the function presents the output of `DynamicalSystems.orbitdiagram`
of the `i`th variable of the `ds`, and allows interactively zooming into it.

Keywords control the name of the parameter, the initial state (used for _any_ parameter)
or whether to add a title above the orbit diagram.

## Interaction

The application is separated in the "OD plot" (left) and the "control panel" (right).
On the OD plot you can interactively click
and drag with the left mouse button to select a region in the OD. This region is then
**re-computed** at a higher resolution.

The options at the control panel are straight-forward, with
* `n` amount of steps recorded for the orbit diagram (not all are in the zoomed region!)
* `t` transient steps before starting to record steps
* `d` density of x-axis (the parameter axis)
* `α` alpha value for the plotted points.

Notice that at each update `n*t*d` steps are taken.
You have to press `update` after changing these parameters.
Press `reset` to bring the OD in the original
state (and variable). Pressing `back` will go back through the history of your exploration
History is stored when the "update" button is pressed or a region is zoomed in.

You can even decide which variable to get the OD for
by choosing one of the variables from the wheel!
Because the y-axis limits can't be known when changing variable, they reset to the size
of the selected variable.

## Accessing the data

What is plotted on the application window is a _true_ orbit diagram, not a plotting
shorthand. This means that all data are obtainable and usable directly.
Internally we always scale the orbit diagram to [0,1]² (to allow `Float64` precision
even though plotting is `Float32`-based). This however means that it is
necessary to transform the data in real scale. This is done through the function
[`scaleod`](@ref) which accepts the 5 arguments returned from the current function:
```julia
figure, oddata = interactive_orbitdiagram(...)
ps, us = scaleod(oddata)
```
"""
function interactive_orbitdiagram(ds, p_index, p_min, p_max, i0::Int = 1;
        u0 = nothing, parname = "p", title = ""
    )

    figure = Figure(resolution = (1200, 600), backgroundcolor = DEFAULT_BG)
    display(figure)
    odax = figure[1,1] = Axis(figure; title)
    for z in (:xpanlock, :ypanlock, :xzoomlock, :yzoomlock)
        setproperty!(odax, z, true)
    end
    controllayout = figure[1, 2] = GridLayout()
    display(figure)

    controllayout.tellheight[] = false
    odax.tellheight = true

    nslider, Tslider, dslider, n, Ttr, d, α, i,
    ▢update, ▢back, ▢reset, ⬜p₋, ⬜p₊, ⬜u₋, ⬜u₊ =
    add_controls!(controllayout, figure, DynamicalSystems.dimension(ds), parname, i0)

    # Initial Orbit diagram data
    DynamicalSystems.reinit!(ds, u0)
    p₋, p₊ = p_min, p_max
    odinit, xmin, xmax = minimal_normalized_od(ds, i[], p_index, p₋, p₊, d[], n[], Ttr[], u0)
    od_obs = Observable(odinit)

    # History stores the variable index and true diagram limits
    history = [(i[], p₋, p₊, xmin, xmax, n[], Ttr[], d[])]
    update_controls!(history[end], i, n, Ttr, d,  ⬜p₋, ⬜p₊, ⬜u₋, ⬜u₊)

    color = lift(a -> RGBAf(0,0,0,a), α)
    scatter!(odax, od_obs; marker = MARKER, markersize = 1px, color = color, strokewidth = 0.0)

    xlims!(odax, 0, 1)
    ylims!(odax, 0, 1)
    odax.xticks = (0:0.25:1, ["$(parname)₋", " ", " ", " ", "$(parname)₊"])
    odax.yticks = (0:0.25:1, ["u₋", " ", " ", " ", "u₊"])
    odax.ylabel = "u"*subscript(i[])
    on(i) do o; odax.ylabel = "u"*subscript(o); end
    odax.xlabel = parname
    MakieLayout.deactivate_interaction!(odax, :rectanglezoom)
    rect = select_rectangle(odax.scene)

    # # Uppon interactively selecting a rectangle, with value `r` (in [0,1]²)
    on(rect) do r
        sp₋, sxmin = r.origin
        sp₊, sxmax = r.origin + r.widths
        # Convert p,x to true values
        j, pp₋, pp₊, pxmin, pxmax = history[end]
        pdif = pp₊ - pp₋; xdif = pxmax - pxmin
        p₋ = sp₋*pdif + pp₋
        p₊ = sp₊*pdif + pp₋
        xmin = sxmin*xdif + pxmin
        xmax = sxmax*xdif + pxmin

        od_obs[] = minimal_normalized_od(
            ds, j,  p_index, p₋, p₊,
            d[], n[], Ttr[], u0, xmin, xmax
        )
        # update history and controls
        push!(history, (j, p₋, p₊, xmin, xmax, n[], Ttr[], d[]))
        update_controls!(history[end], i, n, Ttr, d,  ⬜p₋, ⬜p₊, ⬜u₋, ⬜u₊)
    end

    # Upon hitting the update button
    on(▢update) do clicks
        j, p₋, p₊, xmin, xmax, m, T, dens = history[end]
        # Check if there was any change:
        if !(⬜p₋[] == p₋ && ⬜p₊[] == p₊ && ⬜u₋[] == xmin && ⬜u₊[] == xmax &&
            j == i[] && m == n[] && T == Ttr[] && dens == d[])

            p₋, p₊, xmin, xmax = ⬜p₋[], ⬜p₊[], ⬜u₋[], ⬜u₊[]
            newj, m, T, dens = i[], n[], Ttr[], d[]

            if newj ≠ j # a new variable is selected, so x limits must be recomputed
                od_obs[], xmin, xmax = minimal_normalized_od(
                    ds, newj, p_index, p₋, p₊, d[], n[], Ttr[], u0
                )
            else
                od_obs[] = minimal_normalized_od(
                    ds, newj, p_index, p₋, p₊, d[], n[], Ttr[], u0, xmin, xmax
                )
            end
            # Update history and controls
            push!(history, (newj, p₋, p₊, xmin, xmax, m, T, dens))
            update_controls!(history[end], i, n, Ttr, d,  ⬜p₋, ⬜p₊, ⬜u₋, ⬜u₊)
        end
    end

    # Upon hitting the "reset" button
    on(▢reset) do clicks
        if length(history) > 1
            deleteat!(history, 2:length(history))
            od_obs[] = odinit
            update_controls!(history[end], i, n, Ttr, d,  ⬜p₋, ⬜p₊, ⬜u₋, ⬜u₊)
        end
    end

    # Upon hitting the "back" button
    on(▢back) do clicks
        if length(history) > 1
            pop!(history)
            j, p₋, p₊, xmin, xmax, m, T, dens = history[end]
            od_obs[] = minimal_normalized_od(
                ds, j, p_index, p₋, p₊, dens, m, T, u0, xmin, xmax
            )
            update_controls!(history[end], i, n, Ttr, d,  ⬜p₋, ⬜p₊, ⬜u₋, ⬜u₊)
        end
    end

    # for the following two buttons the slider values must be updated
    onany(▢reset, ▢back) do val1, val2
        j, p₋, p₊, xmin, xmax, m, T, dens = history[end]
        set_close_to!(nslider, m)
        set_close_to!(Tslider, T)
        set_close_to!(dslider, dens)
    end

    return figure, (od_obs, ⬜p₋, ⬜p₊, ⬜u₋, ⬜u₊)
end


"""
    minimal_normalized_od(ds, i, p_index, p₋, p₊, d, n, Ttr, u0) → od, xmin, xmax
    minimal_normalized_od(ds, i, p_index, p₋, p₊, d, n, Ttr, u0, xmin, xmax) → od

Compute and return a minimal and normalized orbit diagram (OD).

All points are stored in a single vector of `Point2f` to ensure fastest possible
plotting. In addition all numbers are scaled to [0, 1]. This allows us to have
64-bit precision while display is only 32-bit!

The version with `xmin, xmax` only keeps points with limits between the
real `xmin, xmax` (in the normal units of the dynamical system).
"""
function minimal_normalized_od(ds, i, p_index, p₋, p₊,
                              d::Int, n::Int, Ttr::Int, u0)

    pvalues = range(p₋, stop = p₊, length = d)
    pdif = p₊ - p₋
    od = Vector{Point2f}() # make this pre-allocated
    xmin = eltype(DynamicalSystems.current_state(ds))(Inf)
    xmax = eltype(DynamicalSystems.current_state(ds))(-Inf)
    @inbounds for p in pvalues
        pp = (p - p₋)/pdif # p to plot, in [0, 1]
        DynamicalSystems.set_parameter!(ds, p_index, p)
        DynamicalSystems.reinit!(ds, u0) # by default `u0` is nothing, so this does nothing
        Ttr > 0 && DynamicalSystems.step!(ds, Ttr)
        for _ in 1:n
            DynamicalSystems.step!(ds)
            x = DynamicalSystems.current_state(ds)[i]
            push!(od, Point2f(pp, x))
            # update limits
            if x < xmin
                xmin = x
            elseif x > xmax
                xmax = x
            end
        end
    end
    # normalize x values to [0, 1]
    xdif = xmax - xmin
    @inbounds for j in eachindex(od)
        x = od[j][2]; p = od[j][1]
        od[j] = Point2f(p, (x - xmin)/xdif)
    end
    return od, xmin, xmax
end

function minimal_normalized_od(ds, i, p_index, p₋, p₊,
                              d::Int, n::Int, Ttr::Int, u0, xmin, xmax)

    pvalues = range(p₋, stop = p₊, length = d)
    pdif = p₊ - p₋; xdif = xmax - xmin
    od = Vector{Point2f}()
    @inbounds for p in pvalues
        pp = (p - p₋)/pdif # p to plot, in [0, 1]
        DynamicalSystems.set_parameter!(ds, p_index, p)
        DynamicalSystems.reinit!(ds, u0)
        Ttr > 0 && DynamicalSystems.step!(ds, Ttr)
        for _ in 1:n
            DynamicalSystems.step!(ds)
            x = DynamicalSystems.current_state(ds)[i]
            if xmin ≤ x ≤ xmax
                push!(od, Point2f(pp, (DynamicalSystems.current_state(ds)[i] - xmin)/xdif))
            end
        end
    end
    return od
end

function  update_controls!(h, i, n, Ttr, d, ⬜p₋, ⬜p₊, ⬜u₋, ⬜u₊)
    j, p₋, p₊, xmin, xmax, m, T, dens = h
    i[] = j; n[] = m; Ttr[] = T; d[] = dens
    ⬜p₋[] = p₋; ⬜p₊[] = p₊
    ⬜u₋[] = xmin; ⬜u₊[] = xmax
    return
end

"""
    scaleod(oddata) -> ps, us
Given the return values of [`interactive_orbitdiagram`](@ref), produce
orbit diagram data scaled correctly in data units. Return the data as a vector of
parameter values and a vector of corresponding variable values.
"""
scaleod(r::Tuple) = scaleod(r...)
function scaleod(od, p₋, p₊, u₋, u₊)
    oddata = od[]; L = length(oddata);
    T = promote_type(typeof(u₋[]), Float32)
    ps = zeros(T, L); us = copy(ps)
    udif = u₊[] - u₋[]; um = u₋[]
    pdif = p₊[] - p₋[]; pm = p₋[]
    @inbounds for i ∈ 1:length(oddata)
        p, u = oddata[i]
        ps[i] = pm + pdif*p; us[i] = um + udif*u
    end
    return ps, us
end

function add_controls!(controllayout, figure, D, parname, i0)
    # Sliders
    sg = SliderGrid(controllayout[1,1],
        (label = "n", range = 1000:100:100000, startvalue = 1000),
        (label = "t", range = 1000:1000:10000, startvalue = 1000),
        (label = "d", range = 100:100:10000, startvalue = 1000),
        (label = "α", range = 0.001:0.001:1, startvalue = 0.1),
    )
    nslider, Tslider, dslider, αslider =  [s for s in sg.sliders]
    # Buttons (incl. variable chooser)
    ▢update = Button(figure, label = "update")
    ▢back = Button(figure, label = "← back")
    ▢reset = Button(figure, label = "reset")
    imenu = Menu(figure, options = [string(j) for j in 1:D], width = 60)
    imenu.i_selected = i0
    buttonslayout = controllayout[2, 1] = GridLayout()
    buttonslayout[1, 1:5] = [▢update, ▢back, ▢reset, Label(figure, "variable:"), imenu]
    # Limit boxes
    # TODO: Make these text input boxes
    ⬜p₋, ⬜p₊, ⬜u₋, ⬜u₊ = Observable.((0.0, 1.0, 0.0, 1.0))
    tsize = 16
    text_p₋ = Label(figure, lift(o -> "$(parname)₋ = $(o)", ⬜p₋),
        halign = :left, width = Auto(false), fontsize = tsize)
    text_p₊ = Label(figure, lift(o -> "$(parname)₊ = $(o)", ⬜p₊),
        halign = :left, width = Auto(false), fontsize = tsize)
    text_u₋ = Label(figure, lift(o -> "u₋ = $(o)", ⬜u₋),
        halign = :left, width = Auto(false), fontsize = tsize)
    text_u₊ = Label(figure, lift(o -> "u₊ = $(o)", ⬜u₊),
        halign = :left, width = Auto(false), fontsize = tsize)
    controllayout[3, 1] = grid!([text_p₋ text_p₊ ; text_u₋ text_u₊])
    ⬜p₋[], ⬜p₊[], ⬜u₋[], ⬜u₊[] = rand(4)
    return nslider, Tslider, dslider,
           nslider.value, Tslider.value, dslider.value, αslider.value,
           imenu.i_selected, ▢update.clicks, ▢back.clicks, ▢reset.clicks,
           ⬜p₋, ⬜p₊, ⬜u₋, ⬜u₊
end

function od_sliders!(figure, controllayout)
    sliders = []
    for (i, (l, vals)) in enumerate(zip(("n =", "t =", "d =", "α ="),
                         (1000:100:100000, 1000:1000:100000, 100:100:10000, 0.001:0.001:1)))

        startvalue = l[1] == 'α' ? 0.1 : l[1] == 'd' ? 1000 : vals[1]
        sll = labelslider!(figure, l, vals; sliderkw = Dict(:startvalue => startvalue))
        push!(sliders, sll.slider)
        controllayout[i, :] = sll.layout
    end
    return sliders
end
