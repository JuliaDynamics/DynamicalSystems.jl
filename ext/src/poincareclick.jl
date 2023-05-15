export interactive_poincaresos
ChaosTools = DynamicalSystems.ChaosTools

"""
    interactive_poincaresos(cds, plane, idxs, complete; kwargs...)
Launch an interactive application for exploring a Poincaré surface of section (PSOS)
of the continuous dynamical system `cds`. Requires `DynamicalSystems`.

The `plane` can only be the `Tuple` type accepted by `DynamicalSystems.poincaresos`,
i.e. `(i, r)` for the `i`th variable crossing the value `r`. `idxs` gives the two
indices of the variables to be displayed, since the PSOS plot is always a 2D scatterplot.
I.e. `idxs = (1, 2)` will plot the 1st versus 2nd variable of the PSOS. It follows
that `plane[1] ∉ idxs` must be true.

`complete` is a three-argument **function** that completes the new initial state
during interactive use, see below.

The function returns: `figure, laststate` with the latter being
an observable containing the latest initial `state`.

## Keyword Arguments
* `direction, rootkw` : Same use as in `DynamicalSystems.poincaresos`.
* `tfinal = (1000.0, 10.0^4)` : A 2-element tuple for the range of values
  for the total integration time (chosen interactively).
* `color` : A **function** of the system's initial condition, that returns a color to
  plot the new points with. The color must be `RGBf/RGBAf`.
   A random color is chosen by default.
* `labels = ("u₁" , "u₂")` : Scatter plot labels.
* `scatterkwargs = ()`: Named tuple of keywords passed to `scatter`.
* `diffeq = NamedTuple()` : Any extra keyword arguments are passed into `init` of DiffEq.

## Interaction
The application is a standard scatterplot, which shows the PSOS of the system,
initially using the system's `u0`. Two sliders control the total evolution time
and the size of the marker points (which is always in pixels).

Upon clicking within the bounds of the scatter plot your click is transformed into
a new initial condition, which is further evolved and its PSOS is computed and then
plotted into the scatter plot.

Your click is transformed into a full `D`-dimensional initial condition through
the function `complete`. The first two arguments of the function are the positions
of the click on the PSOS. The third argument is the value of the variable the PSOS
is defined on. To be more exact, this is how the function is called:
```julia
x, y = mouseclick; z = plane[2]
newstate = complete(x, y, z)
```
The `complete` function can throw an error for ill-conditioned `x, y, z`.
This will be properly handled instead of breaking the application.
This `newstate` is also given to the function `color` that
gets a new color for the new points.
"""
function interactive_poincaresos(ds, plane, idxs, complete;
        # PSOS kwargs:
        direction = -1,
        tfinal = (1000.0, 10.0^4),
        rootkw = (xrtol = 1e-6, atol = 1e-6),
        # Makie kwargs:
        color = randomcolor,
        scatterkwargs = (),
        labels = ("u₁" , "u₂"),
        # DiffEq kwargs:
        diffeq = NamedTuple()
    )

    error("this function has not yet been updated to DynamicalSystems.jl v3.0. PR welcomed!")

    @assert typeof(plane) <: Tuple
    @assert length(idxs) == 2
    @assert eltype(idxs) == Int
    @assert plane[1] ∉ idxs
    u0 = DynamicalSystems.get_state(ds)

    # This is the low-level call of poincaresos:
    DynamicalSystems.DynamicalSystemsBase.check_hyperplane_match(plane, DynamicalSystems.dimension(ds))
    integ = DynamicalSystems.integrator(ds, u0; diffeq)
    planecrossing = ChaosTools.PlaneCrossing(plane, direction > 0)
    i = DynamicalSystems.SVector{2, Int}(idxs)

    figure = Figure(resolution = (1000, 800), backgroundcolor = DEFAULT_BG)

    T_slider, m_slider = _add_psos_controls!(figure, tfinal)
    ax = figure[0, :] = Axis(figure)

    # Initial Section
    f = (t) -> planecrossing(integ(t))
    data = ChaosTools._poincaresos(integ, f, planecrossing, T_slider[], i, rootkw)
    length(data) == 0 && error(ChaosTools.PSOS_ERROR)

    positions_node = Observable(data)
    colors = (c = color(u0); [c for i in 1:length(data)])
    colors_node = Observable(colors)
    scatter!(
        ax, positions_node, color = colors_node,
        markersize = lift(o -> o*px, m_slider), marker = MARKER, scatterkwargs...
    )

    ax.xlabel, ax.ylabel = labels
    laststate = Observable(u0)

    # Interactive clicking on the psos:
    MakieLayout.deactivate_interaction!(ax, :rectanglezoom)
    spoint = select_point(ax.scene)
    on(spoint) do pos
        x, y = pos; z = plane[2] # third variable comes from plane
        newstate = try
           complete(x, y, z)
        catch err
           @error "Could not get state, got error: " exception=err
           return
        end

        DynamicalSystems.reinit!(integ, newstate)
        data = ChaosTools._poincaresos(integ, f, planecrossing, T_slider[], i, rootkw)
        positions = positions_node[]; colors = colors_node[]
        append!(positions, data)
        c = color(newstate)
        append!(colors, fill(c, length(data)))
        # update all the observables with Array as value:
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
