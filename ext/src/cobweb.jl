function DynamicalSystems.interactive_cobweb(
    ds, prange, O::Int = 3;
    fkwargs = [(linewidth = 4.0, color = randomcolor()) for i in 1:O],
    trajcolor = :black,
    pname = "p",
    xmin = 0.0,
    xmax = 1.0,
    Tmax = 1000,
    pindex = 1,
    x0s = range(xmin, xmax; length = 101),
    xs = range(xmin, xmax; length = 1000),
)

@assert O ≥ 1
@assert dimension(ds) == 1
@assert !DynamicalSystems.isinplace(ds)

figure = Figure(resolution = (1000, 800))
axts = figure[1, :] = Axis(figure)
axmap = figure[2, :] = Axis(figure)

sg = SliderGrid(figure[3, :],
    (label = pname, range = prange),
    (label = "n", range = 1:Tmax, startvalue = 20),
    (label = "x₀", range = x0s, startvalue = initial_state(ds)[1]),
)
r_observable, L, x0 = [sg.sliders[i].value for i in 1:3]

# Timeseries plot
function seriespoints(x)
    n = 0:length(x)+1
    return [Point2f(n[i], x[i]) for i in 1:length(x)]
end

function get_float_trajectory(ds, L, x0)
    y = x0 isa Real ? SVector(x0) : x0
    z = trajectory(ds, L, y)[1]
    w = reinterpret(Float64, vec(z))
    return Float32.(w)
end

x = Observable(get_float_trajectory(ds, L[], x0[]))
xn = lift(a -> seriespoints(a), x)
lines!(axts, xn; color = trajcolor, lw = 2.0)
scatter!(axts, xn; color = trajcolor, markersize = 10)
xlims!(axts, 0, 20) # this is better than autolimits
ylims!(axts, xmin, xmax)

# Cobweb diagram
set_parameter!(ds, pindex, prange[1])

fobs = Any[Observable(
    map(x -> dynamic_rule(ds)(x, current_parameters(ds), 0)[1], xs)
)]
for order in 2:O
    push!(fobs, Observable(
        map(x -> dynamic_rule(ds)(x, current_parameters(ds), 0)[1], fobs[end][])
        # dynamic_rule(ds).(fobs[end][], Ref(current_parameters(ds)), 0)
    ))
end

# plot diagonal and fⁿ
lines!(axmap, [xmin,xmax], [xmin,xmax]; linewidth = 2, color = :gray, linestyle=:dash)
fcurves = Any[]
for i in 1:O
    _c = lines!(axmap, xs, fobs[i]; fkwargs[i]...)
    push!(fcurves, _c)
end

cobs = lift(a -> cobweb(a), x)
ccurve = lines!(axmap, cobs; color = (trajcolor, 0.5))
# cscatter = scatter!(axmap, cobs; color = trajcolor, markersize = 2)

xlims!(axmap, xmin, xmax)
ylims!(axmap, xmin, xmax)

# On trigger r-slider update all plots:
on(r_observable) do r
    set_parameter!(ds, pindex, r)
    x[] = get_float_trajectory(ds, L[], x0[])
    fobs[1][] = map(x -> dynamic_rule(ds)(x, current_parameters(ds), 0)[1], xs)
    for order in 2:O
        fobs[order][] = map(
            x -> dynamic_rule(ds)(x, current_parameters(ds), 0)[1], fobs[order-1][]
        )
    end
end

on(L) do l
    res = get_float_trajectory(ds, l, x0[])
    x[] = res
    xlims!(axts, 0, l) # this is better than autolimits
end

on(x0) do u
    x[] = get_float_trajectory(ds, L[], u)
end

# Finally add buttons to hide/show elements of the plot
cbutton = Button(figure; label = "cobweb")
fbuttons = Any[]
for i in 1:O
    _b = Button(figure; label = "f$(superscript(i))")
    push!(fbuttons, _b)
end
figure[4, :] = buttonlayout = GridLayout(tellwidth = false)
buttonlayout[:, 1:O+1] = [cbutton, fbuttons...]

# And add triggering for buttons
on(cbutton.clicks) do click
    ccurve.attributes.visible = !(ccurve.attributes.visible[])
    # cscatter.attributes.visible = !(cscatter.attributes.visible[])
end

for i in 1:O
    on(fbuttons[i].clicks) do click
        fcurves[i].attributes.visible = !(fcurves[i].attributes.visible[])
    end
end

display(figure)
return figure
end

function cobweb(t) # transform timeseries x into cobweb (point2D)
    # TODO: can be optimized to become in-place instead of pushing
    c = Point2f[]
    sizehint!(c, 2length(t))
    push!(c, Point2f(t[1], 0))
    for i ∈ 1:length(t)-1
        push!(c, Point2f(t[i], t[i]))
        push!(c, Point2f(t[i], t[i+1]))
    end
    return c
end