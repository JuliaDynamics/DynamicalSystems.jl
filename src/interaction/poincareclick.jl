using DynamicalSystems, Makie
using DynamicalSystemsBase: CDS

function interactivepsos(ds::CDS{IIP, S, D}, plane, tf, idxs, complete;
                         # PSOS kwargs:
                         u0 = get_state(ds), direction = -1, Ttr::Real = 0.0,
                         warning = true, rootkw = (xrtol = 1e-6, atol = 1e-6),
                         # Makie kwargs:
                         color = _randomcolor, resolution = (750, 750),
                         makiekwargs = (markersize = 0.005,),
                         # DiffEq kwargs:
                         diffeq...) where {IIP, S, D}

    @assert typeof(plane) <: Tuple
    @assert length(idxs) == 2
    @assert eltype(idxs) == Int
    @assert !(plane[1] in idxs)

    # This is the internal code of poincaresos. We use the integrator directly!
    ChaosTools._check_plane(plane, D)
    integ = integrator(ds, u0; diffeq...)
    planecrossing = ChaosTools.PlaneCrossing{D}(plane, direction > 0 )
    f = (t) -> planecrossing(integ(t))

    i = SVector{2, Int}(idxs)
    data = ChaosTools._initialize_output(get_state(ds), i)

    ChaosTools.poincare_cross!(data, integ, f, planecrossing, tf, Ttr, i, rootkw)
    warning && length(data) == 0 && @warn ChaosTools.PSOS_ERROR

    # Create the first trajectory on the section:
    ui, ms = AbstractPlotting.textslider(10 .^ range(-6, stop=1, length=1000),
    "markersize", start=0.01)
    scene = Makie.Scene(resolution = (1500, 1000))
    positions_node = Node(data)
    colors = (c = color(u0); [c for i in 1:length(data)])
    colors_node = Node(colors)

    scplot = Makie.scatter(positions_node, color = colors_node, markersize = ms)

    to_screen(scene, mpos) = Point2f0(mpos) .- Point2f0(minimum(pixelarea(scene)[]))
    mouseposition(scene) = to_world(scene, to_screen(scene, events(scene).mouseposition[]))
    # Interactive part:
    on(scplot.events.mousebuttons) do buttons
        if (ispressed(scplot, Mouse.left) && !ispressed(scplot, Keyboard.space) &&
            AbstractPlotting.is_mouseinside(scplot))

            pos = mouseposition(scplot)

            x, y = pos; z = plane[2] # third variable comes from plane

            newstate = try
               complete(x, y, z)
            catch err
               @error "Could not get state, got error:" exception=err
               return
            end

            reinit!(integ, newstate)

            data = ChaosTools._initialize_output(integ.u, i)
            @time ChaosTools.poincare_cross!(
                data, integ, f, planecrossing, tf, Ttr, i, rootkw
            )

            positions = positions_node[]; colors = colors_node[]
            append!(positions, data)
            append!(colors, (c = color(newstate); [c for i in 1:length(data)]))

            # Notify the signals
            positions_node[] = positions; colors_node[] = colors

            # Makie.scatter!(scplot, data; makiekwargs..., color = color(newstate))
        end
        # display(scene)
        # return scene
    end
    hbox(ui, scplot, parent=scene)
    display(scene)
    return scene
end

ds = Systems.henonheiles()


potential(x, y) = 0.5(x^2 + y^2) + (x^2*y - (y^3)/3)
energy(x,y,px,py) = 0.5(px^2 + py^2) + potential(x,y)
const E = energy(get_state(ds)...)
function complete(y, py, x)
    V = potential(x, y)
    Ky = 0.5*(py^2)
    Ky + V ≥ E && error("Point has more energy!")
    px = sqrt(2(E - V - Ky))
    ic = [x, y, px, py]
    return ic
end

_randomcolor(args...) = RGBf0(rand(3)...)

chaotic = get_state(ds)
stable = [0., 0.1, 0.5, 0.]

plane = (1, 0.0)
tf = 20000.0

psos = interactivepsos(ds, plane, tf, (2, 4), complete; makiekwargs = (markersize = 0.01,))
