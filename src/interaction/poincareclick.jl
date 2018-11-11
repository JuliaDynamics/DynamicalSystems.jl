using DynamicalSystems, Makie
using DynamicalSystemsBase: CDS

function interactivepsos(ds::CDS{IIP, S, D}, plane, tf, idxs, complete;
                         # PSOS kwargs:
                         u0 = get_state(ds), direction = +1, Ttr::Real = 0.0,
                         warning = true, rootkw = (xrtol = 1e-6, atol = 1e-6),
                         # Makie kwargs:
                         color = _randomcolor, resolution = (750, 750),
                         makiekwargs = (markersize = 0.005,),
                         # DiffEq kwargs:
                         diffeq...) where {IIP, S, D}

    @assert typeof(plane) <: Tuple
    @assert length(idxs) == 2
    @assert eltype(idxs) == Int

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
    scene = Makie.Scene(resolution = (750, 750))
    Makie.scatter!(scene, data; makiekwargs..., color = color(u0))

    # Interactive part:
    on(scene.events.mousebuttons) do buttons
        if ispressed(scene, Mouse.left)
            pos = to_world(scene, Point2f0(scene.events.mouseposition[]))

            x, y = pos; z = plane[2] # third variable comes from plane

            newstate = try
               complete(x, y, z)
            catch err
               @error "Could not set state, got error:" exception=err
               return
            end

            println("new state:")
            println(newstate)

            reinit!(integ, newstate)

            data = ChaosTools._initialize_output(integ.u, i)
            ChaosTools.poincare_cross!(
                data, integ, f, planecrossing, tf, Ttr, i, rootkw
            )

            Makie.scatter!(scene, data; makiekwargs..., color = color(newstate))
        end
        display(scene)
        return scene
    end
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
    return ic = [x, y, px, -py]
end

_randomcolor(args...) = RGBf0(rand(3)...)

chaotic = get_state(ds)
stable = [0., 0.1, 0.5, 0.]

plane = (1, 0.0)
tf = 20000.0


psos = interactivepsos(ds, plane, tf, (2, 4), complete)


#
#
# scene = Makie.Scene(resolution = (750, 750))
# psos = poincaresos(ds, plane, tf, force_dtmin=true, u0 = [0.0, 0.236638, 0.443869, 0.0762546])
# Makie.scatter!(scene, psos[:, 2], psos[:, 4], markersize = 0.01, color=RGBf0(rand(3)...))
#
# on(scene.events.mousebuttons) do buttons
#    if ispressed(scene, Mouse.left)
#        pos = to_world(scene, Point2f0(scene.events.mouseposition[]))
#
#         x, y = pos; z = plane[2] # third variable comes from plane
#
#         state = try
#            complete(x, y, z)
#        catch
#            println("Could not set state")
#            return
#        end
#
#        psos = poincaresos(ds, plane, tf; u0 = state, warning = false, force_dtmin=true)
#        Makie.scatter!(scene ,psos[:, 2], psos[:, 4], markersize = 0.01, color=RGBf0(rand(3)...))
#
#    end
#    display(scene)
# end
