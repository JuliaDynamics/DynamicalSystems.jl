using DynamicalSystems, Makie, LinearAlgebra
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

            reinit!(integ, newstate)

            data = ChaosTools._initialize_output(integ.u, i)
            @time ChaosTools.poincare_cross!(
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

@inline function ż(z, p, t)
    @inbounds begin
        A = p.A; B=p.B; D=p.D
        p₀, p₂ = z[3], z[4]
        q₀, q₂ = z[1], z[2]

        return SVector{4}(
            A * p₀,
            A * p₂,
            -A * q₀ - 3 * B / √2 * (q₂^2 - q₀^2) - D * q₀ * (q₀^2 + q₂^2),
            -q₂ * (A + 3 * √2 * B * q₀ + D * (q₀^2 + q₂^2))
        )
    end
end

z0 = SVector{4}(0.0, -2.723, 2.349, -9.801)
params = (A=1, B=0.55, D=0.4)
ds = ContinuousDynamicalSystem(ż, z0, params)

function T(p, A)
    A / 2 * norm(p)^2
end
function potential(q, params)
    @unpack A, B, D = params
    A / 2 * (q[1]^2 + q[2]^2) + B / √2 * q[1] * (3 * q[2]^2 - q[1]^2) + D / 4 * (q[1]^2 + q[2]^2)^2
end

H(p, q, params=(A=1, B=0.55, D=0.4)) = T(p, params.A) + potential(q, params)

const E = H([ds.u0[3], ds.u0[4]],[ds.u0[1], ds.u0[2]], params)

function complete(y, py, x)
    V = potential((x, y), params)
    Ky = 0.5*(py^2)
    Ky + V ≥ E && error("Point has more energy!")
    px = sqrt(2(E - V - Ky))
    ic = [x, y, px, py]
    println("new state:")
    println(ic)
    return ic
end

_randomcolor(args...) = RGBf0(rand(3)...)

chaotic = get_state(ds)
stable = [0., 0.1, 0.5, 0.]

plane = (1, 0.0)
tf = 5000.0

psos = interactivepsos(ds, plane, tf, (2, 4), complete; makiekwargs = (markersize = 0.1,))
