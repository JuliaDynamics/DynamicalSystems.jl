ChaosTools = DynamicalSystems.ChaosTools

function DynamicalSystems.interactive_poincaresos(ds, plane, idxs, complete;
        # PSOS kwargs:
        direction = -1,
        tfinal = (1000.0, 10.0^4),
        rootkw = (xrtol = 1e-6, atol = 1e-6),
        # Makie kwargs:
        color = randomcolor,
        scatterkwargs = (),
        labels = ("u₁", "u₂")
    )

    # Basic sanity checks on the method arguments
    @assert typeof(plane) <: Tuple
    @assert length(idxs) == 2
    @assert eltype(idxs) == Int
    @assert plane[1] ∉ idxs
    u0 = DynamicalSystems.get_state(ds)

    i = DynamicalSystems.SVector{2, Int}(idxs)

    # Construct a new `PoincareMap` structure with the given parameters
    pmap = DynamicalSystems.DynamicalSystemsBase.PoincareMap(ds, plane;
        direction, u0, rootkw, Tmax = tfinal[2])

    z = plane[2] # third variable comes from plane
    return interactive_clicker(pmap;
        tfinal = tfinal,
        complete = (x, y) -> complete(x, y, z),
        project = tr -> tr[:, i],
        color = color,
        scatterkwargs = scatterkwargs,
        labels = labels
    )
end

