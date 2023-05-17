module DynamicalSystemsVisualizations

using DynamicalSystems, Makie
using DynamicalSystemsBase.SciMLBase

# TODO: Somehow extract this from an online repo...?
COLORS = [
    "#7143E0",
    "#191E44",
    "#0A9A84",
    "#AF9327",
    "#791457",
    "#6C768C",
]

include("src/trajectory_observable.jl")

include("src/cobweb.jl")
include("src/orbitdiagram.jl")
include("src/poincareclick.jl")
include("src/brainscan.jl")

end