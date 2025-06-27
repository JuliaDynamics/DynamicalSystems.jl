module DynamicalSystemsVisualizations

using DynamicalSystems, Makie
using DynamicalSystemsBase.SciMLBase

include("src/utils.jl")
include("src/dynamicalsystemobservable.jl")
include("src/interactive_trajectory.jl")
include("src/cobweb.jl")
include("src/orbitdiagram.jl")
include("src/poincareclick.jl")
include("src/brainscan.jl")
include("src/clicker.jl")

subscript = DynamicalSystemsVisualizations.subscript

end
