module DynamicalSystems

export DynamicalSystem

abstract type DynamicalSystem end

include("continuous.jl")
include("discrete.jl")
include("famous_systems.jl")
end # module
