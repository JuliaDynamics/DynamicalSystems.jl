module DynamicalSystems

abstract type DynamicalSystem end

include("continuous.jl")
#include("discrete.jl")
include("famous_systems.jl")
end # module
