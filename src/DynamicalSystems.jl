module DynamicalSystems

export DynamicalSystem

abstract type DynamicalSystem end
include("helpful_functions.jl")
include("discrete.jl")
#include("continuous.jl")
include("lyapunovs.jl")
include("famous_systems.jl")



end # module


using DynamicalSystems
ds = DynamicalSystems.Systems.towel()
