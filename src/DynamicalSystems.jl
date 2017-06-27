"""
A Julia package for the exploration of continuous
and discrete dynamical systems.
"""
module DynamicalSystems

"""
    DynamicalSystem
Abstract type representing a dynamical system. Has the following concrete sub-types:
* `DiscreteDS`
* `DiscreteDS1D`
* `ContinuousDS`
"""
abstract type DynamicalSystem end

export DynamicalSystem, Systems

include("helpful_functions.jl")
include("discrete.jl")
include("continuous.jl")
include("lyapunovs.jl")
include("famous_systems.jl")

end # module
