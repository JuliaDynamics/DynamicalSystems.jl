__precompile__()

"""
A Julia package for the exploration of continuous
and discrete dynamical systems.
"""
module DynamicalSystems

"""
    DynamicalSystem
Abstract type representing a dynamical system. Has the following concrete sub-types:
* `DiscreteDS1D`
* `DiscreteDS`
* `ContinuousDS`
"""
abstract type DynamicalSystem end

export DynamicalSystem, Systems
# System definition and evolution:
include("discrete.jl")
include("continuous.jl")
include("famous_systems.jl")

# Lyapunovs:
include("lyapunovs.jl")

# Entropies and Dimension Estimation:
include(joinpath("dimensions", "entropies.jl"))
include(joinpath("dimensions", "dims.jl"))

end # module
