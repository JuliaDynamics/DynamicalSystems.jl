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


using DynamicalSystems, BenchmarkTools, StaticArrays, OrdinaryDiffEq
# lor = Systems.lorenz()
# sol = timeseries(lor, 1.0)
