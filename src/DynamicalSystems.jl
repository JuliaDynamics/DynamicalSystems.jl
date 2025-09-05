module DynamicalSystems

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end DynamicalSystems

using Reexport

# Core
@reexport using StateSpaceSets
@reexport using DynamicalSystemsBase
# observed/measured data
@reexport using RecurrenceAnalysis
@reexport using FractalDimensions
@reexport using DelayEmbeddings
@reexport using ComplexityMeasures
@reexport using TimeseriesSurrogates
@reexport using SignalDecomposition
# dynamical systems
@reexport using PredefinedDynamicalSystems
@reexport using Attractors
@reexport using ChaosTools
@reexport using PeriodicOrbits
# visualizations (singleton methods for package extension)
using DataStructures
include("visualizations.jl")

end # module
