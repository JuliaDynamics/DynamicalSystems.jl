"""
A Julia suite for chaos and nonlinear dynamics
"""
module DynamicalSystems

using StaticArrays

export SVector, SMatrix, @SVector, @SMatrix, Size

using Reexport

@reexport using DynamicalSystemsBase
@reexport using ChaosTools

# Visualization routines:
# using Requires
# @require PyPlot include("visualizations.jl")

end # module
