__precompile__()

"""
A Julia suite for chaos and nonlinear dynamics
"""
module DynamicalSystems

using StaticArrays

export SVector, SMatrix, @SVector, @SMatrix

using Reexport

@reexport using DynamicalSystemsBase
@reexport using ChaosTools

end # module
