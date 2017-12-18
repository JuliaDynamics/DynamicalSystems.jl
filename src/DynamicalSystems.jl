__precompile__()

"""
A Julia suite for chaos and nonlinear dynamics
"""
module DynamicalSystems

using Reexport

@reexport using DynamicalSystemsBase
@reexport using ChaosTools

end # module
