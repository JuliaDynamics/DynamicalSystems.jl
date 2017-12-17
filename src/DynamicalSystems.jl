__precompile__()

"""
A Julia suite for chaos and nonlinear dynamics
"""
module DynamicalSystems

using Reexport

@reexport using DynamicalSystemsDef
@reexport using ChaosTools

end # module
