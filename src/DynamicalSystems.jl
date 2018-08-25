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

update_name = "update_v1.0.0"

if !isfile(joinpath(@__DIR__, update_name))
printstyled(stdout,
"""
\nUpdate message: DynamicalSystems v1.0

First major release of DynamicalSystems is out!

Minor breaking changes involve renaming `Reconstruction`
to `reconstruct`. Many internals have also been reworked
for extra clarity! Please see the "News" page of the
official documentation!\n
"""; color = :light_magenta)
touch(joinpath(@__DIR__, update_name))
end

end # module
