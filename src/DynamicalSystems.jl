"""
A Julia suite for chaos and nonlinear dynamics
"""
module DynamicalSystems

using StaticArrays

export SVector, SMatrix, @SVector, @SMatrix, Size

using Reexport

@reexport using DelayEmbeddings
@reexport using DynamicalSystemsBase
@reexport using ChaosTools
@reexport using RecurrenceAnalysis

display_update = true
update_name = "update_v1.1.0"

if display_update
if !isfile(joinpath(@__DIR__, update_name))
printstyled(stdout,
"""
\nUpdate message: DynamicalSystems v1.1

A new package has joined DynamicalSystems: RecurrenceAnalysis !

This new package offers tools to compute and analyze
recurrences in your timeseries!

Check out the documentation for more!\n
"""; color = :light_magenta)
touch(joinpath(@__DIR__, update_name))
end
end

end # module
