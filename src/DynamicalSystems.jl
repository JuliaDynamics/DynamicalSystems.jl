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
update_name = "update_v1.3.0"

if display_update
if !isfile(joinpath(@__DIR__, update_name))
printstyled(stdout,
"""
\nUpdate message: DynamicalSystems v1.3

A method that estimates the predictability properties of a
dynamical system has been implemented, following the work of:

Wernecke, H., Sándor, B. & Gros, C.
*How to test for partially predictable chaos*.

See the function `predictability`.\n
"""; color = :light_magenta)
touch(joinpath(@__DIR__, update_name))
end
end

end # module
