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
update_name = "update_v1.2.0"

if display_update
if !isfile(joinpath(@__DIR__, update_name))
printstyled(stdout,
"""
\nUpdate message: DynamicalSystems v1.2

The default solver for continuous systems has changed
from `Vern9` to `SimpleATsit5` from SimpleDiffEq.jl.

This drops the OrdinaryDiffEq dependency, reduces precompilation times
and significantly improves first run times!

Even though not an API change, the numeric results you obtain
(in case you used the default solver) will change slightly.
Please be aware of this!

You can do `using OrdinaryDiffEq` to access the previous solver
and pass keyword `alg = Vern9()` to use it.\n
"""; color = :light_magenta)
touch(joinpath(@__DIR__, update_name))
end
end

end # module
