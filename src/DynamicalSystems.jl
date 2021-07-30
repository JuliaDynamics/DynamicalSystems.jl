"""
A Julia suite for chaos and nonlinear dynamics
"""
module DynamicalSystems

using StaticArrays

export SVector, SMatrix, @SVector, @SMatrix, Size

using Reexport

@reexport using DelayEmbeddings
@reexport using DynamicalSystemsBase
@reexport using Entropies
@reexport using ChaosTools
@reexport using RecurrenceAnalysis

# Update messages:
using Scratch
display_update = true
version_number = "2.0"
update_name = "update_v$(version_number)"

if display_update
    # Get scratch space for this package
    versions_dir = @get_scratch!("versions")
    if !isfile(joinpath(versions_dir, update_name))
        printstyled(
            stdout,
            """
            \nUpdate message: DynamicalSystems v$(version_number)
            Welcome to this new (and breaking) release of DynamicalSystems.jl! 
            All features, old and new, can be found in the online docs.
            The *breaking* changes in this release are:
            1. Keyword propagation to DifferentialEquations.jl does not happen anymore.
               E.g.: `trajectory(ds, 100; reltol = 1e-9, alg = Tsit5())` is invalid.
               Now keyword propagation to DiffEq happens always via the keyword
               `diffeq`, whose value is a `NamedTuple/Dict` (the keyword-value pairs).
               E.g.: `trajectory(ds, 100; diffeq = (reltol = 1e-9, alg = Tsit5())`.
               Notice that this change affects a large list of fuctions of DynamicalSystems.jl!
            2. If `A` is a `Dataset` then `A[range_of_integers]` now returns a `Dataset`.
               Before it used to return `Vector{SVector}`.
            """;
            color = :light_magenta,
        )
        touch(joinpath(versions_dir, update_name))
    end
end


end # module
