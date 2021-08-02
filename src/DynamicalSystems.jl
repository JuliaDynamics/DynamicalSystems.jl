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
version_number = "2.0.1"
update_name = "update_v$(version_number)"

if display_update
    # Get scratch space for this package
    versions_dir = @get_scratch!("versions")
    if !isfile(joinpath(versions_dir, update_name))
        printstyled(
            stdout,
            """
            \nUpdate message: DynamicalSystems v$(version_number)\n
            Welcome to this new (and slightly breaking) release of DynamicalSystems.jl!
            You can find an announcement of the new version on the Julia discourse:

            https://discourse.julialang.org/t/nonlinear-dynamics-textbook-dynamicalsystems-jl-2-0/65665

            The *breaking* changes in this release are:
            1. The keyword `dt` of many functions has been renamed to `Δt`.
               This keyword had conflicts with the options of DifferentialEquations.jl.
               No warning can be thrown for this change, and users still using `dt` will
               have it silently propagated as keyword to the diffeq solvers.
               Functions affected: `trajectory, lyapunov, lyapunovspectrum, gali, expansionentropy, orbitdiagram`
            2. If `A` is a `Dataset` then `A[range_of_integers]` now returns a `Dataset`.
               Before it used to return `Vector{SVector}`.
            """;
            color = :light_magenta,
        )
        touch(joinpath(versions_dir, update_name))
    end
end


end # module
