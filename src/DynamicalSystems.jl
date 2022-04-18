"""
A Julia suite for chaos and nonlinear dynamics
"""
module DynamicalSystems

using Reexport

@reexport using DelayEmbeddings
@reexport using DynamicalSystemsBase
@reexport using Entropies
@reexport using ChaosTools
@reexport using RecurrenceAnalysis

# Also export some static array stuff
using DelayEmbeddings.StaticArrays
export SVector, SMatrix, @SVector, @SMatrix, Size

# Update messages:
using Scratch
display_update = true
version_number = "2.3.0"
update_name = "update_v$(version_number)"

if display_update
    # Get scratch space for this package
    versions_dir = @get_scratch!("versions")
    if !isfile(joinpath(versions_dir, update_name))
        printstyled(
            stdout,
            """
            \nUpdate message: DynamicalSystems v$(version_number)\n
            Interactive GUI for exploring dynamical systems are now in the documentation.
            Made with Makie.jl + InteractiveDynamics.jl!
            """;
            color = :light_magenta,
        )
        touch(joinpath(versions_dir, update_name))
    end
end


end # module
