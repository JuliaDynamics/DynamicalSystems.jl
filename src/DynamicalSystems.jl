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
version_number = "2.2.1"
update_name = "update_v$(version_number)"

if display_update
    # Get scratch space for this package
    versions_dir = @get_scratch!("versions")
    if !isfile(joinpath(versions_dir, update_name))
        printstyled(
            stdout,
            """
            \nUpdate message: DynamicalSystems v$(version_number)\n
            Brand new `AttractorMapper` infastructure that provides a generic framework
            for mapping initial conditions to attractors, and hence calculating basins of
            attraction. See documentation online for more!
            """;
            color = :light_magenta,
        )
        touch(joinpath(versions_dir, update_name))
    end
end


end # module
