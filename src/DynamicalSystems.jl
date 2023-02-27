module DynamicalSystems

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end DynamicalSystems

using Reexport

# Core
@reexport using StateSpaceSets
@reexport using DynamicalSystemsBase
# observed/measured data
@reexport using RecurrenceAnalysis
@reexport using FractalDimensions
@reexport using DelayEmbeddings
@reexport using ComplexityMeasures
@reexport using TimeseriesSurrogates
# dynamical systems
@reexport using PredefinedDynamicalSystems
@reexport using Attractors
@reexport using ChaosTools

# Update messages:
using Scratch
display_update = true
version_number = "3.0.0"
update_name = "update_v$(version_number)"

if display_update
    # Get scratch space for this package
    versions_dir = @get_scratch!("versions")
    if !isfile(joinpath(versions_dir, update_name))
        printstyled(
            stdout,
            """
            \nUpdate message: DynamicalSystems v$(version_number)\n
            This is an exciting new MAJOR release with a plethora of improvements,
            a near-completely rewrite of many library internals,
            and some unfortunate breaking changes!

            Please see online changelog for more!
            """;
            color = :light_magenta,
        )
        touch(joinpath(versions_dir, update_name))
    end
end


end # module
