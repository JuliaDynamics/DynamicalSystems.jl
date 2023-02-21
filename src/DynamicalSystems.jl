"""
**DynamicalSystems.jl** is an award-winning Julia software library for nonlinear
dynamics and nonlinear timeseries analysis. The _module_ `DynamicalSystems`
builds the documentation and exports all packages composing DynamicalSystems.jl.

To install it, run `import Pkg; Pkg.add("DynamicalSystems")`.

All further information is provided in the documentation,
which you can either find [online](https://juliadynamics.github.io/DynamicalSystems.jl/dev/)
or build locally by running the `docs/make.jl` file.
"""
module DynamicalSystems

using Reexport

# Core
@reexport using StateSpaceSets
@reexport using DynamicalSystemsBase
# observed/measured data
@reexport using RecurrenceAnalysis
@reexport using FractalDimensions
@reexport using DelayEmbeddings
@reexport using ComplexityMeasures
# dynamical systems
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
            Interactive GUI for exploring dynamical systems are now in the documentation.
            Made with Makie.jl + InteractiveDynamics.jl!
            """;
            color = :light_magenta,
        )
        touch(joinpath(versions_dir, update_name))
    end
end


end # module
