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
version_number = "2.2.0"
update_name = "update_v$(version_number)"

if display_update
    # Get scratch space for this package
    versions_dir = @get_scratch!("versions")
    if !isfile(joinpath(versions_dir, update_name))
        printstyled(
            stdout,
            """
            \nUpdate message: DynamicalSystems v$(version_number)\n
            Direct keyword acceptance of DifferentialEquations.jl related keywords is
            no longer supported in DynamicalSystems.jl. This means that this is deprecated:
            ```
            lyapunov(ds, T; alg = Tsit5(), abstol = 1e-9)
            ```
            Instead, use the explicit keyword `diffeq`, and give it a named tuple with
            the keyword arguments, i.e.:
            ```
            lyapunov(ds, T; diffeq = (alg = Tsit5(), abstol = 1e-9))
            ```

            Also, the online documentation now uses Makie.jl,
            references InteractiveDynamics.jl and showcases some of its functions.
            """;
            color = :light_magenta,
        )
        touch(joinpath(versions_dir, update_name))
    end
end


end # module
