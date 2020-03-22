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

display_update = false
update_name = "update_v1.4.0"

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
@doc """
    numericallyapunov(R::Dataset, ks;  refstates, w, distance, ntype)
Return `E = [E(k) for k ∈ ks]`, where `E(k)` is the average logarithmic distance
between states of a [`neighborhood`](@ref)
that are evolved in time for `k` steps (`k` must be integer).
## Keyword Arguments
* `refstates = 1:(length(R) - ks[end])` : Vector of indices
  that notes which
  states of the reconstruction should be used as "reference states", which means
  that the algorithm is applied for all state indices contained in `refstates`.
* `w::Int = 1` : The Theiler window, which determines
  whether points are separated enough in time to be considered separate trajectories
  (see [^1] and [`neighborhood`](@ref)).
* `ntype::AbstractNeighborhood = FixedMassNeighborhood(1)` : The method to
  be used when evaluating the neighborhood of each reference state. See
  [`AbstractNeighborhood`](@ref) or [`neighborhood`](@ref) for more info.
* `distance::Metric = Cityblock()` : The distance function used in the
  logarithmic distance of nearby states. The allowed distances are `Cityblock()`
  and `Euclidean()`. See below for more info.
## Description
If the dataset/reconstruction
exhibits exponential divergence of nearby states, then it should clearly hold
```math
E(k) \\approx \\lambda\\Delta t k + E(0)
```
for a *well defined region* in the `k` axis, where ``\\lambda`` is
the approximated
maximum Lyapunov exponent. ``\\Delta t`` is the time between samples in the
original timeseries.
You can use [`linear_region`](@ref) with arguments `(ks .* Δt, E)` to
identify the slope
(= ``\\lambda``)
immediatelly, assuming you
have choosen sufficiently good `ks` such that the linear scaling region is bigger
than the saturated region.
The algorithm used in this function is due to Parlitz [^1], which itself
expands upon Kantz [^2]. In sort, for
each reference state a neighborhood is evaluated. Then, for each point in this
neighborhood, the logarithmic distance between reference state and neighborhood
state is
calculated as the "time" index `k` increases. The average of the above over
all neighborhood states over all reference states is the returned result.
If the `Metric` is `Euclidean()` then use the Euclidean distance of the
full `D`-dimensional points (distance ``d_E`` in ref. [^1]).
If however the `Metric` is `Cityblock()`, calculate
the absolute distance of *only the first elements* of the `m+k` and `n+k` points
of the reconstruction `R` (distance ``d_F`` in ref. [^1]).
## References
[^1] : Skokos, C. H. *et al.*, *Chaos Detection and Predictability* - Chapter 1
(section 1.3.2), Lecture Notes in Physics **915**, Springer (2016)
[^2] : Kantz, H., Phys. Lett. A **185**, pp 77–87 (1994)
""" ChaosTools.numericallyapunov

end # module
