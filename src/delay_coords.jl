using NearestNeighbors, Requires, StaticArrays
using LsqFit: curve_fit
using StatsBase: autocor
using Distances: Metric, Cityblock, Euclidean
import NearestNeighbors: KDTree

export reconstruct, Cityblock, Euclidean, AbstractNeighborhood
export FixedMassNeighborhood, FixedSizeNeighborhood, numericallyapunov
export estimate_delay, neighborhood, Reconstruction, KDTree
#####################################################################################
#                            Reconstruction Object                                  #
#####################################################################################
"""
    Reconstruction{D, T, τ} <: AbstractDataset{D, T}
`D`-dimensional delay-coordinates reconstruction object with delay `τ`,
created from a timeseries `s` with `T` type numbers.

Use `Reconstruction(s::AbstractVector{T}, D, τ)` to create an instance.

## Description
The ``n``th row of a `Reconstruction` is the `D`-dimensional vector
```math
(s(n), s(n+\\tau), s(n+2\\tau), \\dots, s(n+(D-1)\\tau))
```

The reconstruction object `R` can have same
invariant quantities (like e.g. lyapunov exponents) with the original system
that the timeseries were recorded from, for proper `D` and `τ` [1, 2].

`R` can be accessed similarly to a [`Dataset`](@ref):
```julia
s = rand(1e6)
R = Reconstruction(s, 4, 1) # dimension 4 and delay 1
R[3] # third point of reconstruction, ≡ (s[3], s[4], s[5], s[6])
R[1, 2] # Second element of first point of reconstruction, ≡ s[2]
```
and can also be given to all functions that accept a `Dataset`,
like [`generalized_dim`](@ref).

The functions `dimension(R)` and `delay(R)` return `D` and `τ` respectively.

## References

[1] : F. Takens, *Detecting Strange Attractors in Turbulence — Dynamical
Systems and Turbulence*, Lecture Notes in Mathematics **366**, Springer (1981)

[2] : T. Sauer *et al.*, J. Stat. Phys. **65**, pp 579 (1991)
"""
type Reconstruction{D, T<:Number, τ} <: AbstractDataset{D, T}
    data::Vector{SVector{D,T}}
end

Reconstruction(s::AbstractVector{T}, D, τ) where {T} =
Reconstruction{D, T, τ}(reconstruct(s, Val{D}(), τ))

@inline Base.eltype(::Reconstruction{D, T, t}) where {T,D,t} = T
@inline delay(::Reconstruction{D, T, t}) where {T,D,t} = t

function reconstruct_impl(::Type{Val{D}}) where D
    gens = [:(s[i + $k*τ]) for k=0:D-1]

    quote
        L = length(s) - ($(D-1))*τ;
        T = eltype(s)
        data = Vector{SVector{$D, T}}(L)
        for i in 1:L
            data[i] = SVector{$D,T}($(gens...))
        end
        V = typeof(s)
        T = eltype(s)
        data
    end
end

@generated function reconstruct(s::AbstractVector{T}, ::Val{D}, τ) where {D, T}
    reconstruct_impl(Val{D})
end


# Pretty print:
matname(d::Reconstruction{D, T, τ}) where {D, T, τ} =
"$(D)-dimensional delay coordinates reconstruction with delay τ=$(τ)"


#####################################################################################
#                      Estimate Reconstruction Parameters                           #
#####################################################################################
"""
    localextrema(y) -> max_ind, min_ind
Find the local extrema of given array `y`, by scanning point-by-point. Return the
indices of the maxima (`max_ind`) and the indices of the minima (`min_ind`).
"""
function localextrema end
@inbounds function localextrema(y)
    l = length(y)
    i = 1
    maxargs = Int[]
    minargs = Int[]
    if y[1] > y[2]
        push!(maxargs, 1)
    elseif y[1] < y[2]
        push!(minargs, 1)
    end

    for i in 2:l-1
        left = i-1
        right = i+1
        if  y[left] < y[i] > y[right]
            push!(maxargs, i)
        elseif y[left] > y[i] < y[right]
            push!(minargs, i)
        end
    end

    if y[l] > y[l-1]
        push!(maxargs, l)
    elseif y[l] < y[l-1]
        push!(minargs, l)
    end
    return maxargs, minargs
end


function exponential_decay_extrema(c::AbstractVector)
    ac = abs.(c)
    ma, mi = localextrema(ac)
    # ma start from 1 but correlation is expected to start from x=0
    ydat = ac[ma]; xdat = ma .- 1
    # Do curve fit from LsqFit
    model(x, p) = @. exp(-x/p[1])
    decay = curve_fit(model, xdat, ydat, [1.0]).param[1]
    return decay
end

function exponential_decay(c::AbstractVector)
    # Do curve fit from LsqFit
    model(x, p) = @. exp(-x/p[1])
    decay = curve_fit(model, 0:length(c)-1, abs.(c), [1.0]).param[1]
    return decay
end

"""
    estimate_delay(s) -> τ
Estimate an optimal delay to be used in [`reconstruct`](@ref),
by performing an exponential fit to
the `abs.(c)` with `c` the auto-correlation function of `s`.
Return the exponential decay time `τ` rounded to an integer.
"""
function estimate_delay(x::AbstractVector)
    c = autocor(x, 0:length(x)÷10)
    i = 1
    # Find 0 crossing:
    while c[i] > 0
        i+= 1
        i == length(c) && break
    end
    # Find exponential fit:
    τ = exponential_decay(c)
    # Is there a method to deduce which one of the 2 is the better approach?
    return round(Int, τ)
end

function estimate_dimension(s::AbstractVector)
  # Estimate number of “false nearest neighbors” due to
  # projection into a too low dimension reconstruction space
end

#####################################################################################
#                    Numerical Lyapunov (from Reconstruction)                       #
#####################################################################################
# Everything in this section is based on Ulrich Parlitz [1]




# Neighborhoods:
"""
    AbstractNeighborhood
Supertype of methods for deciding the neighborhood of points for a given point.

Concrete subtypes:
* `FixedMassNeighborhood(K::Int)`  : The neighborhood of a point consists of the `K`
  nearest neighbors of the point.
* `FixedSizeNeighborhood(ε::Real)` : The neighborhood of a point consists of all
  neighbors that have distance < `ε` from the point.

Notice that these distances are always computed using the `Euclidean()` distance
in `D`-dimensional space.

See also [`neighborhood`](@ref) or [`numericallyapunov`](@ref).
"""
abstract type AbstractNeighborhood end
struct FixedMassNeighborhood <: AbstractNeighborhood
    K::Int
end
FixedMassNeighborhood() = FixedMassNeighborhood(1)
struct FixedSizeNeighborhood <: AbstractNeighborhood
    ε::Float64
end
FixedSizeNeighborhood() = FixedSizeNeighborhood(0.001)

"""
    neighborhood(n, point, tree::KDTree, method::AbstractNeighborhood)
Return a vector of indices which are the neighborhood of `point`, whose index
in the original data is `n`.

If the original data is `data <: AbstractDataset`, then
use `tree = KDTree(data)` to obtain the `tree` instance (which also
contains a copy of the data).
Both `point` and `n` must be provided because the
`tree` has indices in different sorting.

The `method` can be a subtype of [`AbstractNeighborhood`](@ref).

`neighborhood` works for *any* subtype of `AbstractDataset`, for example
```julia
R = some_dataset
tree = KDTree(R)
neigh = neighborhood(n, R[n], tree, method)
```

## References

`neighborhood` simply interfaces the functions
`knn` and `inrange` from
[NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) by using
the last argument, `method`.
"""
function neighborhood(
    n, point, tree::KDTree, method::FixedMassNeighborhood)
    idxs, = knn(tree, point, method.K, false, i -> i==n)
    return idxs
end
function neighborhood(
    n, point, tree::KDTree, method::FixedSizeNeighborhood)
    idxs = inrange(tree, point, method.ε)
    deleteat!(idxs, findin(idxs, n)) # unfortunately this has to be done...
    return idxs
end
neighborhood(n, point, tree::KDTree) =
neighborhood(n, point, tree, FixedMassNeighborhood(1))

KDTree(D::AbstractDataset) = KDTree(D.data, Euclidean())



"""
```julia
numericallyapunov(R::Reconstruction, ks;  refstates, distance, method)
```
Return `E = [E(k) for k ∈ ks]`, where `E(k)` is the average logarithmic distance for
nearby states that are evolved in time for `k` steps (`k` must be integer).

## Keyword Arguments

* `refstates::AbstractVector{Int} = 1:(length(R) - ks[end])` : Vector of indices
  that notes which
  states of the reconstruction should be used as "reference states", which means
  that the algorithm is applied for all state indices contained in `refstates`.
* `method::AbstractNeighborhood = FixedMassNeighborhood(1)` : The method to
  be used when evaluating
  the neighborhood of each reference state. See
  [`AbstractNeighborhood`](@ref) or [`neighborhood`](@ref) for more info.
* `distance::Metric = Cityblock()` : The distance function used in the
  logarithmic distance of nearby states. The allowed distances are `Cityblock()`
  and `Euclidean()`. See below for more info.


## Description
If the reconstruction
exhibits exponential divergence of nearby states, then it should clearly hold:
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

The algorithm used in this function is due to Parlitz [1], which itself
expands upon Kantz [2]. In sort, for
each reference state a neighborhood is evaluated. Then, for each point in this
neighborhood, the logarithmic distance between reference state and neighborhood
state is
calculated as the "time" index `k` increases. The average of the above over
all neighborhood states over all reference states is the returned result.

If the `Metric` is `Euclidean()` then use the Euclidean distance of the
full `D`-dimensional points (distance ``d_E`` in ref. [1]).
If however the `Metric` is `Cityblock()`, calculate
the absolute distance of *only the first elements* of the `m+k` and `n+k` points
of the
reconstruction `R`(distance
``d_F`` in
ref. [1]). Notice that
the distances used are defined in the package
[Distances.jl](https://github.com/JuliaStats/Distances.jl), but are re-exported here
for ease-of-use.

This function assumes that the Theiler window (see [1]) is the same as the delay time:
``w  = \\tau``.

## References

[1] : Skokos, C. H. *et al.*, *Chaos Detection and Predictability* - Chapter 1
(section 1.3.2), Lecture Notes in Physics **915**, Springer (2016)

[2] : Kantz, H., Phys. Lett. A **185**, pp 77–87 (1994)
"""
function numericallyapunov(R::Reconstruction, ks;
                           refstates = 1:(length(R) - ks[end]),
                           distance = Cityblock(),
                           method = FixedMassNeighborhood(1))
    Ek = numericallyapunov(R, ks, refstates, distance, method)
end

function numericallyapunov(R::Reconstruction{D, T, τ},
                           ks::AbstractVector{Int},
                           ℜ::AbstractVector{Int},
                           distance::Metric,
                           method::AbstractNeighborhood) where {D, T, τ}

    # ℜ = \Re<tab> = set of indices that have the points that one finds neighbors.
    # n belongs in ℜ and R[n] is the "reference state".
    # Thus, ℜ contains all the reference states the algorithm will iterate over.
    # ℜ is not estimated. it is given by the user. Most common choice:
    # ℜ = 1:(length(R) - ks[end])

    # ⩅(n) = \Cup<tab> = neighborhood of reference state n
    # which is evaluated for each n and for the given neighborhood method

    # Initialize:
    timethres = length(R) - ks[end]
    if maximum(ℜ) > timethres
        erstr = "Maximum index of reference states is > length(R) - ks[end] "
        erstr*= "and the algorithm cannot be performed on it. You have to choose "
        erstr*= "reference state indices of at most up to length(R) - ks[end]."
        throw(ArgumentError(erstr))
    end
    E = zeros(T, length(ks))
    E_n = copy(E); E_m = copy(E)
    data = R.data #tree data
    tree = KDTree(data, Euclidean()) # this creates a copy of `data`
    skippedm = 0; skippedn = 0

    for n in ℜ
        # The ⋓(n) can be evaluated on the spot instead of being pre-calculated
        # for all reference states. (it would take too much memory)
        # Since U[n] doesn't depend on `k` one can then interchange the loops:
        # Instead of k being the outermost loop, it becomes the innermost loop!
        point = data[n]
        ⋓ = neighborhood(n, point, tree, method)
        for m in ⋓
            # If `m` is nearer to the end of the timeseries than k allows
            # is it completely skipped (and length(⋓) reduced).
            # If m is closer to n than the Theiler window allows, also skip.
            # It is assumed that w = τ (the Theiler window is the delay time)
            if m > timethres || abs(m - n) <= τ
                skippedm += 1
                continue
            end
            for (j, k) in enumerate(ks) #ks should be small (of order 10 to 100 MAX)
                E_m[j] = log(delay_distance(distance, R, m, n, k))
            end
            E_n .+= E_m # no need to reset E_m
        end
        if skippedm >= length(⋓)
            skippedn += 1
            skippedm = 0
            continue # be sure to continue if no valid point!
        end
        E .+= E_n ./ (length(⋓) - skippedm)
        skippedm = 0
        fill!(E_n, zero(T)) #reset distances for n reference state
    end
    #plot E[k] versus k and boom, you got lyapunov in the linear scaling region.
    if skippedn >= length(ℜ)
        ers = "skippedn == length(ℜ)\n"
        ers*= "Could happen because all the neighbors fall within the Theiler "
        ers*= "window. Fix: increase neighborhood size."
        error(ers)
    end
    E ./= length(ℜ) - skippedn
end

@inline @inbounds function delay_distance(di::Cityblock, R, m, n, k)
    abs(R[m+k][1] - R[n+k][1])
end

@inline @inbounds function delay_distance(di::Euclidean,
    R::Reconstruction{D, T, τ}, m, n, k) where {D, T, τ}
    return norm(R[m+k] - R[n+k])
end
