using NearestNeighbors, Requires, StaticArrays
using LsqFit: curve_fit
using StatsBase: autocor
using Distances: Metric, Cityblock, Euclidean
import NearestNeighbors: KDTree

export reconstruct, Cityblock, Euclidean, AbstractNeighborhood
export FixedMassNeighborhood, FixedSizeNeighborhood, numericallyapunov
export estimate_delay, neighborhood, Reconstruction, KDTree, Reconstruction
#####################################################################################
#                            Reconstruction Object                                  #
#####################################################################################
"""
    Reconstruction{V, T, D, τ}
`D`-dimensional reconstruction object with delay `τ`, created from a vector `s::V`.
`s` is the only field of `Reconstruction`.

The `n`th row of this object is formally the `D`-dimensional vector:
```math
(s(n), s(n+\\tau), s(n+2\\tau), \\dots, s(n+(D-1)\\tau))
```
which is created only on demand.

See [`reconstruct`](@ref) for usage.
"""
type Reconstruction{V<:AbstractVector, T<:Real, D, τ}
    s::V
end

# Type information:
@inline Base.eltype(R::Reconstruction{V, T, D, τ}) where {V, T, D, τ} = T
@inline Base.length(R::Reconstruction{V, T, D, τ}) where {V, T, D, τ} =
length(R.s) - (D-1)*τ
@inline dimension(R::Reconstruction{V, T, D, τ}) where {V, T, D, τ} = D
@inline delay(R::Reconstruction{V, T, D, τ}) where {V, T, D, τ} = τ
@inline Base.size(R::Reconstruction{V, T, D, τ}) where {V, T, D, τ} = (length(R), D)
# Indexing:
@inline function Base.getindex(
    R::Reconstruction{V, T, D, τ}, n::Int) where {V, T, D, τ}
    gen = (R.s[n+k*τ] for k in 0:(D-1))
    SVector{D,T}(gen...)
end
@inline function Base.getindex(
    R::Reconstruction{V, T, D, τ}, n::Int, j::Int) where {V, T, D, τ}
    return R.s[n + (j-1)*τ]
end
@inline function Base.getindex(
    R::Reconstruction{V, T, D, τ}, n::Int, j::Colon) where {V, T, D, τ}
    return R[n]
end
@inline function Base.getindex(
    R::Reconstruction{V, T, D, τ}, n::Colon, j::Int) where {V, T, D, τ}
    return R.s[1+(j-1)*τ:length(R)+(j-1)*τ]
end
# Conversions:
@inbounds function Dataset(R::Reconstruction{V, T, D, τ}) where {V, T, D, τ}
    data = SVector{D, T}[]
    for i in 1:length(R)
        push!(data, R[i])
    end
    return Dataset(data)
end
# Extend methods:
non0hist(ε, R::Reconstruction) = non0hist(ε, Dataset(R))
genentropy(α, ε, R::Reconstruction) = genentropy(α, ε, Dataset(R))
generalized_dim(α, R::Reconstruction) = generalized_dim(α, Dataset(R))

# Pretty-print
import Base.show
function Base.show(io::IO, R::Reconstruction{V, T, D, τ}) where {V, T, D, τ}
    print(io, "$(D)-dimensional Reconstruction with delay τ=$(τ)")
end

@require Juno begin
  function Juno.render(
    i::Juno.Inline,  R::Reconstruction{V, T, D, τ}) where {V, T, D, τ}
    t = Juno.render(i, Juno.defaultrepr(R))
    t[:head] = Juno.render(i,
    Text("$(D)-dimensional Reconstruction with delay τ=$(τ) from vector s"))
    return t
  end
end

"""
    reconstruct(s::AbstractVector, D::Int, τ::Int) -> R::Reconstruction
Create and return an efficient [`Reconstruction`](@ref) data structure that serves as
the
delay coordinates embedding [1, 2] reconstruction of the signal `s`.
The reconstuction has
dimension `D` and delay `τ` (measured in indices). This object can have same
invariant quantities (like e.g. lyapunov exponents) with the original system
that the timeseries were recorded from, for proper `D` and `τ` [1, 2].

`R` interfaces `s` and can be accessed similarly to a [`Dataset`](@ref):
```julia
R = reconstruct(s, 4, 1) # delay coords. reconstruction of dimension 4 and delay 1
R[3] # third point of reconstruction, ≡ (s[3], s[4], s[5], s[6])
R[1, 2] # Second element of first point of reconstruction, ≡ s[2]
```
*(this is only smart indexing, no Vectors or SVectors are created during the
construction of `R`)*

`R` can also be given to all functions that accept a `Dataset`, but
it is first converted to a `Dataset` at each call.
This means that you should first convert it yourself, using `Dataset(R)` if you
call functions like [`generalized_dim`](@ref) multiple times.

The functions `dimension(R)` and `delay(R)` return `D` and `τ` respectively. Notice
that `length(R) = length(s) - (D-1)*τ` (i.e. the amount of D-dimensional points
"contained" in `R`) but `size(R) = (length(R), D)`.

[1] : F. Takens, *Detecting Strange Attractors in Turbulence — Dynamical
Systems and Turbulence*, Lecture Notes in Mathematics **366**, Springer (1981)

[2] : T. Sauer *et al.*, J. Stat. Phys. **65**, pp 579 (1991)
"""
function reconstruct(s::AbstractVector, D::Int,  τ::Int)
    V = typeof(s)
    T = eltype(s)
    return Reconstruction{V, T, D, τ}(s)
end

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
function estimate_delay(s::AbstractVector)
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
# along with some clever indexing and dispatching that allows for minimal
# computations of distances.
# [1] : Skokos, C. H. *et al.*, *Chaos Detection and Predictability* - Chapter 1
# (section 1.3.2), Lecture Notes in Physics **915**, Springer (2016)



# Trees:
# Crazy hacking cheating tree that contains 1-dimensional points LOLZIES
function treedata(R::Reconstruction{V, T, D, τ}, di::Cityblock) where {V, T, D, τ}
    SX = reinterpret(SVector{1, T}, R.s, (length(R.s),))
end
# Tree based on the fact that the Euclidean distance will be used
function treedata(R::Reconstruction{V, T, D, τ}, di::Euclidean) where {V, T, D, τ}
    Dataset(R).data
end



# Neighborhoods:
"""
Supertype of methods for deciding the neighborhood of points for a given point.

Concrete subtypes:
* `FixedMassNeighborhood(K::Int)`  : The neighborhood of a point consists of the `K`
  nearest neighbors of the point.
* `FixedSizeNeighborhood(ϵ::Real)` : The neighborhood of a point consists of all
  neighbors that have distance < `ϵ` from the point.

Notice that these distances are always computed using the `Euclidean()` distance
in `D`-dimensional space, irrespectively of the `distance` used in the
function [`numericallyapunov`](@ref).

See also [`neighborhood`](@ref) or [`numericallyapunov`](@ref).
"""
abstract type AbstractNeighborhood end
struct FixedMassNeighborhood <: AbstractNeighborhood
    K::Int
end
FixedMassNeighborhood() = FixedMassNeighborhood(1)
struct FixedSizeNeighborhood <: AbstractNeighborhood
    ϵ::Float64
end
FixedSizeNeighborhood() = FixedSizeNeighborhood(0.001)

"""
    neighborhood(n, point, tree::KDTree, method::AbstractNeighborhood)
Return a vector of indices which are the neighborhood of `point`, whose index
in the original data is `n`. Both `point` and `n` must be provided because the
`tree` has indices in different sorting (thus making `tree.data[n]` incorrect).
The `method` can be a subtype of [`AbstractNeighborhood`](@ref) (see its documentation
string for more).

`neighborhood` can be used for *any* dataset. Just do:
```julia
R = some_dataset
tree = KDTree(R)
neigh = neighborhood(n, R[n], tree, method)
```
where `R` can be *either* a [`Dataset`](@ref) or a [`Reconstruction`](@ref).

Notice that the distances in the trees are always computed using the `Euclidean()`
distance in `D`-dimensional space, irrespectively of the `distance` used in the
`numericallyapunov` function.

`neighborhood` **simply interfaces** the functions
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
    idxs = inrange(tree, point, method.ϵ)
    deleteat!(idxs, findin(idxs, n)) # unfortunately this has to be done...
    return idxs
end
neighborhood(n, point, tree::KDTree) =
neighborhood(n, point, tree, FixedMassNeighborhood(1))

KDTree(D::Dataset) = KDTree(D.data, Euclidean())
KDTree(R::Reconstruction) = KDTree(Dataset(R).data, Euclidean())



"""
```julia
numericallyapunov(R::Reconstruction, ks;  refstates, distance, method)
```
Return `E = [E(k) for k ∈ ks]`, where `E(k)` is the average logarithmic distance for
nearby states that are evolved in time for `k` steps (`k` must be integer).
If the reconstruction
exhibits exponential divergence of nearby states, then it should clearly hold:
```math
E(k) \\approx \\lambda\\Delta t k + E(0)
```
for a **well defined region** in the `k` axis, where ``\\lambda`` is
the approximated
maximum Lyapunov exponent. `Δt` is the time between samples in the
original timeseries.
You can use [`linear_region`](@ref) with arguments `(ks, E)` to identify the slope
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

The following keywords tune the algorithm behavior:

* `refstates::AbstractVector{Int} = 1:(length(R) - ks[end])` : Vector of indices
  that notes which
  states of the reconstruction should be used as "reference states", which means
  that the algorithm is applied for all state indices contained in `refstates`.
* `distance::Metric = Cityblock()` : The distance function used in the
  logarithmic distance of nearby states. The allowed distances are `Cityblock()`
  and `Euclidean()` from the
  package `Distances.jl` (re-exported here). See below for more info on this
  choice which has fundamental impact on the algorithm.
* `method::AbstractNeighborhood = FixedMassNeighborhood(1)` : The method to
  be used when evaluating
  the neighborhood of each reference state. See
  `AbstractNeighborhood` for more info on this choice.

If the `Metric` is `Euclidean()` then calculate the Euclidean distance of the
full `D`-dimensional points (distance ``d_E`` in ref. [1]).
If however the `Metric` is `Cityblock()`, calculate
the absolute distance of **only the first elements** of the `m+k` and `n+k` points
of the
reconstruction `R`, which are the `m+k` and `n+k` elements of vector `R.s` (distance
``d_F`` in
ref. [1]). Notice that
the distances used are defined in the package `Distances.jl`, but are re-exported
in `DynamicalSystems.jl` for ease-of-use (the
distances are used for dispatch purposes *only*).

It is shown in [1] that the second type of distance
function makes little difference in the lyapunov estimation versus e.g. the
Euclidean distance, but it is **much cheaper to evaluate**, you only need to
perform one subtraction and one absolute value!

This function assumes that the Theiler window (see [1]) is the same as the delay time:
``w  = \\tau``.

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

function numericallyapunov(R::Reconstruction{V, T, D, τ},
                           ks::AbstractVector{Int},
                           ℜ::AbstractVector{Int},
                           distance::Metric,
                           method::AbstractNeighborhood) where {V, T, D, τ}

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
    td = Dataset(R).data #tree data
    tree = KDTree(td, Euclidean()) # this creates a copy of `td`
    skippedm = 0; skippedn = 0

    for n in ℜ
        # The ⋓(n) can be evaluated on the spot instead of being pre-calculated
        # for all reference states. (it would take too much memory)
        # Since U[n] doesn't depend on `k` one can then interchange the loops:
        # Instead of k being the outermost loop, it becomes the innermost loop!
        point = td[n]
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
    abs(R.s[m+k] - R.s[n+k])
end

@inline @inbounds function delay_distance(di::Euclidean,
    R::Reconstruction{V, T, D, τ}, m, n, k) where {V, T, D, τ}
    suma = zero(T)
    for j in 0:D-1
        suma += (R.s[m+k+j*τ] - R.s[n+k+j*τ])^2
    end
    sqrt(suma)
end
