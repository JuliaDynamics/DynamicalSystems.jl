export non0hist, renyi
using StaticArrays

"""
```julia
d2v(dataset) -> vectors
```
Convert a dataset to a bundle of vectors, using views.
"""
function d2v(dataset)
  D = size(dataset)[2] #dimension of the system (dynamic variables)
  vectors = typeof(view(dataset, :, 1))[] #initialize array that will store vectors
  for k in 1:D
    push!(vectors, view(dataset, :, k))
  end
  return vectors
end

@inbounds function perform_non0hist(vectors, ranges, ε)
  D = length(vectors)
  L = length(vectors[1])
  d = Dict{NTuple{D, Int}, Int}()
  minima = [minimum(v) for v in vectors]
  for j in 1:length(vectors[1])
    ind = NTuple{D, Int}(Int(floor((vectors[k][j] - minima[k])/ε)) + 1 for k in 1:D)
    haskey(d, ind) || (d[ind] = 0)
    d[ind] += 1
  end
  return [val/L for val in values(d)]
end

function perform_non0hist_statsbase(vectors, ranges, ε)
  pks = fit(Histogram, (vectors...), (ranges...), closed=:left).weights/L
  return T[x for x in pks if x != 0]
end

"""
```julia
non0hist(ε, dataset)
non0hist(ε, vectors...)
```
Partition a data-set into tabulated intervals (boxes) of
size `ε` and return the *sum-normalized* histogram in an **unordered 1D form**,
discarding all non-zero elements.

Use the `fit(Histogram, ...)` from package `StatsBase` if you
wish to keep information about the edges of the binning as well
as the zero elements.
"""
function non0hist end
@inbounds function non0hist(ε::Real, vectors::Vararg{AbstractVector{T}}) where {T<:Real}
  # Initialize:
  D = length(vectors); L = length(vectors[1])
  for i in 1:D
    length(vectors[i]) == L ||
    throw(ArgumentError("All vectors given to `non0hist` must be of equal length!"))
  end

  if T == BigFloat || T == BigInt
    ranges = StepRangeLen{T, T, T}[]
  else
    ranges = StepRangeLen{T, Base.TwicePrecision{T},Base.TwicePrecision{T}}[]
  end
  # Create ranges:
  for v in vectors
    push!(ranges, minimum(v):ε:maximum(v)+ε)  #be sure to have that +ε at the end!
  end
  # Perform Histogram:
  perform_non0hist(vectors, ranges, ε)
end
@inbounds non0hist(ε::Real, dataset::AbstractMatrix{<:Real}) =
non0hist(ε, d2v(dataset)...)


"""
```julia
renyi(α, ε, dataset)
renyi(α, ε, vectors...)
```
Compute the `α` order generalized Rényi entropy [1] of a dataset,
by first partitioning it into boxes of length `ε` (log base-e is used).

```julia
renyi(α, p::AbstractArray)
```
Compute the entropy of an Array `p` directly, assuming that `p` is
sum-normalized (log base-e is used).

The Rényi entropy `R_α(p) = 1/(1-α) * sum_i(p_i^α)` generalizes other known entropies,
like e.g. the Shannon entropy
(α = 1, see *the* Shannon paper [2]), the Hartley entropy (α = 0, also known as
maximum entropy), or the Collision entropy (α = 2, also known as correlation entropy).

[1] : A. Rényi, Proceedings of the fourth Berkeley Symposium on Mathematics,
Statistics and Probability, pp 547 (1960)

[2] : C. E. Shannon, Bell System Technical Journal **27**, pp 379 (1948)
"""
function renyi(α::Real, ε::Real, vectors::Vararg{AbstractVector{T}}) where {T<:Real}
  p = non0hist(ε, vectors...)
  renyi(α, p)
end
renyi(α::Real, ε::Real, dataset::AbstractMatrix{<:Real}) =
renyi(α, ε, d2v(dataset)...)

function renyi{T<:Real}(α::Real, p::AbstractArray{T})
  α < 0 && throw(ArgumentError("Order of Rényi entropy must be ≥ 0."))

  if α ≈ 0
    return log(length(p)) #Hartley entropy, max-entropy
  elseif α ≈ 1
    return -sum( x*log(x) for x in p ) #Shannon entropy, information to locate with ε
  elseif (isinf(α))
    return -log(maximum(p)) #Min entropy
  else
    return (1/(1-α))*log( sum(x^α for x in p) )
  end
end
