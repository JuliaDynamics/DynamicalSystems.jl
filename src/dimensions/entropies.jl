export non0hist, genentropy, renyi, d2v
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
  D = length(vectors); L = length(vectors[1])
  # `d` is a dictionary that contains all the histogram information
  # the keys are the bin edges indeces and the values are the amount of
  # datapoints in each bin
  d = Dict{NTuple{D, Int}, Int}()
  minima = [minimum(v) for v in vectors]
  for j in 1:L
    ind = NTuple{D, Int}( #index of datapoint in the ranges space
    Int(floor((vectors[k][j] - minima[k])/ε)) + 1 for k in 1:D)

    haskey(d, ind) || (d[ind] = 0) #check if you need to create key (= bin)
    d[ind] += 1 #add 1 to an existing bin
  end
  return [val/L for val in values(d)]
end

function perform_non0hist_statsbase(vectors, ranges, ε)
  # StatsBase version
  pks = fit(Histogram, (vectors...), (ranges...), closed=:left).weights/L
  return T[x for x in pks if x != 0]
end

"""
```julia
non0hist(ε, dataset), non0hist(ε, vectors...)
```
Partition a data-set into tabulated intervals (boxes) of
size `ε` and return the *sum-normalized* histogram in an **unordered 1D form**,
**discarding all zero** elements.

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
  # Create ranges (bin edges):
  for v in vectors
    push!(ranges, minimum(v):ε:maximum(v)+ε)  #be sure to have that +ε at the end!
  end
  # Perform Histogram:
  perform_non0hist(vectors, ranges, ε)
end
non0hist(ε::Real, dataset::AbstractMatrix{<:Real}) = non0hist(ε, d2v(dataset)...)


"""
```julia
genentropy(α, ε, dataset), genentropy(α, ε, vectors...)
```
Compute the `α` order generalized (Rényi) entropy [1] of a dataset,
by first partitioning it into boxes of length `ε` (log base-e is used).

Other aliases: `renyi`.

```julia
genentropy(α, p::AbstractArray)
```
Compute the entropy of an Array `p` directly, assuming that `p` is
sum-normalized (log base-e is used).

The Rényi entropy:
```math
R_\\alpha(P) = \\frac{1}{1-\\alpha}\sum_{p \\in P}(p^\\alpha)
```
generalizes other known entropies,
like e.g. the information entropy
(α = 1, see *the* Shannon paper [2]), the maximum entropy (α = 0, also known as
Hartley entropy), or the correlation entropy (α = 2, also known as collision entropy).

[1] : A. Rényi, *Proceedings of the fourth Berkeley Symposium on Mathematics,
Statistics and Probability*, pp 547 (1960)

[2] : C. E. Shannon, Bell System Technical Journal **27**, pp 379 (1948)
"""
function genentropy(α::Real, ε::Real, vectors::Vararg{AbstractVector{T}}) where {T<:Real}
  p = non0hist(ε, vectors...)
  genentropy(α, p)
end
genentropy(α::Real, ε::Real, dataset::AbstractMatrix{<:Real}) =
genentropy(α, ε, d2v(dataset)...)

function genentropy{T<:Real}(α::Real, p::AbstractArray{T})
  α < 0 && throw(ArgumentError("Order of Rényi entropy must be ≥ 0."))

  if α ≈ 0
    return log(length(p)) #Hartley entropy, max-entropy
  elseif α ≈ 1
    return -sum( x*log(x) for x in p ) #Shannon entropy, information to locate with ε
  elseif (isinf(α))
    return -log(maximum(p)) #Min entropy
  else
    return (1/(1-α))*log( sum(x^α for x in p) ) #genentropy entropy
  end
end

renyi = genentropy
