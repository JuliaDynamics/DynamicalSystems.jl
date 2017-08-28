export non0hist, genentropy, renyi, shannon, hartley

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
  return [x for x in pks if x != 0]
end

"""
```julia
non0hist(ε, dataset)
```
Partition a data-set into tabulated intervals (boxes) of
size `ε` and return the *sum-normalized* histogram in an **unordered 1D form**,
**discarding all zero** elements. This method is extremely effecient in both memory
and speed, because it uses a dictionary to collect the information of bins with
elements, while it completely disregards empty bins.

Use e.g. `fit(Histogram, ...)` from `StatsBase` if you
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
  ranges = typeof(minimum(vectors[1]):ε:maximum(vectors[1])+ε)[]
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
genentropy(α, ε, dataset)
```
Compute the `α` order generalized (Rényi) entropy [1] of a dataset,
by first partitioning it into boxes of length `ε`.
```julia
genentropy(α, p::AbstractArray)
```
Compute the entropy of an Array `p` directly, assuming that `p` is
sum-normalized.

log base-e is used in both cases, i.e. units of "nat".

The Rényi entropy `R_α(p) = (1/1-α)*sum(pi^α for pi ∈ p)`
generalizes other known entropies,
like e.g. the information entropy
(α = 1, see *the* Shannon paper [2]), the maximum entropy (α = 0, also known as
Hartley entropy), or the correlation entropy (α = 2, also known as collision entropy).

The following aliases are provided:

* `renyi = genentropy`
* `shannon(args...) = genentropy(1, args...)`
* `hartley(args...) = genentropy(0, args...)`

[1] : A. Rényi, *Proceedings of the fourth Berkeley Symposium on Mathematics,
Statistics and Probability*, pp 547 (1960)

[2] : C. E. Shannon, Bell Systems Technical Journal **27**, pp 379 (1948)
"""
function genentropy(α::Real, ε::Real, vectors::Vararg{AbstractVector{T}}) where {T<:Real}
  ε < 0 && throw(ArgumentError("Box-size for entropy calculation must be ≥ 0."))
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
# Aliases:
renyi = genentropy

"shannon(args...) = genentropy(1, args...)"
shannon(args...) = genentropy(1, args...)

"hartley(args...) = genentropy(0, args...)"
hartley(args...) = genentropy(0, args...)
