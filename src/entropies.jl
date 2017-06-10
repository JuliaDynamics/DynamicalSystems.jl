"""
    non0hist(ε, vectors...)
Partition a data-set contained in the `vectors` into tabulated intervals (boxes) of
size `ε` and return the *normalized* histogram in a 1D form, discarding all
non-zero elements. Use the `fit(Histogram, ...)` from `StatsBase.jl` if you wish to
keep information about the edges of the binning as well as the zero elements.

It is assumed that from the given vectors, each one describes one
dimension of the data-set, e.g. x, y, etc.
"""
@inbounds function non0hist(ε, vectors::Vararg{AbstractVector{T}}) where {T<:Real}
  # Initialize:
  D = length(vectors); L = length(vectors[1])
  for i in 1:D
    @assert length(v[i]) == L "All vectors given to `non0hist` be of equal length!"
  end

  if T == BigFloat || T == BigInt
    ranges = StepRangeLen{T, T, T}[]
  else
    ranges = StepRangeLen{T, Base.TwicePrecision{T},Base.TwicePrecision{T}}[]
  end
  # Create ranges:
  @simd for v in vectors
    push!(ranges, minimum(v):ε:maximum(v)+ε)  #be sure to have that +ε at the end!
  end
  # Perform Histogram:
  pks = fit(Histogram, vectors, (ranges...), closed=:left).weights/L
  return Float64[x for x in pks if x != 0]
end



"""
```julia
renyi(α, ε, vectors...)
```
Compute the `α` order generalized Rényi entropy of a data-set composed of equally
sized `vectors`, by first partitioning them into boxes of of size `ε`. (log base-e is used)

```julia
renyi(α, p::AbstractArray)
```
Compute the `α` order generalized Rényi entropy of an Array `p` directly,
assuming that `p` is normalized.
"""
function renyi(α, ε, vectors::Vararg{AbstractVector{T}}) where {T<:Real}
  p = non0hist(ε, vectors...)
  renyi(α, p)
end
function renyi{T<:Real}(α::Real, p::AbstractArray{T})
  α < 0 && throw(ArgumentError("Order of Rényi entropy must be ≥0."))
  s = sum(p)
  if !(s ≈ 1.0)
    p ./= s
  end

  if q ≈ 0
    return log(length(p)) #Hartley entropy, max-entropy
  elseif q ≈ 1
    return -sum( x*log(x) for x in p ) #Shannon entropy, information to locate with ε
  elseif (isinf(α))
    return -log(maximum(p)) #Min entropy
  else
    return (1/(1-q))*log( sum(x^α for x in p) )
  end
end
