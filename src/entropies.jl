export non0hist, renyi
using StatsBase

"""
```julia
non0hist(ε, dataset)
non0hist(ε, vectors...)
```
Partition a data-set into tabulated intervals (boxes) of
size `ε` and return the *sum-normalized* histogram in a 1D form, discarding all
non-zero elements. Use the `fit(Histogram, ...)` from package `StatsBase` if you
wish to keep information about the edges of the binning as well as the zero elements.
"""
function non0hist end
@inbounds function non0hist(ε, vectors::Vararg{AbstractVector{T}}) where {T<:Real}
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
  # Write my OWN Version, which would be many times over faster,
  # since I don't need the 0 elements.
  pks = fit(Histogram, vectors, (ranges...), closed=:left).weights/L
  return T[x for x in pks if x != 0]
end



"""
```julia
renyi(α, ε, dataset)
renyi(α, ε, vectors...)
```
Compute the `α` order generalized Rényi entropy [1] of a dataset,
by first partitioning them into boxes of of size `ε` (log base-e is used).

```julia
renyi(α, p::AbstractArray)
```
Compute the `α` order generalized Rényi entropy of an Array `p` directly,
assuming that `p` is sum-normalized.

The Rényi entropy `R_α(p) = 1/(1-α) * sum_i(p_i^α)` generalizes other known entropies,
like e.g. the Shannon entropy
(α = 1, see *the* Shannon paper [2]), the Hartley entropy (α = 0, also known as
maximum entropy), or the Collision (or correlation) entropy (α = 2).

[1] : A. Rényi, Proceedings of the fourth Berkeley Symposium on Mathematics,
Statistics and Probability, pp 547 (1960)

[2] : C. E. Shannon, Bell System Technical Journal **27**, pp 379 (1948)
"""
function renyi(α, ε, vectors::Vararg{AbstractVector{T}}) where {T<:Real}
  p = non0hist(ε, vectors...)
  renyi(α, p)
end
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



# # Linear fit on dimension:
# path = "C:\\Users\\datseris\\PowerFolders\\MyScience\\Lectures\\Introduction to Physics of Complex Systems\\Exercises 2016\\10 - Generilized Dimensions And Synchronization\\data1.txt"
# f = readdlm(path)
# f1 = view(f, :, 1); f2 = view(f, :, 2)
#
# es = logspace(-1, -3, 6)
# D = zeros(es)
# for (i, ε) in enumerate(es)
#   println("ε = $ε")
#   D[i] = renyi(1, ε, f1, f2)
# end
#
# using LsqFit
# model(x, p) = p[1].*x .+ p[2]
#
# cf = curve_fit(model, -ln.(es), D, rand(2))
# plot(lnes, -cf.param[1] .* ln.(es) .+ cf.param[2])
