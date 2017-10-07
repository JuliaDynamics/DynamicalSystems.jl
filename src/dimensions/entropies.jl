export non0hist, genentropy, renyi, shannon, hartley

@inbounds function perform_non0hist{D, T<:Real, V}(data::Dataset{D,T, V}, ranges, ε)
    L = length(data)
    # `d` is a dictionary that contains all the histogram information
    # the keys are the bin edges indices and the values are the amount of
    # datapoints in each bin
    d = Dict{SVector{D, Int}, Int}()
    mini = minima(data)
    for point in data
        # index of datapoint in the ranges space:
        # It is not necessary to convert floor to Int (dunno why)
        ind::SVector{D, Int} = SVector{3}(
        (Int(floor( (point[i] - mini[i])/ε) ) for i in 1:D)...)
        # Add 1 to the bin that contains the datapoint:
        haskey(d, ind) || (d[ind] = 0) #check if you need to create key (= bin)
        d[ind] += 1
    end
    return [val/L for val in values(d)]
end


"""
```julia
non0hist(ε, dataset)
```
Partition a dataset into tabulated intervals (boxes) of
size `ε` and return the *sum-normalized* histogram in an **unordered 1D form**,
**discarding all zero** elements. This method is extremely effecient in both memory
and speed, because it uses a dictionary to collect the information of bins with
elements, while it completely disregards empty bins.

Use e.g. `fit(Histogram, ...)` from `StatsBase` if you
wish to keep information about the edges of the binning as well
as the zero elements.
"""
function non0hist end
@inbounds function non0hist{D, T<:Real, V}(ε::Real, data::Dataset{D, T, V})
    # Initialize:
    mini = minima(data); maxi = maxima(data)
    ranges = [mini[i]:ε:maxi[i]+ε for i in 1:D]
    # Perform Histogram:
    perform_non0hist(data, ranges, ε)
end
non0hist(ε::Real, matrix::AbstractMatrix) =
non0hist(ε, convert(Dataset, matrix))


"""
```julia
genentropy(α, ε, dataset)
```
Compute the `α` order generalized (Rényi) entropy [1] of a dataset,
by first partitioning it into boxes of length `ε` using [`non0hist`](@ref).
```julia
genentropy(α, p::AbstractArray)
```
Compute the entropy of an Array `p` directly, assuming that `p` is
sum-normalized. *log base-e is used in both cases, i.e. units of "nat".*

The Rényi entropy
```math
R_\\alpha(p) = \\frac{1}{1-\\alpha}\\sum_i p_i^\\alpha
```
generalizes other known entropies,
like e.g. the information entropy
(``\\alpha = 1``, see *the* Shannon paper [2]), the maximum entropy (``\\alpha=0``,
also known as Hartley entropy), or the correlation entropy
(``\\alpha = 2``, also known as collision entropy).

The following aliases are provided:

* `renyi = genentropy`
* `shannon(args...) = genentropy(1, args...)`
* `hartley(args...) = genentropy(0, args...)`

[1] : A. Rényi, *Proceedings of the fourth Berkeley Symposium on Mathematics,
Statistics and Probability*, pp 547 (1960)

[2] : C. E. Shannon, Bell Systems Technical Journal **27**, pp 379 (1948)
"""
function genentropy(α::Real, ε::Real, data::Dataset)
    ε < 0 && throw(ArgumentError("Box-size for entropy calculation must be ≥ 0."))
    p = non0hist(ε, data)
    return genentropy(α, p)
end
genentropy(α::Real, ε::Real, matrix::AbstractMatrix) =
genentropy(α, ε, convert(Dataset, matrix))

function genentropy{T<:Real}(α::Real, p::AbstractArray{T})
  α < 0 && throw(ArgumentError("Order of Rényi entropy must be ≥ 0."))

  if α ≈ 0
    return log(length(p)) #Hartley entropy, max-entropy
  elseif α ≈ 1
    return -sum( x*log(x) for x in p ) #Shannon entropy, information to locate with ε
  elseif isinf(α)
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
