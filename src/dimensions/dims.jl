export boxcounting_dim, capacity_dim, generalized_dim,
information_dim, correlation_dim, collision_dim, estimate_ε,
kaplanyorke_dim

magnitude(x::Real) = round(Int, log10(x))

"""
```julia
estimate_ε(dataset; m = 10, k::Int = 12, n::Int = 4)
```
Return `logspace(magnitude(x), magnitude(x) - n, k)` where `x` is
the `minimum( maximum(abs.(v)) - minimum(abs.(v)) for v in vectors )/m`.

In essense, get a `k`-element `logspace` with maximum being the `1/m` of the
relative order of magnitude of the vectors,
and the minimum being `n` orders of magnitude less.
"""
function estimate_ε(vectors::Vararg{AbstractVector{<:Real}};
  m = 10, k::Int = 12, n::Int = 4)
  # maximum ε is 1/m of maximum - minimum
  maxε = Inf
  for v in vectors
    vv = abs.(v)
    ma = maximum(vv)
    mi = minimum(vv)
    d = (ma - mi)/10
    if d < maxε
      maxε = d
    end
  end
  logspace(magnitude(maxε), magnitude(maxε)-n, k)
end
estimate_ε(dataset::AbstractMatrix{<:Real}) = estimate_ε(d2v(dataset)...)

"""
    generalized_dim(α, dataset)
Return the `α` order generalized dimension that corresponds to the given dataset.
This quantity corresponds to the
power law exponent of the scaling of the `genentropy` versus the box size `ε`.

**WARNING** - This call performs a lot of automated steps:

  1. A vector of box sizes is decided by calling `es = estimate_ε(dataset)`.
  2. For each `ε ∈ es` the appropriate entropy is
     calculated, through `d[i] = genentropy(α, es[i], dataset)`. Let `x = -log.(es)`.
  3. The curve d(x) is decomposed into linear regions, using `linear_regions(x, d)`.
  4. The biggest linear region is chosen, and a fit for the slope of that
     region is performed using the package `LsqFit` (see `linear_region`).
  5. This fitted slope is returned.

By doing these steps one by one yourself, you can adjust the keyword arguments
given to each of these function calls, refining the accuracy of the result.

The following aliases are provided:

  * α = 0 : `boxcounting_dim`, `capacity_dim`
  * α = 1 : `information_dim`
  * α = 2 : `correlation_dim`, `collision_dim`
"""
function generalized_dim(α, vectors::Vararg{AbstractVector{<:Real}})
  es = estimate_ε(vectors...)
  dd = zeros(es)
  for i in 1:length(es)
    dd[i] = genentropy(α, es[i], vectors...)
  end
  return linear_region(-log.(es), dd)[2]
end
generalized_dim(α, dataset::AbstractMatrix{<:Real}) =
generalized_dim(α, d2v(dataset)...)
# Aliases
"correlation_dim(args...) = generalized_dim(2, args...)"
correlation_dim(args...) = generalized_dim(2, args...)
collision_dim = correlation_dim

"capacity_dim(args...) = generalized_dim(0, args...)"
capacity_dim(args...) = generalized_dim(0, args...)
boxcounting_dim = capacity_dim

"information_dim(args...) = generalized_dim(1, args...)"
information_dim(args...) = generalized_dim(1, args...)

"""
```julia
kaplanyorke_dim(lyapunovs::AbstractVector)
```
Calculate the Kaplan-Yorke dimension [1]. This simply is the point where
`cumsum(lyapunovs)` becomes zero (interpolated). Returns the dimension of the system
if the sum of the exponents never becomes negative.

[1] :  J. Kaplan & J. Yorke,
*Chaotic behavior of multidimensional difference equations*,
Lecture Notes in Mathematics vol. **730**, Springer (1979)
"""
function kaplanyorke_dim(v::AbstractVector)
  issorted(v, rev = true) || throw(ArgumentError(
  "The lyapunov vector must be sorted from most positive to most negative"))

  s = cumsum(v); k = length(v)
  # Find k such that sum(λ_i for i in 1:k) is still possitive
  for i in eachindex(s)
    if s[i] < 0
      k = i-1
      break
    end
  end

  if k == 0
    return zero(v[1])
  elseif k < length(v)
    return k + s[k]/abs(v[k+1])
  else
    return typeof(v[1])(length(v))
  end
end
