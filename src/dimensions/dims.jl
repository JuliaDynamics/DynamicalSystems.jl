export boxcounting_dim, capacity_dim, generalized_dim,
information_dim, correlation_dim

magnitude(x::Real) = round(Int, log10(x))

"""
```julia
estimate_ε(dataset; m = 10, k::Int = 9)
```
Return `logspace(magnitude(x), magnitude(x) - 4, k)` where `x` is
the `minimum( maximum(abs.(v)) - minimum(abs.(v)) for v in vectors )/m`.

In essense, get a `k`-element `logspace` with maximum being the 1/m of the
relative order of magnitude of the vectors,
and the minimum being 4 orders of magnitude less.
"""
function estimate_ε(vectors::Vararg{AbstractVector{<:Real}}; m = 10, k::Int = 9)
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
  logspace(magnitude(maxε), magnitude(maxε)-4, k)
end
estimate_ε(dataset::AbstractMatrix{<:Real}) = estimate_ε(d2v(dataset)...)

"""
    generalized_dim(α, dataset)
Return the capacity (a.k.a. box counting) dimension
that corresponds to the given dataset.

**WARNING** : This call performs a lot of automated steps:
1. A vector of box sizes is decided by calling `es = estimate_e(dataset)`.
2. For each `ε ∈ es` the appropriate entropy is
  calculated, through `d[i] = genentropy(α, es[i], dataset)`. Let `x = -log.(es)`.
3. The curve d(x) is decomposed into linear regions, using `linear_regions(x, d)`.
4. The biggest linear region is chosen, and a fit for the slope of that
  region is performed
  using the package `LsqFit` (see `linear_region`).
5. This fitted slope is returned.

By doing these steps one by one yourself, you can adjust the keyword arguments
given to each of these function calls, refining the accuracy of the result.
"""
function generalized_dim(α, vectors::Vararg{AbstractVector{<:Real}})
  es = estimate_ε(vectors...)
  dd = zeros(es)
  for i in 1:length(es)
    dd[i] = genentropy(0, es[i], vectors...)
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
