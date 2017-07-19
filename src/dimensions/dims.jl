export boxcounting_dim, capacity_dim, generalized_dim

estimate_ε(dataset) = estimate_ε(d2v(dataset)...)
magnitude(x::Real) = round(Int, log10(x))

"""
    estimate_ε(vectors..., m=10)
Return `logspace(magnitude(x), magnitude(x) - 4, 9)` where `x` is
the `minimum( maximum(abs.(v)) - minimum(abs.(v)) for v in vectors )/m`.

In essense, get a `logspace` with maximum being the 1/m of the relative order of
magnitude of the vectors, and the minimum being 4 orders of magnitude less.
"""
function estimate_ε(vectors..., m=10)
  # maximum ε is 1/100 of maximum - minimum
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
  logspace(magnitude(maxε), magnitude(maxε)-4, 9)
end




# boxcounting_dim(dataset)
# const capacity_dim = boxcounting_dim
