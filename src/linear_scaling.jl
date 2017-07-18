using LsqFit: curve_fit
export linear_region, linear_regions
#######################################################################################
# Functions and methods to deduce linear scaling regions
#######################################################################################
"""
    isevenly(a::AbstractVector)
Check if `a` is evenly spaced.
"""
function isevenly(a::AbstractVector)
  test = a[2] - a[1]
  for i in 2:length(a)-1
    if !(a[i+1] - a[i] ≈ test)
      throw(ArgumentError("x-axis is not evenly spaced!"))
    end
  end
  true
end

"""
```julia
linear_region(x, y; dxi::Int = 1, tol = 0.1) -> ([ind1, ind2], slope)
```
Call `linear_regions`, identify the largest linear region and approximate
the slope of this region using least squares fit. Return the indeces where
the region starts and stops (`x[ind1:ind2]`) as well as the approximated `slope`.
"""
function linear_region(x::AbstractVector, y::AbstractVector,
  dxi::Int = 1, tol::Real = 0.1)

  # Find biggest linear region:
  reg_ind = max_linear_region(linear_regions(x,y; dxi=dxi, tol=tol)...)
  # Prepare least squares fit:
  xfit = view(x, reg_ind[1]:reg_ind[2])
  yfit = view(y, reg_ind[1]:reg_ind[2])
  p0 = [1.0, 1.0]
  model(x, p) = p[1].*x .+ p[2]
  # Find fit of tangent:
  fit = curve_fit(model, xfit, yfit, p0)
  approx_tang = fit.param[1]
  return reg_ind, approx_tang
end

"""
    max_linear_region(lrs::Vector{Int}, tangents::Vector{Float64})
Find the biggest linear region and return it.
"""
function max_linear_region(lrs::Vector{Int}, tangents::Vector{Float64})
  dis = 0
  tagind = 0
  for i in 1:length(lrs)-1
    if lrs[i+1] - lrs[i] > dis
      dis = lrs[i+1] - lrs[i]
      tagind = i
    end
  end
  return [lrs[tagind], lrs[tagind+1]]
end

"""
```julia
linear_regions(x, y; dxi::Int = 1, tol = 0.1) -> (lrs, tangents)
```
Identify regions where the curve `y(x)` is linear, by scanning the
`x`-axis every `Dx` indeces (e.g. at `x[1]:x[5], x[5]:x[10], x[10]:x[15], ...`
if `Dx=5`).

If the slope (calculated using LsqFit) of a region of width `Dx` is
approximate to the previous region,
within tolerance `tol`,
then these two regions belong to the same linear region.

Return the indeces of `x` that correspond to linear regions, `lrs`,
and the approximated `tangents` at each region. `lrs` is a vector of `Int`.

Example:
```julia
lrs, tangents = linear_regions(xdata, ydata)
x[lrs[1]:lrs[2]] #this is the first linear region
tangents[1] #this is the tangent approximating the first linear region
```
"""
function linear_regions(x::AbstractVector, y::AbstractVector;
  dxi::Int = 1, tol::Real = 0.1)

  maxit = length(x) ÷ dxi

  tangents = Float64[slope(view(x, 1:max(dxi,2)), view(y, 1:max(dxi, 2)))]

  prevtang = tangents[1]
  lrs = Int[1] #start of first linear region is always 1
  lastk = 1


  for k in 1:maxit-1
    #tang = (y[(k+1)*dxi] - y[(k)*dxi]) / (x[(k+1)*dxi] - x[(k)*dxi])
    tang = slope(view(x, k*dxi:(k+1)*dxi), view(y, k*dxi:(k+1)*dxi))
    if isapprox(tang, prevtang, rtol=tol)
      # Tanget is similar with initial previous one (based on tolerance)
      continue
    else
      # Tangent is not similar.
      # Push new tangent for a new linear region
      push!(tangents, tang)
      # Set the START of a new linear region
      # which is also the END of the previous linear region
      push!(lrs, k*dxi)
      lastk = k
    end

    # Set new previous tangent (only if it was not the same as current)
    prevtang = tang
  end
  push!(lrs, length(x))
  return lrs, tangents
end

function slope(xfit, yfit)
  p0 = [1.0, 1.0]
  model(x, p) = p[1].*x .+ p[2]
  # Find fit of tangent:
  curve_fit(model, xfit, yfit, p0).param[1]
end
