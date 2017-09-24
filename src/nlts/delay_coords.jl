using RecursiveArrayTools, Requires, StatsBase
using LsqFit: curve_fit

export reconstruct, Reconstruction
#####################################################################################
#                            Reconstruction Object                                  #
#####################################################################################
"""
    Reconstruction{T, N, A} <: AbstractVectorOfArray{T, N}
See `reconstruct`.
"""
type Reconstruction{T, N, A} <: AbstractVectorOfArray{T, N}
  u::A # A <: AbstractVector{<: AbstractArray{T, N - 1}}
  τ::Int
  d::Int
end
# struct Reconstruction{
#     D, T<:Number, V::SVector{D,T}} <: AbstractDataset{D}
#     data::Vector{V}
#     τ::Int
# end

# Assume that the first element is representative for all other elements
Reconstruction{T, N}(vec::AbstractVector{T}, dims::NTuple{N}, τ, d) =
Reconstruction{eltype(T), N, typeof(vec)}(vec, τ, d)

Reconstruction(vec::AbstractVector, τ, d) =
Reconstruction(vec, (size(vec[1])..., length(vec)), τ, d)

# Pretty-print
import Base.show
function Base.show{T, N, A}(io::IO, R::Reconstruction{T, N, A})
  print(io, "$(R.d)-dimensional Reconstruction{$T} with delay τ=$(R.τ)")
end

@require Juno begin
  function Juno.render{T, N, A}(i::Juno.Inline, R::Reconstruction{T, N, A})
    t = Juno.render(i, Juno.defaultrepr(R))
    t[:head] = Juno.render(i,
    Text("$(R.d)-dimensional Reconstruction{$T} with delay τ=$(R.τ)"))
    pop!(t[:children]); pop!(t[:children])

    t[:children][1][:child][:head][:contents][1][:contents][1] = "SubArray"
    t[:children][1][:label][:contents][1] = "u (columns) → "
    return t
  end
end

"""
    reconstruct(s::AbstractVector, τ::Int, d::Int) -> R
Create and return an efficient `Reconstruction` data structure that serves as the
delay-coordinates reconstruction of the signal `s`. The reconstuction has
dimension `d` and delay `τ` (measured in indeces). This object can have same
invariant quantities (like e.g. lyapunov exponents) with the original system
that the timeseries were recorded from [1, 2].

The returned `R` is a `VectorOfArrays` from `RecursiveArrayTools.jl` and
stores all the information about the embedding without allocating new arrays,
using `view`. It can however be used as a normal matrix:
```julia
R[:, 2] # get the second column the reconstructed matrix
R[5, 1] # get the 5th element of the first column of the matrix
```

[1] : F. Takens, *Detecting Strange Attractors in Turbulence— Dynamical
Systems and Turbulence*, Lecture Notes in Mathematics **366** (1981)

[2] : T. Sauer *et al.*, J. Stat. Phys. **65**, pp 579 (1991)
"""
function reconstruct(s::AbstractVector, τ::Int, d::Int)
  N = length(s)
  u = typeof(view(s, 1:2))[]

  for i in 0:d-1
    push!(u, view(s, (1 + (i*τ)):(N - (d-i-1)*τ) ) )
  end
  return Reconstruction(u, τ, d)
end

#####################################################################################
#                      Estimate Reconstruction Parameters                           #
#####################################################################################
"""
    localextrema(y) -> max_ind, min_ind
Find the local extrema of given array `y`, by scanning point-by-point. Return the
indeces of the maxima (`max_ind`) and the indeces of the minima (`min_ind`).
"""
function localextrema(y)
    l = length(y)
    i = 1
    maxargs = Int[]
    minargs = Int[]
    if y[1] > y[2]
        push!(maxargs, 1)
    elseif y[1] < y[2]
        push!(minargs, 1)
    end

    for i in 2:l-1
        left = i-1
        right = i+1
        if  y[left] < y[i] > y[right]
            push!(maxargs, i)
        elseif y[left] > y[i] < y[right]
            push!(minargs, i)
        end
    end

    if y[l] > y[l-1]
        push!(maxargs, l)
    elseif y[l] < y[l-1]
        push!(minargs, l)
    end
    return maxargs, minargs
end


function exponential_decay_extrema(c::AbstractVector)
    ac = abs.(c)
    ma, mi = localextrema(ac)
    # ma start from 1 but correlation is expected to start from x=0
    ydat = ac[ma]; xdat = ma .- 1
    # Do curve fit from LsqFit
    model(x, p) = @. exp(-x/p[1])
    decay = curve_fit(model, xdat, ydat, [1.0]).param[1]
    return decay
end

function exponential_decay(c::AbstractVector)
    # Do curve fit from LsqFit
    model(x, p) = @. exp(-x/p[1])
    decay = curve_fit(model, 0:length(c)-1, abs.(c), [1.0]).param[1]
    return decay
end

"""
    estimate_delay(s) -> τ
Estimate an optimal delay by performing an exponential fit to
the `abs.(c)` with `c` the auto-correlation function of `s`.
Return the exponential decay time `τ` rounded to an integer.
"""
function estimate_delay(s::AbstractVector)
    c = autocor(x, 0:length(x)÷10)
    i = 1
    # Find 0 crossing:
    while c[i] > 0
        i+= 1
        i == length(c) && break
    end
    # Find exponential fit:
    τ = exponential_decay(c)
    # Is there a method to deduce which one of the 2 is the better approach?
    return round(Int, τ)
end

function estimate_d(s::AbstractVector)
  # Estimate number of “false nearest neighbors” due to
  # projection into a too low dimension reconstruction space
end

#####################################################################################
#                    Numerical Lyapunov (from Reconstruction)                       #
#####################################################################################
d_F(x1, x2) = abs(x1[1] - x2[1])
function lyapunov(R::Reconstruction)
    τ = R.τ; D = R.d
end
# ℜ = set of indeces that have the points that one finds neighbors for .
# n belongs in ℜ
# ⩅ = indeces of nearest neighbors of point n (so it is a vector of vectors)
# or a dictionary.

# R[n] is the "reference state" (vector not number)
# but R[n][1] is the first number of the point. This is used in the
# distance calulation. Notice that R[n][i] coincides with s[i], with
# s being the generator of the Reconstruction.

# Calculate E(k)
E = zeros(100)
for k in 1:100 #how to decide the range of k???
    dis = 0.0
    for n in ℜ
        for m in ⋓[n]
            dis += ln(d_F(R[m][k], R[n][k]))
        end
        dis /= length(⋓[n])
    end
    E[k] = dis/length(ℜ)
end
#plot E[k] versus k and boom, you got lyapunov.

# Method to estimate set ℜ :
ℜ = nothing

# Method to estimate neighborhoods for each n : keyword: neighbors="fixed_size"
# or neighbors="fixed_mass". The latter uses keyword neighbor_mass::Int = 1
# while the first uses keyword neighbor_size = 0.01
⋓ = Dict{Int, Vector{Int}}()
for n in ℜ
    ⋓[n] = K nearest neighbors using KDTree BUT that have |m - n| > τ
    # OR
    ⋓[n] = all points `m` that have distance from `n` < d and that have |m - n| > τ
end
