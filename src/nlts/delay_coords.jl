using RecursiveArrayTools, Requires, StatsBase
export reconstruct
#####################################################################################
#                            Reconstruction Object                                  #
#####################################################################################
type Reconstruction{T, N, A} <: AbstractVectorOfArray{T, N}
  u::A # A <: AbstractVector{<: AbstractArray{T, N - 1}}
  τ::Int
  d::Int
end

# Assume that the first element is representative for all other elements
Reconstruction{T, N}(vec::AbstractVector{T}, dims::NTuple{N}, τ, d) =
Reconstruction{eltype(T), N, typeof(vec)}(vec, τ, d)

Reconstruction(vec::AbstractVector, τ, d) =
Reconstruction(vec, (size(vec[1])..., length(vec)), τ, d)

d2v(R::Reconstruction) = R.u

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
function estimate_τ(s::AbstractVector)
  c = autocor(x, 0:length(x)÷10)
  # First approach: find zero
  # Second approach: if all positive, perform exponential fit

end

function estimate_d(s::AbstractVector)
  # Estimate number of “false nearest neighbors” due to
  # projection into a too low dimension reconstruction space
end
