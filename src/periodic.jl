import Base: *

# Method to find unstable periodic orbits of any map:
_plot_dataset(d::Dataset; kwargs...) = _plot_dataset(d.data; kwargs...)
function _plot_dataset(S::Vector{SVector{2,T}}; kwargs...) where {T}
    mat = Matrix{T}(length(S), D)
    for i in 1:length(S)
      mat[i,:] .= S[i]
    end
    PyPlot.scatter(view(mat, :, 1), view(mat, :, 1), kwargs...)
end
function _phasespace(ds; xs = linspace(0, 2π, 20), ys = linspace(0, 2π, 20),
    maxiters = 1000)
    f = ds.eom
    dataset = timeseries(ds, maxiters)
    for x in xs
        for y in ys
            ds.state = SVector{2}(x, y)
            append!(dataset, timeseries(ds, maxiters))
        end
    end
    m = Matrix(dataset)
    PyPlot.scatter(view(m, :, 1), view(m, :, 2), s= 1, color = "black")
    PyPlot.xlim(xs[1], xs[end])
    PyPlot.ylim(ys[1], ys[end])
end




"""
    LambdaMatrix{T<:Real}
Type representing the matrix ``\\mathbf{\\Lambda}_k`` used to create a new
dynamical system with some unstable fixed points turned to stable. See eq. (1)
of [1]. `D` is the dimension of the matrix (and thus of the dynamical system).

# Fields
* `λ::T`
* `ind::Vector{Int}`
* `sings::Vector{Int8}`

`λ` is the multiplier of the ``C_k`` matrix, with `0<λ<1`.
`ind` is a vector of integers.
The ith entry of this vector gives the *row* of the nonzero element of the ith
column of `Ck`. This element is +1 if `signs[i] > 0` and -1 otherwise. Each element
of `ind` **must be unique** such that the resulting matrix is orthogonal
**and** represents the group of special reflections and permutations.

Besides the default constructor, the following is also provided:

    LambdaMatrix(λ, D::Integer, random::Bool = true)
If `random == true` construct a `LambdaMatrix` of dimension `D` with a
"random" ``C_k`` matrix, by using `randperm(D)` as well as random `signs`. Else,
use `ind = collect(1:D)` and `sings = ones(Int8, D)`.

All possible combinations for `ind` and `sings` can be obtained by:
```julia
using Combinatorics
indperms = collect(permutations([1:D;], D))
p = falses(D); singsperm = [p[:]]
for i = 1:D
    p[i] = true
    append!(singsperm, multiset_permutations(p, D))
end
```

[1] : P. Schmelcher & F. Diakonos, Phys. Rev. Lett **78**, pp 4733 (1997)
"""
struct LambdaMatrix{T<:Real}
    λ::T
    ind::Vector{Int}
    sings::Vector{Int8}
end

dimension(L::LambdaMatrix) = size(L.Ck)[1]

function LambdaMatrix(λ::T, D::Integer, random::Bool = true) where {T<:Real}
    if random
        positions = randperm(D)
        signs = [rand(Bool) ? Int8(1) : Int8(-1) for i in 1:D]
    else
        positions = collect(1:D)
        signs = ones(Int8, D)
    end
    return LambdaMatrix(λ, positions, signs)
end

function *(Λ::LambdaMatrix, s::SVector{D,T}) where {D, T}
    gen = (Λ.λ*Λ.sings[i]*s[Λ.ind[i]] for i in 1:D)
    SVector{D,T}(gen...)
end

# All possible `ind`` (without the sign part)
# allCk = [nthperm([1:D;], n) for n in 1:factorial(D)]

ds = Systems.standardmap()
xs = linspace(0, 2π, 20); ys = linspace(0, 2π, 20)
o = order = 4
D = dimension(ds); T = eltype(ds.state)

maxiter = 100000
roundtol = 4
disttol = 1e-12
FP = SVector{D,T}[]

f = ds.eom

for i in 1:4
    Λ = LambdaMatrix(0.005, D, true)

    Sk = (state) -> state + Λ*(iterate(state, f, o) - state)


    for x in xs
        for y in ys
            st = SVector{2}(x, y)

            for i in 1:maxiter
                prevst = st
                st = Sk(prevst)
                if norm(prevst - st) < 1e-10
                    unist = round.(st, roundtol)
                    !(unist in FP) && push!(FP, unist)
                    break
                end
            end
        end
    end
end

length(FP) == 0 && error("no point converged")
FP = Dataset(FP)
using PyPlot
_phasespace(ds)
_plot_dataset(FP)
