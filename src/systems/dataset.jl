using StaticArrays, Requires
using IterTools: chain
export Dataset

abstract type AbstractDataset{D} end

# Size:
@inline Base.length(d::AbstractDataset) = length(d.data)
@inline Base.size(d::AbstractDataset{D}) where {D} = (length(d.data), D)
@inline Base.size(d::AbstractDataset, i::Int) = size(d)[i]
@inline Base.iteratorsize(d::AbstractDataset) = Base.HasLength()
# 1D indexing  over the container elements:
@inline Base.getindex(d::AbstractDataset, i::Int) = d.data[i]
@inline Base.endof(d::AbstractDataset) = endof(d.data)
# 2D indexing exactly like if the dataset was a matrix
# with each column a dynamic variable
@inline Base.getindex(d::AbstractDataset, i::Int, j::Int) = d.data[i][j]
@inline Base.getindex(d::AbstractDataset, i::Colon, j::Int) =
[d.data[k][j] for k in 1:length(d)]
@inline Base.getindex(d::AbstractDataset, i::Int, j::Colon) = d.data[i]
@inline Base.getindex(d::AbstractDataset, r::Range) = d.data[r]
# Itereting interface:
@inline Base.eachindex(D::AbstractDataset) = Base.OneTo(length(D.data))
@inline Base.start(d::AbstractDataset) = 1
@inline Base.next(d::AbstractDataset, state) = (d[state], state + 1)
@inline Base.done(d::AbstractDataset, state) = state >= length(d.data) + 1
# Append
Base.append!(d1::AbstractDataset, d2::AbstractDataset) = append!(d1.data, d2.data)

# Other commonly used functions:
@inline Base.push!(d::AbstractDataset, new_item) = push!(d.data, new_item)
@inline dimension(::AbstractDataset{D}) where {D} = D

"""
    Dataset{D, T} <: AbstractDataset{D}
A dedicated interface for datasets, i.e. vectors of vectors.
It contains **equally-sized datapoints** of length `D`,
represented by `SVector{D, T}`, containing numbers of type `T`.

`Dataset` has methods for iterating, `convert`, `push!`, `append!` and other
basic functions.

The internal data representation is more efficient than having a `Matrix` and
also leads
to faster numerical computation of other quantities (like e.g. entropies). However,
it can be used exactly like a matrix that has each of the columns be the
timeseries of each of the dynamic variables. [`trajectory`](@ref) always returns
a `Dataset`.

For example,
```julia
data = trajectory(ds, 100.0) #this gives a dataset that behaves like a matrix
data[:, 2] # this is the second variable timeseries
data[1] == data[1, :] # this is the first datapoint of the dataset (D-dimensional)
data[5, 3] # this is the value of the third variable, at the 5th timepoint
```

Use `convert(Matrix, dataset)` to create a `Matrix`, and `convert(Dataset, matrix)`
to create a `Dataset` from a matrix. **Notice:**  `convert(Dataset, matrix)` assumes
that each column of the matrix represents one dynamic variable. If instead each
column of the matrix represents a datapoint, use `reinterpret(Dataset, matrix)`.

If you have various timeseries vectors `x, y, z, ...` pass them like
`Dataset(x, y, z, ...)`.
"""
struct Dataset{D, T<:Number} <: AbstractDataset{D}
    data::Vector{SVector{D,T}}
end

function Dataset(v::Vector{Vector{T}}) where {T<:Number}
    D = length(v[1])
    L = length(v)
    data = Vector{SVector{3,Float64}}(L)
    for i in 1:length(v)
        D != length(v[i]) && throw(ArgumentError(
        "All data-points in a Dataset must have same size"
        ))
        data[i] = SVector{D,T}(v[i])
    end
    return Dataset{D, T}(data)
end

function Dataset(vecs::Vararg{Vector{T}}) where {T<:Real}
    length(vecs) == 1 && throw(ArgumentError("Give more than 1 vectors to Dataset."))
    L = length(vecs[1])
    D = length(vecs)
    for i in 2:length(vecs)
        length(vecs[i]) != L && throw(argumentError(
        "All vectors that make the dataset must have equal length"))
    end
    data = SVector{D, T}[]
    for i in 1:L
        t = (vecs[j][i] for j in 1:D) #datapoint generator
        push!(data, SVector{D, T}(t...))
    end
    return Dataset(data)
end

@inline Base.eltype(d::Dataset{D, T}) where {D,T} = T

# Conversions:
function Base.convert(::Type{Matrix}, d::Dataset{D,T}) where {D, T}
  mat = Matrix{T}(length(d), D)
  for i in 1:length(d)
    mat[i,:] .= d.data[i]
  end
  mat
end

function Base.convert(::Type{Dataset}, mat::AbstractMatrix)
    D = size(mat, 2); T = eltype(mat)
    d = SVector{D, T}[]
    for i in 1:size(mat, 1)
        push!(d, SVector{D, T}(view(mat, i, :)))
    end
    return Dataset(d)
end

function Base.reinterpret(::Type{Dataset}, mat::Array{T,2}) where {T<:Real}
    s = size(mat)
    D = s[1]; N = s[2]
    Dataset(reinterpret(SVector{D, T}, mat, (N,)))
end

function Base.convert(::Type{Dataset}, y::Vector{T}) where {T}
    data = reinterpret(SVector{1,T}, y, (length(y),))
    return Dataset(data)
end



### Pretty printing
function matstring(d::Dataset)
    N = length(d); D = dimension(d)
    if N > 50
        mat = zeros(eltype(d), 50, D)
        for (i, a) in enumerate(chain(1:25, N-24:N))
            mat[i, :] .= d[a]
        end
    else
        mat = convert(Matrix, d)
    end
    s = sprint(io -> show(IOContext(io, limit=true), MIME"text/plain"(), mat))
    s = join(split(s, '\n')[2:end], '\n')
    tos = "$D-dimensional Dataset with $N points:"*"\n"*s
    return tos
end

@require Juno begin
    function Juno.render(i::Juno.Inline, d::Dataset{D, T}) where {D, T}
    N = length(d)
    tos = matstring(d)
    Juno.render(Juno.Tree(Text(tos), []))
    end
end

function Base.show(io::IO, d::Dataset{D, T}) where {D, T}
    mat = convert(Matrix, d)
    n = length(d)
    s = sprint(io -> show(IOContext(io, limit=true), MIME"text/plain"(), mat))
    s = join(split(s, '\n')[2:end], '\n')
    println(io, "$D-dimensional Dataset with $n points:")
    print(io, Text(s))
end


#####################################################################################
#                                 Minima and Maxima                                 #
#####################################################################################
function minima(data::Dataset{D, T}) where {D, T<:Real}
    m = zeros(T, D) .+ T(Inf)
    for point in data
        for i in 1:D
            if point[i] < m[i]
                m[i] = point[i]
            end
        end
    end
    return SVector{D,T}(m)
end

function maxima(data::Dataset{D, T}) where {D, T<:Real}
    m = zeros(T, D) .+ T(-Inf)
    for point in data
        for i in 1:D
            if point[i] > m[i]
                m[i] = point[i]
            end
        end
    end
    return SVector{D, T}(m)
end
