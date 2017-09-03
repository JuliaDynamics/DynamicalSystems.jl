using StaticArrays, Requires
using IterTools: chain
export Dataset

"""
    Dataset{D, T, V}
A `Dataset` is an interface for vectors of vectors, inspired by
RecursiveArrayTools.jl. It contains **equally-sized datapoints** of length `D`,
represented by vectors of type `V`, containing numbers of type `T`.

This data representation is most of the time better than having a `Matrix`, but can be
used exactly like a matrix that has each of the columns be the timeseries of each of
the dynamic variables. For example,
```julia
data = timeseries(ds, 100.0) #this gives a dataset
data[:, 2] # this is the timeseries vector of the second variable
data[1] == data[1, :] # this is the first datapoint of the timeseries
data[5, 3] # this is the value of the third variable, at the 5th timepoint
```

Use `convert(Matrix, dataset)` to create a `Matrix`, and `convert(Dataset, matrix)`
to create a `Dataset` from a matrix. **Notice:**  `convert(Dataset, matrix)` assumes
that each column of the matrix represents one dynamic variable. If instead each
column of the matrix represents a datapoint, use `reinterpret(Dataset, matrix)`.
"""
struct Dataset{
  D, T<:Number, V<:Union{Vector{T}, SVector{D,T}}}
  data::Vector{V}
end
Dataset(v::Vector{SVector{D,T}}) where {D, T<:Number} = Dataset{D, T, SVector{D,T}}(v)
function Dataset(v::Vector{Vector{T}}) where {T<:Number}
  D = length(v[1])
  for i in 1:length(v)
    D != length(v[i]) && throw(ArgumentError(
    "All data-points in a Dataset must have same size"
    ))
  end
  return Dataset{D, T, Vector{T}}(v)
end

@inline dimension{D, T, V}(::Dataset{D, T, V}) = D
@inline Base.eltype{D, T, V}(d::Dataset{D, T, V}) = T
@inline vectype{D, T, V}(::Dataset{D, T, V}) = V

# Size:
@inline Base.length(d::Dataset) = length(d.data)
@inline Base.size(d::Dataset{D, T, V}) where {D, T, V} = (length(d.data), D)
@inline Base.size(d::Dataset, i::Int) = size(d)[i]
@inline Base.iteratorsize(d::Dataset) = Base.HasLength()
# 1D indexing  over the container elements:
@inline Base.getindex(d::Dataset, i::Int) = d.data[i]
@inline Base.endof(d::Dataset) = endof(d.data)
# 2D indexing exactly like if the dataset was a matrix
# with each column a dynamic variable
@inline Base.getindex(d::Dataset, i::Int, j::Int) = d.data[i][j]
@inline Base.getindex(d::Dataset, i::Colon, j::Int) = [d.data[k][j] for k in 1:length(d)]
@inline Base.getindex(d::Dataset, i::Int, j::Colon) = d.data[i]
@inline Base.getindex(d::Dataset, r::Range) = d.data[r]
# Itereting interface:
@inline Base.eachindex(D::Dataset) = Base.OneTo(length(D.data))
Base.start(d::Dataset) = 1
Base.next(d::Dataset, state) = (d[state], state + 1)
Base.done(d::Dataset, state) = state >= length(d.data) + 1

# Other commonly used functions:
function Base.push!(d::Dataset{D, T, V}, new_item::S) where {D, T, V, S}
  if V != S
    throw(ArgumentError("Can only push vectors of type $V into this Dataset"))
  end
  push!(d.data, new_item)
end

# Conversions:
function Base.convert{D, T}(::Type{Matrix}, d::Dataset{D,T})
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

### Pretty printing
function vectorname(d::Dataset{D, T, V}) where {D, T, V}
  if V <: Vector
    s = "Vector{$T}"
  elseif V <: SVector
    s = "SVector{$D,$T}"
  end
  return s
end

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
  tos = "$D-dimensional Dataset with $N points ($(vectorname(d))) :"*"\n"*s
  return tos
end

@require Juno begin
  function Juno.render(i::Juno.Inline, d::Dataset{D, T, V}) where {D, T, V}
    N = length(d)
    tos = matstring(d)
    Juno.render(
      Juno.Tree(Text(tos), []))
  end
end

function Base.show(io::IO, d::Dataset{D, T, V}) where {D, T, V}
  mat = convert(Matrix, d)
  n = length(d)
  s = sprint(io -> show(IOContext(io, limit=true), MIME"text/plain"(), mat))
  s = join(split(s, '\n')[2:end], '\n')
  println(io, "$D-dimensional Dataset with $n points ($(vectorname(d))) :")
  print(io, Text(s))
end
