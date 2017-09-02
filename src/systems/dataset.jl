using StaticArrays, Requires

"""
    Dataset{D, T, V}
A `Dataset` is an interface for `Vector`s of `Vector`s, inspired by
RecursiveArrayTools.jl. It contains numbers of type `T` and represents datapoints in
`D` dimensions, represented by vectors type `V`.

This data representation is most of the time better than having a `Matrix`, but can be
used exactly like a matrix that has each of the columns be the timeseries of each of
the dynamic variables.
"""
mutable struct Dataset{
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

### Pretty printing
function vectorname(d::Dataset{D, T, V}) where {D, T, V}
  if V <: Vector
    s = "Vector{$T}"
  elseif V <: SVector
    s = "SVector{$D,$T}"
  end
  return s
end

@require Juno begin
  function Juno.render(i::Juno.Inline, d::Dataset{D, T, V}) where {D, T, V}
    N = length(d)
    mat = convert(Matrix, d)
    s = sprint(io -> show(IOContext(io, limit=true), MIME"text/plain"(), mat))
    s = join(split(s, '\n')[2:end], '\n')
    Juno.render(
      Juno.Tree(Text("$D-dimensional Dataset with $N points ($(vectorname(d))) :"),
        [Text(s)])
    )
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

### Karlesson comments
# # A dataset will be a vector of vectors
# # you can always move quickly from a
# # vector of static arrays to a matrix using reinterpret though:
# << v = [rand(SVector{3}) for i in 1:1000];
#
# << @time V = reinterpret(Float64, v, (3, 1000));
#   0.000009 seconds (7 allocations: 288 bytes)
#
# << typeof(V)
# >> Array{Float64,2}
#
# # and then back,
# << reinterpret(SVector{3, Float64}, V, (1000,))
# >> 1000-element Array{SVector{3,Float64},1}:
#  [0.511705, 0.880423, 0.502126]
# ...
