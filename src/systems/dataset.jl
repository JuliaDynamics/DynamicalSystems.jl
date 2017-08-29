using StaticArrays, Requires
import Base: size, getindex, convert
# Dataset:
const Dataset{D, T} = Vector{SVector{D, T}}
v = [rand(SVector{3}) for i in 1:1000]

@inline Base.getindex(d::Dataset, i::Int, j::Int) = d[i][j]
@inline Base.getindex(d::Dataset, i::Colon, j::Int) = [d[k][j] for k in 1:length(d)]
@inline Base.getindex(d::Dataset, i::Int, j::Colon) = d[i]
@inline dimension{D,T}(::Dataset{D,T}) = D

function convert{D, T}(::Type{Matrix}, d::Dataset{D,T})
  mat = Matrix{T}(length(d), D)
  for i in 1:length(d)
    mat[i,:] .= d[i]
  end
  mat
end

function convert(::Type{Dataset}, mat::AbstractMatrix)
  D = size(mat, 2); T = eltype(mat)
  d = SVector{D, T}[]
  for i in 1:size(mat, 1)
    push!(d, SVector{D, T}(view(mat, i, :)))
  end
  return d
end

### Pretty printing
@require Juno begin
  function Juno.render(i::Juno.Inline, d::Dataset{D,T}) where {D,T}
    N = length(d)
    Juno.render(
      Juno.Tree(Text("$D-dimensional Dataset{$T} with $N points"),
        [Juno.Text(stringrep(d))])
    )
  end
end
function stringrep(d::Dataset, nshow::Int=21, ndigits::Int=5)
  iseven(nshow) && throw(ArgumentError("nshow must be odd"))
  n = min(nshow, length(d))
  s = ""
  if n < nshow
    for i in 1:n
      s *= string(v[i], "\n")
    end
  else
    for i in 1:(n÷2)
      s *= string(v[i], "\n")
    end
    padlen = length(string(v[n÷2]))÷2
    s *= lpad("⋮\n", padlen)
    for i in (length(v)-(n÷2)+1):length(v)
      s *= string(v[i], "\n")
    end
  end
  return s
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
