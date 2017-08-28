# Mathematical functions that do some stuff
# very fast or very conveniently.
# Also includes helper functions
#####################################################################################
#                                   Conversions                                     #
#####################################################################################
"""
```julia
d2v(matrix) -> vectors
```
Convert a dataset (Matrix) to a tuple of vectors, using views.
"""
function d2v(matrix)
  D = size(matrix)[2] #dimension of the system (dynamic variables)
  vectors = typeof(view(matrix, :, 1))[] #initialize array that will store vectors
  for k in 1:D
    push!(vectors, view(matrix, :, k))
  end
  return vectors
end

v2d(vectors...) = hcat(vectors...)
v2d(vectors) = hcat(vectors...)

#####################################################################################
#                                Pairwse Distance                                   #
#####################################################################################
using NearestNeighbors, StaticArrays
export min_pairwise_distance

function min_pairwise_distance(cts::AbstractMatrix)
    if size(cts, 1) > size(cts, 2)
        error("Points must be close (transpose the Matrix)")
    end
    tree = KDTree(cts)
    min_d = Inf
    min_pair = (0, 0)
    for p in 1:size(cts, 2)
        inds, dists = knn(tree, view(cts, :, p), 1, false, i -> i == p)
        ind, dist = inds[1], dists[1]
        if dist < min_d
            min_d = dist
            min_pair = (p, ind)
        end
    end
    return min_pair, min_d
end

function min_pairwise_distance(pts::Vector{SVector{N, T}}) where {N, T}
    tree = KDTree(pts)
    min_d = Inf
    min_pair = (0, 0)
    for p in 1:length(pts)
        inds, dists = knn(tree, pts[p], 1, false, i -> i == p)
        ind, dist = inds[1], dists[1]
        if dist < min_d
            min_d = dist
            min_pair = (p, ind)
        end
    end
    return min_pair, min_d
end
