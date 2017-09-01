# Mathematical functions that do some stuff
# very fast or very conveniently.
# Also includes helper functions

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

function min_pairwise_distance{D, T<:Real}(pts::Dataset{D, T})
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

#####################################################################################
#                                 Minima and Maxima                                 #
#####################################################################################
function minima{D, T<:Real}(data::Dataset{D, T})
    m = zeros(T, D).*T(Inf)
    for point in data
        for i in 1:D
            if point[i] < m[i]
                m[i] = point[i]
            end
        end
    end
    return SVector{D,T}(m)
end

function maxima{D, T<:Real}(data::Dataset{D, T})
    m = zeros(T, D).*T(-Inf)
    for point in data
        for i in 1:D
            if point[i] > m[i]
                m[i] = point[i]
            end
        end
    end
    return SVector{D, T}(m)
end
