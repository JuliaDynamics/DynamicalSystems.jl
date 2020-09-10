# Numerical Data

!!! info "Trajectory and Timeseries"
    The word "timeseries" can be confusing, because it can mean a univariate (also called scalar or one-dimensional)
    timeseries or a multivariate (also called multi-dimensional) timeseries. To resolve this confusion, in
    **DynamicalSystems.jl** we have the following convention: **"timeseries"** always
    refers to a one-dimensional vector of numbers, which exists with respect to
    some other one-dimensional vector of numbers that corresponds to a time-vector.
    On the other hand,
    the word **"trajectory"** is used to refer to a *multi-dimensional* timeseries,
    which is of course simply a group/set of one-dimensional timeseries.
    A trajectory is represented by a [`Dataset`](@ref) and is the return type of [`trajectory`](@ref).


Trajectories (and in general sets in state space) in **DynamicalSystems.jl** are represented by a structure called `Dataset`
(while timeseries are standard Julia `Vector`s).
```@docs
Dataset
```

In essence a `Dataset` is simply a wrapper for a `Vector` of `SVector`s.
However, it is visually represented as a matrix, similarly to how numerical data would be printed on a spreadsheet (with time being the *column* direction).
It also offers a lot more functionality than just pretty-printing.
Besides the examples in the documentation string, you can also do:
```julia
using DynamicalSystems
hen = Systems.henon()
data = trajectory(hen, 10000) # this returns a dataset
for point in data
    # do stuff with each data-point
    # (vector with as many elements as system dimension)
end
```

Most functions from **DynamicalSystems.jl** that manipulate and use data are expecting a `Dataset`.
This allows us to define efficient methods that coordinate well with each other, like e.g. [`embed`](@ref).

## Dataset Functions
Functions that operate on datasets.
```@docs
minima
maxima
minmaxima
columns
```

## Dataset I/O
Input/output functionality for an `AbstractDataset` is already achieved using base Julia, specifically `writedlm` and `readdlm`.
To write and read a dataset, simply do:

```julia
using DelimitedFiles

data = Dataset(rand(1000, 2))

# I will write and read using delimiter ','
writedlm("data.txt", data, ',')

# Don't forget to convert the matrix to a Dataset when reading
data = Dataset(readdlm("data.txt", ',', Float64))
```
---

## Neighborhoods
Combining the excellent performance of [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) with the `AbstractDataset` allows us to define a function that calculates a "neighborhood" of a given point, i.e. finds other points near it. The different "types" of the neighborhoods are subtypes of `AbstractNeighborhood`.
```@docs
neighborhood
AbstractNeighborhood
```
