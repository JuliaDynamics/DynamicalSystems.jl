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


## Datasets

Trajectories (and in general sets in state space) in **DynamicalSystems.jl** are represented by a structure called `Dataset`
(while timeseries are standard Julia `Vector`s).
```@docs
Dataset
```

In essence a `Dataset` is simply a wrapper for a `Vector` of `SVector`s.
However, it is visually represented as a matrix, similarly to how numerical data would be printed on a spreadsheet (with time being the *column* direction).
It also offers a lot more functionality than just pretty-printing.
Besides the examples in the documentation string, you can e.g. iterate over data points
```julia
using DynamicalSystems
hen = Systems.henon()
data = trajectory(hen, 10000) # this returns a dataset
for point in data
    # stuff
end
```

Most functions from **DynamicalSystems.jl** that manipulate and use multidimensional data are expecting a `Dataset`.
This allows us to define efficient methods that coordinate well with each other, like e.g. [`embed`](@ref).

## Dataset Functions
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

## Neighborhoods
Neighborhoods refer to the common act of finding points in a dataset that are nearby a given point (which typically belongs in the dataset).
**DynamicalSystems.jl** bases this interface on [Neighborhood.jl](https://julianeighbors.github.io/Neighborhood.jl/dev/).
You can go to its documentation if you are interested in finding neighbors in a dataset for e.g. a custom algorithm implementation.

For **DynamicalSystems.jl**, what is relevant are the two types of neighborhoods that exist:
```@docs
NeighborNumber
WithinRange
```


## Theiler window
The Theiler window is a concept that is useful when finding neighbors in a dataset that is coming from the sampling of a continuous dynamical system.
As demonstrated in the figure below, it tries to eliminate spurious "correlations" (wrongly counted neighbors) due to a potentially dense sampling of the trajectory (e.g. by giving small sampling time in [`trajectory`](@ref)).

The figure below demonstrates a typical `WithinRange` search around the black point with index `i`. Black, red and green points are found neighbors, but points within indices `j` that satisfy `|i-j| ≤ w` should *not* be counted as "true" neighbors.
These neighbors are typically the same around _any_ state space point, and thus wrongly bias calculations by providing a non-zero baseline of neighbors.
For the sketch below, `w=3` would have been used.

Typically a good choice for `w` coincides with the choice an optimal delay time, see [`estimate_delay`](@ref), for any of the timeseries of the dataset.

![](theiler.png)
