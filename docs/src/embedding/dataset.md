# Numerical Data
Numerical data in **DynamicalSystems.jl** is most often represented by a structure called `Dataset`
```@docs
Dataset
```
---
In essence a `Dataset` is simply a container for a `Vector` of `SVector`s.
However, it
is visually represented as a matrix, similarly to how numerical data would be printed
on a spreadsheet (with time being the *column* direction). It also offers a lot more
functionality than just pretty-printing.
Besides the examples in the documentation string,
you can also do:
```julia
using DynamicalSystems
hen = Systems.henon()
data = trajectory(hen, 10000) # this returns a dataset
for point in data
# do stuff with each datapoint
# (vector with as many elements as system dimension)
end
```

Most functions from **DynamicalSystems.jl** that manipulate and use data are expecting an `AbstractDataset` subtype. This allows us to define efficient methods that coordinate well with each other, like e.g. [`neighborhood`](@ref).

If given a matrix, we first convert to `Dataset`. This means that you should *first
convert* your data to a `Dataset` if you want to call functions more than once, to avoid
constantly converting.

## Dataset Functions
Functions that operate on datasets.
```@docs
minima
maxima
minmaxima
columns
```
---
## Dataset I/O
Input/output functionality for an `AbstractDataset` is already achieved using base Julia, specifically `writedlm` and `readdlm`.

The thing to note is that all data of an `AbstractDataset` is contained within its field `data`.

To write and read a dataset, simply do:

```julia
using DelimitedFiles

data = Dataset(rand(1000, 2))

# I will write and read using delimiter ','
writedlm("data.txt", data.data, ',')

# Don't forget to convert the matrix to a Dataset when reading
data = Dataset(readdlm("data.txt", ',', Float64))
```
---

## Neighborhoods in a dataset
Combining the excellent performance of [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) with the `AbstractDataset` allows us to define a function that calculates a "neighborhood" of a given point, i.e. finds other points near it. The different "types" of the neighborhoods are subtypes of `AbstractNeighborhood`.
```@docs
neighborhood
AbstractNeighborhood
```
---
