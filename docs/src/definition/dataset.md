# Numerical Data
Numerical data in **DynamicalSystems.jl** is represented by a structure called
`Dataset`
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
hen == Systems.henon()
data = trajectory(hen, 10000) # this returns a dataset
for point in data
# do stuff with each datapoint (vector with as many elements as system dimension)
end
```

All functions of our ecosystem that manipulate and use data are expecting an `AbstractDataset` subtype. This allows us to define efficient methods that coordinate
well with other packages, like e.g. [`neighborhood`](@ref).

If given a matrix, we first convert to `Dataset`. This means that you should *first
convert* your data to a `Dataset` if you want to call functions more than once, to avoid
constantly converting.

## Extrema
Functions that find extrema of a dataset for each timeseries:
```@docs
minima
maxima
minmaxima
```

## Dataset IO
In addition to the above, we also offer (very basic) functions that read/write a
`Dataset` from/to a delimited text file:
```@docs
read_dataset
write_dataset
```
---
For example
```julia
using DynamicalSystems

ds = Systems.towel()
data = trajectory(ds, 1000)

# Write comma-delimited file:
write_dataset("test.csv", data, ',')
# Read comma-delimited file:
read_dataset("test.csv", Dataset{2, Float64}, ',')
```
