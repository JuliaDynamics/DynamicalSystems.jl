# Estimating Delay Embedding Parameters
The following functions can estimate good values that can be used in
[`reconstruct`](@ref) for either the delay time or the
number of temporal neighbors.

## Delay Time
```@docs
estimate_delay
```

## Number of Temporal neighbors
This also corresponds to the "embedding dimension" but only for univariate timeseries.
```@docs
estimate_dimension
stochastic_indicator
```
---


## Mutual Information
WIP
