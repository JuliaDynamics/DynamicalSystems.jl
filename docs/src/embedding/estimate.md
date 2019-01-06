# Estimating Delay Embedding Parameters
The following functions can estimate good values that can be used in
[`reconstruct`](@ref) for either the delay time or the
number of temporal neighbors.

## Delay Time
```@docs
estimate_delay
```
### Mutual Information
```@docs
mutualinformation
```

---

Besides the above method, there also exists code that computes mutual information in two other ways. Both ways are in the file `DelayEmbedding\src\old_mutual_info.jl`. The first way is the original algorithm of Fraser, while the second is the algorithm of Kraskov. Both of these implementations are inferior to the one exposed here (performance-wise).

## Embedding Dimension
```@docs
estimate_dimension
DelayEmbeddings.fnn
DelayEmbeddings.afnn
DelayEmbeddings.f1nn
stochastic_indicator
```
---
