# Estimating Delay Embedding Parameters
The following functions can estimate good values that can be used in
[`reconstruct`](@ref) for either the delay time or the
number of temporal neighbors.

## Delay Time
```@docs
estimate_delay
exponential_decay_fit
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
```

### Example
```@example estimateD
using DynamicalSystems, PyPlot

ds = Systems.roessler()
tr = trajectory(ds, 1000.0; dt = 0.05)

τ = estimate_delay(tr[:, 1], "mi_min") # first minimum of mutual information

figure();
for method in ["afnn", "fnn", "f1nn"]
    Ds = estimate_dimension(tr[:, 1], τ, 1:6, method)
    plot(1:6, Ds ./ maximum(Ds), label = method, marker = "o")
end
legend(); xlabel("\$\\gamma\$ (temporal neighbors)")
tight_layout()
savefig("estimateD.png"); nothing # hide
```
![](estimateD.png)

### Functions
```@docs
DelayEmbeddings.fnn
DelayEmbeddings.afnn
DelayEmbeddings.f1nn
stochastic_indicator
```
---
