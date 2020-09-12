# Optimal DCE Parameters
This page discusses and provides algorithms for estimating optimal parameters to do Delay Coordinates Embedding (DCE) with.

The approaches can be grouped into two schools:
1. **Independent**, where one tries to independently find the best value for a delay time `τ` and an embedding dimension `d`.
2. **Unified**, where at the same time an optimal combination of `τ, d` is found.

The independent approach is something "old school", while recent scientific research has shifted almost exclusively to unified approaches.

In addition, the unified approaches are the only ones that can accommodate multi-variate inputs. This means that if you have multiple measured input timeseries, you should be able to take advantage of all of them for the best possible embedding of the dynamical system's set.

## Independent delay time
```@docs
estimate_delay
exponential_decay_fit
```
### Mutual Information
```@docs
mutualinformation
```

## Independent embedding dimension
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
DelayEmbeddings.stochastic_indicator
```
---

## Unified approach
Several algorithms have been created to implement a unified approach to delay coordinates embedding. You can find some implementations below:
```@docs
pecora
uzal_cost
garcia_almeida_embedding
mdop_embedding
```

### Low-level functions of unified approach
```@docs
DelayEmbeddings.n_statistic
DelayEmbeddings.beta_statistic
DelayEmbeddings.mdop_maximum_delay
```
