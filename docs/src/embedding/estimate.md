# Optimal DCE Parameters
This page discusses and provides algorithms for estimating optimal parameters to do Delay Coordinates Embedding (DCE) with.

The approaches can be grouped into two schools:
1. **Independent** (also called **traditional**), where one tries to independently find the best value for a delay time `τ` and an embedding dimension `d`.
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
optimal_traditional_de
delay_afnn
delay_ifnn
delay_fnn
delay_f1nn
```

### Example
```@example MAIN
using DynamicalSystems, PyPlot

ds = Systems.roessler()
# This trajectory is a chaotic attractor with fractal dim ≈ 2
# therefore the set needs at least embedding dimension of 3
tr = trajectory(ds, 1000.0; dt = 0.05)
x = tr[:, 1]

dmax = 7
fig = figure()
for (i, method) in enumerate(["afnn", "fnn", "f1nn", "ifnn"])
    # Plot statistic used to estimate optimal embedding
    # as well as the automated output embedding
    𝒟, τ, E = optimal_traditional_de(x, method; dmax)
    plot(1:dmax, E; label = method, marker = "o", ms = 5, color = "C$(i-1)")
    optimal_d = size(𝒟, 2)
    scatter(optimal_d, E[optimal_d]; marker = "s", s = 100, color = "C$(i-1)")
end
legend(); xlabel("embedding dimension")
ylabel("estimator")
tight_layout()
fig.savefig("estimateD.png"); nothing # hide
```
![](estimateD.png)

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
