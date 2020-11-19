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
pecuzal_embedding
```

### Example
In following we illustrate the functionality of the PECUZAL method on three
examples. We start with a univariate case, i.e. we only feed in one time series,
here the *x*-component of the Lorenz system.  
```@example MAIN
using DelayEmbeddings

lo = Systems.lorenz([1.0, 1.0, 50.0])
tr = trajectory(lo, 100; dt = 0.01, Ttr = 10)

s = vec(tr[:, 1]) # input timeseries = x component of lorenz
theiler = estimate_delay(s, "mi_min") # estimate a Theiler window
Tmax = 100 # maximum possible delay

Y, τ_vals, ts_vals, Ls , εs = pecuzal_embedding(s; τs = 0:Tmax , w = theiler)

println(τ_vals)
println(ts_vals)
println(Ls)
```
The output reveals that PECUZAL suggests a 3-dimensional embedding out of the
un-lagged time series as the 1st component of the reconstruction, the time
series lagged by 18 samples as the 2nd component and the time series lagged by
9 samples as the 3rd component. The minimum obtained *L*-value in the 3rd
embedding cycle has been ~-2.63, after which the algorithm breaks.
```@example MAIN
using PyPlot

figure(figsize=(14., 8.))
subplot(1,2,1, projection="3d")
plot3D(Y[:,1], Y[:,2], Y[:,3],"gray")
title("PECUZAL reconstructed x-component of Lorenz System")
xlabel("x(t+$(τ_vals[1]))")
ylabel("x(t+$(τ_vals[2]))")
zlabel("x(t+$(τ_vals[3]))")
grid()

subplot(1,2,2, projection="3d")
plot3D(tr[:,1], tr[:,2], tr[:,3],"gray")
title("Original Lorenz System")
xlabel("x(t)")
ylabel("y(t)")
zlabel("z(t)")
grid()

tight_layout()
savefig("pecuzal_uni.png"); nothing # hide
```
![](pecuzal_uni.png)

We can also look at the output of the low-level function leading to the results,
here the *continuity statistic*.
```@example MAIN
using PyPlot

figure(figsize=(8., 5.))
plot(εs[:,1], label="1st embedding cycle")
scatter([τ_vals[2]], [εs[τ_vals[2],1]])
plot(εs[:,2], label="2nd embedding cycle")
scatter([τ_vals[3]], [εs[τ_vals[3],2]])
plot(εs[:,3], label="3rd embedding cycle")
title("Continuity statistics for PECUZAL embedding of Lorenz x-component")
xlabel("delay τ")
ylabel("⟨ε⋆⟩")
legend(loc="upper left")
grid()
savefig("continuity_uni.png"); nothing # hide
```
![](continuity_uni.png)

Similar to the approach in the preceding example, we now highlight the capability
of the PECUZAL embedding method for a multivariate input. The idea is now to feed
in all three time series to the algorithm, even though this is a very
far-from-reality example. We already have an adequate representation of the
system we want to reconstruct, namely the three time series from the numerical
integration. But let us see what PECUZAL suggests for a reconstruction.
```@example MAIN
# compute Theiler window
w1 = estimate_delay(tr[:,1], "mi_min")
w2 = estimate_delay(tr[:,2], "mi_min")
w3 = estimate_delay(tr[:,3], "mi_min")
w = maximum(hcat(w1,w2,w3))
Y_m, τ_vals_m, ts_vals_m, Ls_m , εs_m = pecuzal_embedding(tr; τs = 0:Tmax , w = theiler)

println(τ_vals_m)
println(ts_vals_m)
println(Ls_m)
```
PECUZAL offers a 3-dimensional embedding using the un-lagged *z*- and *x*-component
as 1st and 3rd component of the reconstruction vectors, as well as the *x*-component
lagged by 12 samples.

```@example MAIN

ts_str = ["x", "y", "z"]

figure(figsize=(14., 8.))
subplot(1,2,1, projection="3d")
plot3D(Y_m[:,1], Y_m[:,2], Y_m[:,3],"gray")
title("PECUZAL reconstructed Lorenz System")
xlabel("$(ts_str[ts_vals_m[1]])(t+$(τ_vals_m[1]))")
ylabel("$(ts_str[ts_vals_m[2]])(t+$(τ_vals_m[2]))")
zlabel("$(ts_str[ts_vals_m[3]])(t+$(τ_vals_m[3]))")
grid()

subplot(1,2,2, projection="3d")
plot3D(tr[:,1], tr[:,2], tr[:,3],"gray")
title("Original Lorenz System")
xlabel("x(t)")
ylabel("y(t)")
zlabel("z(t)")
grid()

tight_layout()
savefig("pecuzal_multi.png"); nothing # hide
```
![](pecuzal_multi.png)

Finally we show what PECUZAL does with a non-deterministic source:

```@example MAIN
using Random

# Dummy input
d1 = randn(1000)
d2 = rand(1000)
Tmax = 100
dummy_set = Dataset(hcat(d1,d2))

w1 = estimate_delay(d1, "mi_min")
w2 = estimate_delay(d2, "mi_min")
theiler = minimum(hcat(w1,w2))

Y_d, τ_vals_d, ts_vals_d, Ls_d , ε★_d = pecuzal_embedding(dummy_set; τs = 0:Tmax , w = theiler)
```

### Low-level functions of unified approach
```@docs
DelayEmbeddings.n_statistic
DelayEmbeddings.beta_statistic
DelayEmbeddings.mdop_maximum_delay
```
