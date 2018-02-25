A timeseries recorded in some manner from a dynamical system can be used to gain information about the dynamics of the entire phase-space of the system. This can be done by reconstructing a new phase-space from the timeseries. One method that can do this is what is known as [delay coordinates embedding](https://en.wikipedia.org/wiki/Takens%27_theorem) or delay coordinates *reconstruction*.

This is done through the `Reconstruction` interface:
```@docs
Reconstruction
```
---
Here are some examples of `Reconstruction`s of a 3D continuous chaotic system:
```julia
using DynamicalSystems, PyPlot

ds = Systems.gissinger()
data = trajectory(ds, 1000.0)

xyz = columns(data)

figure(figsize = (12,10))
k = 1
for i in 1:3
    for τ in [5, 30, 100]
        R = Reconstruction(xyz[i], 2, τ)
        ax = subplot(3,3,k)
        plot(R[:, 1], R[:, 2], color = "C$(k-1)", lw = 0.8)
        title("var = $i, τ = $τ")
        k+=1
    end
end

tight_layout()
suptitle("2D Reconstructions")
subplots_adjust(top=0.9)
```
![Example reconstructions](https://i.imgur.com/OZDBvu5.png)

A `Reconstruction` can also be made from a trajectory (i.e. multidimensional timeseries). For this to be possible, the number of timeseries must be known by Type:
```julia
a = rand(1000, 3) # my trajectory
R = Reconstruction(a, 2, 2) # errorino

A = Size(1000, 3)(a) # create array with the size as Type information
R = Reconstruction(A, 2, 2) #aaaall good

ds = Systems.towel(); tr = trajectory(ds, 10000)
R = Reconstruction(tr, 2, 2) # Dataset size is also known by Type!
```
```
(D=6, τ=2) - delay coordinates Reconstruction
 0.085     -0.121       0.075     0.76827   -0.038933    0.672094
 0.285813  -0.0675286   0.238038  0.681871   0.0508933   0.825263
 ⋮                                                       ⋮       
```

## Estimating Reconstruction Parameters
The following functions can (sometimes) estimate good values that can be used in
[`Reconstruction`](@ref). There are no guarantees though!
```@docs
estimate_delay
```
---
