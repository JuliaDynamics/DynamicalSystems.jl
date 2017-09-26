# Nonlinear Timeseries Analysis
## Delay Coordinates Reconstruction
A mono-dimensional timeseries recorded in some manner from a dynamcal system can be used to gain information about the dynamics of the entire phase-space of the system. This can be done by reconstructing a new phase-space from the timeseries. One method that can do this is
what is known as [delay coordinates embedding](https://en.wikipedia.org/wiki/Takens%27_theorem).

In `DynamicalSystems.jl` this is done through the `reconstruct` interface:
```@docs
reconstruct
```
The returned `Reconstruction` object can be immediately passed into e.g. a method that calculates the
attractor dimension, or simply be used as any matrix. For example:
```julia
using DynamicalSystems
he = Systems.henon()
ts = timeseries(he, 20000)
D1 = information_dim(ts) # around 1.20
x = ts[:, 1] # some "recorded" timeseries
R = reconstruct(x, 2, 1) # delay coords. reconstruction of dimension 2 and delay 1
R[1] # first point of reconstruction, ≡ (x[1], x[2])
R[:, 2] # Second COLUMN of the reconstruction, ≡ x[2:end] since τ=1
D2 = information_dim(R) #around 1.20
```
The 2 numbers `D1` and `D2` are **very close**, but of course I knew before-hand good parameter
values for `D` and `τ`.

## Numerical Lyapunov Estimation
Given anytimeseries, one can first `reconstruct` it, and then calculate a maximum
lyapunov exponent for it, provided that the system the timeseries was recorded
from actually exhibits exponential separation of nearby trajectories. This is done
with
```@docs
numericallyapunov
```
This function has a total of 4 different approaches for the algorithmic process.
The methods for the neighborhoods are
```@docs
AbstractNeighborhood
neighborhood
```
Example:
```julia
using DynamicalSystems, PyPlot

ds = Systems.henon()
data = timeseries(ds, 100000)
x = data[:, 1]

ks = 1:20
ℜ = 1:10000
for di in [Euclidean(), Cityblock()]
  method = FixedMassNeighborhood(2)

  figure()
  title("Distance: $(di)")
  for D in [2, 4, 7]
    R = reconstruct(x, D, 1)
    E = numericallyapunov(R, ks, ℜ, di, method)
    λ = linear_region(ks, E)[2]
    plot(ks-1, E-E[1], label = "D=$D, λ=$(round(λ, 3))")
  end
  legend()
  tight_layout()
end
```
