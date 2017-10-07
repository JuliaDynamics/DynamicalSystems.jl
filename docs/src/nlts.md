# Nonlinear Timeseries Analysis
## Delay Coordinates Reconstruction
A mono-dimensional timeseries recorded in some manner from a dynamcal system can be used to gain information about the dynamics of the entire phase-space of the system. This can be done by reconstructing a new phase-space from the timeseries. One method that can do this is
what is known as [delay coordinates embedding](https://en.wikipedia.org/wiki/Takens%27_theorem).

In `DynamicalSystems.jl` this is done through the `reconstruct` interface:
```@docs
Reconstruction
reconstruct
```
---
As an example, let's pass a `Reconstruction` into e.g. a method that calculates the
attractor dimension:
```julia
using DynamicalSystems
he = Systems.henon()
ts = timeseries(he, 100000)
D1 = information_dim(ts) # around 1.20
x = ts[:, 1] # some "recorded" timeseries
R = reconstruct(x, 2, 1) # delay coords. reconstruction of dimension 2 and delay 1
R[1] # first point of reconstruction, ≡ (x[1], x[2])
R[:, 2] # Second COLUMN of the reconstruction, ≡ x[2:end] since τ=1
D2 = information_dim(R) #around 1.20
println("D2 - D1 = $(abs(D2- D1))")
# Prints something like 0.0032000971757613073
```
The 2 numbers `D1` and `D2` are **very close**, but of course I knew before-hand good parameter values for `D` and `τ` (I cheated, huhu!).

## Estimating Reconstruction Parameters
The following functions are provided estimate good values that can be used in
[`reconstruct`](@ref):
```
estimate_delay
```
More coming soon!

## Numerical Lyapunov Estimation
Given any timeseries, one can first `reconstruct` it, and then calculate a maximum
lyapunov exponent for it, provided that the system the timeseries was recorded
from actually exhibits exponential separation of nearby trajectories. This is done
with
```@docs
numericallyapunov
```
---
### Neighborhoods
The function `numericallyapunov` has a total of 4 different approaches for the algorithmic process, by
combining 2 types of distances with 2 types of neighborhoods.

The methods for the neighborhoods are subtypes of `AbstractNeighborhood`, and offer
a convenient way to find neighboring points to a given point in a dataset.
```@docs
AbstractNeighborhood
neighborhood
```
---
As you can see, the function `neighborhood` is generally applicable!

### Example of Numerical Lyapunov computation
```julia
using DynamicalSystems, PyPlot

ds = Systems.henon()
data = timeseries(ds, 100000)
x = data[:, 1]

ks = 1:20
ℜ = 1:10000
fig = figure(figsize=(10,6))
i = 1

for (i, di) in enumerate([Euclidean(), Cityblock()])
  subplot(1, 2, i)
  i+=1
  method = FixedMassNeighborhood(2)

  title("Distance: $(di)", size = 18)
  for D in [2, 4, 7]
    R = reconstruct(x, D, 1)
    E = numericallyapunov(R, ks;
    refstates = ℜ, distance = di, method = method)
    # The following operation:
    Δt = 1
    λ = linear_region(ks.*Δt, E)[2]
    # gives the linear slope, i.e. the Lyapunov exponent
    plot(ks-1, E-E[1], label = "D=$D, λ=$(round(λ, 3))")
    legend()
    tight_layout()
  end


end
```
which gives the result
![Lyapunov exponent of a timeseries](https://i.imgur.com/vbKo1vV.png)

### Bad Time-axis (`ks`) length

!!! danger "Large `ks`"
    This simply cannot be stressed enough! It is just too easy to overshoot
    the range at which the exponential expansion region is valid!

Let's revisit the example of the previous section:
```julia
ds = Systems.henon()
data = timeseries(ds, 100000)
x = data[:, 1]
```
The timeseries of length 100000 could be considered big. A time length of 100 seems
very small. Yet it turns out it is way too big! The following
```julia
ks = 1:100
R = reconstruct(x, 2, 1)
E = numericallyapunov(R, ks, method = FixedMassNeighborhood(2))
figure()
plot(ks-1, E-E[1])
println("Lyappunov: ", linear_region(ks, E)[2])
```
gives this plot:
![Bad time-vector example](https://i.imgur.com/wbpuBis.png)
and prints "Lyapunov: 0.4161...".

Notice that even though this value
for the Lyapunov exponent is correct, it happened to be correct simply due to the
jitter of the saturated region. Since the saturated region is much bigger
than the linear scaling region, if it wasn't that jittery the function
[`linear_region`](@ref) would not give the scaling of the linear region, but instead
a slope near 0! (or if you were to give bigger tolerance as a keyword argument)

### Case of a Continuous system
The process for continuous systems works identically with discrete, but one must be
a bit more thoughtful when choosing parameters. The following example
has comments to help the users get familiar with the process:
```julia
using DynamicalSystems, PyPlot

ds = Systems.lorenz() # Max lyapunov is around 0.90
# create a timeseries of 1 dimension
dt = 0.05
x = timeseries(ds, 1000.0; dt = dt)[:, 1]

τ1 = estimate_delay(x) #gives 7

# Reconstruct it (see section on estimating parameters)
figure()
for D in [4, 8], τ in [τ1, 15]
  R = reconstruct(x, D, τ)

  # I now know that I have to use much bigger ks than 1:20, as in the
  # continuous case! (See reference given in `numericallyapunovs`)
  ks1 = 0:200
  # I also know that I do not need that dense computations, since 1 increment
  # in k means increment of 0.05 real time
  ks2 = 0:4:200

  # Calculate lyapunovs:
  method = FixedMassNeighborhood(5) #5 nearest neighbors of each state

  # E1 = numericallyapunov(R, ks1; method = method)
  # λ1 = linear_region(ks1 .* dt, E1)[2]
  E2 = numericallyapunov(R, ks2; method = method)
  λ2 = linear_region(ks2 .* dt, E2)[2]


  # plot(ks1,E1.-E1[1], label = "dense, D=$(D), τ=$(τ), λ=$(round(λ1, 3))")
  plot(ks2,E2.-E2[1], label = "D=$(D), τ=$(τ), λ=$(round(λ2, 3))")
end

legend()
xlabel("k (0.05×t)")
ylabel("E - E(0)")
title("Continuous Reconstruction Lyapunov")
tight_layout()
```
As you can see, using ``\\tau = 15`` makes almost no sense! The estimates with
`τ = 7` though are very good (the actual value is around `λ = 0.89...`).
