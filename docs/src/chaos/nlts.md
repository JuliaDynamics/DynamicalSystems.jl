# Nonlinear Timeseries Analysis

## Estimating Reconstruction Parameters
The following functions can estimate good values that can be used in
[`Reconstruction`](@ref) for either the delay time or the
dimension.
```@docs
estimate_delay
estimate_dimension
stochastic_indicator
```
---

## Numerical Lyapunov Exponent
Given any timeseries, one can first obtain a [`Reconstruction`](@ref) from it using
delay coordinates, and then calculate a maximum
Lyapunov exponent for it. This is done
with
```@docs
numericallyapunov
```
---
The function `numericallyapunov` has a total of 4 different approaches for the algorithmic process, by
combining 2 types of distances with 2 types of neighborhoods.

### Example of Numerical Lyapunov computation
```julia
using DynamicalSystems, PyPlot

ds = Systems.henon()
data = trajectory(ds, 100000)
x = data[:, 1] #fake measurements for the win!

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
    R = Reconstruction(x, D, 1)
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
data = trajectory(ds, 100000)
x = data[:, 1]
```
The timeseries of length 100000 could be considered big. A time length of 100 seems
very small. Yet it turns out it is way too big! The following
```julia
ks = 1:100
R = Reconstruction(x, 2, 1)
E = numericallyapunov(R, ks, method = FixedMassNeighborhood(2))
figure()
plot(ks-1, E-E[1])
println("Lyappunov: ", linear_region(ks, E)[2])
```
gives this plot:
![Bad time-vector example](https://i.imgur.com/wbpuBis.png)
and prints
```
Lyapunov: 0.4161...
```

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
x = trajectory(ds, 1000.0; dt = dt)[:, 1]

# Reconstruct it
figure()
for D in [4, 8], τ in [τ1, 15]
    R = Reconstruction(x, D, τ)

    # I now know that I have to use much bigger ks than 1:20, because this is a
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
which produces:
![Continuous Reconstruction exaple](https://i.imgur.com/lgyGLDv.png)
As you can see, using `τ = 15` makes almost no sense! The estimates with
`τ = 7` though are very good (the actual value is around `λ ≈ 0.89...`).

## Broomhead-King Coordinates
```@docs
broomhead_king
```
---
This alternative/improvement of the traditional delay coordinates can be a very
powerful tool. An example where it shines is noisy data where there is the effect
of superficial dimensions due to noise.

Take the following example where we produce noisy data from a system and then use
Broomhead-King coordinates as an alternative to "vanilla" delay coordinates:

```julia
using DynamicalSystems, PyPlot
using Distributions # for random numbers

ds = Systems.gissinger()
data = trajectory(ds, 1000.0, dt = 0.05)
x = data[:, 1]

L = length(x)
distrib = Normal(0, 0.1)
s = x .+ rand(distrib, L)

U, S = broomhead_king(s, 40)

figure(figsize= (10,6))
subplot(1,2,1)
plot(U[:, 1], U[:, 2])
title("Broomhead-King of s")

subplot(1,2,2)
R = Reconstruction(s, 2, 30)
plot(columns(R)...; color = "C3")
title("2D Reconstruction of s")

tight_layout()
```
![Broomhead-King example](https://i.imgur.com/xVWDjuh.png)

we have used the same system as in the [delay coordinates reconstruction](/definition/reconstruction) example, and picked the optimal
delay time of `τ = 30` (for same `dt = 0.05`). Regardless, the vanilla delay coordinates fail spectacularly
when compared with the Broomhead-King coordinates.
