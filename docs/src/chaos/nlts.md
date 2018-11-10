# Nonlinear Timeseries Analysis

## Numerical Lyapunov Exponent
Given any timeseries, one can first [`reconstruct`](@ref) it using
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
```@example entropy
using DynamicalSystems, PyPlot

ds = Systems.henon()
data = trajectory(ds, 100000)
x = data[:, 1] #fake measurements for the win!

ks = 1:20
ℜ = 1:10000
fig = figure(figsize=(10,6))

for (i, di) in enumerate([Euclidean(), Cityblock()])
    subplot(1, 2, i)
    ntype = FixedMassNeighborhood(2)
    title("Distance: $(di)", size = 18)
    for D in [1, 3, 6]
        R = reconstruct(x, D, 1)
        E = numericallyapunov(R, ks;
        refstates = ℜ, distance = di, ntype = ntype)
        Δt = 1
        λ = linear_region(ks.*Δt, E)[2]
        # gives the linear slope, i.e. the Lyapunov exponent
        plot(ks .- 1, E .- E[1], label = "D=$D, λ=$(round(λ, digits = 3))")
        legend()
        tight_layout()
    end
end
savefig("numerlyap.png"); nothing # hide
```
![](numerlyap.png)


### Bad Time-axis (`ks`) length

!!! danger "Large `ks`"
    This simply cannot be stressed enough! It is just too easy to overshoot
    the range at which the exponential expansion region is valid!

Let's revisit the example of the previous section:
```@example entropy
ds = Systems.henon()
data = trajectory(ds, 100000)
x = data[:, 1]
length(x)
```
The timeseries of such length could be considered big. A time length of 100 seems
very small. Yet it turns out it is way too big! The following
```@example entropy
ks = 1:100
R = reconstruct(x, 1, 1)
E = numericallyapunov(R, ks, ntype = FixedMassNeighborhood(2))
figure()
plot(ks .- 1, E .- E[1])
title("Lyappunov: $(linear_region(ks, E)[2])")
savefig("badlyap.png"); nothing # hide
```
![](badlyap.png)

Notice that even though this value
for the Lyapunov exponent is correct, it happened to be correct simply due to the
jitter of the saturated region. Since the saturated region is much bigger
than the linear scaling region, if it wasn't that jittery the function
[`linear_region`](@ref) would not give the scaling of the linear region, but instead
a slope near 0! (or if you were to give bigger tolerance as a keyword argument)

### Case of a Continuous system
The process for continuous systems works identically with discrete, but one must be
a bit more thoughtful when choosing parameters. The following example helps the users get familiar with the process:
```@example entropy
using DynamicalSystems, PyPlot

ntype = FixedMassNeighborhood(5) #5 nearest neighbors of each state

ds = Systems.lorenz()
# create a timeseries of 1 dimension
dt = 0.05
x = trajectory(ds, 1000.0; dt = dt)[:, 1]
```

We know that we have to use much bigger `ks` than `1:20`, because this is a continuous case! (See reference given in `numericallyapunovs`)
```@example entropy
ks1 = 0:200
```
and in fact it is even better to not increment the `ks` one by one but instead do
```@example entropy
ks2 = 0:4:200
```
Now we plot some example computations
```@example entropy
figure()
for D in [3, 7], τ in [7, 15]
    r = reconstruct(x, D, τ)

    # E1 = numericallyapunov(r, ks1; ntype = ntype)
    # λ1 = linear_region(ks1 .* dt, E1)[2]
    E2 = numericallyapunov(r, ks2; ntype = ntype)
    λ2 = linear_region(ks2 .* dt, E2)[2]

    # plot(ks1,E1.-E1[1], label = "dense, D=$(D), τ=$(τ), λ=$(round(λ1, 3))")
    plot(ks2,E2.-E2[1], label = "D=$(D), τ=$(τ), λ=$(round(λ2, digits = 3))")
end

legend()
xlabel("k (0.05×t)")
ylabel("E - E(0)")
title("Continuous Reconstruction Lyapunov")
tight_layout()
savefig("continuousnumlyap.png"); nothing # hide
```
![](continuousnumlyap.png)

As you can see, using `τ = 15` is not a great choice! The estimates with
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

```@example entropy
using DynamicalSystems, PyPlot

ds = Systems.gissinger()
data = trajectory(ds, 1000.0, dt = 0.05)
x = data[:, 1]

L = length(x)
s = x .+ 0.5rand(L) #add noise

U, S = broomhead_king(s, 40)
summary(U)
```

Now let's simply compare the above result with the one you get from doing a "standard" call to [`reconstruct`](@ref):
```@example entropy
figure(figsize= (10,6))
subplot(1,2,1)
plot(U[:, 1], U[:, 2])
title("Broomhead-King of s")

subplot(1,2,2)
R = reconstruct(s, 1, 30)
plot(columns(R)...; color = "C3")
title("2D reconstruction of s")
tight_layout()
savefig("broomhead_king.png"); nothing # hide
```
![](broomhead_king.png)

we have used the same system as in the [delay coordinates reconstruction](/definition/reconstruction) example, and picked the optimal
delay time of `τ = 30` (for same `dt = 0.05`). Regardless, the vanilla delay coordinates is much worse than the Broomhead-King coordinates.
