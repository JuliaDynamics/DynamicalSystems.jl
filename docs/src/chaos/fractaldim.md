# Fractal Dimension

There are numerous methods that one can use to calculate a so-called "dimension" of a dataset which in the context of dynamical systems is called the [Fractal dimension](https://en.wikipedia.org/wiki/Fractal_dimension).
Several variants of a computationally feasible fractal dimension exist and a simple usage example is shown in the [Fractal dimension example](@ref) subsection.

## Generalized dimension
Based on the definition of the Generalized entropy ([`genentropy`](@ref)), one can calculate an appropriate dimension, called *generalized dimension*:
```@docs
generalized_dim
molteno_dim
molteno_boxing
```

!!! danger "Be wary when using `generalized_dim`"
    As stated clearly by the documentation string, calling `generalized_dim` performs a lot of automated steps by calling other functions (see below)
    with default arguments. It is actually more like a convenient bundle than
    an actual function and therefore you should be careful
    when considering the validity of the returned number.

## Fractal dimension example
For an example of using entropies to compute the dimension of an attractor let's use everyone's favorite system:
```@example MAIN
using DynamicalSystems, PyPlot
lor = Systems.lorenz()
```

Our goal is to compute entropies for many different partition sizes `ε`, so let's get down to it:
```@example MAIN
tr = trajectory(lor, 100.0; Ttr = 10.0)

ες = ℯ .^ (-3.5:0.5:3.5) # semi-random guess
Hs = genentropy.(Ref(tr), ες; q = 1)
```

```@example MAIN
xs = @. -log(ες)
fig = figure()
plot(xs, Hs)
ylabel("\$H_1\$")
xlabel("\$-\\log (\\epsilon)\$");
fig.tight_layout(pad=0.3); fig
```

The slope of the linear scaling region of the above plot is the generalized dimension (of order q = 1) for the attractor of the Lorenz system.

Given that we _see_ the plot, we can estimate where the linear scaling region starts and ends. However, we can use the function [`linear_region`](@ref) to get an estimate of the result as well. First let's visualize what it does:

```@example MAIN
lrs, slopes = linear_regions(xs, Hs, tol = 0.25)

fig = figure()
for i in 1:length(lrs)-1
    plot(xs[lrs[i]:lrs[i+1]], Hs[lrs[i]:lrs[i+1]], marker = "o")
end
ylabel("\$H_1\$")
xlabel("\$-\\log (\\epsilon)\$");
fig.tight_layout(pad=0.3); fig
```

The [`linear_region`](@ref) function  computes the slope of the largest region:

```@example MAIN
Δ = linear_region(xs, Hs)[2]
```
This result is an approximation of the information dimension (because we used `q = 1`) of the Lorenz attractor.

---

The above pipeline is bundled in [`generalized_dim`](@ref).
For example, the dimension of the strange attractor of the
[`Systems.henon`](@ref) map, following the above approach but taking automated steps, is:
```@example MAIN
using DynamicalSystems
hen = Systems.henon()
tr = trajectory(hen, 200000)
D_hen = generalized_dim(tr; q = 1)
```

As a side note, be sure that you have enough data points, otherwise the values you will get will never be correct, as is demonstrated by
J.-P. Eckmann and D. Ruelle (see Physica D **56**, pp 185-187 (1992)).


## Linear scaling regions
And other utilities, especially [`linreg`](@ref), used in both [`generalized_dim`] and [`grassberger_dim`](@ref).
```@docs
linear_regions
linear_region
linreg
estimate_boxsizes
```

## Correlation sum based dimension
```@docs
correlationsum
grassberger_dim
boxed_correlationsum
data_boxing
estimate_r0_buenoorovio
estimate_r0_theiler
correlationsum_fixedmass
```

## Takens' best estimate
```@docs 
takens_best_estimate
```

## Kaplan-Yorke Dimension
```@docs
kaplanyorke_dim
```
Notice that calling this function requires you to pass the Lyapunov exponents in an ordered vector form (largest to smallest). Example:
```@example MAIN
using DynamicalSystems
hen = Systems.henon()
D_kp = kaplanyorke_dim(lyapunovspectrum(hen, 100000))
```
