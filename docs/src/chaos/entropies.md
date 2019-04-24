# Entropies and Dimensions

## Generalized Entropy
In the study of dynamical systems there are many quantities that identify as "entropy".
Notice that these quantities are not the more commonly known
[thermodynamic ones](https://en.wikipedia.org/wiki/Entropy), used in Statistical Physics. Rather, they are more like the to the entropies of [information theory](https://en.wikipedia.org/wiki/Entropy_(information_theory)), which represents
information contained within a dataset, or information about the dimensional
scaling of a dataset.

---

The main way of computing entropies in **DynamicalSystems.jl** is the "generalized entropy":
```@docs
genentropy
```
---
Basically, given a [dataset](system_definition/#numerical-data) you can partition it into boxes to calculate an entropy. See below for a detailed example.

!!! tip "Worried about memory overflow? Don't be!"
    Partitioning the dataset (i.e. doing a *histogram*) is in general a costly
    operation that depends exponentially on the number of dimensions of the data
    and algebraically to the box size `ε`.

    However, in this specific case the partition process has some special aspects
    that can be taken advantage
    of, reducing tremendously the memory allocation and spent time!


The function used internally by `genentropy` is `non0hist`:
```@docs
non0hist
```
---

## Attractor Dimension Estimation
There are numerous methods that one can use to calculate a so-called "dimension" of a
dataset, like for example the [Fractal dimension](https://en.wikipedia.org/wiki/Fractal_dimension). This real number can offer
a lot of information about the object that the dataset represents.

Based on the definition of the [generalized entropy](#ChaosTools.genentropy), one can calculate an appropriate
dimension, called *generalized dimension*:
```@docs
generalized_dim
```
---
!!! danger "Be wary when using `generalized_dim`"
    As stated clearly by the documentation string, calling `generalized_dim` performs a lot of automated steps by calling other functions (see below)
    with default arguments. It is actually more like a convenient bundle than
    an actual function and therefore you should be careful
    when considering the validity of the returned number.

```@docs
estimate_boxsizes
linear_regions
linear_region
```
---

## Example
For an example of using entropies to compute the dimension of an attractor let's use everyone's favorite system:
```@example entropy
using DynamicalSystems, PyPlot
lor = Systems.lorenz()
```

Our goal is to compute entropies for many different partition sizes `ε`, so let's get down to it:
```@example entropy
tr = trajectory(lor, 100.0; Ttr = 10.0)

ες = ℯ .^ (-3.5:0.5:3.5) # semi-random guess
Hs = genentropy.(1, ες, Ref(tr))
```

```@example entropy
xs = @. -log(ες)
figure()
plot(xs, Hs)
ylabel("\$H_1\$")
xlabel("\$-\\log (\\epsilon)\$");
savefig("genentropy1.png"); nothing # hide
```
![](genentropy1.png)

The slope of the linear scaling region of the above plot is the generalized dimension (of order α = 2) for the attractor of the Lorenz system.

Given that we _see_ the plot, we can estimate where the linear scaling region starts and ends. However, we can use the function [`linear_region`](@ref) to get an estimate of the result as well. First let's visualize what it does:

```@example entropy
lrs, slopes = linear_regions(xs, Hs, tol = 0.25)

figure()
for i in 1:length(lrs)-1
    plot(xs[lrs[i]:lrs[i+1]], Hs[lrs[i]:lrs[i+1]], marker = "o")
end
ylabel("\$H_1\$")
xlabel("\$-\\log (\\epsilon)\$");
savefig("genentropy2.png"); nothing # hide
```
![](genentropy2.png)

The [`linear_region`](@ref) function  computes the slope of the largest region:

```@example entropy
linear_region(xs, Hs)[2]
```
This result is an approximation of the information dimension (because we used `α = 1`) of the Lorenz attractor.

---

The above pipeline is bundled in [`generalized_dim`](@ref).
For example, the dimension of the strange attractor of the
[Hénon map](system_definition/#DynamicalSystems.Systems.henon), following the above approach but taking automated steps, is:
```@example entropy
using DynamicalSystems
hen = Systems.henon()
ts = trajectory(hen, 200000)
D_hen = generalized_dim(1, ts)
```

As a side note, be sure that you have enough data points, otherwise the values you will
get will never be correct, as is demonstrated by
J.-P. Eckmann and D. Ruelle (see Physica D **56**, pp 185-187 (1992)).

---


## Other Entropies and Dimensions

### Permutation Entropy
The permutation entropy is introduced by C. Bandt and B. Pompe as a
"A Natural Complexity Measure for Timeseries", which directly applies to arbitrary real-world data and is particularly useful in the presence of dynamical or observational noise.

```@docs
permentropy
```

For example, we will compute and compare the [`lyapunov`](@ref) exponent of the logistic
map with the order-6 permutation entropy, like in the original paper.
```@example entropy
using DynamicalSystems, PyPlot
ds = Systems.logistic()
rs = 3.5:0.001:4
ls = Float64[]; hs = Float64[]
for r in rs
    ds.p[1] = r
    push!(ls, lyapunov(ds, 100000))
    # For 1D systems `trajectory` returns a vector
    push!(hs, permentropy(trajectory(ds, 10000), 6))
end

f = figure(figsize = (10,6))
a1 = subplot(211)
plot(rs, ls); ylim(-2, log(2)); ylabel("\$\\lambda\$")
a1.axes.get_xaxis().set_ticklabels([])
xlim(rs[1], rs[end]);

a2 = subplot(212)
plot(rs, hs; color = "C1"); ylabel("\$h_6\$")
xlim(rs[1], rs[end]); xlabel("\$r\$")
tight_layout()
savefig("permentropy.png"); nothing # hide
```
![](permentropy.png)


!!! info "Permutation Entropy performance"
    Even though the current implementation is fine and runs reasonably fast for
    moderate orders, it can get slow for high orders. Issue [ChaosTools.jl#22](https://github.com/JuliaDynamics/ChaosTools.jl/issues/22)
    keeps track of this, and contains information on how to improve performance.


### Kaplan-Yorke Dimension
```@docs
kaplanyorke_dim
```
---
Notice that calling this function requires you to pass the Lyapunov exponents in an
ordered vector form (largest to smallest). Example:
```@example lyap
using DynamicalSystems
hen = Systems.henon()
D_kp = kaplanyorke_dim(lyapunovs(hen, 100000))
```
