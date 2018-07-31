# Entropies and Dimensions

## Entropies
In the study of dynamical systems there are many quantities that identify as "entropy".
Notice that these quantities are not the more commonly known
[thermodynamic ones](https://en.wikipedia.org/wiki/Entropy), used in Statistical Physics. Rather, they are more like the to the entropies of [information theory](https://en.wikipedia.org/wiki/Entropy_(information_theory)), which represents
information contained within a dataset, or information about the dimensional
scaling of a dataset.

### Generalized Entropy
```@docs
genentropy
```
---
Basically, given a [dataset](system_definition/#numerical-data) you can
partition it into boxes to calculate an entropy.

!!! tip "Worried about memory overflow? Don't be!"
    Partitioning the dataset (i.e. doing a *histogram*) is in general a costly
    operation that depends exponentially on the number of dimensions of the data
    and algebraically to the box size `ε`.
    However, in this specific case the partition process has some special aspects
    that can be taken advantage
    of, reducing tremendously the memory allocation and spent time!

    In fact, there is an upper bound to the memory allocated by `non0hist`: A constant
    multiplied by the length of the array, `N = length(p)`. No matter how small `ε` or how many dimensions the data has, the method can at most assign `N` dictionary entries.


The function used internally by `genentropy` is `non0hist`:
```@docs
non0hist
```
---
For example, the Shannon entropy of a coin-flip process should be one bit,
[by definition](https://en.wikipedia.org/wiki/Shannon_(unit)). Let's see...
```@example entropy
using DynamicalSystems
y = Float64.(rand(Bool, 1000000)) # just some coin tosses
sh = genentropy(1, 0.1, y)  # this is the shannon entropy
isapprox(sh, log(2),  rtol = 1e-4)
```
Because all entropies are by default calculated on base-$e$, the unit of measurement is "nat" and one bit is $\log(2)\times$nat.

### Permutation Entropy
The permutation entropy is introduced by C. Bandt and B. Pompe as a
"A Natural Complexity Measure for Timeseries", which directly applies to arbitrary real-world data and is particularly useful in the presence of dynamical or observational noise.

```@docs
permentropy
```

For example, we will compute and compare the [`lyapunov`](@ref) exponent of the logistic
map with the order-6 permutation entropy, like in the original paper.
```julia
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
a1[:axes][:get_xaxis]()[:set_ticklabels]([])
xlim(rs[1], rs[end]);

a2 = subplot(212)
plot(rs, hs; color = "C1"); ylabel("\$h_6\$")
xlim(rs[1], rs[end]); xlabel("\$r\$")
tight_layout()
```
![Permutation Entropy](https://i.imgur.com/tsqSA7a.png)

!!! info "Permutation Entropy performance"
    Even though the current implementation is fine and runs reasonably fast for
    moderate orders, it can get slow for high orders. Issue [ChaosTools.jl#22](https://github.com/JuliaDynamics/ChaosTools.jl/issues/22)
    keeps track of this, and contains information on how to improve performance.

## Attractor Dimension Estimation
There are numerous methods that one can use to calculate a so-called "dimension" of a
dataset, like for example the [Fractal dimension](https://en.wikipedia.org/wiki/Fractal_dimension). This real number can offer
a lot of information about the object that the dataset represents.

### Generalized Dimensions
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

#### Example
For example, the dimension of the strange attractor of the
[Hénon map](system_definition/#DynamicalSystems.Systems.henon) is:
```julia
using DynamicalSystems
hen = Systems.henon()
ts = trajectory(hen, 1000000)
D_hen = information_dim(ts)
```
```
1.2279316105815665
```

As a side note, be sure that you have enough data points, otherwise the values you will
get will never be correct, as is demonstrated by
J.-P. Eckmann and D. Ruelle (see Physica D **56**, pp 185-187 (1992)).

---

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
