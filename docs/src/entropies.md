# Entropies and Dimensions

## Entropies
In the study of dynamical systems there are many quantities that identify as "entropy".
Notice that these quantities are not the more commonly known
[thermodynamic ones](https://en.wikipedia.org/wiki/Entropy), used in Statistical Physics. Rather, they are more like the to the entropies of [information theory](https://en.wikipedia.org/wiki/Entropy_(information_theory)), which represents
information contained within a dataset, or information about the dimensional
scaling of a dataset.

`DynamicalSystems.jl` defines a lot entropies, summarized in the following sections.

### Generalized Entropy & Co.
The generalized entropy is a concept mainly attributed to Rényi (see below).
```@docs
genentropy
```
---
Basically, given a [dataset](system_definition/#numerical-data) you can
partition it into boxes to calculate an entropy.

!!! tip "Worried about memory overflow? Don't be!"
    Partitioning the dataset (i.e. doing a *histogram*) is in general a costly
    operation that depends exponentially on the number of dimensions of the system.
    However, in this specific case the partition process has some special aspects
    that can be taken advantage
    of, reducing tremendously the memory allocation and spent time!

The function used internally is `non0hist`:
```@docs
non0hist
```
---
It typically outperforms traditional histograms
by **several orders of magnitude** in both memory and speed. You can compare
`DynamicalSystems.perform_non0hist` with `fit(Histogram, ...)` of [`StatsBase`](http://juliastats.github.io/StatsBase.jl/stable/)
for specific numbers on your machine.

For example, the Shannon entropy of a coin-flip process should be one bit,
[by definition](https://en.wikipedia.org/wiki/Shannon_(unit)). Let's see...
```julia
using DynamicalSystems
y = Int.(rand(Bool, 10000)) # just some coin tosses
sh = shannon(0.01, y)       # ≡ genentropy(1, 0.01, y)
isapprox(sh, log(2),  rtol = 1e-3) # true!
```
Because all entropies are calculated on base-e, the unit of measurement is "nat", and
one bit is log(2)×nat.

### Kolmogorov-Sinai Entropy
TBA.

## Attractor Dimension Estimation
There are numerous methods that one can use to calculate a so-called "dimension" of a
dataset, like for example the [Fractal dimension](https://en.wikipedia.org/wiki/Fractal_dimension). This real number can offer
a lot of information about the object that the dataset represents.

### Generalized Dimensions & Co.
Based on the definition of the [generalized entropy](entropies/#DynamicalSystems.genentropy), one can calculate an appropriate
dimension, called *generalized dimension*:
```@docs
generalized_dim
```
---
As stated clearly, this call performs a lot of automated steps. One is always better
better of performing the steps one by one to attain maximum control.

For example, we will calculate the dimensions of the strange attractors of the
[Hénon map](system_definition/#DynamicalSystems.Systems.henon) and the [Lorenz system](system_definition/#DynamicalSystems.Systems.lorenz):
```julia
using DynamicalSystems
hen = Systems.henon(-rand(2))
ts = timeseries(hen, 200000)
D_hen = information_dim(ts)

lor = Systems.lorenz(rand(3))
ts = timeseries(lor, 5000, dt = 0.05)
D_lor = capacity_dim(ts)
```
You will find that `D_hen` is around `1.2` and `D_lor` is around `1.95`, both of which
[are correct values](http://www.dt.fee.unicamp.br/~tiago/courses/dinamica_caotica/Lyapunov.pdf). As
a side note, be sure that you have enough data points, otherwise the values you will
get will never be correct, as is demonstrated by
J.-P. Eckmann and D. Ruelle (see Physica D **56**, pp 185-187 (1992)).

### Kaplan-Yorke Dimension
Also known as Lyapunov Dimension:
```@docs
kaplanyorke_dim
```
---
Notice that calling this function requires you to pass the lyapunov exponents in an
ordered vector form (largest to smallest). Example:
```julia
hen = Systems.henon()
D_kp = kaplanyorke_dim(lyapunovs(hen, 200000))
```
