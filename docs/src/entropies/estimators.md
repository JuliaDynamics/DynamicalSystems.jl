# Probabilities Estimators
In this page we list the various estimators (and further functions) that can be used to obtain probabilities representing a given dataset, or entropies directly. See [Entropies & Probabilities](@ref) for more.

## Visitation frequency (binning)

```@docs
VisitationFrequency
```

### Specifying binning/boxes

```@docs
RectangularBinning
```


## CountOccurrences (counting)

```@docs
CountOccurrences
```

## Kernel density

```@docs
NaiveKernel
```

### Distance evaluation methods

```@docs
TreeDistance
DirectDistance
```

### Example

Here, we draw some random points from a 2D normal distribution. Then, we use kernel
density estimation to associate a probability to each point `p`, measured by how many
points are within radius `1.5` of `p`. Plotting the actual points, along with their
associated probabilities estimated by the KDE procedure, we get the following surface
plot.

```@example
using Distributions, PyPlot, DelayEmbeddings, Entropies
𝒩 = MvNormal([1, -4], 2)
N = 500
D = Dataset(sort([rand(𝒩) for i = 1:N]))
x, y = columns(D)
p = probabilities(D, NaiveKernel(1.5))
figure()
surf(x, y, p.p)
xlabel("x"); ylabel("y")
savefig("kernel_surface.png"); nothing # hide
```

![](kernel_surface.png)

## Nearest neighbor estimators

### Kraskov

```@docs
Kraskov
```

### Kozachenko-Leonenko

```@docs
KozachenkoLeonenko
```

### Example

This example reproduces Figure in Charzyńska & Gambin (2016)[^Charzyńska2016]. Both
estimators nicely converge to the true entropy with increasing time series length.
For a uniform 1D distribution ``U(0, 1)``, the true entropy is `0` (red line).

```@example
using DynamicalSystems, PyPlot
import Distributions: Uniform, Normal

Ns = [100:100:500; 1000:1000:10000]
Ekl = Vector{Vector{Float64}}(undef, 0)
Ekr = Vector{Vector{Float64}}(undef, 0)

est_nn = KozachenkoLeonenko(w = 0)
# with k = 1, Kraskov is virtually identical to KozachenkoLeonenko, so pick a higher
# number of neighbors
est_knn = Kraskov(w = 0, k = 3)

nreps = 50
for N in Ns
    kl = Float64[]
    kr = Float64[]
    for i = 1:nreps
        pts = Dataset([rand(Uniform(0, 1), 1) for i = 1:N]);
        push!(kl, genentropy(pts, est_nn))
         # with k = 1 almost identical
        push!(kr, genentropy(pts, est_knn))
    end
    push!(Ekl, kl)
    push!(Ekr, kr)
end

# Plot
using PyPlot, StatsBase
f = figure(figsize = (5,6))
ax = subplot(211)
px = PyPlot.plot(Ns, mean.(Ekl); color = "C1", label = "KozachenkoLeonenko");
PyPlot.plot(Ns, mean.(Ekl) .+ StatsBase.std.(Ekl); color = "C1", label = "");
PyPlot.plot(Ns, mean.(Ekl) .- StatsBase.std.(Ekl); color = "C1", label = "");

xlabel("Time step"); ylabel("Entropy (nats)"); legend()
ay = subplot(212)
py = PyPlot.plot(Ns, mean.(Ekr); color = "C2", label = "Kraskov");
PyPlot.plot(Ns, mean.(Ekr) .+ StatsBase.std.(Ekr); color = "C2", label = "");
PyPlot.plot(Ns, mean.(Ekr) .- StatsBase.std.(Ekr); color = "C2", label = "");

xlabel("Time step"); ylabel("Entropy (nats)"); legend()
tight_layout()
PyPlot.savefig("nn_entropy_example.png"); nothing # hide
```

![](nn_entropy_example.png)

[^Charzyńska2016]: Charzyńska, A., & Gambin, A. (2016). Improvement of the k-NN entropy estimator with applications in systems biology. Entropy, 18(1), 13.

## Permutation (symbolic)

```@docs
SymbolicPermutation
```

### Example

This example reproduces an example from Bandt and Pompe (2002), where the permutation
entropy is compared with the largest Lyapunov exponents from time series of the chaotic
logistic map. Entropy estimates using [`SymbolicWeightedPermutation`](@ref)
and [`SymbolicAmplitudeAwarePermutation`](@ref) are added here for comparison.

```@example
using DynamicalSystems, PyPlot

ds = Systems.logistic()
rs = 3.4:0.001:4
N_lyap, N_ent = 100000, 10000
m, τ = 6, 1 # Symbol size/dimension and embedding lag

# Generate one time series for each value of the logistic parameter r
lyaps = Float64[]
hs_perm = Float64[]
hs_wtperm = Float64[]
hs_ampperm = Float64[]

base = Base.MathConstants.e
for r in rs
    ds.p[1] = r
    push!(lyaps, lyapunov(ds, N_lyap))

    x = trajectory(ds, N_ent) # time series
    hperm = Entropies.genentropy(x, SymbolicPermutation(m = m, τ = τ), base = base)
    hwtperm = Entropies.genentropy(x, SymbolicWeightedPermutation(m = m, τ = τ), base = base)
    hampperm = Entropies.genentropy(x, SymbolicAmplitudeAwarePermutation(m = m, τ = τ), base = base)

    push!(hs_perm, hperm); push!(hs_wtperm, hwtperm); push!(hs_ampperm, hampperm)
end

f = figure(figsize = (6, 8))
a1 = subplot(411)
plot(rs, lyaps); ylim(-2, log(2)); ylabel("\$\\lambda\$")
a1.axes.get_xaxis().set_ticklabels([])
xlim(rs[1], rs[end]);

a2 = subplot(412)
plot(rs, hs_perm; color = "C2"); xlim(rs[1], rs[end]);
xlabel(""); ylabel("\$h_6 (SP)\$")

a3 = subplot(413)
plot(rs, hs_wtperm; color = "C3"); xlim(rs[1], rs[end]);
xlabel(""); ylabel("\$h_6 (SWP)\$")

a4 = subplot(414)
plot(rs, hs_ampperm; color = "C4"); xlim(rs[1], rs[end]);
xlabel("\$r\$"); ylabel("\$h_6 (SAAP)\$")
tight_layout()
savefig("permentropy.png"); nothing # hide
```

![](permentropy.png)

## Time-scale (wavelet)

```@docs
TimeScaleMODWT
```

### Example

The scale-resolved wavelet entropy should be lower for very regular signals (most of the
energy is contained at one scale) and higher for very irregular signals (energy spread
more out across scales).

```@example
using DynamicalSystems, PyPlot
N, a = 1000, 10
t = LinRange(0, 2*a*π, N)

x = sin.(t);
y = sin.(t .+  cos.(t/0.5));
z = sin.(rand(1:15, N) ./ rand(1:10, N))

est = TimeScaleMODWT()
h_x = genentropy(x, est)
h_y = genentropy(y, est)
h_z = genentropy(z, est)

f = figure(figsize = (10,6))
ax = subplot(311)
px = plot(t, x; color = "C1", label = "h=$(h=round(h_x, sigdigits = 5))");
ylabel("x"); legend()
ay = subplot(312)
py = plot(t, y; color = "C2", label = "h=$(h=round(h_y, sigdigits = 5))");
ylabel("y"); legend()
az = subplot(313)
pz = plot(t, z; color = "C3", label = "h=$(h=round(h_z, sigdigits = 5))");
ylabel("z"); xlabel("Time"); legend()
tight_layout()
savefig("waveletentropy.png"); nothing # hide
```

![](waveletentropy.png)

## Utility methods

Some convenience functions for symbolization are provided.

```@docs
Entropies.encode_as_bin
Entropies.joint_visits
Entropies.marginal_visits
Entropies.symbolize
Entropies.encode_motif
```
