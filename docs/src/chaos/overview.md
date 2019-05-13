# Features Overview
The features offered in this documentation section come from the package [ChaosTools.jl](https://github.com/JuliaDynamics/ChaosTools.jl). If you are encountering an issue with some of the methods, you can report/open a new issue at the GitHub Issues page.

### [Orbit Diagrams](orbitdiagram)

1. Orbit diagrams (aka bifurcation diagrams) of maps: [`orbitdiagram`](@ref).
2. Poincaré surfaces of section for continuous systems: [`poincaresos`](@ref).
3. Automated production of orbit diagrams for continuous systems: [`produce_orbitdiagram`](@ref).

### [Lyapunov Exponents](lyapunovs)

The following treat systems where the equations of motion are known:

1. Maximum Lyapunov exponent for both discrete and continuous systems: [`lyapunov`](@ref).
2. Lyapunov *spectrum* for both discrete and continuous systems: [`lyapunovs`](@ref).

### [Categorizing Chaos](chaos_detection)

1. The Generalized Alignment Index: $\text{GALI}_k$ : [`gali`](@ref).
    * Implemented for both discrete and continuous systems.
2. A test to categorize strong chaos, partially predictable chaos and regular behavior: [`predictability`](@ref).
    * Implemented for both discrete and continuous systems.

### [Entropies and Dimensions](entropies)

1. Generalized (Renyi) entropy: [`genentropy`](@ref).
2. Permutation entropy: [`permentropy`](@ref).
2. Fast and cheap (memory-wise) method for computing entropies of large datasets.
3. Generalized dimensions (e.g. capacity dimension, information dimension, etc.): [`generalized_dim`](@ref).
3. Kaplan-Yorke dimension: [`kaplanyorke_dim`](@ref).
4. Automated detection of best algorithmic parameters for calculating attractor dimensions.

And, in order to automatically deduce dimensions, we also offer methods for:

* Partitioning a function $y(x)$ vs. $x$ into regions where it is approximated by a straight line, using a flexible algorithm with a lot of control over the outcome. See [`linear_regions`](@ref).
* Detection of largest linear region of a function $y(x)$ vs. $x$ and extraction of the slope of this region.

### [Nonlinear Timeseries Analysis](nlts)

1. Broomhead-King coordinates: [`broomhead_king`](@ref).
3. Numerically determining the maximum Lyapunov exponent of a (e.g. experimentally) measured timeseries: [`numericallyapunov`](@ref).

### [Periodicity](periodicity)

1. Numerical method to find unstable and stable fixed points of *any order* $n$ of a discrete map (of any dimensionality): [`periodicorbits`](@ref).
    * Convenience functions for defining and realizing all possible combinations of $\mathbf{\Lambda}_k$ matrices required in the above method.
