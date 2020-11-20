# Entropies & Probabilities

Here we discuss the main API for obtaining probabilities from a given dataset (that typically represents a trajectory or set in the state space of a dynamical system).
The data are expected in the form discussed in [Numerical Data](@ref).


The main **API** for this is contained in two functions:

* [`probabilities`](@ref) which computes probability distributions of given datasets
* [`genentropy`](@ref) which uses the output of [`probabilities`](@ref), or a set of
    pre-computed [`Probabilities`](@ref), to calculate entropies.

These functions dispatch on subtypes of [`ProbabilitiesEstimator`](@ref), which are summarized in the [Probabilities Estimators](@ref) page.

## Probabilities

```@docs
Probabilities
probabilities
probabilities!
ProbabilitiesEstimator
```

## Fast histograms

```@docs
binhist
```

## Generalized entropy

```@docs
genentropy
```
