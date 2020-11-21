# Entropies & Probabilities

Here we discuss obtaining probabilities and entropies from a given dataset (that typically represents a trajectory or set in the state space of a dynamical system).
The data are expected in the form discussed in [Numerical Data](@ref).

The main **API** for this is contained in two functions:

* [`probabilities`](@ref) which computes probability distributions of given datasets
* [`genentropy`](@ref) which uses the output of [`probabilities`](@ref), or a set of pre-computed [`Probabilities`](@ref), to calculate entropies.

These functions dispatch on subtypes of `ProbabilitiesEstimator`, which are summarized in the [Probabilities Estimators](@ref) page.

## Probabilities

```@docs
Probabilities
probabilities
probabilities!
```

## Fast histograms

```@docs
binhist
```

## Generalized entropy

In the study of dynamical systems there are many quantities that identify as "entropy".
Notice that these quantities are not the [thermodynamic ones](https://en.wikipedia.org/wiki/Entropy), used in Statistical Physics.
Rather, they are more like the to the entropies of [information theory](https://en.wikipedia.org/wiki/Entropy_(information_theory)).

All of the entropy-related quantities boil down to one thing: first extracting probabilities from a dataset and then applying the generalized entropy formula using [`genentropy`](@ref).

```@docs
genentropy
permentropy
```
