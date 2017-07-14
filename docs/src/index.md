`DynamicalSystems.jl` is a Julia package for the exploration of continuous and discrete dynamical systems. it aims to be a useful and powerful companion for students and scientists treading
on the field of Chaos, nonlinear dynamics and dynamical systems in general. The package
treats discrete and continuous systems of the forms:
```math
\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}) \quad \text{or}\quad \vec{x}_{n+1} = \vec{f}(\vec{x}_n)
```
as well as numerical data.

One of a major goals of this package is to be completely transparent as to what is
going on "under the hood". In scientific research, you never want to use *black boxes*,
e.g. functions that give a result without telling you how it was calculated. To this end,
**almost every documentation string** gives
direct references to the original papers where the algorithm is taken from, in case you don't understand (or simply don't want to read) the source code. For example,
the documentation string of [`?lyapunovs`](https://datseris.github.io/DynamicalSystems.jl/latest/lyapunovs/#DynamicalSystems.lyapunovs) will cite:
```
[1] : A.M. Lyapunov, The General Problem of the Stability of Motion, Taylor & Francis (1992)

[2] : K. Geist et al, Progr. Theor. Phys. 83, pp 875 (1990)
```
where the first reference introduces the concept of Lyapunov exponents and the second
contains a plethora of methods for their computation (including the method used in
the aforementioned function).

## Contents
This is the (non-final) list of what this package aims to offer:

1. Intuitive, consistent APIs for the definition of dynamical systems, see [System Definition](system_definition).
2. Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user, see [System Definition](system_definition).
3. Lyapunov exponent estimation, see [Lyapunov Exponents](lyapunovs).
4. Entropy estimation, see [Entropies and Dimensions](entropies).
5. Attractor dimension estimation, see [Entropies and Dimensions](entropies).
6. Entropy/Attractor dimension/Lyapunov exponents for *numerical data*, see [Lyapunov Exponents](lyapunovs) and [Entropies and Dimensions](entropies).
6. Attractor reconstruction, embedding and all that jazz.
7. Numeric Computation of Kolmogorov-Sinai entropy.
8. Definition of chaos, by Ott.
7. Chaos control, TBA.
8. Other stuff I have not yet decided upon, since this is like a pre-alpha version.
8. Suggest or Contribute more stuff! (see [Contributor Guide](contributors_guide)).
