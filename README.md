![DynamicalSystems.jl logo: The Double Pendulum](https://i.imgur.com/nFQFdB0.gif)

A Julia package for the exploration of continuous and discrete dynamical systems.

| **Documentation**   |  **Travis**     | **AppVeyor** | Gitter |
|:--------:|:-------------------:|:-----:|:-----:|
|[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaDynamics.github.io/DynamicalSystems.jl/latest) | [![Build Status](https://travis-ci.org/JuliaDynamics/DynamicalSystems.jl.svg?branch=master)](https://travis-ci.org/JuliaDynamics/DynamicalSystems.jl) | [![Build status](https://ci.appveyor.com/api/projects/status/iqss550vb2dik7b2?svg=true)](https://ci.appveyor.com/project/JuliaDynamics/dynamicalsystems-jl) | [![Gitter](https://img.shields.io/gitter/room/nwjs/nw.js.svg)](https://gitter.im/JuliaDynamics/Lobby)

`DynamicalSystems.jl` aims to be a useful companion for students and scientists treading
on the field of Chaos, nonlinear dynamics and dynamical systems in general.
See the [documentation](https://JuliaDynamics.github.io/DynamicalSystems.jl/latest) for more.

## Installation
This package is registered. Simply use `Pkg.add("DynamicalSystems")` to install it.

Bug-fixes and upgrades are constantly fed to the master branch and for this it is
also advised to use `Pkg.checkout("DynamicalSystems")` after installing. This ensures
that you get all the bug-fixes without having to wait for a formal release tag.

## Contents
This is the (non-final) list of what this package offers:

1. Intuitive, consistent APIs for the definition of general dynamical systems.
2. Automatic "completion" of the dynamics of the system with numerically computed Jacobians, in case they are not provided by the user.
3. Lyapunov exponent estimation.
4. Entropy and Attractor Dimension estimation.
6. Attractor reconstruction, delay coordinates embedding and all that jazz.
6. Entropy/Attractor dimension/Lyapunov exponents for *numerical data*.
7. Finding unstable periodic orbits of any period of Discrete maps.

The following are examples of "wanted features", that we would like to
have in our package and may (*or may not*) be implemented soon.

* Numeric Computation of Kolmogorov-Sinai entropy.
* Definition of chaos, by Ott.
* GALI method for estimating chaotic behavior, TBA.

For more, see the The [wanted features GitHub page](https://github.com/JuliaDynamics/DynamicalSystems.jl/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Awanted_feature) GitHub page.

## Double pendulum video:
Checkout this *amazing* video by Cormullion, featuring the `Systems.double_pendulum()`
from our package:

[![Double Pendulums Video](http://img.youtube.com/vi/vLDpLxU2fEg/0.jpg)](
https://www.youtube.com/watch?v=vLDpLxU2fEg)
