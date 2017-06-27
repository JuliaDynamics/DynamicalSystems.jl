`DynamicalSystems.jl` is a Julia package for the exploration of continuous and discrete dynamical systems. it aims to be a useful and powerful companion for students and scientists treading
on the field of Chaos, nonlinear dynamics and dynamical systems in general. The package
treats discrete and continuous systems of the forms:
```math
\frac{d\vec{u}}{dt} = \vec{f}(\vec{u}) \quad \text{or}\quad \vec{x}_{n+1} = \vec{f}(\vec{x}_n)
```
as well as numerical data. More information is contained in the
[Systems Definition](system_definition) page.




This is the (non-final) list of what this package aims to offer:

1. Intuitive, consistent APIs for the definition of dynamical systems.
2. Automatic "completion" of the dynamics of the system with numerically computed
  Jacobians, in case they are not provided by the user.
3. Lyapunov exponent estimation.
4. Entropy estimation.
5. Attractor dimension estimation.
6. Entropy/Attractor dimension/Lyapunov exponents for *numerical data*.
7. Chaos control.
8. Other stuff I have not yet decided upon, since this is like a pre-alpha version.
8. Suggest or Contribute more stuff! (see contributors guide).
