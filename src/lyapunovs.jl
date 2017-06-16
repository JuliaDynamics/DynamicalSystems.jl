export lyapunovs, lyapunov
#######################################################################################
#                                      Discrete                                       #
#######################################################################################
"""
```julia
lyapunovs(ds::DynamicalSystem, N::Int; kwargs...)
```
Calculate the spectrum of lyapunov exponents of the system `ds` by applying the
QR-decomposition method [1] for `N` steps. Returns a vector with the *final*
values of the lyapunov exponents.
# Keywords:
* `Ttr` : Extra "transient" time to evolve the system if one is not sure it is
  already on the attractor. Should be `Int` for discrete systems. Defaults are
  system type dependent.
* `dt = 0.1` : (only for continuous) Length of time intervals between each application
  of the QR-decomposition algorithm.

[1] : Add reference here.
"""
function lyapunovs(ds::DiscreteDS, N::Real; Ttr::Int= 100)
  N = convert(Int, N)
  u = deepcopy(ds.state)
  dim = length(u)
  eom = ds.eom
  jac = ds.jacob
  # Transient iterations
  for i in 1:Ttr
    u = eom(u)
  end

  # Initialization
  λ = zeros(eltype(u), dim)
  J = jac(u)
  Q = @MMatrix eye(eltype(J), dim);
  # Main algorithm
  for i in 1:N
    u = eom(u)
    J = jac(u)
    D = J*Q

    F = qrfact(D)
    Q .= F[:Q]
    A_mul_signB!(Q, F[:R]) #ensure QR gives positive diagonal to R
    for i in 1:dim
      λ[i] += log(abs(F[:R][i,i]))
    end
  end
  λ./N
end# this works. add tests and good to go.

"""
```julia
lyapunov(ds::DynamicalSystem, Τ; kwargs...)
```
Calculate the maximum lyapunov exponent using a method due to Benettin [1], which simply
evolves two neighboring trajectories while constantly rescaling one of the two. `T` is
the total time of evolution (should be `Int` for discrete systems).

Of course, this method only works if one assumes that the system has at least
one positive lyapunov exponent. Use `lyapunovs` otherwise.

# Keywords:
* `Ttr` : Extra "transient" time to evolve the system if one is not sure it is
  already on the attractor. Should be `Int` for discrete systems. Defaults are
  system type dependent.
* `d0 = 1e-7` : Initial/rescaling distance between two neighboring trajectories.
* `threshold = 10^3*d0` : Threshold to rescale the test trajectory.
* `et = 1.0` : (only for continuous) Time of individual evolutions
  between sucessive checks for rescaling.

[1] : Benettin *et al.*, Phys. Rev. A **14**, pp 2338 (1976)
"""
function lyapunov(ds::DiscreteDS, N = 100000; Ttr::Int = 100,
  d0=1e-7*one(eltype(ds.state)), threshold=10^3*d0)

  N = convert(Int, N)
  threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))
  eom = ds.eom
  st1 = deepcopy(ds.state)

  # transient
  for i in 1:Ttr
    st1 = eom(st1)
  end

  st2 = st1 + d0
  dist = d0
  λ = zero(eltype(st1))
  i = 0
  while i < N
    #evolve until rescaling:
    while dist < threshold
      st1 = eom(st1)
      st2 = eom(st2)
      # if hasnan
      dist = norm(st1 - st2)
      i+=1
    end
    a = dist/d0
    λ += log(a)
    #rescale:
    st2 = st1 + (st2 - st1)/a
    dist = d0
  end
  λ /= i
end

function lyapunovs(ds::DiscreteDS1D, N=10000; Ttr = 100)

  N = convert(Int, N)
  eom = ds.eom
  der = ds.deriv
  x = deepcopy(ds.state)

  #transient
  for i in 1:Ttr
    x = eom(x)
  end

  λ = log(abs(der(x)))
  for i in 1:N
    x = eom(x)
    λ += log(abs(der(x)))
  end
  λ/N
end
lyapunov(ds::DiscreteDS1D, N::Int=10000; Ttr::Int = 100) = lyapunovs(ds, N, Ttr=Ttr)

#######################################################################################
#                                    Continuous                                       #
#######################################################################################
function lyapunov(ds::ContinuousDS, T = 10.0; Ttr::Int = 100,
  d0=1e-7*one(eltype(ds.state)), threshold=10^3*d0)

  N = convert(Int, N)
  threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))
  eom = ds.eom
  st1 = deepcopy(ds.state)

  # transient
  for i in 1:Ttr
    st1 = eom(st1)
  end

  st2 = st1 + d0
  dist = d0
  λ = zero(eltype(st1))
  i = 0
  while i < N
    #evolve until rescaling:
    while dist < threshold
      st1 = eom(st1)
      st2 = eom(st2)
      # if hasnan
      dist = norm(st1 - st2)
      i+=1
    end
    a = dist/d0
    λ += log(a)
    #rescale:
    st2 = st1 + (st2 - st1)/a
    dist = d0
  end
  λ /= i
end
