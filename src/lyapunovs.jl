export lyapunovs, lyapunov
#####################################################################################
#                                    Discrete                                       #
#####################################################################################
"""
```julia
lyapunovs(ds::DynamicalSystem, N; kwargs...) -> [λ1, λ2, ..., λD]
```
Calculate the spectrum of lyapunov [1] exponents of `ds` by applying the
QR-decomposition method `N` times (see method "H2" of [2], or directly the original
paper(s) [3]).
Returns a vector with the *final*
values of the lyapunov exponents in descending order.
## Keyword Arguments:
* `Ttr` : Extra "transient" time to evolve the system before application of the
  algorithm. Should be `Int` for discrete systems. Defaults are
  system type dependent.
* `dt = 1.0` : (only for continuous) Time of individual evolutions
  between sucessive orthonormalization steps.
* `diff_eq_kwargs = Dict()` : (only for continuous)
  Keyword arguments passed into the solvers of the
  `DifferentialEquations` package (see `evolve` or `timeseries` for more info).

[1] : A. M. Lyapunov, *The General Problem of the Stability of Motion*,
Taylor & Francis (1992)

[2] : K. Geist *et al.*, Progr. Theor. Phys. **83**, pp 875 (1990)

[3] : G. Benettin *et al.*, Meccanica **15**, pp 9-20 & 21-30 (1980)
"""
function lyapunovs(ds::DiscreteDS, N::Real; Ttr::Real = 100)

  u = deepcopy(ds.state)
  D = length(u)
  eom = ds.eom
  jac = ds.jacob
  # Transient iterations
  for i in 1:Ttr
    u = eom(u)
  end

  # Initialization
  λ = zeros(eltype(u), D)
  Q = eye(eltype(u), D)
  K = copy(Q)
  # Main algorithm
  for i in 1:N
    u = eom(u)
    A_mul_B!(K, jac(u), Q)

    Q, R = qr(K)
    for i in 1:D
      λ[i] += log(abs(R[i,i]))
    end
  end
  λ./N
end# this works. add tests and good to go.

"""
```julia
lyapunov(ds::DynamicalSystem, Τ; kwargs...) -> λ
```
Calculate the maximum lyapunov exponent `λ` using a method due to Benettin [1],
which simply
evolves two neighboring trajectories while constantly rescaling one of the two.
`T`  denotes the total time of evolution (should be `Int` for discrete systems).

## Keyword Arguments:
* `Ttr` : Extra "transient" time to evolve the system before application of the
  algorithm. Should be `Int` for discrete systems. Defaults are
  system type dependent.
* `d0 = 1e-9` : Initial & rescaling distance between two neighboring trajectories.
* `threshold = 10^3*d0` : Threshold to rescale the test trajectory.
* `diff_eq_kwargs = Dict()` : (only for continuous)
  Keyword arguments passed into the solvers of the
  `DifferentialEquations` package (see `evolve` or `timeseries` for more info).
* `dt = 0.1` : (only for continuous) Time of evolution between each check of
  distance exceeding the `threshold`.

*Warning*: Default values have been choosen to give accurate & fast results for
maximum lyapunov exponent expected between 0.1 to 1.0. Be sure to adjust
them properly for your system. Specifically for the continuous systems,
be sure that `exp(λ*dt) < threshold/d0`.

[1] : G. Benettin *et al.*, Phys. Rev. A **14**, pp 2338 (1976)
"""
function lyapunov(ds::DiscreteDS, N::Real = 100000; Ttr::Int = 100,
  d0=1e-9*one(eltype(ds.state)), threshold=10^4*d0)

  threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))
  eom = ds.eom
  st1 = deepcopy(ds.state)

  # transient system evolution
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
      dist = norm(st1 - st2)
      i+=1
      i>=N && break # this line is nessesary for safety! (if systems never go apart)
    end
    # local lyapunov exponent is simply the relative distance of the trajectories
    a = dist/d0
    λ += log(a)
    #rescale:
    st2 = st1 + (st2 - st1)/a
    dist = d0
  end
  λ /= i
end

function lyapunovs(ds::DiscreteDS1D, N::Real = 10000; Ttr = 100)

  eom = ds.eom
  der = ds.deriv
  x = deepcopy(ds.state)

  #transient system evolution
  for i in 1:Ttr
    x = eom(x)
  end

  # The case for 1D systems is trivial: you add log(abs(der(x))) at each step
  λ = log(abs(der(x)))
  for i in 1:N
    x = eom(x)
    λ += log(abs(der(x)))
  end
  λ/N
end
lyapunov(ds::DiscreteDS1D, N::Int=10000; Ttr::Int = 100) = lyapunovs(ds, N, Ttr=Ttr)

#####################################################################################
#                                  Continuous                                       #
#####################################################################################
function lyapunov(ds::ContinuousDS, T::Real = 10000.0; Ttr = 10.0,
  d0=1e-9*one(eltype(ds.state)), threshold=10^4*d0, dt = 0.1,
  diff_eq_kwargs = Dict())

  const dict = Dict()

  threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))

  # Transient system evolution
  if Ttr != 0
    evolve!(ds, Ttr; diff_eq_kwargs = diff_eq_kwargs)
  end
  # initialize:
  st1 = ds.state
  st2 = st1 + d0
  prob1 = ODEProblem(ds, dt)
  prob2 = ODEProblem(ds, dt)
  prob2.u0 = st2
  dist = d0
  λ = zero(eltype(st1))
  t = zero(T)
  i = 0
  # start evolution and rescaling:
  while t < T
    #evolve until rescaling:
    while dist < threshold
      # evolve
      if diff_eq_kwargs == dict
        st1 = evolve!(prob1)
        st2 = evolve!(prob2)
      else
        st1 = evolve!(prob1, diff_eq_kwargs)
        st2 = evolve!(prob2, diff_eq_kwargs)
      end
      dist = norm(st1 - st2)
      t += dt
      i += 1
      t ≥ T && break #once again, necessary for systems that don't go apart
    end
    # add computed scale to accumulator (scale = local lyaponov exponent):
    a = dist/d0
    # Warning message for bad decision of `thershold` or `d0`:
    if a > 1e^4 && i <= 1
      warnstr = "Distance between test and original trajectory exceeded threshold "
      warnstr*= "after just 1 evolution step. "
      warnstr*= "Please decrease `dt`, increase `threshold` or decrease `d0`."
      warn(warnstr)
      errorstr = "Parameters choosen for `lyapunov` with "
      errorstr*= "`ContinuousDS` are not fitting the algorithm."
      throw(ArgumentError(errorstr))
    end
    λ += log(a)
    # Rescale and reset everything:
    prob2.u0 = st1 + (st2 - st1)/a
    dist = d0; i = 0
  end
  λ /= t
end


function lyapunovs(ds::ContinuousDS, N::Real;
  Ttr::Real = 1.0, diff_eq_kwargs::Dict = Dict(), dt::Real = 1.0)

  D = dimension(ds)
  T = eltype(ds.state)
  jac = ds.jacob
  const dict = Dict()
  # Transient
  if Ttr != 0
    evolve!(ds, Ttr; diff_eq_kwargs = diff_eq_kwargs)
  end
  # Initialization (uses tangent dynamics)
  tbprob = tangentbundle_setup(ds, dt)
  λ = zeros(T, D)
  Q = eye(eltype(ds.state), D)

  # Main algorithm
  for i in 1:N
    tbprob.u0[:, 2:end] .= Q # update tangent dynamics state
    # evolve system and tangent dynamics:
    if diff_eq_kwargs == dict
      evolve!(tbprob)
    else
      evolve!(tbprob, diff_eq_kwargs)
    end

    # Perform QR (on the tangent flow):
    Q, R = qr(view(tbprob.u0, :, 2:D+1))
    # Add correct (positive) numbers to Lyapunov spectrum
    for j in 1:D
      λ[j] += log(abs(R[j,j]))
    end
  end
  λ./(N*dt) #return spectrum
end
