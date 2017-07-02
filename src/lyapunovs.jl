export lyapunovs, lyapunov
#######################################################################################
#                                      Discrete                                       #
#######################################################################################
"""
```julia
lyapunovs(ds::DynamicalSystem, N; kwargs...) -> [λ1, λ2, ..., λD]
```
Calculate the spectrum of lyapunov exponents of the system `ds` by applying the
QR-decomposition method [1] `N` times. Returns a vector with the *final*
values of the lyapunov exponents.
# Keyword Arguments:
* `Ttr` : Extra "transient" time to evolve the system before application of the
  algorithm. Should be `Int` for discrete systems. Defaults are
  system type dependent.
* `dt = 1.0` : (only for continuous) Time of individual evolutions
  between sucessive orthonormalization steps.
* `diff_eq_kwargs = Dict()` : (only for continuous)
  Keyword arguments passed into the solvers of the
  `DifferentialEquations` package (see `evolve` or `timeseries` for more info).

[1] : K. Geist *et al*, Progr. Theor. Phys. **83**, pp 875 (1990)
"""
function lyapunovs(ds::DiscreteDS, N::Real; Ttr::Int= 100)
  N = convert(Int, N)
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
Calculate the maximum lyapunov exponent using a method due to Benettin [1], which simply
evolves two neighboring trajectories while constantly rescaling one of the two. `T` is
an optional argument which denotes
the total time of evolution (should be `Int` for discrete systems).

# Keyword Arguments:
* `Ttr` : Extra "transient" time to evolve the system before application of the
  algorithm. Should be `Int` for discrete systems. Defaults are
  system type dependent.
* `d0 = 1e-7` : Initial & rescaling distance between two neighboring trajectories.
* `threshold = 10^3*d0` : Threshold to rescale the test trajectory.
* `diff_eq_kwargs = Dict()` : (only for continuous)
  Keyword arguments passed into the solvers of the
  `DifferentialEquations` package (see `evolve` or `timeseries` for more info).
* `dt = 10.0` : (only for continuous) Time between each check of distance exceeding
  the `threshold`.

*Warning*: Default values have been choosen to give accurate & fast results for
maximum lyapunov exponent expected between 0.1 to 1.0. Be sure to adjust
them properly for your system.

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
      i>=N && break
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
function lyapunov(ds::ContinuousDS, T = 10000.0; Ttr = 10.0,
  d0=1e-7*one(eltype(ds.state)), threshold=10^3*d0, dt = 10.0,
  diff_eq_kwargs = Dict())

  threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))

  # transient:
  ds = evolve(ds, Ttr, diff_eq_kwargs)
  # initialize:
  st1 = ds.state
  st2 = st1 + d0
  prob1 = ODEProblem(ds, dt)
  prob2 = ODEProblem(setu(st2, ds), dt)
  dist = d0
  λ = zero(eltype(st1))
  t = zero(T)
  # start evolution and rescaling:
  while t < T
    #evolve until rescaling:
    while dist < threshold
      # evolve
      st1 = evolve!(prob1, diff_eq_kwargs)
      st2 = evolve!(prob2, diff_eq_kwargs)
      dist = norm(st1 - st2)
      t += dt
      t >= T && break
    end
    a = dist/d0
    λ += log(a)
    #rescale:
    st2 = st1 + (st2 - st1)/a
    prob2.u0 = st2
    dist = d0
  end
  λ /= t
end


function lyapunovs(ds::ContinuousDS, N::Int;
  Ttr::Real = 1.0, diff_eq_kwargs::Dict = Dict(), dt::Real = 1.0)

# Ttr = 1.0
# N = 1000
# et = 1.0
# diff_eq_kwargs = Dict()
# ds = Systems.lorenz()


D = dimension(ds)
T = eltype(ds.state)
# Transient iterations
ds = evolve(ds, Ttr, diff_eq_kwargs)

# Initialization
λ::Vector{T} = zeros(T, D)
tbprob = tangentbundle_setup(ds, dt)
Q = @MMatrix eye(eltype(ds.state), D)
S = tbprob.u0
# Main algorithm
for i in 1:N
  tbprob.u0[:, 1] .= view(S, :, 1)
  tbprob.u0[:, 2:D+1] .= Q

  S = tangentbundle_evolve(tbprob, diff_eq_kwargs)
  G = view(S, :, 2:D+1)*Q
  F = qrfact(G)
  Q .= F[:Q]
  A_mul_signB!(Q, F[:R]) #ensure QR gives positive diagonal to R
  for j in 1:D
    λ[j] += log(abs(F[:R][j,j]))
  end
  # Set initial conditions for next iteration:

end
λ./(N*dt)

end# this works. add tests and good to go.
