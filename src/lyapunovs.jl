export lyapunovs, lyapunov, lyapunov_full
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
* `Ttr = 0` : Extra "transient" time to evolve the system before application of the
  algorithm. Should be `Int` for discrete systems.
* `dt = 1.0` : (only for continuous) Time of individual evolutions
  between sucessive orthonormalization steps.
* `diff_eq_kwargs = Dict()` : (only for continuous)
  Keyword arguments passed into the solvers of the
  `DifferentialEquations` package (see `timeseries` for more info).

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
  Q = @SMatrix eye(eltype(u), D)
  K = copy(Q)
  # Main algorithm
  for i in 1:N
    u = eom(u)
    K = jac(u)*Q

    Q, R = qr_sq(K)
    for i in 1:D
      λ[i] += log(abs(R[i, i]))
    end
  end
  λ./N
end

"""
```julia
lyapunov(ds::DynamicalSystem, Τ; kwargs...) -> λ
```
Calculate the maximum lyapunov exponent `λ` using a method due to Benettin [1],
which simply
evolves two neighboring trajectories while constantly rescaling one of the two.
`T`  denotes the total time of evolution (should be `Int` for discrete systems).

## Keyword Arguments:
* `Ttr = 0` : Extra "transient" time to evolve the system before application of the
  algorithm. Should be `Int` for discrete systems.
* `d0 = 1e-9` : Initial & rescaling distance between two neighboring trajectories.
* `threshold = 10^3*d0` : Threshold to rescale the test trajectory.
* `diff_eq_kwargs = Dict()` : (only for continuous)
  Keyword arguments passed into the solvers of the
  `DifferentialEquations` package (`timeseries` for more info).
* `dt = 0.1` : (only for continuous) Time of evolution between each check of
  distance exceeding the `threshold`.
* `displacement = .+` : (only for continuous) Function of the form `f(state, d0)`
  which given a `state` and `d0` produces a new state that is (approximately) `d0`
  far away from the given `state` but in the way the user wants. To be used for
  Hamiltonian systems, where the test trajectory must be choosen such that it
  would have the same energy as the reference trajectory.

This function returns only the final value of the computation of the maximum lyapunov
exponent. Use the function `lyapunov_full` to get a vector of the convergence
of the computation versus time.

[1] : G. Benettin *et al.*, Phys. Rev. A **14**, pp 2338 (1976)
"""
function lyapunov(ds::DiscreteDS, N::Real = 100000; Ttr::Int = 100,
  d0=1e-9, threshold=10^3*d0)

  threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))
  eom = ds.eom
  st1 = deepcopy(ds.state)

  # transient system evolution
  for i in 1:Ttr
    st1 = eom(st1)
  end

  st2 = st1 + d0
  dist = d0*one(eltype(ds.state))
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


"""
```julia
lyapunov_full(ds::DynamicalSystem, Τ; kwargs...) -> λ_ts, ts
```
Do the same as `lyapunov` but return the vector of convergence of the computation
`λ_ts` versus time `ts`.
"""
function lyapunov_full(ds::DiscreteDS, N::Real = 100000; Ttr::Int = 100,
  d0=1e-9, threshold=10^3*d0)

  threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))
  eom = ds.eom
  st1 = deepcopy(ds.state)

  # transient system evolution
  for i in 1:Ttr
    st1 = eom(st1)
  end

  st2 = st1 + d0
  dist = d0*one(eltype(ds.state))
  λ = zero(eltype(st1))
  λ_ts = Vector{eltype(st1)}(0)   # the timeseries for the Lyapunov exponent
  ts = Vector{Int}(0)
  i::Int = 0
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
    push!(λ_ts, λ/i)
    push!(ts, i)
    #rescale:
    st2 = st1 + (st2 - st1)/a
    dist = d0
  end
  return λ_ts, ts
end

function lyapunovs(ds::DiscreteDS1D, N::Real = 10000; Ttr::Int = 100)

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
#                              Lyapunov Helpers                                     #
#####################################################################################
function tangentbundle_setup_integrator(ds::ContinuousDynamicalSystem, t_final;
  diff_eq_kwargs=Dict())

  D = dimension(ds)
  f! = ds.eom!
  jac = ds.jacob

  # the equations of motion `tbeom!` evolve the system and the tangent dynamics
  # The e.o.m. for the system is f!(t, u , du).
  # The e.o.m. for the tangent dynamics is simply:
  # dY/dt = J(u) ⋅ Y
  # with J the Jacobian of the system (NOT the flow), at the current state
  tbeom! = (t, u, du) -> begin
    f!(view(du, :, 1), u)
    A_mul_B!(
    view(du, :, 2:D+1),
    jac(view(u, :, 1)),
    view(u, :, 2:D+1)
    )
  end

  # S is the matrix that keeps the system state in the first column
  # and tangent dynamics (Jacobian of the Flow) in the rest of the columns
  S = [ds.state eye(eltype(ds.state), D)]

  tbprob = ODEProblem(tbeom!, S, (zero(t_final), t_final))
  if haskey(diff_eq_kwargs, :solver)
    solver = diff_eq_kwargs[:solver]
    pop!(diff_eq_kwargs, :solver)
    tb_integ = init(tbprob, solver; diff_eq_kwargs..., save_everystep=false)
  else
    tb_integ = init(tbprob, Tsit5(); diff_eq_kwargs..., save_everystep=false)
  end
  return tb_integ
end

function check_tolerances(d0, dek)
  defatol = 1e-6; defrtol = 1e-3
  atol = haskey(dek, :abstol) ? dek[:abstol] : defatol
  rtol = haskey(dek, :reltol) ? dek[:reltol] : defrtol
  if atol > 10d0
    warn("Absolute tolerance (abstol) of integration is much bigger than `d0`.")
  end
  if rtol > 10d0
    warn("Relative tolerance (reltol) of integration is much bigger than `d0`.")
  end
end

#####################################################################################
#                            Continuous Lyapunovs                                   #
#####################################################################################
function lyapunov(ds::ContinuousDynamicalSystem, T = 10000.0; Ttr = 0.0,
  d0=1e-9, threshold=10^3*d0, dt = 0.1,
  diff_eq_kwargs = Dict(:abstol=>d0, :reltol=>d0), displacement = .+)

  check_tolerances(d0, diff_eq_kwargs)
  D = dimension(ds)
  T = convert(eltype(ds.state), T)
  threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))

  # Transient system evolution
  Ttr != 0 && evolve!(ds, Ttr; diff_eq_kwargs = diff_eq_kwargs)

  # initialize:
  st1 = ds.state
  integ1 = ODEIntegrator(ds, T; diff_eq_kwargs=diff_eq_kwargs)
  integ1.opts.advance_to_tstop=true
  ds.state .= displacement(st1, d0)
  integ2 = ODEIntegrator(ds, T; diff_eq_kwargs=diff_eq_kwargs)
  integ2.opts.advance_to_tstop=true
  ds.state = st1
  return lyapunov(integ1, integ2, T; d0=d0, threshold=threshold, dt=dt)
end

function lyapunov(integ1::ODEIntegrator, integ2::ODEIntegrator, T::Real;
  d0=1e-9, threshold=10^3*d0, dt = 0.1)

  dist = d0*one(eltype(integ1.u))
  λ = zero(eltype(integ1.u))
  i = 0;
  tvector = dt:dt:T

  # start evolution and rescaling:
  for τ in tvector
    # evolve until rescaling:
    push!(integ1.opts.tstops, τ); step!(integ1)
    push!(integ2.opts.tstops, τ); step!(integ2)
    dist = norm(integ1.u .- integ2.u)
    i += 1
    # Rescale:
    if dist ≥ threshold
      # add computed scale to accumulator (scale = local lyaponov exponent):
      a = dist/d0
      # Warning message for bad decision of `thershold` or `d0`:
      if a > threshold/d0 && i ≤ 1
        warnstr = "Distance between test and original trajectory exceeded threshold "
        warnstr*= "after just 1 evolution step. "
        warnstr*= "Please decrease `dt`, increase `threshold` or decrease `d0`."
        warn(warnstr)
        errorstr = "Parameters choosen for `lyapunov` with "
        errorstr*= "`ContinuousDynamicalSystem` are not fitted for the algorithm."
        throw(ArgumentError(errorstr))
      end
      λ += log(a)
      # Rescale and reset everything:
      integ2.u = integ1.u .+ (integ2.u .- integ1.u)./a
      u_modified!(integ2, true)
      set_proposed_dt!(integ2, integ1)
      dist = d0; i = 0
    end
  end
  λ /= T
end

"""
Compute the timeseries of the maximum Lyapunov exponent. This function
should be used to check the convergence of the series.
"""
function lyapunov_full(ds::ContinuousDynamicalSystem, T = 10000.0; Ttr = 0.0,
  d0=1e-9, threshold=10^3*d0, dt = 0.1, displacement = .+,
  diff_eq_kwargs = Dict(:abstol=>d0, :reltol=>d0))

  check_tolerances(d0, diff_eq_kwargs)

  T = convert(eltype(ds.state), T)
  threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))

  # Transient system evolution
  Ttr != 0 && evolve!(ds, Ttr; diff_eq_kwargs = diff_eq_kwargs)

  # initialize:
  st1 = ds.state
  integ1 = ODEIntegrator(ds, T; diff_eq_kwargs=diff_eq_kwargs)
  integ1.opts.advance_to_tstop=true
  ds.state .= displacement(st1, d0)
  integ2 = ODEIntegrator(ds, T; diff_eq_kwargs=diff_eq_kwargs)
  integ2.opts.advance_to_tstop=true
  ds.state = st1
  return lyapunov_full(integ1, integ2, T; d0=d0, threshold=threshold, dt=dt)
end

function lyapunov_full(integ1::ODEIntegrator, integ2::ODEIntegrator, T::Real;
  d0=1e-9, threshold=10^3*d0, dt = 0.1)


  dist = d0*one(eltype(integ1.u))
  λ = zero(eltype(integ1.u))
  λ_ts = Array{eltype(integ1.u)}(0)   # the timeseries for the Lyapunov exponent
  ts = Array{eltype(T)}(0)            # the time points of the timeseries
  i = 0;
  tvector = dt:dt:T

  # start evolution and rescaling:
  for τ in tvector
    # evolve until rescaling:
    push!(integ1.opts.tstops, τ); step!(integ1)
    push!(integ2.opts.tstops, τ); step!(integ2)
    dist = norm(integ1.u .- integ2.u)
    i += 1
    # Rescale:
    if dist ≥ threshold
      # add computed scale to accumulator (scale = local lyaponov exponent):
      a = dist/d0
      # Warning message for bad decision of `thershold` or `d0`:
      if a > threshold/d0 && i ≤ 1
        warnstr = "Distance between test and original trajectory exceeded threshold "
        warnstr*= "after just 1 evolution step. "
        warnstr*= "Please decrease `dt`, increase `threshold` or decrease `d0`."
        warn(warnstr)
        errorstr = "Parameters choosen for `lyapunov` with "
        errorstr*= "`ContinuousDynamicalSystem` are not fitted for the algorithm."
        throw(ArgumentError(errorstr))
      end
      λ += log(a)
      push!(λ_ts, λ/τ)
      push!(ts, τ)
      # Rescale and reset everything:
      integ2.u = @. integ1.u + (integ2.u - integ1.u)/a
      u_modified!(integ2, true)
      set_proposed_dt!(integ2, integ1)
      dist = d0; i = 0
    end
  end
  λ_ts, ts
end


function lyapunovs(ds::ContinuousDynamicalSystem, N::Real=1000;
  Ttr::Real = 0.0, diff_eq_kwargs::Dict = Dict(), dt::Real = 0.1)

  tstops = dt:dt:N*dt
  D = dimension(ds)
  λ = zeros(eltype(ds.state), D)
  Q = eye(eltype(ds.state), D)

  # Transient evolution:
  Ttr != 0 && evolve!(ds, Ttr; diff_eq_kwargs = diff_eq_kwargs)

  # Create integrator for dynamics and tangent space:
  integ = tangentbundle_setup_integrator(
  ds, tstops[end]; diff_eq_kwargs = diff_eq_kwargs)
  integ.opts.advance_to_tstop=true

  # Main algorithm
  for τ in tstops
    integ.u[:, 2:end] .= Q # update tangent dynamics state (super important!)
    push!(integ.opts.tstops, τ)
    step!(integ)

    # Perform QR (on the tangent flow):
    Q, R = qr_sq(view(integ.u, :, 2:D+1))
    # Add correct (positive) numbers to Lyapunov spectrum
    for j in 1:D
      λ[j] += log(abs(R[j,j]))
    end
  end
  λ./(N*dt) #return spectrum
end
