export lyapunovs, lyapunov
#####################################################################################
#                                    Discrete                                       #
#####################################################################################
"""
```julia
lyapunovs(ds::DynamicalSystem, N; kwargs...) -> [λ1, λ2, ..., λD]
```
Calculate the spectrum of Lyapunov exponents [1] of `ds` by applying the
QR-decomposition method `N` times (see method "H2" of [2], or directly the original
paper(s) [3]).
Returns a vector with the *final*
values of the lyapunov exponents in descending order.
### Keyword Arguments:
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
lyapunov(ds::DynamicalSystem, Τ, ret_conv = Val{false}; kwargs...)
```
Calculate the maximum Lyapunov exponent `λ` using a method due to Benettin [1],
which simply
evolves two neighboring trajectories (one called "given" and one called "test")
while constantly rescaling the test one.
`T`  denotes the total time of evolution (should be `Int` for discrete systems). The
Lyapunov exponent is the average of the time-local Lyapunov exponents,
meaning
```math
\\lambda = \\frac{1}{t_{n}}\\sum_{i=1}^{n}
\\ln\\left( a_i \\right),\\quad a_i = \\frac{d(t_{i})}{d_0}
```
with ``d(t_i)`` the distance between test and given trajectory at
time ``t_i``, which is also the time of the ``i``-th rescaling.

If `ret_conv` is `Val{true}` return the convergence timeseries of the Lyapunov
exponent
`λts` as well as the corresponding time vector `ts`. If `ret_conv` is `Val{false}`
(default) return the converged Lyapunov value `λts[end]` instead.

### Keyword Arguments:

  * `Ttr = 0` : Extra "transient" time to evolve the system before application of the
    algorithm. Should be `Int` for discrete systems.
  * `d0 = 1e-9` : Initial & rescaling distance between two neighboring trajectories.
  * `threshold = 1e-5` : Whenever the distance between given and test trajectory
    exceeds `threshold`, the test
    trajectory is rescaled back to having distance `d0` from the given trajectory.
  * `diff_eq_kwargs = Dict(:abstol=>d0, :reltol=>d0)` : (only for continuous)
    Keyword arguments passed into the solvers of the
    `DifferentialEquations` package (see `timeseries` for more info).
  * `dt = 0.1` : (only for continuous) Time of evolution between each check of
    distance exceeding the `threshold`.

  * `inittest = (st1, d0) -> st1 .+ d0/sqrt(D)` :
    A function that given
    `(st1, d0)` initializes the test state with distance
    `d0` from the given state `st1` (`D` is the dimension
    of the system). This function can be used when you want to avoid
    the test state appearing in a region of the phase-space where it would have
    e.g. different energy or escape to infinity.

For the continuous case, the algorithm becomes faster with increasing `dt`, since
integration is interrupted less frequenty. For the fastest performance you want to
fine-tune `dt, d0, threshold` such that you have the minimum amount of rescalings
but you want to still be well within the linearized dynamics region.

[1] : G. Benettin *et al.*, Phys. Rev. A **14**, pp 2338 (1976)
"""
function lyapunov(ds::DiscreteDS,
                  N::Int,
                  return_convergence::Type{Val{B}} = Val{false};
                  Ttr::Int = 0,
                  d0=1e-9,
                  threshold=1e-5,
                  inittest = inittest_default(dimension(ds))
                  ) where {B}

    threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))
    eom = ds.eom
    st1 = deepcopy(ds.state)
    st2 = inittest(st1, d0)

    # transient system evolution
    for i in 1:Ttr
        st1 = eom(st1)
    end

    if B
        λts, ts = _lyapunov_full(eom, st1, st2, N, d0, threshold)
        return λts, ts
    else
        λ = _lyapunov_final(eom, st1, st2, N, d0, threshold)
        return λ
    end
end

function _lyapunov_full(eom, st1, st2, N, d0, threshold)
    dist = d0
    λ = zero(eltype(st1))
    λs = eltype(st1)[]
    ts = Int[]
    i = 0
    while i < N
        #evolve until rescaling:
        while dist < threshold
            st1 = eom(st1)
            st2 = eom(st2)
            dist = norm(st1 - st2)
            i+=1
            i>=N && break
        end
        # local lyapunov exponent is simply the relative distance of the trajectories
        a = dist/d0
        λ += log(a)
        push!(λs, λ/i)
        push!(ts, i)
        i>=N && break
        #rescale:
        st2 = st1 + (st2 - st1)/a #must rescale in direction of difference
        dist = d0
    end
    return λs, ts
end

function _lyapunov_final(eom, st1, st2, N, d0, threshold)
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
            i>=N && break
        end
        # local lyapunov exponent is simply the relative distance of the trajectories
        a = dist/d0
        λ += log(a)
        i>=N && break
        #rescale:
        st2 = st1 + (st2 - st1)/a #must rescale in direction of difference
        dist = d0
    end
    return λ/i
end

inittest_default(D) = (state1, d0) -> state1 .+ d0/sqrt(D)

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
#                          Continuous Lyapunov Helpers                              #
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
        warnstr = "Absolute tolerance (abstol) of integration is much bigger than "
        warnstr*= "`d0`! It is highly suggested to increase it using `diff_eq_kwargs`."
        warn(warnstr)
    end
    if rtol > 10d0
        warnstr = "Relative tolerance (reltol) of integration is much bigger than "
        warnstr*= "`d0`! It is highly suggested to increase it using `diff_eq_kwargs`."
        warn(warnstr)
    end
end



#####################################################################################
#                            Continuous Lyapunovs                                   #
#####################################################################################
function lyapunov(ds::ContinuousDynamicalSystem,
                  T::Real,
                  return_convergence::Type{Val{B}} = Val{false};
                  Ttr = 0.0,
                  d0=1e-9,
                  threshold=1e-5,
                  dt = 1.0,
                  diff_eq_kwargs = Dict(:abstol=>d0, :reltol=>d0),
                  inittest = inittest_default(dimension(ds)),
                  ) where {B}

    check_tolerances(d0, diff_eq_kwargs)
    S = eltype(ds)

    T = convert(eltype(ds.state), T)
    threshold <= d0 && throw(ArgumentError("Threshold must be bigger than d0!"))

    # Transient system evolution
    Ttr != 0 && evolve!(ds, Ttr; diff_eq_kwargs = diff_eq_kwargs)

    # Create a copy integrator with different state
    # (workaround for https://github.com/JuliaDiffEq/DiffEqBase.jl/issues/58)
    # initialize:
    st1 = copy(ds.state)
    integ1 = ODEIntegrator(ds, T; diff_eq_kwargs=diff_eq_kwargs)
    integ1.opts.advance_to_tstop=true
    ds.state .= inittest(st1, d0)
    integ2 = ODEIntegrator(ds, T; diff_eq_kwargs=diff_eq_kwargs)
    integ2.opts.advance_to_tstop=true
    ds.state .= st1

    if B
        λts::Vector{S}, ts::Vector{S} =
        lyapunov_full(integ1, integ2, T, d0, threshold, dt, diff_eq_kwargs)
        return λts, ts
    else
        λ::S =
        lyapunov_final(integ1, integ2, T, d0, threshold, dt, diff_eq_kwargs)
        return λ
    end
end

function lyapunov_full(integ1::ODEIntegrator,
                  integ2::ODEIntegrator,
                  T::Real,
                  d0=1e-9,
                  threshold=1e-5,
                  dt = 0.1,
                  diff_eq_kwargs = Dict(:abstol=>d0, :reltol=>d0),
                  )

    S = eltype(integ1.u)
    dist = d0
    λ::S = 0.0
    λ_ts::Vector{S} = [0.0]  # the timeseries for the Lyapunov exponent
    ts::Vector{typeof(dt)} = [0.0]    # the time points of the timeseries
    i = 0;
    tvector = dt:dt:T

    # start evolution and rescaling:
    for τ in tvector
        # evolve until rescaling:
        push!(integ1.opts.tstops, τ)
        step!(integ1)
        push!(integ2.opts.tstops, τ)
        step!(integ2)
        dist = norm(integ1.u .- integ2.u)
        # Rescale:
        if dist ≥ threshold
            # add computed scale to accumulator (scale = local lyaponov exponent):
            a = dist/d0
            λ += log(a)
            push!(λ_ts, λ/τ)
            push!(ts, τ)
            # Rescale and reset everything:
            # Must rescale towards difference direction:
            @. integ2.u = integ1.u + (integ2.u - integ1.u)/a
            u_modified!(integ2, true)
            set_proposed_dt!(integ2, integ1)
            dist = d0
        end
    end
    λ_ts, ts
end

function lyapunov_final(integ1::ODEIntegrator,
                  integ2::ODEIntegrator,
                  T::Real,
                  d0=1e-9,
                  threshold=10^5*d0,
                  dt = 0.1,
                  diff_eq_kwargs = Dict(:abstol=>d0, :reltol=>d0),
                  )

    dist = d0
    λ::eltype(iteg1.u) = 0.0
    i = 0;
    tvector = dt:dt:T
    finalτ = dt

    # start evolution and rescaling:
    for τ in tvector
        # evolve until rescaling:
        push!(integ1.opts.tstops, τ)
        step!(integ1)
        push!(integ2.opts.tstops, τ)
        step!(integ2)
        dist = norm(integ1.u .- integ2.u)
        # Rescale:
        if dist ≥ threshold
            # add computed scale to accumulator (scale = local lyaponov exponent):
            a = dist/d0
            λ += log(a)
            finalτ = τ
            # Rescale and reset everything:
            # Must rescale towards difference direction:
            @. integ2.u = integ1.u + (integ2.u - integ1.u)/a
            u_modified!(integ2, true)
            set_proposed_dt!(integ2, integ1)
            dist = d0
        end
    end
    return λ/finalτ
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
