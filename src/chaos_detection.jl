export gali
#####################################################################################
#                               Continuous GALI                                     #
#####################################################################################
function variational_eom_gali(ds::ContinuousDS, k::Int)
    f! = ds.eom!
    jac! = ds.jacob!
    J = ds.J
    # the equations of motion `veom!` evolve the system and the
    # deviation vectors
    # The e.o.m. for the system is f!(t, u , du).
    # The e.o.m. for the deviation vectors are tricky;
    # u[:, i] is the i deviation vector, for i≥2
    veom! = (t, u, du) -> begin
        f!(view(du, :, 1), u)
        jac!(J, u)
        # @views A_mul_B!(du[:, 2:end], ds.J, u[:, 2:end])
        for i in 1:k
            du[:, i+1] .= J*view(u, :, i+1)
        end
    end
    return veom!
end



"""
    gali(ds::DynamicalSystem, k::Int, tmax [, ws]; kwargs...) -> GALI_k, t
Compute ``\\text{GALI}_k`` [1] for a given `k` up to time `tmax`.
Return
``\\text{GALI}_k(t)`` and time vector ``t``.

`ws` is an optional argument
containing the deviation vectors ``w_i`` for ``i \\in [1,k]``, expected either
as a matrix with each column a deviation vector, or as a vector of vectors.
If not given,
random orthonormal vectors are chosen.

## Keyword Arguments
* `threshold` : If `GALI_k` falls below the `threshold` iteration is terminated.
  Default value is `1e-12`.
* `dt=1.0` : Time step between variational vector normalizations for continuous systems.
* `diff_eq_kwargs` : See [`trajectory`](@ref).

## Description
The Generalized Alignment Index,
``\\text{GALI}_k``, is an efficient (and very fast) indicator of chaotic or regular
behavior type in ``D``-dimensional *Hamiltonian* systems
(``D`` is number of variables). The behavior of
``\\text{GALI}_k(t)`` depends critically of
the type of orbit resulting
from the initial condition `ds.state`. If it is a chaotic orbit, then
```math
\\text{GALI}_k(t) \\sim
\\exp\\left[\\sum_{j=1}^k (\\lambda_1 - \\lambda_j)t \\right]
```
with ``\\lambda_1`` being the maximum [`lyapunov`](@ref) exponent.
If on the other hand the orbit is regular, corresponding
to movement in ``d``-dimensional torus with `` 1 \\le d \\le D/2``
then it holds
```math
\\text{GALI}_k(t) \\sim
    \\begin{cases}
      \\text{const.}, & \\text{if} \\;\\; 2 \\le k \\le d  \\; \\; \\text{and}
      \\; \\;d > 1 \\\\
      t^{-(k - d)}, & \\text{if} \\;\\;  d < k \\le D - d \\\\
      t^{-(2k - D)}, & \\text{if} \\;\\;  D - d < k \\le D
    \\end{cases}
```

Traditionally, if ``\\text{GALI}_k(t)`` does not become less than
the `threshold` until `tmax`
the given orbit is said to be chaotic, otherwise it is regular.

The entirety of our implementation is not based on the original paper, but rather in
the method described in [2], which uses the product of the singular values of ``A``,
a matrix that has as *columns* the deviation vectors.

## Performance Notes
If you want to do repeated evaluations of `gali` for many initial conditions and for
continuous systems, you can take advantage of the function:

    gali(integrator, k, W, tmax, dt, threshold)

in conjuction with `reinit!(integrator, W)` for various `W=cat(2, state, ws)`.
See the source code on how to
set-up the `integrator` and `W` for the first time.

## References

[1] : Skokos, C. H. *et al.*, Physica D **231**, pp 30–54 (2007)

[2] : Skokos, C. H. *et al.*, *Chaos Detection and Predictability* - Chapter 5
(section 5.3.1 and ref. [85] therein), Lecture Notes in Physics **915**,
Springer (2016)
"""
function gali(ds::ContinuousDS, k::Int, tmax::Real, ws::Matrix;
    threshold = 1e-12, dt = 1.0, diff_eq_kwargs = Dict())

    veom! = variational_eom_gali(ds, k)
    W = cat(2, ds.state, ws)
    prob = ODEProblem(veom!, W, (zero(dt), oftype(dt, tmax)))

    if haskey(diff_eq_kwargs, :saveat)
        pop!(diff_eq_kwargs, :saveat)
    end
    solver = get_solver!(diff_eq_kwargs)
    integrator = init(prob, solver; diff_eq_kwargs...,
    save_everystep=false, dense=false)

    return gali(integrator, k, W, tmax, dt, threshold)
end

function gali(ds::ContinuousDS, k::Int, tmax::Real;
    threshold = 1e-12, dt = 1.0, diff_eq_kwargs = Dict())
    D = dimension(ds)
    ws = qr(rand(D, D))[1][:, 1:k]
    gali(ds, k, tmax, ws;
    threshold = threshold, dt = dt, diff_eq_kwargs = diff_eq_kwargs)
end

function gali(ds::ContinuousDS, k::Int, tmax::Real, ws::AbstractVector;
    threshold = 1e-12, dt = 1.0,  diff_eq_kwargs = Dict())
    WS = cat(2, ws...)
    gali(ds, k, tmax, WS;
    threshold = threshold, dt = dt, diff_eq_kwargs = diff_eq_kwargs)
end

@inbounds function gali(integrator, k, W, tmax, dt, threshold)

    # warn("GALI has *not* been tested with periodic orbits of continuous systems!")
    rett = 0:dt:tmax
    gali_k = ones(eltype(W), length(rett))

    ti=1

    for ti in 2:length(rett)
        τ = rett[ti]
        # Evolve:
        while integrator.t < τ
            step!(integrator)
        end
        # Interpolate:
        integrator(W, τ)
        # Normalize
        for j in 1:k
            normalize!(view(W, :, j+1))
        end
        # Calculate singular values:
        zs = svdfact(view(W, :, 2:k+1))[:S]
        gali_k[ti] = prod(zs)
        if gali_k[ti] < threshold
            break
        end
        for j in 1:k
            normalize!(view(integrator.u, :, j+1))
        end
        u_modified!(integrator, true)
    end

    return gali_k[1:ti], rett[1:ti]

end

#####################################################################################
#                                 Discrete GALI                                     #
#####################################################################################
function gali(ds::DiscreteDynamicalSystem, k::Int, tmax;
    threshold = 1e-12)

    D = dimension(ds)
    Ws = qr(rand(D, D))[1][:, 1:k]
    return gali(ds, k, tmax, Ws; threshold = threshold)
end

function gali(ds::DiscreteDS{D, S, F, JJ}, k::Int, tmax, Ws::Matrix;
    threshold = 1e-12) where {D,S,F,JJ}

    ws = Vector{SVector{D, S}}(k)
    for i in 1:k
        ws[i] = SVector{D, S}(Ws[:, i])
    end
    return gali(ds, k, tmax, ws; threshold = threshold)
end

@inbounds function gali(ds::DiscreteDS{D, S, F, JJ}, k::Int,
    tmax, ws::Vector{SVector{D,S}};
    threshold = 1e-12) where {D,S,F,JJ}

    f = ds.eom
    J = ds.jacob
    x = ds.state

    rett = 0:Int(tmax)
    gali_k = ones(S, length(rett))

    ti=1

    for ti in 2:length(rett)
        # evolve state:
        x = f(x)
        # evolve all deviation vectors:
        jac = J(x)
        for i in 1:k
            ws[i] = normalize(jac*ws[i]) #gotta normalize bro!!!
        end
        # Calculate singular values:
        At = cat(2, ws...) # transpose of "A" in the paper, ref [2].
        zs = svdfact(At)[:S]
        gali_k[ti] =  prod(zs)
        if gali_k[ti] < threshold
            break
        end
    end

    return gali_k[1:ti], rett[1:ti]
end

function gali(ds::BigDiscreteDS, k::Int, tmax, Ws::AbstractVector;
    threshold = 1e-12)

    ws = cat(2, Ws...)
    return gali(ds, k, tmax, ws; threshold = threshold)
end
@inbounds function gali(ds::BigDiscreteDS, k::Int,
    tmax, ws::Matrix; threshold = 1e-12)

    f! = ds.eom!
    jacob! = ds.jacob!
    J = ds.J
    x = copy(ds.state)
    xprev = copy(x)
    wsdummy = copy(ws)

    rett = 0:Int(tmax)
    gali_k = ones(eltype(x), length(rett))

    ti=1
    zs = zeros(k)

    for ti in 2:length(rett)
        # evolve state:
        xprev .= x
        f!(x, xprev)
        jacob!(J, x)
        # evolve all deviation vectors:
        for i in 1:k
            A_mul_B!(view(ws, :, i), J, view(wsdummy, :, i))
            normalize!(@view ws[:, i]) #gotta normalize bro!!!
        end
        wsdummy .= ws
        # SVD fact for gali:
        zs .= svdfact(ws)[:S]
        # println("zs = $zs")

        gali_k[ti] =  prod(zs)
        if gali_k[ti] < threshold
            break
        end
    end

    return gali_k[1:ti], rett[1:ti]
end
