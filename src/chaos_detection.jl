#####################################################################################
#                            Variational Equations                                  #
#####################################################################################
function variational_eom(ds::DiscreteDS)
    @inline veom(u) = ds.jacob(u)*u
    return veom
end


# variational equations:
#
# Maps:
# Jacobian × State
# Flows:
# Jacobians × Vector INTEGRATE? see lecture notes

#####################################################################################
#                                        GALI                                       #
#####################################################################################
"""
    gali(ds::DynamicalSystem, k::Int, tmax [, Ws]; kwargs...) -> GALI_k, t
Compute ``\\text{GALI}_k`` [1] for a given `k` up to time `tmax`. Return
``\\text{GALI}_k(t)`` and time vector ``t``.

`Ws` is an optional argument
containing the deviation vectors ``w_i`` for ``i \\in [1,k]``. If not given,
random orthonormal vectors are chosen using `qr`.

### Keywords
* `threshold` : If `GALI_k` reaches the threshold iteration is terminated.
  Default values are `1e-15` for discrete and `1e-12` for continuous systems.
* `dt=0.1` : Time step of integration for continuous systems.

### Description
The Generalized Alignment Index,
``\\text{GALI}_k``, is an efficient (and very fast) indicator of chaotic or regular
behavior type in ``D``-dimensional *Hamiltonian* systems (``D`` is number of variables).
``\\text{GALI}_k`` depends critically of
the type of orbit resulting
from the initial condition `ds.state`. If it is a chaotic orbit, then
```math
\\text{GALI}_k(t) \\approx
\\exp\\left[\\sum_{j=1}^k (\\lambda_1 - \\lambda_j)t \\right]
```
with ``\\lambda_1`` being the maximum [`lyapunov`](@ref) exponent.
If on the other hand the orbit is regular (movement in (D/2)-dimensional tori)
then it holds
```math
\\text{GALI}_k(t) \\approx
    \\begin{cases}
      \\text{const.}, & \\text{if} \\;\\; 2 \\le k \\le D/2  \\\\
      t^{-(2k - D)}, & \\text{if} \\;\\;  D/2 < k \\le D
    \\end{cases}
```
Traditionally, if ``\\text{GALI}_k(t)`` does not exceed the `threshold` until `tmax`,
the given orbit is said to be chaotic.

The entirety of our implementation is not based on the original paper, but rather in
the method described in [2], which uses the product of the singular values of ``A``,
a matrix that has as *columns* the deviation vectors.

### References

[1] : Skokos, C. H. *et al.*, Physica D **231**, pp 30–54 (2007)

[2] : Skokos, C. H. *et al.*, *Chaos Detection and Predictability* - Chapter 5
(section 5.3.1 and ref. [85] therein), Lecture Notes in Physics **915**,
Springer (2016)
"""
function gali(ds::DiscreteDS{D, S, F, J}, k::Int, tmax;
    threshold = 1e-15) where {D, S, F, J}

    Ws = qr(rand(D, D))[1]
    ws = Vector{SVector{D, S}}(k)
    for i in 1:k
        ws[i] = SVector{D, S}(Ws[:, i])
    end
    return gali(ds, k, tmax, ws; threshold = threshold)
end

function gali(ds::DiscreteDS{D, S, F, JJ}, k::Int, tmax, ws::Vector{SVector{D,S}};
    threshold = 1e-15) where {D,S,F,JJ}


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
