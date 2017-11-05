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
    gali(ds::DynamicalSystem, k::Int, t [, Ws]; kwargs...) -> GALI_k, rett
Compute ``\\text{GALI}_k`` [1] for a given `k` and time vector return it with
the time vector `rett`, which may be shorter than `t`. `Ws` is an optional argument
containing the deviation vectors ``w_i`` for ``i \\in [1,k]``.

The entirety of our implementation is not based on the origianl paper, but rather in
the method described in [2], which uses the product of the singular values of ``A``,
a matrix that has as *columns* the deviation vectors:
```math
A(t) = \\text{cat}(1, w_1, w_2, \\ldots, w_k) \\\\
Z(t)  =\\text{SVD}(A(t))[2]  \\\\
\\text{GALI}_k(t) = \\prod_1^k Z(t)  \\\\
```

### Keywords:
* `threshold` : If `GALI_k` reaches the threshold iteration is terminated.

### References:

[1] : Skokos, C. H. *et al.*, Physica D **231**, pp 30–54 (2007)

[2] : Skokos, C. H. *et al.*, *Chaos Detection and Predictability* - Chapter 5
(section 5.3.1 and ref. [85] therein), Lecture Notes in Physics **915**,
Springer (2016)
"""
function gali(ds, k::Int, A = qr(rand(dimension(ds), dimension(ds)))[1])
end

ds = Systems.towel()
t = 1:1000
k = 2
threshold = 1e-16
D = dimension(ds)
S = eltype(ds.state)
A = qr(rand(D, D))[1]
ws = Vector{SVector{D, S}}(k)
for i in 1:k
    ws[i] = SVector{D, S}(A[:, i])
end

f = ds.eom
J = ds.jacob

x = ds.state

gali_k = S[1.0]
rett = eltype(t)[0]

for τ in t

    # evolve all deviation vectors:
    for i in 1:k
        ws[i] = normalize(J(x)*ws[i]) #gotta normalize bro!!!
    end
    # evolve state:
    x = f(x)
    # Calculate singular values:
    At = cat(2, ws...) # transpose of "A" in the paper
    zs = svdfact(At)[:S]
    push!(gali_k, prod(zs))
    push!(rett, τ)
    if gali_k[end] < threshold
        break
    end

end
