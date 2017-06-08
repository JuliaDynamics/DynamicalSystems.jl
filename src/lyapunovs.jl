export λspectrum, λmax
#######################################################################################
#                                      Discrete                                       #
#######################################################################################
""" spectrum of lyapunovs """
function λspectrum(ds::DiscreteDS, N::Real; Ntr::Int= 100)
  N = convert(Int, N)
  u = deepcopy(ds.state)
  dim = length(u)
  eom = ds.eom
  jac = ds.jacob
  # Transient iterations
  for i in 1:Ntr
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
    R = F[:R]
    A_mul_signB!(Q, R) #ensure QR gives positive diagonal to R
    for i in 1:dim
      λ[i] += log(abs(R[i,i]))
    end
  end
  λ./N
end# this works. add tests and good to go.

"""
```julia
λmax(ds::DynamicalSystem, args...)
```
Calculate the maximum lyapunov exponent using a method due to Benettin [1], which simply
evolves two neighboring trajectories while constantly rescaling one of the two.

[1] : Benettin *et al.*, Phys. Rev. A **14**, pp 2338 (1976)
"""
function λmax(ds::DiscreteDS, N = 100000;
  d0=1e-7*one(eltype(ds.state)), threshold=10^3*d0)

  eom = ds.eom
  st1 = deepcopy(ds.state)
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

function λmax(ds::DiscreteDS1D, N=100000)

  eom = ds.eom
  der = ds.deriv
  x = deepcopy(ds.state)
  λ = log(abs(der(x)))
  for i in 1:N
    x = eom(x)
    λ += log(abs(der(x)))
  end
  λ/N
end
