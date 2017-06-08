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
function λmax(ds::DynamicalSystem)
  d0 = 1e-9
  ds2 = setu
end
