module Systems
using DynamicalSystems, StaticArrays
#######################################################################################
#                                    Continuous                                       #
#######################################################################################
"""
    lorenz(u0=[0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8/3)
The famous three dimensional system due to Lorenz [1], shown to exhibit
so-called "deterministic nonperiodic flow". It was originally invented to study a
simplified form of atmospheric convection.

Currently, it is most famous for its strange attractor (occuring at the default
parameters), which resembles a butterfly. For the same reason it is
also associated with the term "butterfly effect" (a term which Lorenz himself disliked)
even though the effect applies generally to dynamical systems.

Default values are the ones used in the original paper.

[1] E. N. Lorenz, J. atmos. Sci. **20**, pp 130 (1963)
"""
function lorenz(u0=[0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8/3)
  @inline eom_lorenz(u) =
  SVector{3}(σ*(u[2]-u[1]), u[1]*(ρ-u[3]) - u[2], u[1]*u[2] - β*u[3])
  @inline function jacob_lorenz(u)
    i = one(eltype(u))
    @SMatrix [-σ*i           σ*i    zero(i);
              (ρ*i - u[3])   (-i)   (-u[1]);
              u[2]           u[1]   (-β*i) ]
  end# should give exponents [0.9056, 0, -14.5723]
  return ContinuousDS(u0, eom_lorenz, jacob_lorenz)
end

#######################################################################################
#                                     Discrete                                        #
#######################################################################################
"""
    towel(u0=[0.085, -0.121, 0.075])

The folded-towel map is a hyperchaotic map due to Rössler [1]. It is famous
for being a mapping that has the smallest possible dimensions necessary for hyperchaos,
having two positive and one negative lyapunov exponent.

The name comes from the fact that when plotted looks like a folded towel, in every
projection.

Default values are the ones used in the original paper.

[1] : O. E. Rössler, Phys. Lett. A, **71A**, pp 155 (1979).
"""
function towel(u0=[0.085, -0.121, 0.075])
  @inline function eom_towel(x)
    x1, x2, x3 = x[1], x[2], x[3]
    SVector( 3.8*x1*(1-x1) - 0.05*(x2+0.35)*(1-2*x3),
    0.1*( (x2+0.35)*(1-2*x3) - 1 )*(1 - 1.9*x1),
    3.78*x3*(1-x3)+0.2*x2 )
  end

  @inline function jacob_towel(x)
    @SMatrix [3.8*(1 - 2x[1]) -0.05*(1-2x[3]) 0.1*(x[2] + 0.35);
    -0.19((x[2] + 0.35)*(1-2x[3]) - 1)  0.1*(1-2x[3])*(1-1.9x[1])  -0.2*(x[2] + 0.35)*(1-1.9x[1]);
    0.0  0.2  3.78(1-2x[3]) ]
  end
  return DiscreteDS(u0, eom_towel, jacob_towel)
end# should result in lyapunovs: [0.432207,0.378834,-3.74638]

function logistic(x0=rand(); r=4.0)
  @inline eom_logistic(x) = r*x*(1-x)
  @inline deriv_logistic(x) = r*(1-2x)
  return DiscreteDS1D(x0, eom_logistic, deriv_logistic)
end

"""
    henon(u0=zeros(2); a = 1.4, b = 0.3)
The Hénon map is a two-dimensional mapping due to Hénon [1] that can display a strange
attractor (at the default parameters). In addition, it also displays many other aspects
of chaos, like period doubling or intermittency, for other parameters.

According to the author, it is a system displaying all the properties of the
Lorentz system (1963) while being
as simple as possible.

Default values are the ones used in the original paper.

[1] : M. Hénon, Commun.Math. Phys. **50**, pp 69 (1976)
"""
function henon(u0=zeros(2); a = 1.4, b = 0.3)
  @inline eom_henon(x) = SVector{2}(1.0 - a*x[1]^2 + x[2], b*x[1])
  @inline jacob_henon(x) = @SMatrix [-2*a*x[1] 1.0; b 0.0]
  # should give exponents 0.4189 -1.6229
  return DiscreteDS(u0, eom_henon, jacob_henon)
end


























end# Systems module
