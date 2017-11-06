"""
Sub-module of the module `DynamicalSystems`, which contains pre-defined
famous systems.
"""
module Systems
using DynamicalSystems, StaticArrays
#######################################################################################
#                                    Continuous                                       #
#######################################################################################
"""
```julia
lorenz(u0=[0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8/3)
```
```math
\\begin{aligned}
\\dot{X} &= \\sigma(Y-X) \\\\
\\dot{Y} &= -XZ + \\rho X -Y \\\\
\\dot{Z} &= XY - \\beta Z
\\end{aligned}
```
The famous three dimensional system due to Lorenz [1], shown to exhibit
so-called "deterministic nonperiodic flow". It was originally invented to study a
simplified form of atmospheric convection.

Currently, it is most famous for its strange attractor (occuring at the default
parameters), which resembles a butterfly. For the same reason it is
also associated with the term "butterfly effect" (a term which Lorenz himself disliked)
even though the effect applies generally to dynamical systems.
Default values are the ones used in the original paper.

[1] : E. N. Lorenz, J. atmos. Sci. **20**, pp 130 (1963)
"""
function lorenz(u0=[0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8/3)
    @inline @inbounds function eom_lorenz!(du, u)
        du[1] = σ*(u[2]-u[1])
        du[2] = u[1]*(ρ-u[3]) - u[2]
        du[3] = u[1]*u[2] - β*u[3]
    end
    @inline @inbounds function jacob_lorenz(u)
        i = one(eltype(u))
        o = zero(eltype(u))
        @SMatrix [-σ*i           σ*i    zero(i);
                  (ρ*i - u[3])   (-i)   (-u[1]);
                  u[2]           u[1]   (-β*i) ]
    end# should give exponents [0.9056, 0, -14.5723]
    return ContinuousDS(u0, eom_lorenz!, jacob_lorenz)
end


"""
```julia
roessler(u0=rand(3); a = 0.2, b = 0.2, c = 5.7)
```
```math
\\begin{aligned}
\\dot{x} &= -y-z \\\\
\\dot{y} &= x+ay \\\\
\\dot{z} &= -b + z(x-c)
\\end{aligned}
```
This three-dimensional continuous system is due to Rössler [1].
It is a system that by design behaves similarly
to the `lorenz` system and displays a (fractal)
strange attractor.

However, it is easier to analyze qualitatively, as for example
the attractor is composed of a single manifold.
Default values are the same as the original paper.

[1] : O. E. Rössler, Phys. Lett. **57A**, pp 397 (1976)
"""
function roessler(u0=rand(3); a = 0.2, b = 0.2, c = 5.7)
    @inline @inbounds function eom_roessler!(du, u)
        du[1] = -u[2]-u[3]
        du[2] = u[1] + a*u[2]
        du[3] = b + u[3]*(u[1] - c)
    end
    @inline @inbounds function jacob_roessler(u)
        i = one(eltype(u))
        o = zero(eltype(u))
        @SMatrix [o     -i      -i;
                  i      a       o;
                  u[3]   o       u[1] - c]
    end
  return ContinuousDS(u0, eom_roessler!, jacob_roessler)
end

"""
    double_pendulum(u0=rand(4); G=10.0, L1 = 1.0, L2 = 1.0, M1 = 1.0, M2 = 1.0)
Famous chaotic double pendulum system (also used for our logo!). Keywords
are gravity (G), lengths of each rod and mass of each ball (all assumed SI units).

The variables order is [θ1, dθ1/dt, θ2, dθ2/dt].

Jacobian is not created! So no `lyapunovs` for you!

(please contribute the Jacobian and the e.o.m. in LaTeX :smile:)
"""
function double_pendulum(u0=rand(4); G=10.0, L1 = 1.0, L2 = 1.0, M1 = 1.0, M2 = 1.0)
    @inline @inbounds function eom_dp!(du, state)
        du[1] = state[2]

        del_ = state[3] - state[1]
        den1 = (M1 + M2)*L1 - M2*L1*cos(del_)*cos(del_)
        du[2] = (M2*L1*state[2]*state[2]*sin(del_)*cos(del_) +
                   M2*G*sin(state[3])*cos(del_) +
                   M2*L2*state[4]*state[4]*sin(del_) -
                   (M1 + M2)*G*sin(state[1]))/den1

        du[3] = state[4]

        den2 = (L2/L1)*den1
        du[4] = (-M2*L2*state[4]*state[4]*sin(del_)*cos(del_) +
                   (M1 + M2)*G*sin(state[1])*cos(del_) -
                   (M1 + M2)*L1*state[2]*state[2]*sin(del_) -
                   (M1 + M2)*G*sin(state[3]))/den2
    end
    return ContinuousDS(u0, eom_dp!, nothing)
end

"""
    henonhelies(u0=[0, -0.25, 0.42081,0]; λ = 1)
```math
\\begin{aligned}
\\dot{x} &= p_x \\\\
\\dot{y} &= p_y \\\\
\\dot{p}_x &= -x -2\\lambda xy \\\\
\\dot{p}_y &= -y -\\lambda (x^2 - y^2)
\\end{aligned}
```

The Hénon–Heiles system [1] was introduced as a simplification of the motion
of a star around a galactic center. It was originally intended to study the
existence of a "third integral of motion" (which would make this 4D system integrable).
In that search, the authors encountered chaos, as the third integral existed
for only but a few initial conditions.

The default initial condition is a typical chaotic orbit.

[1] : Hénon, M. & Heiles, C., The Astronomical Journal **69**, pp 73–79 (1964)
"""
function henonhelies(u0=[0, -0.25, 0.42081, 0]; λ = 1)
    i = one(eltype(u0))
    o = zero(eltype(u0))
    @inline @inbounds function eom_hh!(du, u)
        du[1] = u[3]
        du[2] = u[4]
        du[3] = -u[1] - 2λ*u[1]*u[2]
        du[4] = -u[2] -λ*(u[1]^2 - u[2]^2)
    end
    @inline @inbounds function jacob_hh(u)
        @SMatrix [o    o     i    o;
                  o    o     o    i;
                  -i   -2λ*u[2]   o    o;
                  -2λ*u[1]  -1-2λ*u[2]  o   o]
    end
    return ContinuousDS(u0, eom_hh!, jacob_hh)
end

#######################################################################################
#                                     Discrete                                        #
#######################################################################################
"""
```julia
towel(u0 = [0.085, -0.121, 0.075])
```
```math
\\begin{aligned}
x_{n+1} &= a x_n (1-x_n) -0.05 (y_n +0.35) (1-2z_n) \\\\
y_{n+1} &= 0.1 \\left( \\left( y_n +0.35 \\right)\\left( 1+2z_n\\right) -1 \\right)
\\left( 1 -1.9 x_n \\right) \\\\
z_{n+1} &= 3.78 z_n (1-z_n) + b y_n
\\end{aligned}
```
The folded-towel map is a hyperchaotic mapping due to Rössler [1]. It is famous
for being a mapping that has the smallest possible dimensions necessary for hyperchaos,
having two positive and one negative lyapunov exponent.

The name comes from the fact that when plotted looks like a folded towel, in every
projection.

Default values are the ones used in the original paper.

[1] : O. E. Rössler, Phys. Lett. **71A**, pp 155 (1979)
"""
function towel(u0=[0.085, -0.121, 0.075])
    @inline @inbounds function eom_towel(x)
    x1, x2, x3 = x[1], x[2], x[3]
    SVector( 3.8*x1*(1-x1) - 0.05*(x2+0.35)*(1-2*x3),
    0.1*( (x2+0.35)*(1-2*x3) - 1 )*(1 - 1.9*x1),
    3.78*x3*(1-x3)+0.2*x2 )
    end

    @inline @inbounds function jacob_towel(x)
        @SMatrix [3.8*(1 - 2x[1]) -0.05*(1-2x[3]) 0.1*(x[2] + 0.35);
        -0.19((x[2] + 0.35)*(1-2x[3]) - 1)  0.1*(1-2x[3])*(1-1.9x[1])  -0.2*(x[2] + 0.35)*(1-1.9x[1]);
        0.0  0.2  3.78(1-2x[3]) ]
        end
    return DiscreteDS(u0, eom_towel, jacob_towel)
end# should result in lyapunovs: [0.432207,0.378834,-3.74638]

"""
```julia
standardmap(u0=0.001rand(2); k = 0.971635)
```
```math
\\begin{aligned}
\\theta_{n+1} &= \\theta_n + p_{n+1} \\\\
p_{n+1} &= p_n + k\\sin(\\theta_n)
\\end{aligned}
```
The standard map (also known as Chirikov standard map) is a two dimensional,
area-preserving chaotic mapping due to Chirikov [1]. It is one of the most studied
chaotic systems and by far the most studied Hamiltonian (area-preserving) mapping.

The map corresponds to the  Poincaré's surface of section of the kicked rotor system.
Changing the non-linearity parameter `k` transitions the system from completely periodic
motion, to quasi-periodic, to local chaos (mixed phase-space) and finally to global
chaos.

The default parameter `k` is the critical parameter where the golden-ratio torus is
destroyed, as was calculated by Greene [2]. The e.o.m. considers the angle variable
`θ` to be the first, and the angular momentum `p` to be the second, while
both variables
are always taken modulo 2π (the mapping is on the [0,2π)² torus).

[1] : B. V. Chirikov, Preprint N. **267**, Institute of
Nuclear Physics, Novosibirsk (1969)

[2] : J. M. Greene, J. Math. Phys. **20**, pp 1183 (1979)
"""
function standardmap(u0=0.001rand(2); k = 0.971635)
    const twopi = 2π
    @inline @inbounds function eom_standard(x)
        theta = x[1]; p = x[2]
        p+=k*sin(theta)
        theta += p
        while theta >= twopi; theta -= twopi; end
        while theta < 0; theta += twopi; end
        while p >= twopi; p -= twopi; end
        while p < 0; p += twopi; end
        return SVector(theta, p)
    end
    @inline @inbounds jacob_standard(x) =
    @SMatrix [1 + k*cos(x[1])    1;
              k*cos(x[1])        1]
  return DiscreteDS(u0, eom_standard, jacob_standard)
end

function coupledstandardmaps(M::Int, u0 = 0.001rand(2M);
    ks = ones(M), γ = 1.0)
    T = eltype(u0)
    const twopi = 2π
    thetaidx = (1:M)
    pidx = (M+1:2M)
    @inbounds function eom_coupledsm!(x)
        θs = @view x[thetaidx]
        ps = @view x[pidx]
        @. ps = mod2pi(
                ps + ks*sin(θs)
                -γ*(sin(circshift(θs, 1) - θs) + sin(circshift(θs, -1) - θs)))
        @. θs = mod2pi(θs + ps)
    end
    @inbounds function jacob_coupledsm!(J, x)
    end


    return BigDiscreteDS
end





"""
```julia
henon(u0=zeros(2); a = 1.4, b = 0.3)
```
```math
\\begin{aligned}
x_{n+1} &= 1 - ax^2_n+y_n \\\\
y_{n+1} & = bx_n
\\end{aligned}
```
The Hénon map is a two-dimensional mapping due to Hénon [1] that can display a strange
attractor (at the default parameters). In addition, it also displays many other aspects
of chaos, like period doubling or intermittency, for other parameters.

According to the author, it is a system displaying all the properties of the
Lorentz system (1963) while being as simple as possible.

Default values are the ones used in the original paper.

[1] : M. Hénon, Commun.Math. Phys. **50**, pp 69 (1976)
"""
function henon(u0=zeros(2); a = 1.4, b = 0.3)
    @inline @inbounds eom_henon(x) = SVector{2}(1.0 - a*x[1]^2 + x[2], b*x[1])
    @inline @inbounds jacob_henon(x) = @SMatrix [-2*a*x[1] 1.0; b 0.0]
    return DiscreteDS(u0, eom_henon, jacob_henon)
end # should give lyapunov exponents [0.4189, -1.6229]

"""
```julia
logistic(x0 = rand(); r = 4.0)
```
```math
x_{n+1} = rx_n(1-x_n)
```
The logistic map is an one dimensional unimodal mapping due to May [1] and is used by
many as the archetypal example of how chaos can arise from very simple equations.

Originally intentend to be a discretized model of polulation dynamics, it is now famous
for its bifurcation diagram, an immensly complex graph that that was shown
be universal by Feigenbaum [2].

[1] : R. M. May, Nature **261**, pp 459 (1976)

[2] : M. J. Feigenbaum, J. Stat. Phys. **19**, pp 25 (1978)
"""
function logistic(x0=rand(); r = 4.0)
    @inline eom_logistic(x) = r*x*(1-x)
    @inline deriv_logistic(x) = r*(1-2x)
    return DiscreteDS1D(x0, eom_logistic, deriv_logistic)
end




end# Systems module
