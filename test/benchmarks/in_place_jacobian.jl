using OrdinaryDiffEq, BenchmarkTools, StaticArrays
abstract type ContinuousDynamicalSystem end

function lorenz1(u0=[0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8//3)
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
    return ContinuousDS1(u0, eom_lorenz!, jacob_lorenz)
end
mutable struct ContinuousDS1{T<:AbstractVector, F, J} <: ContinuousDynamicalSystem
    state::T
    eom!::F
    jacob::J
end




function lorenz2(u0=[0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8//3)
    @inline @inbounds function eom_lorenz!(du, u)
        du[1] = σ*(u[2]-u[1])
        du[2] = u[1]*(ρ-u[3]) - u[2]
        du[3] = u[1]*u[2] - β*u[3]
    end
    i = one(eltype(u0))
    o = zero(eltype(u0))
    J = zeros(eltype(u0), 3, 3)
    J[1,:] .= [-σ*i,        σ*i,    o]
    J[2,:] .= [ρ*i - u0[3],   -i,   -u0[1]]
    J[3,:] .= [u0[2],        u0[1],   (-β*i)]

    @inline @inbounds function jacob_lorenz!(J, u)
        J[2,1] = ρ - u[3]
        J[2,3] = -u[1]
        J[3,1] = u[2]; J[3,2] = u[1]
    end
    return ContinuousDS2(u0, eom_lorenz!, jacob_lorenz!, J)
end
mutable struct ContinuousDS2{T<:Real, F, J} <: ContinuousDynamicalSystem
    state::Vector{T}
    eom!::F
    jacob::J
    J::Matrix{T}
end

### Test jacobian function speed
function bjac()
    ds1 = lorenz1()
    b1 = @benchmark $ds1.jacob($ds1.state)
    println("Jacobian of 1: ")
    Juno.render(ds1.jacob(ds1.state))
    Juno.render(b1)


    ds2 = lorenz2()
    b2 = @benchmark $ds2.jacob($ds2.J, $ds2.state)
    println("Jacobian of 2: ")
    Juno.render(ds2.J)
    Juno.render(b2)
end
bjac()


### TANGENT BUNDLE!!!
function tangentbundle_setup_integrator(ds::ContinuousDS1, t_final;
  diff_eq_kwargs=Dict())

    D = dimension(ds)
    f! = ds.eom!
    jac = ds.jacob

    tbeom! = (t, u, du) -> begin
        f!(view(du, :, 1), u)
        A_mul_B!(
            view(du, :, 2:D+1),
            jac(view(u, :, 1)),
            view(u, :, 2:D+1)
        )
    end

    S = [ds.state eye(eltype(ds.state), D)]
    tbprob = ODEProblem(tbeom!, S, (zero(t_final), t_final))
    tb_integ = init(tbprob, Tsit5(); diff_eq_kwargs..., save_everystep=false)
    return tb_integ
end


function tangentbundle_setup_integrator(ds::ContinuousDS2, t_final;
  diff_eq_kwargs=Dict())

    D = dimension(ds)
    f! = ds.eom!

    tbeom! = (t, u, du) -> begin
        us = view(u, :, 1)
        f!(view(du, :, 1), us)
        ds.jacob(ds.J, us)
        A_mul_B!(
            view(du, :, 2:D+1),
            ds.J,
            view(u, :, 2:D+1)
        )
    end

    S = [ds.state eye(eltype(ds.state), D)]
    tbprob = ODEProblem(tbeom!, S, (zero(t_final), t_final))
    tb_integ = init(tbprob, Tsit5(); diff_eq_kwargs..., save_everystep=false)
    return tb_integ
end



# Benchmark it


function lyapunovs1(ds::ContinuousDynamicalSystem, N::Real=1000;
    Ttr::Real = 0.0, diff_eq_kwargs::Dict = Dict(), dt::Real = 0.1)
    # Initialize
    tstops = dt:dt:N*dt
    D = dimension(ds)
    λ = zeros(eltype(ds.state), D)
    Q = eye(eltype(ds.state), D)
    # Transient evolution:
    Ttr != 0 && evolve!(ds, Ttr; diff_eq_kwargs = diff_eq_kwargs)
    # Create integrator for dynamics and tangent space:
    integ = tangentbundle_setup_integrator(
    ds, tstops[end]; diff_eq_kwargs = diff_eq_kwargs)

    # Main algorithm
    for τ in tstops
        integ.u[:, 2:end] .= Q # update tangent dynamics state (super important!)
        u_modified!(integ, true)
        # Integrate
        while integ.t < τ
            step!(integ)
        end
        # Perform QR (on the tangent flow):
        Q, R = qr(view(integ.u, :, 2:D+1))
        # Add correct (positive) numbers to Lyapunov spectrum
        for j in 1:D
            λ[j] += log(abs(R[j,j]))
        end
    end
    λ./(integ.t) #return spectrum
end

function lyapunovs1(ds::ContinuousDynamicalSystem, N::Real=1000;
    Ttr::Real = 0.0, diff_eq_kwargs::Dict = Dict(), dt::Real = 0.1)
    # Initialize
    tstops = dt:dt:N*dt
    D = dimension(ds)
    λ = zeros(eltype(ds.state), D)
    Q = eye(eltype(ds.state), D)
    # Transient evolution:
    Ttr != 0 && evolve!(ds, Ttr; diff_eq_kwargs = diff_eq_kwargs)
    # Create integrator for dynamics and tangent space:
    integ = tangentbundle_setup_integrator(
    ds, tstops[end]; diff_eq_kwargs = diff_eq_kwargs)

    # Main algorithm
    for τ in tstops
        integ.u[:, 2:end] .= Q # update tangent dynamics state (super important!)
        u_modified!(integ, true)
        # Integrate
        while integ.t < τ
            step!(integ)
        end
        # Perform QR (on the tangent flow):
        Q, R = qr(view(integ.u, :, 2:D+1))
        # Add correct (positive) numbers to Lyapunov spectrum
        for j in 1:D
            λ[j] += log(abs(R[j,j]))
        end
    end
    λ./(integ.t) #return spectrum
end
dimension(ds) = 3

function binteg()
    A = Dict(:abstol=> 1e-6, :reltol=>1e-6)
    ds1 = lorenz1()
    ls1 = lyapunovs1(ds1, 5000.0, diff_eq_kwargs=A)
    b1 = @benchmark lyapunovs1($ds1, 500.0, diff_eq_kwargs=$A)
    println("Lyapunovs 1: $(ls1)")
    Juno.render(b1)

    ds2 = lorenz2()
    ls2 = lyapunovs1(ds2, 5000.0, diff_eq_kwargs=A)
    b2 = @benchmark lyapunovs1($ds2, 500.0, diff_eq_kwargs=$A)
    println("Lyapunovs 2: $(ls2)")

    Juno.render(b2)

end
binteg()



#
# b1 = @benchmark ds1.jacob(ds1.state)
# println("Jacobian of 1: ")
# Juno.render(b1)
#
#
# ds2 = lorenz2()
# b2 = @benchmark ds2.jacob(ds2.J, ds2.state)
# println("Jacobian of 2: ")
# Juno.render(b2)
