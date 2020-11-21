using DynamicalSystems, BenchmarkTools, StaticArrays

ds = Systems.towel()

function lyapunovspectrum_normal(ds::DiscreteDS, N::Real; Ttr::Real = 100)

    u = evolve(ds, Ttr)
    D = length(u)
    eom = ds.eom
    jac = ds.jacob

    # Initialization
    λ = zeros(eltype(u), D)
    Q = @SMatrix eye(eltype(u), D)
    K = copy(Q)
    # Main algorithm
    for i in 1:N
        u = eom(u)
        K = jac(u)*Q

        Q, R = DynamicalSystems.qr_sq(K)
        for i in 1:D
            λ[i] += log(abs(R[i, i]))
        end
    end
    λ./N
end

function lyapunovspectrum_static(ds::DiscreteDS, N::Real; Ttr::Real = 100)

    u = evolve(ds, Ttr)
    D = length(u)
    eom = ds.eom
    jac = ds.jacob

    # Initialization
    λ = zeros(eltype(u), D)
    Q = @SMatrix eye(eltype(u), D)
    K = copy(Q)
    # Main algorithm
    for i in 1:N
        u = eom(u)
        K = jac(u)*Q

        Q, R = qr(K)
        for i in 1:D
            λ[i] += log(abs(R[i, i]))
        end
    end
    λ./N
end

N = 100000
lyapunovspectrum_normal(ds, 100000)
lyapunovspectrum_static(ds, 100000)

a = @btime lyapunovspectrum_normal(ds, 100000);
b = @btime lyapunovspectrum_static(ds, 100000);

# qr static is 3 times faster and does 0 allocations.
