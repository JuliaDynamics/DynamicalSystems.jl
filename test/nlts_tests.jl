if current_module() != DynamicalSystems
  using DynamicalSystems
end
using Base.Test
using Distances: Cityblock, Euclidean

test_value = (val, vmin, vmax) -> @test vmin <= val <= vmax

println("\nTesting delay coordinates reconstruction...")
@testset "Henon Reconstruction" begin
    ds = Systems.henon()
    data = trajectory(ds, 100000)
    x = data[:, 1] # some "recorded" timeseries
    @testset "Sizes" begin
        for τ in [1, 2, 7]
            for D in [2, 3, 6]
                R = reconstruct(x, D, τ)
                @test length(size(R)) == 2
                @test length(R) == length(x) - (D-1)*τ
                @test length(R[1]) == D
            end
        end
    end
    @testset "Dimension" begin
        τ = 1; D = 2
        R = reconstruct(x, D, τ)
        D2 = information_dim(R)
        test_value(D2, 1.1, 1.3)
    end
    ks = 1:20
    @testset "Numerical Lyapunov" begin
        @testset "meth = $meth" for meth in
            [FixedMassNeighborhood(1), FixedMassNeighborhood(4),
            FixedSizeNeighborhood(0.01)]
            @testset "distance = $di" for di in [Euclidean(), Cityblock()]
                for D in [2, 4]
                    R = reconstruct(x, D, 1)
                    E = numericallyapunov(R, ks,
                    refstates = 1:1000, distance=di, method=meth)
                    λ = linear_region(ks, E)[2]
                    test_value(λ, 0.3, 0.5)
                end
            end
        end
    end
end
