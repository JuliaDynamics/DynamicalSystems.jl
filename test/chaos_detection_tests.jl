using DynamicalSystems, Base.Test
using StaticArrays
using LsqFit: curve_fit

test_value = (val, vmin, vmax) -> @test vmin <= val <= vmax

println("\nTesting chaos detection algorithms...")

@testset "GALI discrete" begin

    @testset "Chaotic - towel map" begin
        ds = Systems.towel()
        model(x,p)= @. exp(-p[1]*x)
        ls = lyapunovs(ds, 10000)
        threshold = 1e-16
        @testset "k=$k" for k in [2,3]
            ex = sum(ls[1] - ls[j] for j in 2:k)
            g, t = DynamicalSystems.gali(ds, k, 1000; threshold = threshold)
            fite = curve_fit(model, t, g, [ex]).param[1]
            @test g[end] < threshold
            @test t[end] < 1000
            @test isapprox(fite, ex, rtol=1)
        end
    end

    @testset "2D map: Standard Map" begin
        ds = Systems.standardmap()
        k = 2
        @testset "chaotic" begin
            g, t = DynamicalSystems.gali(ds, k, 1000)
            @test g[end] ≤ 1e-12
            @test t[end] < 1000
        end
        @testset "regular" begin
            ds.state = [π, rand()]
            g, t = DynamicalSystems.gali(ds, k, 1000)
            @test t[end] == 1000
            @test g[end] > 1/1000^2
        end
    end
end
