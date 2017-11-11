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
            g, t = gali(ds, k, 1000; threshold = threshold)
            fite = curve_fit(model, t, g, [ex]).param[1]
            @test g[end] < threshold
            @test t[end] < 1000
            @test isapprox(fite, ex, rtol=1)
        end
    end
end
