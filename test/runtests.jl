using DynamicalSystems, Tests, GLMakie

@testset "DynamicalSystems" begin
    ps = Dict(
    1 => 1:0.1:30,
    2 => 10:0.1:50,
    3 => 1:0.01:10.0)
    pnames = Dict(1 => "σ", 2 => "ρ", 3 => "β")

    lims = ((-30, 30), (-30, 30), (0, 100))

    ds = PredefinedDynamicalSystems.lorenz()

    @test interactive_orbitdiagram(ds, 1, 1.0, 3.0)
end

