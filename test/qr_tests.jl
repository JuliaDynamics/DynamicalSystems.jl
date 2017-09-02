using DynamicalSystems, Base.Test

@testset "QR-decomposition" begin
    @testset "Base.Matrix" begin
        A = [1. 1 0;
             1  0 1;
             0  1 1]

        QA, RA = DynamicalSystems.qr_sq(A)
        QtA =  [1/√2   1/√6  -1/√3;
                1/√2  -1/√6   1/√3;
                0      2/√6   1/√3]
        RtA =  [2/√2,  3/√6,  2/√3]
        for i in length(QA)
            @test isapprox(QA[i], QtA[i], rtol=1e-6)
        end
        for i in 1:3
            @test RA[i] ≈ RtA[i]
        end
    end
end
