using Interesso

@testset "dae" begin
    function F(du, u, p, t)
        res1 = du[1] + 1.0 * u[1];
        res2 = du[2] + 2.0 * u[2];
        return [res1, res2]
    end
    prob = DAEProblem(F, [-1.0, -4.0], [1.0, 2.0], (0.0, 10.0));

    alg = interesso(stages = 5, stateApproximationDegree = 3,
        residualsQuadratureOrder = 3
    );

    M = Mesh(prob, alg);

    @test M.xSize == 2
    @test M.fSize == 2

    @test M.stages == 5
    @test M.tauScale == ones(M.stages)
    @test M.tauShift == [2.0 * i - 1.0 for i in 1:5]

    @test M.residualsQuadratureOrder == 3
    @test M.residualsQuadraturePoints ≈ [-0.7745966692414834, 0.0, 0.7745966692414834] atol=1e-12
    @test M.residualsQuadratureWeights ≈ [0.5555555555555556, 0.8888888888888888, 0.5555555555555556] atol=1e-12
    @test M.residualsQuadratureTimes[5] ≈ [8.225403330758517, 9.0, 9.774596669241483] atol=1e-12

    @test length(M.stateInterpolationPoints) == 4
    @test M.stateVectorSize == 8
    @test M.stageVectorSize == 8

    @test length(M.d0) == 40
end