using Interesso

@testset "exponentials" begin
    function F(du, u, p, t)
        res1 = du[1] + 1.0 * u[1];
        res2 = du[2] + 2.0 * u[2];
        return [res1, res2]
    end
    prob = DAEProblem(F, [-1.0, -4.0], [1.0, 2.0], (0.0, 10.0));
    
    alg = interesso(stages = 5, stateApproximationDegree = 3,
        residualsQuadratureOrder = 3
    );

    sol = solve(prob, alg);

    @test sol.u[1,:] ≈ [1.0, 2.0] atol=1e-8
end