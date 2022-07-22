using Interesso: flexPoints, flexWeights, InteressoQuadrature

@testset "flexPoints" begin
    normalisedPoints = collect(range(-1.0, stop=1.0, length=5));
    tA = 25.0;
    tB = 35;
    flexedPoints = flexPoints(normalisedPoints, tA, tB);
    @test flexedPoints[1] == tA;
    @test flexedPoints[3] == (tA+tB)/2;
    @test flexedPoints[5] == tB;
end

@testset "flexWeights" begin
    normalisedWeights = collect(range(-1.0, stop=1.0, length=5));
    tA = -10.0;
    tB = 10;
    flexedWeights = flexWeights(normalisedWeights, tA, tB);
    @test flexedWeights[1] == -10.0;
    @test flexedWeights[2] == -5.0;
    @test flexedWeights[3] == 0.0;
    @test flexedWeights[4] == 5.0;
    @test flexedWeights[5] == 10.0;
end

@testset "Quadrature" begin
    @testset "GaussLegendre" begin
        stagePoints = collect(range(0.0, stop=15.0, length=6));
        Q = InteressoQuadrature((:GaussLegendre, 3), stagePoints);
        @test Q.points[1][2] == 1.5;
        @test Q.points[5][2] == 13.5;
        @test Q.weights[1][1] == 5/9 * 3/2;
        @test Q.weights[2][2] == 8/9 * 3/2;
    end
end

@testset "numericalIntegration" begin
    stagePoints = collect(range(0.0, stop=15.0, length=6));
    Q = InteressoQuadrature((:GaussLegendre, 4), stagePoints);
    f(x) = 0.2 * x^3 - 0.4 * x^2 + 2.0 * x + 5.0;
    integration = sum(sum(Q.weights[i][j] * f(Q.points[i][j]) for j in 1:4) for i in 1:5);
    @test integration ≈ 2381.25;
end