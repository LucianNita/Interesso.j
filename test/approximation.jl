using Interesso: interpolationPoints, interpolationMatrices

@testset "points" begin
    approximationDegree = 3;
    points = interpolationPoints(approximationDegree + 1);
    @test points isa Vector
    @test length(points) == approximationDegree + 1
end

@testset "residuals" begin
    variableSize = 3;
    approximationDegree = 3;
    evaluationPoints = [-0.7745966692414834, 0.0, 0.7745966692414834];

    (A, ADot) = interpolationMatrices(variableSize, approximationDegree, evaluationPoints);
    approximationOrder = approximationDegree + 1;

    @test size(A, 1) == variableSize * length(evaluationPoints)
    @test size(ADot, 1) == variableSize * length(evaluationPoints)
    @test size(A, 2) == variableSize * approximationOrder
    @test size(ADot, 2) == variableSize * approximationOrder

    d = ones(variableSize * approximationOrder);
    χ = A * d;
    χDot = ADot * d;

    @test χ ≈ ones(variableSize * length(evaluationPoints)) atol=1e-14
    @test χDot ≈ zeros(variableSize * length(evaluationPoints)) atol=1e-14
end

@testset "lhs" begin
    variableSize = 3;
    approximationDegree = 3;
    evaluationPoints = [-1.0];

    (lhsA, lhsADot) = interpolationMatrices(variableSize, approximationDegree, evaluationPoints);
    approximationOrder = approximationDegree + 1;
    
    @test size(lhsA, 1) == variableSize
    @test size(lhsADot, 1) == variableSize
    @test size(lhsA, 2) == variableSize * approximationOrder
    @test size(lhsADot, 2) == variableSize * approximationOrder

    d = ones(variableSize * approximationOrder);
    lhsχ = lhsA * d;
    lhsχDot = lhsADot * d;

    @test lhsχ ≈ ones(variableSize) atol=1e-14
    @test lhsχDot ≈ zeros(variableSize) atol=1e-14
end

@testset "rhs" begin end