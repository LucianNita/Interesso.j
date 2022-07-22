function barycentricInterpolation(
    differentialParameterizations::Vector,
    algebraicParameterizations::Vector,
    quadraturePoints::Vector)

    dSize = length(differentialParameterizations);
    aSize = length(algebraicParameterizations);

    dInterpolants = [Chebyshev2{d-1}() for d in differentialParameterizations];
    aInterpolants = [Chebyshev1{a-1}() for a in algebraicParameterizations];

    dInterpolationMatrices = [interpolation_matrix(dInterpolants[d], quadraturePoints) for d in 1:dSize];
    dDotInterpolationMatrices = [dInterpolationMatrices[d] * differentiation_matrix(dInterpolants[d]) for d in 1:dSize];
    aInterpolationMatrices = [interpolation_matrix(aInterpolants[a], quadraturePoints) for a in 1:aSize];

    quadratureOrder = length(quadraturePoints);
    A_d = spzeros(quadratureOrder * dSize, sum(differentialParameterizations) + sum(algebraicParameterizations));
    A_d = BlockArray(
        A_d,
        repeat([quadratureOrder], dSize),
        vcat(repeat([1], dSize), differentialParameterizations .- 2, algebraicParameterizations, repeat([1], dSize)));
    A_dDot = deepcopy(A_d);
    for d in 1:dSize
        A_d[Block(d,d)] = reshape(dInterpolationMatrices[d][:,1], quadratureOrder, 1);
        A_dDot[Block(d,d)] = reshape(dDotInterpolationMatrices[d][:,1], quadratureOrder, 1);
        A_d[Block(d,dSize+d)] = dInterpolationMatrices[d][:,2:(end-1)];
        A_dDot[Block(d,dSize+d)] = dDotInterpolationMatrices[d][:,2:(end-1)];
        A_d[Block(d,dSize*2+aSize+d)] = reshape(dInterpolationMatrices[d][:,end], quadratureOrder, 1);
        A_dDot[Block(d,dSize*2+aSize+d)] = reshape(dDotInterpolationMatrices[d][:,end], quadratureOrder, 1);
    end

    A_a = spzeros(quadratureOrder * aSize, sum(differentialParameterizations) + sum(algebraicParameterizations));
    A_a = BlockArray(
        A_a,
        repeat([quadratureOrder], aSize),
        vcat(repeat([1], dSize), differentialParameterizations .- 2, algebraicParameterizations,repeat([1], dSize)));
    for a in 1:aSize
        A_a[Block(a,dSize*2+a)] = aInterpolationMatrices[a];
    end

    A_d = Array(A_d);
    A_dDot = Array(A_dDot);
    A_a = Array(A_a);

    return (A_d, A_dDot, A_a)
end