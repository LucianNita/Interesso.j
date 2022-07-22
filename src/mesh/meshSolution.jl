struct MeshSolution
    intervalPoints::Vector
    dMeshVariables::Vector
    aMeshVariables::Vector
end

function getMeshSolution(o::MOI.AbstractOptimizer, M::FixedMesh, e::InteressoEvaluator{Fixed})
    xStaged = [o.inner.x[xIndices] for (_, xIndices, _) in e.xIndexing];
    iSize = length(M.dMeshVariables);
    
    dParameterizations = [length(M.dMeshVariables[1][d]) for d in 1:e.dSize];
    aParameterizations = [length(M.aMeshVariables[1][a]) for a in 1:e.aSize];
    stageIndices = (
        1:e.dSize,
        e.dSize+1 : sum(dParameterizations.-1),
        sum(dParameterizations.-1)+1 : sum(dParameterizations.-1)+sum(aParameterizations),
        sum(dParameterizations.-1)+sum(aParameterizations)+1 : sum(dParameterizations)+sum(aParameterizations)
    );

    dMeshVariables = [[vcat(
        xStaged[i][stageIndices[1]][d],
        xStaged[i][stageIndices[2]][sum(dParameterizations[1:d].-2)-(dParameterizations[d]-2)+1 : sum(dParameterizations[1:d].-2)],
        xStaged[i][stageIndices[4]][d])
        for d in 1:e.dSize]
        for i in 1:iSize];

    aMeshVariables = [[
        xStaged[i][stageIndices[3]][sum(aParameterizations[1:a])-aParameterizations[a]+1 : sum(aParameterizations[1:a])]           
        for a in 1:e.aSize]
        for i in 1:iSize];
    return MeshSolution(M.intervalPoints, dMeshVariables, aMeshVariables)
end

function getMeshSolution(o::MOI.AbstractOptimizer, M::FlexibleMesh, e::InteressoEvaluator{Flexible})
    xAug = vcat(e.tspan[1], o.inner.x, e.tspan[2]);
    tStages = [xAug[xAugIndices][1] for (_, xAugIndices, _) in e.xIndexing];
    intervalPoints = vcat(tStages, e.tspan[2]);

    xStaged = [xAug[xAugIndices][2:end-1] for (_, xAugIndices, _) in e.xIndexing];
    iSize = length(M.dMeshVariables);
    
    dParameterizations = [length(M.dMeshVariables[1][d]) for d in 1:e.dSize];
    aParameterizations = [length(M.aMeshVariables[1][a]) for a in 1:e.aSize];
    stageIndices = (
        1:e.dSize,
        e.dSize+1 : sum(dParameterizations.-1),
        sum(dParameterizations.-1)+1 : sum(dParameterizations.-1)+sum(aParameterizations),
        sum(dParameterizations.-1)+sum(aParameterizations)+1 : sum(dParameterizations)+sum(aParameterizations)
    );

    dMeshVariables = [[vcat(
        xStaged[i][stageIndices[1]][d],
        xStaged[i][stageIndices[2]][sum(dParameterizations[1:d].-2)-(dParameterizations[d]-2)+1 : sum(dParameterizations[1:d].-2)],
        xStaged[i][stageIndices[4]][d])
        for d in 1:e.dSize]
        for i in 1:iSize];

    aMeshVariables = [[
        xStaged[i][stageIndices[3]][sum(aParameterizations[1:a])-aParameterizations[a]+1 : sum(aParameterizations[1:a])]           
        for a in 1:e.aSize]
        for i in 1:iSize];
    return MeshSolution(intervalPoints, dMeshVariables, aMeshVariables)
end