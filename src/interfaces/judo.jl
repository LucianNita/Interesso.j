mutable struct JuDOSolver
    prob::JuDOBase.AbstractDOProblem
    alg::InteressoAlgorithm
    mesh::InteressoMesh
    optimizer::MOI.AbstractOptimizer
    nlpBlockData::MOI.NLPBlockData
end

function JuDOBase.solve(prob::JuDOBase.AbstractDOProblem, alg::InteressoAlgorithm; kwargs...)
    return solve!(init(prob, alg))
end

function solve!(s::JuDOSolver)
    #while true
    MOI.empty!(s.optimizer)

    addMeshVariables(s.optimizer, s.mesh);

    MOI.set(s.optimizer, MOI.NLPBlock(), s.nlpBlockData);
    MOI.set(s.optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE);

    MOI.optimize!(s.optimizer);

    meshSolution = getMeshSolution(s.optimizer, s.mesh, s.nlpBlockData.evaluator)

    optimizerResults = Dict(
        "objective" => MOI.get(s.optimizer, MOI.ObjectiveValue())
    )
    #judge solution

    #break

    #refineMesh!(s.mesh);
    #end
    return buildSolution(s.prob, s.alg, s.nlpBlockData.evaluator, meshSolution, optimizerResults)
end

function init(prob::JuDOBase.AbstractDOProblem, alg::InteressoAlgorithm)
    mesh = initialiseMesh(prob, alg);

    evaluator = initialiseEvaluator(prob, alg, mesh);

    nlpBlockData = MOI.NLPBlockData([], evaluator, true);

    return JuDOSolver(prob, alg, mesh, alg.optimizer, nlpBlockData)
end

function initialiseMesh(prob::JuDOBase.AbstractDOProblem, alg::InteressoAlgorithm)
    #Uniform Mesh
    intervalPoints = collect(range(prob.tspan[1], stop=prob.tspan[2], length=(alg.stages+1)));

    #Discretise Variables
    dMeshVariables = [[[MeshVariable(d.bounds, d.guess) for _ in 1:d.parameterization]
        for d in prob.differentialVariables] for _ in 1:alg.stages];
    aMeshVariables = [[[MeshVariable(a.bounds, a.guess) for _ in 1:a.parameterization]
        for a in prob.algebraicVariables] for _ in 1:alg.stages];

    #Overwrite initial and terminal conditions
    for d in 1:length(prob.differentialVariables)
        dMeshVariables[1][d][1] = MeshVariable((
            prob.differentialVariables[d].initial,
            prob.differentialVariables[d].initial),
            prob.differentialVariables[d].initial);
        dMeshVariables[end][d][end] = MeshVariable((
            prob.differentialVariables[d].terminal,
            prob.differentialVariables[d].terminal),
            prob.differentialVariables[d].terminal);
    end

    if alg.meshFlexibility isa Nothing
        return FixedMesh(dMeshVariables, aMeshVariables, intervalPoints)
    else
        tMeshVariables = [MeshVariable((-Inf, Inf), p) for p in intervalPoints[2:end-1]];
        #Overwrite initial and final nodes flexibility
        tMeshVariables[1] = MeshVariable(
            (alg.meshFlexibility[1]+prob.tspan[1], alg.meshFlexibility[2]+prob.tspan[1]),
            intervalPoints[2]);
        tMeshVariables[end] = MeshVariable(
            (prob.tspan[2]-alg.meshFlexibility[2], prob.tspan[2]-alg.meshFlexibility[1]),
            intervalPoints[end-1]);
        return FlexibleMesh(dMeshVariables, aMeshVariables, tMeshVariables, alg.meshFlexibility)
    end
end

function initialiseEvaluator(prob::JuDOBase.AbstractDOProblem, alg::InteressoAlgorithm, M::FixedMesh)
    dSize = length(prob.differentialVariables);
    aSize = length(prob.algebraicVariables);

    diffScale = [2.0 / (M.intervalPoints[i+1] - M.intervalPoints[i]) for i in 1:alg.stages];

    xStageSize = +(sum([d.parameterization for d in prob.differentialVariables]),
        sum([a.parameterization for a in prob.algebraicVariables]));  
    xIndexing = preComputexIndexing(alg.stages, xStageSize, Fixed());

    (HBitMatrix, HStructure) = preComputeHStructure(alg.stages, xStageSize, Fixed());

    residualsIntegrand(d, dDot, a, t) = sum(prob.residuals(d, dDot, a, t).^2);

    (tRes, wRes) = gauss(alg.residualsQuadratureOrder);

    (Ares_d, Ares_dDot, Ares_a) = barycentricInterpolation(
        [d.parameterization for d in prob.differentialVariables],
        [a.parameterization for a in prob.algebraicVariables],
        tRes);
        
    (tRes, wRes) = stageQuadrature(alg.stages, M.intervalPoints, tRes, wRes);

    return InteressoEvaluator{Fixed}(prob.tspan, dSize, aSize, diffScale,
        xIndexing, HBitMatrix, HStructure,
        residualsIntegrand,
        alg.residualsQuadratureOrder, tRes, wRes,
        Ares_d, Ares_dDot, Ares_a)
end

function initialiseEvaluator(prob::JuDOBase.AbstractDOProblem, alg::InteressoAlgorithm, ::FlexibleMesh)

    dSize = length(prob.differentialVariables);
    aSize = length(prob.algebraicVariables);

    diffScale = zeros(alg.stages); #not used

    xStageSize = +(sum([d.parameterization for d in prob.differentialVariables]),
        sum([a.parameterization for a in prob.algebraicVariables])); 
    xIndexing = preComputexIndexing(alg.stages, xStageSize, Flexible());

    (HBitMatrix, HStructure) = preComputeHStructure(alg.stages, xStageSize, Flexible());

    residualsIntegrand(d, dDot, a, t) = sum(prob.residuals(d, dDot, a, t).^2);

    (tRes, wRes) = gauss(alg.residualsQuadratureOrder);

    (Ares_d, Ares_dDot, Ares_a) = barycentricInterpolation(
        [d.parameterization for d in prob.differentialVariables],
        [a.parameterization for a in prob.algebraicVariables],
        tRes);

    return InteressoEvaluator{Flexible}(prob.tspan, dSize, aSize, diffScale,
        xIndexing, HBitMatrix, HStructure,
        residualsIntegrand,
        alg.residualsQuadratureOrder, tRes, wRes,
        Ares_d, Ares_dDot, Ares_a)
end

function buildSolution(prob::JuDOBase.AbstractDOProblem, alg::InteressoAlgorithm, e::InteressoEvaluator, s::MeshSolution, optimizerResults::Dict)
    dParameterizations = [d.parameterization for d in prob.differentialVariables];
    dInterpolants = [[Chebyshev2{p-1}(s.intervalPoints[i], s.intervalPoints[i+1]) for p in dParameterizations] for i in 1:alg.stages];
    dSolutions = [JuDOBase.DynamicSolution(
        [(interpolate(dInterpolants[i][d], s.dMeshVariables[i][d]), s.intervalPoints[i], s.intervalPoints[i+1]) for i in 1:alg.stages]
        ) for d in 1:e.dSize];
    aParameterizations = [a.parameterization for a in prob.algebraicVariables];
    aInterpolants = [[Chebyshev1{p-1}(s.intervalPoints[i], s.intervalPoints[i+1]) for p in aParameterizations] for i in 1:alg.stages];
    aSolutions = [JuDOBase.DynamicSolution(
        [(interpolate(aInterpolants[i][a], s.aMeshVariables[i][a]), s.intervalPoints[i], s.intervalPoints[i+1]) for i in 1:alg.stages]
        ) for a in 1:e.aSize];
    return JuDOBase.DOProblemSolution(dSolutions, aSolutions, optimizerResults)
end