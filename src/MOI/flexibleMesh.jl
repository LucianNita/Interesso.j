struct Flexible <: MeshFlexibility end

function integrateResiduals(e::InteressoEvaluator{Flexible}, xAug_i)
    x_i = xAug_i[2:end-1];
    τScale = (xAug_i[end] - xAug_i[1]) * 0.5;
    τShift = (xAug_i[end] + xAug_i[1]) * 0.5;
    diffScale = 1 / τScale;
    
    d_i = reshape(e.A_res_d * x_i, e.residualsQuadratureOrder, e.dSize);
    dDot_i = diffScale .* reshape(e.A_res_dDot * x_i, e.residualsQuadratureOrder, e.dSize);
    a_i = reshape(e.A_res_a * x_i, e.residualsQuadratureOrder, e.aSize);

    return sum(τScale * e.w_res .*
        [e.residualsIntegrand(d_i[j,:], dDot_i[j,:], a_i[j,:], τScale * e.t_res[j] + τShift)
        for j in 1:e.residualsQuadratureOrder]
    )
end

#function integrateCost(e::InteressoEvaluator{Flexible}, xAug_i) end

function MOI.eval_objective(e::InteressoEvaluator{Flexible}, x)
    xAug = vcat(e.tspan[1], x, e.tspan[2]);
    integratedResiduals = sum(
        integrateResiduals(e, xAug[xAugIndices]) for (_, xAugIndices, _) in e.xIndexing
    );
    #cost
    return integratedResiduals #needs to be normalised by stage duration
end

function MOI.eval_objective_gradient(e::InteressoEvaluator{Flexible}, ∇objective, x)
    xAug = vcat(e.tspan[1], x, e.tspan[2]);
    ∇objectiveAug = zeros(length(xAug));

    for (_, xAugIndices, _) in e.xIndexing
        ∇integratedResiduals = gradient(xAug_i -> integrateResiduals(e, xAug_i), xAug[xAugIndices])[1];
        #∇cost
        ∇objectiveAug[xAugIndices] .+= ∇integratedResiduals;
    end
    ∇objective[1:end] = ∇objectiveAug[2:end-1];
    return nothing
end

function MOI.eval_hessian_lagrangian(e::InteressoEvaluator{Flexible}, H, x, σ, μ)
    xAug = vcat(e.tspan[1], x, e.tspan[2]);

    for (stage, xAugIndices, HIndices) in e.xIndexing
        H_integratedResiduals = hessian(xAug_i -> integrateResiduals(e, xAug_i), xAug[xAugIndices]);
        #H_cost
        if stage == 1
            H[HIndices] = H_integratedResiduals[2:end,2:end][e.HBitMatrix[1:end-1,1:end-1]];
        elseif stage == length(e.xIndexing)
            H[HIndices] = H_integratedResiduals[1:end-1,1:end-1][e.HBitMatrix[1:end-1,1:end-1]];
        else
            H[HIndices] = H_integratedResiduals[e.HBitMatrix];
        end
    end
    return nothing
end

function MOI.hessian_lagrangian_structure(e::InteressoEvaluator{Flexible})
    return e.HStructure
end

function preComputexIndexing(stages::Integer, xStageSize::Integer, ::Flexible)
    HFirstSize = Integer((xStageSize+1)*(xStageSize+2)/2);
    xFirstIndexing = (1, 1:(xStageSize+1)+1, 1:HFirstSize);

    HOtherSize = Integer((xStageSize+2)*(xStageSize+3)/2);
    xOtherIndexing = [
        (i, #Stage
        (i-1)*(xStageSize+1)+1 : i*(xStageSize+1)+1, #xAugIndices
        HFirstSize+(i-2)*HOtherSize+1 : HFirstSize + (i-1)*HOtherSize) #HIndices
        for i in 2:(stages-1)];
    xLastIndexing = (
        stages,
        (stages-1)*(xStageSize+1)+1 : stages*(xStageSize+1)+1,
        HFirstSize+(stages-2)*HOtherSize+1 : 2*HFirstSize + (stages-2)*HOtherSize
        );
    return vcat(xFirstIndexing, xOtherIndexing, xLastIndexing);
end

function preComputeHStructure(stages::Integer, xStageSize::Integer, ::Flexible)
    HBitMatrix = tril!(trues(xStageSize+2, xStageSize+2), 0);

    HFirstStructure = [
        (k,j) for j in 1:(xStageSize+1) for k in 1:(xStageSize+1) if j <= k
    ];
    HOtherStructure = [
        ((i-1)*(xStageSize+1)-1+k,(i-1)*(xStageSize+1)-1+j)
        for i in 2:(stages-1) for j in 1:(xStageSize+2) for k in 1:(xStageSize+2) if j <= k
    ];
    HLastStructure = [
        ((stages-1)*(xStageSize+1)-1+k, (stages-1)*(xStageSize+1)-1+j)
        for j in 1:(xStageSize+1) for k in 1:(xStageSize+1) if j <= k
    ];
    return (HBitMatrix, vcat(HFirstStructure, HOtherStructure, HLastStructure))
end