struct Fixed <: MeshFlexibility end

function integrateResiduals(e::InteressoEvaluator{Fixed}, x_i, w_res_i, t_res_i, diffScale_i)
    d_i = reshape(e.A_res_d * x_i, e.residualsQuadratureOrder, e.dSize);
    dDot_i = diffScale_i .* reshape(e.A_res_dDot * x_i, e.residualsQuadratureOrder, e.dSize);
    a_i = reshape(e.A_res_a * x_i, e.residualsQuadratureOrder, e.aSize);

    return sum(w_res_i .*
        [e.residualsIntegrand(d_i[j,:], dDot_i[j,:], a_i[j,:], t_res_i[j])
        for j in 1:e.residualsQuadratureOrder]
    )
end

#function integrateCost(e::InteressoEvaluator{Fixed}, x_i, WCost_i, tCost_i) end

function MOI.eval_objective(e::InteressoEvaluator{Fixed}, x)
    integratedResiduals = sum(
        integrateResiduals(e, x[xIndices], e.w_res[i], e.t_res[i], e.diffScale[i])
        for (i, xIndices, _) in e.xIndexing
    );
    #cost     
    return integratedResiduals #needs to be normalised by stage duration
end

function MOI.eval_objective_gradient(e::InteressoEvaluator{Fixed}, ∇objective, x)
    for (i, xIndices, _) in e.xIndexing
        ∇integratedResiduals = gradient(x_i -> integrateResiduals(e, x_i, e.w_res[i], e.t_res[i], e.diffScale[i]), x[xIndices])[1];
        #∇cost
        ∇objective[xIndices] = ∇integratedResiduals;
    end
    return nothing
end

function MOI.eval_hessian_lagrangian(e::InteressoEvaluator{Fixed}, H, x, σ, μ)
    for (i, xIndices, HIndices) in e.xIndexing
        H_integratedResiduals = hessian(x_i -> integrateResiduals(e, x_i, e.w_res[i], e.t_res[i], e.diffScale[i]), x[xIndices]);
        #H_cost
        H[HIndices] = H_integratedResiduals[e.HBitMatrix];
    end
    return nothing
end

function MOI.hessian_lagrangian_structure(e::InteressoEvaluator{Fixed})
    return e.HStructure
end

function preComputexIndexing(stages::Integer, xStageSize::Integer, ::Fixed)
    HStageSize = Integer(xStageSize*(xStageSize+1)/2);
    return [
        (i, #stage
        (i-1)*xStageSize+1 : i*xStageSize, #xIndices
        (i-1)*HStageSize+1 : i*HStageSize) #HIndices
    for i in 1:stages]
end

function preComputeHStructure(stages::Integer, xStageSize::Integer, ::Fixed)
    HBitMatrix = tril!(trues(xStageSize, xStageSize), 0);
    HStructure = [
        ((i-1)*xStageSize+k,(i-1)*xStageSize+j)
        for i in 1:stages for j in 1:xStageSize for k in 1:xStageSize if j <= k
    ];
    return (HBitMatrix, HStructure)
end