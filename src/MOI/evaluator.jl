abstract type MeshFlexibility end

struct InteressoEvaluator{MF<:MeshFlexibility} <: MOI.AbstractNLPEvaluator
    tspan::Tuple{Real,Real}
    dSize::Integer
    aSize::Integer
    diffScale::Vector
    xIndexing::Vector
    HBitMatrix::BitMatrix
    HStructure::Vector{Tuple}

    residualsIntegrand::Function
    residualsQuadratureOrder::Integer
    t_res::Vector
    w_res::Vector
    A_res_d::SparseMatrixCSC
    A_res_dDot::SparseMatrixCSC
    A_res_a::SparseMatrixCSC
end

function MOI.initialize(e::InteressoEvaluator, requested_features::Vector{Symbol})
    for feature in requested_features
        #check no requested features other than :Grad, :Hess
    end
    return requested_features
end

function MOI.features_available(::InteressoEvaluator)
    return [:Grad, :Hess]
end

function MOI.eval_constraint(::InteressoEvaluator, g, x)
    return nothing
end

function MOI.eval_constraint_jacobian(::InteressoEvaluator, J, x)
    return nothing
end