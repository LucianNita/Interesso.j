struct LeastSquares <: InteressoAlgorithm
    optimizer::MOI.AbstractOptimizer
    stages::Integer
    meshFlexibility::Union{Nothing,Tuple{Real, Real}}
    residualsQuadratureOrder::Integer
    #costQuadratureOrder::Integer
    hessianApproximation::Bool
end

function LeastSquares(;
    optimizer=Ipopt.Optimizer(),
    stages=4,
    meshFlexibility=nothing,
    residualsQuadratureOrder=5,
    #costQuadratureOrder,
    hessianApproximation=false)
    return LeastSquares(
        optimizer,
        stages,
        meshFlexibility,
        residualsQuadratureOrder,
        #costQuadratureOrder,
        hessianApproximation
    )
end