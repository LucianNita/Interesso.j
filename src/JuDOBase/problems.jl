export FixedTimeDOProblem

mutable struct FixedTimeDOProblem <: AbstractDOProblem
    tspan::Tuple{Real, Real}

    differentialVariables::Vector{DynamicVariable}
    algebraicVariables::Vector{DynamicVariable}

    residuals::Function #F(d, dDot, a, t)
end