struct DOSystem#{T <: Nonlinear}
    cost
    equations
    constraints
    differentialVariables::Vector
    algebraicVariables::Vector
    staticVariables::Vector
end