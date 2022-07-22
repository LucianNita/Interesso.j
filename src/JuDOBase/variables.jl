export DynamicVariable, StaticVariable

struct DynamicVariable
    parameterization::Integer
    forceContinuity::Bool
    bounds::Tuple{Real,Real}

    initial::Union{Nothing,Real}
    terminal::Union{Nothing,Real}

    guess
end

struct StaticVariable
    bounds::Tuple{Real,Real}
    guess::Real
end