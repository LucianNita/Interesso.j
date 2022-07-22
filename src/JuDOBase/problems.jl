# Dynamic Optimisation Problem

# Dynamic Feasibility Problem
struct DFProblem{F} <: JuDOProblem{F}
    tspan::Tuple{F,F}
    x::Vector{DynamicVariable{F}}
    u::Vector{DynamicVariable{F}}
    p::Vector{StaticVariable{F}}
    F
    ∇F!
end

# Constructors