struct DynamicSolution
    interpolants::Vector
end

struct DOProblemSolution <: AbstractDOSolution
    differentialVariables::Vector{DynamicSolution}
    algebraicVariables::Vector{DynamicSolution}
    optimizerResults::Dict
end

@recipe function f(d::DynamicSolution)
    seriestype --> :path
    linecolor  --> :black
    legend     --> false
    for s in d.interpolants
        @series begin
            s[1], s[2], s[3]
        end
    end
end