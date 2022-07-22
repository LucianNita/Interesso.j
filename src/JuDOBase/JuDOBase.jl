module JuDOBase

using RecipesBase

# Variables
abstract type JuDOVariable{F<:AbstractFloat} end
include("variables.jl")
export StaticVariable, DynamicVariable

# Problems
abstract type JuDOProblem{F<:AbstractFloat} end
include("problems.jl")
export DFProblem

# Solutions
abstract type JuDOSolution{F<:AbstractFloat} end
include("solutions.jl")

end