module JuDOBase

using Reexport
@reexport using Symbolics
using RecipesBase

@variables(t);
Dot = Differential(t);
export t, Dot

abstract type AbstractDOProblem end

abstract type AbstractDOSolution end

function solve(prob::AbstractDOProblem, alg; kwargs...) end

include("variables.jl")
#include("system.jl")
include("problems.jl")
include("solutions.jl")


end #module