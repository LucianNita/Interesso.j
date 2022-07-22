module Interesso

import Ipopt
using MathOptInterface
const MOI = MathOptInterface;

using SparseArrays
using BlockArrays
using BarycentricInterpolation
using QuadGK
using LinearAlgebra: tril!
using Zygote: gradient
using ReverseDiff: hessian
using Reexport

include("quadrature.jl")
include("interpolation.jl")
include("mesh/mesh.jl")

abstract type InteressoAlgorithm end
include("algorithms/leastSquares.jl")
export LeastSquares

include("MOI/evaluator.jl")
include("MOI/fixedMesh.jl")
include("MOI/flexibleMesh.jl")

include("mesh/meshSolution.jl")

#@reexport using JuDOBase
include("JuDOBase/JuDOBase.jl")
include("interfaces/judo.jl")

#@reexport using DiffEqBase
#inlcude("interfaces/diffEq.jl")


end #module