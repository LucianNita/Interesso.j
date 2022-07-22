using SafeTestsets

@safetestset "Approximation" begin include("approximation.jl") end

@safetestset "Mesh" begin include("mesh.jl") end

@safetestset "DAEs" begin include("dae.jl") end