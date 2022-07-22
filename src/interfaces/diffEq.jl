"""
    solve(prob::DiffEqBase.AbstractDAEProblem, alg::InteressoAlgorithm)

Reexport DiffEqBase `solve` to solve the problem using an Interesso algorithm.

# Arguments
- `prob::DiffEqBase.AbstractDAEProblem`: DAE problem definition as required by DifferentialEquations package.
- `alg::InteressoAlgorithm`: Interesso algorithm definition (for more details see docs for `interesso`)
"""
function DiffEqBase.__solve(prob::DiffEqBase.AbstractDAEProblem, alg::InteressoAlgorithm)

    M = Mesh(prob, alg);

    optimizer = alg.optimizer;

    #while true
    #MOI.empty!(optimizer);

    d = MOI.add_variables(optimizer, M.dSize);
    MOI.set(optimizer, MOI.VariablePrimalStart(), d, M.dStart);

    MOI.set(optimizer, MOI.NLPBlock(), M.nlpBlockData);

    MOI.optimize!(optimizer);

    dOptimal = MOI.get(optimizer, MOI.VariablePrimal(), d);

    #Assemble solution
    #DiffEqBase.build_solution(prob, alg, t, u) #interp=LagrangeInterpolation()
    return
end