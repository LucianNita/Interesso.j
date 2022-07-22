abstract type InteressoMesh end

struct MeshVariable
    bounds::Tuple{Real, Real}
    guess::Real
end

struct FixedMesh <: InteressoMesh
    dMeshVariables::Vector{Vector{Vector{MeshVariable}}}
    aMeshVariables::Vector{Vector{Vector{MeshVariable}}}
    intervalPoints::Vector
end

struct FlexibleMesh <: InteressoMesh
    dMeshVariables::Vector{Vector{Vector{MeshVariable}}}
    aMeshVariables::Vector{Vector{Vector{MeshVariable}}}
    tMeshVariables::Vector{MeshVariable}
    flexibility::Tuple{Real,Real}
end

function addMeshVariables(o::MOI.AbstractOptimizer, M::FixedMesh)
    iSize = length(M.dMeshVariables);
    lhsIndices = addMeshIntervalVariables(o, M.dMeshVariables[1], M.aMeshVariables[1])
    for i in 2:iSize
        lhsIndices = addMeshIntervalVariables(o, M.dMeshVariables[i], M.aMeshVariables[i], lhsIndices);
    end
    return nothing
end

function addMeshVariables(o::MOI.AbstractOptimizer, M::FlexibleMesh)
    iSize = length(M.dMeshVariables);
    lhsIndices = addMeshIntervalVariables(o, M.dMeshVariables[1], M.aMeshVariables[1])
    lhsIndex = addTimeVariable(o, M.tMeshVariables[1])
    for i in 2:(iSize-1)
        lhsIndices = addMeshIntervalVariables(o, M.dMeshVariables[i], M.aMeshVariables[i], lhsIndices);
        lhsIndex = addTimeVariable(o, M.tMeshVariables[i], lhsIndex, M.flexibility);
    end
    addMeshIntervalVariables(o, M.dMeshVariables[iSize], M.aMeshVariables[iSize], lhsIndices);
    return nothing
end

function addMeshIntervalVariables(o::MOI.AbstractOptimizer,
    dMeshVariables::Vector{Vector{MeshVariable}}, aMeshVariables::Vector{Vector{MeshVariable}})

    addMeshVariables(o, [d[1] for d in dMeshVariables]);
    addMeshVariables.(o, [d[2:end-1] for d in dMeshVariables]);
    addMeshVariables.(o, [a[:] for a in aMeshVariables]);
    return addMeshVariables(o, [d[end] for d in dMeshVariables]);
end

function addMeshIntervalVariables(o::MOI.AbstractOptimizer,
    dMeshVariables::Vector{Vector{MeshVariable}}, aMeshVariables::Vector{Vector{MeshVariable}}, lhsIndices::Vector)

    rhsIndices = addMeshVariables(o, [d[1] for d in dMeshVariables]);
    addContinuityConstraint.(o, lhsIndices, rhsIndices);
    addMeshVariables.(o, [d[2:end-1] for d in dMeshVariables]);
    addMeshVariables.(o, [a for a in aMeshVariables]);
    return addMeshVariables(o, [d[end] for d in dMeshVariables]);
end

function addMeshVariables(o::MOI.AbstractOptimizer, variables::Vector{MeshVariable})
    vSize = length(variables);
    indices = MOI.add_variables(o, vSize);
    MOI.set(o, MOI.VariablePrimalStart(), indices, [v.guess for v in variables]);
    for v in 1:vSize
        if variables[v].bounds !== (-Inf, Inf)
            MOI.add_constraint(o, indices[v], MOI.GreaterThan(variables[v].bounds[1]));
            MOI.add_constraint(o, indices[v], MOI.LessThan(variables[v].bounds[2]));
        end
    end
    return indices
end

function addContinuityConstraint(o::MOI.AbstractOptimizer, lhsIndex, rhsIndex)
    terms = [MOI.ScalarAffineTerm(1.0, lhsIndex), MOI.ScalarAffineTerm(-1.0, rhsIndex)];
    constraint = MOI.ScalarAffineFunction(terms, 0.0);
    return MOI.add_constraint(o, constraint, MOI.EqualTo(0.0))
end

function addTimeVariable(o::MOI.AbstractOptimizer, variable::MeshVariable)
    index = MOI.add_variable(o);
    MOI.set(o, MOI.VariablePrimalStart(), index, variable.guess);
    if variable.bounds !== (-Inf, Inf)
        MOI.add_constraint(o, index, MOI.GreaterThan(variable.bounds[1]));
        MOI.add_constraint(o, index, MOI.LessThan(variable.bounds[2]));
    end
    return index
end

function addTimeVariable(o::MOI.AbstractOptimizer, variable::MeshVariable, lhsIndex, flexibility)
    rhsIndex = MOI.add_variable(o);
    MOI.set(o, MOI.VariablePrimalStart(), rhsIndex, variable.guess);
    if variable.bounds !== (-Inf, Inf)
        MOI.add_constraint(o, rhsIndex, MOI.GreaterThan(variable.bounds[1]));
        MOI.add_constraint(o, rhsIndex, MOI.LessThan(variable.bounds[2]));
    end
    addFlexibilityConstraint(o, lhsIndex, rhsIndex, flexibility);
    return rhsIndex
end

function addFlexibilityConstraint(o::MOI.AbstractOptimizer, lhsIndex, rhsIndex, flexibility)
    terms = [MOI.ScalarAffineTerm(1.0, lhsIndex), MOI.ScalarAffineTerm(-1.0, rhsIndex)];
    constraint = MOI.ScalarAffineFunction(terms, 0.0);
    return MOI.add_constraints(o, [constraint, constraint], [MOI.LessThan(-flexibility[1]), MOI.GreaterThan(-flexibility[2])])
end