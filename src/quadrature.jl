τScale(tA, tB) = (tB - tA)/2;
τShift(tA, tB) = (tA + tB)/2;

function stageQuadrature(stages::Integer, stagePoints::Vector, points::Vector, weights::Vector)
    τScales = [τScale(stagePoints[i], stagePoints[i+1]) for i in 1:stages];
    τShifts = [τShift(stagePoints[i], stagePoints[i+1]) for i in 1:stages];
    
    t = [points * τScales[i] .+ τShifts[i] for i in 1:stages];
    w = [weights * τScales[i] for i in 1:stages];
    return (t, w)
end