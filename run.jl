using Interesso
#=using Plots

dParameterization = 3;
q1 = DynamicVariable(dParameterization, true, (-Inf, Inf), 0.0, 0.96592582628, 0.0);
q2 = DynamicVariable(dParameterization, true, (-Inf, Inf), 0.0, 0.0, 0.0);
q3 = DynamicVariable(dParameterization, true, (-Inf, Inf), 0.0, 0.0, 0.0);
q4 = DynamicVariable(dParameterization, true, (-Inf, Inf), 1.0, 0.2588190451, 1.0);
ω1 = DynamicVariable(dParameterization, true, (-Inf, Inf), 0.0, 0.0, 0.0);
ω2 = DynamicVariable(dParameterization, true, (-Inf, Inf), 0.0, 0.0, 0.0);
ω3 = DynamicVariable(dParameterization, true, (-Inf, Inf), 0.0, 0.0, 0.0);

aParameterization = 1;
u1 = DynamicVariable(aParameterization, false, (-50.0, 50.0), nothing, nothing, 50.0);
u2 = DynamicVariable(aParameterization, false, (-50.0, 50.0), nothing, nothing, -50.0);
u3 = DynamicVariable(aParameterization, false, (-50.0, 50.0), nothing, nothing, -50.0);

t0 = 0.0;
tf = 28.630408;
tspan = (t0, tf);

function F(d, dDot, a, t)
    dq = dDot[1:4];
    q  = d[1:4];
    dω = dDot[5:7];
    ω  = d[5:7];
    u  = a;

    I = [5621.0, 4547.0, 2364.0];
    res1 = dq[1] - 0.5 * (ω[1]*q[4] - ω[2]*q[3] + ω[3]*q[2]);
    res2 = dq[2] - 0.5 * (ω[1]*q[3] + ω[2]*q[4] - ω[3]*q[1]);
    res3 = dq[3] - 0.5 * (-ω[1]*q[2] + ω[2]*q[1] + ω[3]*q[4]);
    #res4 = dq[4] - 0.5 * (-ω[1]*q[1] - ω[2]*q[2] - ω[3]*q[3]);
    res4 = q[1]^2 + q[2]^2 + q[3]^2 + q[4]^2 - 1.0;
    res5 = dω[1] - u[1] / I[1] + (I[3] - I[2]) / I[1] * ω[2] * ω[3];
    res6 = dω[2] - u[2] / I[2] + (I[1] - I[3]) / I[2] * ω[1] * ω[3];
    res7 = dω[3] - u[3] / I[3] + (I[2] - I[1]) / I[3] * ω[1] * ω[2];
    return [res1, res2, res3, res4, res5, res6, res7]
end

dop = FixedTimeDOProblem(tspan, [q1, q2, q3, q4, ω1, ω2, ω3], [u1, u2, u3], F);

#Fixed Mesh
s = Interesso.init(dop, Interesso.LeastSquares(
    stages=8,
    meshFlexibility=nothing));

#Flexible Mesh
intervals = 4;
flex = 0.3;
s = Interesso.init(dop, Interesso.LeastSquares(
    meshFlexibility=((tf/intervals)*(1.0-flex), (tf/intervals)*(1.0+flex))
));

solution = Interesso.solve!(s);

## Plotting
pq1 = plot(solution.differentialVariables[1]);
pq2 = plot(solution.differentialVariables[2]);
pq3 = plot(solution.differentialVariables[3]);
pq4 = plot(solution.differentialVariables[4]);
l = @layout [a; b; c; d];
plot(pq1, pq2, pq3, pq4, layout=l)

plot(solution.differentialVariables[5])
plot!(solution.differentialVariables[6])
plot!(solution.differentialVariables[7])

pu1 = plot(solution.algebraicVariables[1]);
pu2 = plot(solution.algebraicVariables[2]);
pu3 = plot(solution.algebraicVariables[3]);
l = @layout [a; b; c];
plot(pu1, pu2, pu3, layout=l)
=#