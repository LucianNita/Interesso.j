 function F(du, u, p, t)
    res1 = du[1] + 1.0 * u[1];
    res2 = du[2] + 2.0 * u[2];
    return [res1, res2]
end
prob = DAEProblem(F, [0.0, 0.0], [2.0, 1.0], (0.0, 10.0));

alg = interesso(stages = 5, stateApproximationDegree = 3,
    residualsQuadratureOrder = 3
);

sol = solve(prob, alg);
plot(sol);

#=
 function F(du, u, p, t)
    res1 = du[1] - 0.25 * u[1] + 0.01*u[1]*u[2];
    res2 = du[2] - 1.0 * u[2] - 0.01*u[1]*u[2];
    return [res1, res2]
end
prob = DAEProblem(F, [10.0, 10.0], [80.0, 30.0], (0.0, 10.0));

alg = interesso(stages = 2, stateApproximationDegree = 15,
    residualsQuadratureOrder = 30
); #

sol = solve(prob, alg);
plot(sol);
=#
#=
function F(du, u, p, t)
   res1 = du[1] - u[2];
   res2 = u[1] - t;
   res3 = du[2];
   return [res1, res2, res3]
end
prob = DAEProblem(F, [1.0, 1.0], [1.0, 1.0], (0.0, 10.0));

alg = interesso(stages = 50, stateApproximationDegree = 3,
   residualsQuadratureOrder = 4
);

sol = solve(prob, alg);
plot(sol);
=#
#=
function F(du, u, p, t)
   res1 = du[1] + 5.0*t*u[1]^2 - 5.0/t + 1/t^2;
   return [res1]
end
prob = DAEProblem(F, [1.0], [1.0], (1.0, 25.0));

alg = interesso(stages = 50, stateApproximationDegree = 3,
   residualsQuadratureOrder = 4
);

sol = solve(prob, alg);
plot(sol);
=#
