using JuMP, Gurobi
solver = GurobiSolver()
m = Model(solver = solver)

@variable(m, x[i=1:3] >= 0, Bin)
@constraint(m, 6*x[1] + 5*x[2] + 5*x[3] <= 10)

@objective(m, Max, 6*x[1] + 4*x[2] + 3*x[3])
solve(m)


model = cutting_planes(m, [1;2;3], 50)
