using JuMP, Gurobi
solver = GurobiSolver()
model = Model(solver = solver)

@variable(model, y >= 0, Int)
@variable(model, x >= 0)

@constraints(model, begin
7y - 2x <= 14
0*y + x <= 3
0*y + 1x >= 0.2
2y - 2x <= 3
1y + 0*x == 2
end)

@objective(model, :Max, 4y - x)

solve(model)


model_F, convergence, z = cutting_planes(model, [1], 10)

println("Status = ", convergence)
println("Z = ", z)
