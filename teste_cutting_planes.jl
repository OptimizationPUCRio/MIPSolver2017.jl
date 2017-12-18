using JuMP, Gurobi
solver = GurobiSolver()
model = Model(solver = solver)

@variable(model, y >= 0, Int)
@variable(model, x >= 0)

@constraints(model, begin
7y - 2x <= 14
x <= 3
2y - 2x <= 3
end)

@objective(model, :Max, 4y - x)

solve(model)


model = cutting_planes(model, [1], 50)
