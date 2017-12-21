workspace()
include("C:/Users/Bianca/Desktop/otimização/Prog Inteira/MIPSolver2017.jl/test/miptests.jl")
include("C:/Users/Bianca/Desktop/otimização/Prog Inteira/MIPSolver2017.jl/src/branch_and_bound.jl")
using Gurobi
using JuMP

test_TSPmip7_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_MIP_medio_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_TSPmip20_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_TSPmip25_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_TSPmip30_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_TSPmip40_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_TSPmip50_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
include("C:/Users/Bianca/Desktop/otimização/Prog Inteira/MIPSolver2017.jl/src/branch_and_bound.jl")
test_MIP_Grande_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
include("C:/Users/Bianca/Desktop/otimização/Prog Inteira/MIPSolver2017.jl/src/branch_and_bound.jl")
testRobustCCUC(solveMIP,GurobiSolver())


include("C:/Users/Bianca/Desktop/otimização/Prog Inteira/MIPSolver2017.jl/src/branch_and_bound.jl")
test_TSPmip7_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_MIP_medio_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_TSPmip20_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_TSPmip25_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_TSPmip30_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_TSPmip40_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_TSPmip50_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_MIP_Grande_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
testRobustCCUC(solveMIP,GurobiSolver())


function test_TSPmip50_Guilherme()
  m = Model(solver = GurobiSolver())
  number_of_nodes = 50
  srand(123)
  C = 1000*rand(number_of_nodes,number_of_nodes)

  @variable(m, X[i=1:number_of_nodes,j=1:number_of_nodes], Bin)
  @variable(m, u[i=1:number_of_nodes])
  @constraint(m, u[1] == 1)
  for i=1:number_of_nodes
      @constraint(m,sum(X[i,j] for j=1:number_of_nodes if j!=i) == 1 )
  end
  for j=1:number_of_nodes
      @constraint(m,sum(X[i,j] for i=1:number_of_nodes if i!=j) == 1 )
  end
  for i=2:number_of_nodes
      @constraint(m, u[i] <= number_of_nodes)
      @constraint(m, u[i] >= 2)
      for j=2:number_of_nodes
          @constraint(m, u[i] - u[j] + 1 <= number_of_nodes*(1 - X[i,j]))
      end
  end
  @objective(m, Min, sum(C[i,j]*X[i,j] for i=1:number_of_nodes, j=1:number_of_nodes))
  return m
end

mod = test_TSPmip50_Guilherme()

solveMIP(mod)
