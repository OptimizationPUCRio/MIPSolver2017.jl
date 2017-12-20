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
test_MIP_Grande_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
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
