include("../../MIPTests.jl/miptests.jl")
include("../../MIPTests.jl/runmiptests.jl")
include("../src\\branch_and_bound.jl")
using Gurobi


test_TSPmip7_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))




                                    test_MIP_medio_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))




test_TSPmip20_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))




                                    test_TSPmip25_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))




test_TSPmip30_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))




                                    test_TSPmip40_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))



test_TSPmip50_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))
test_MIP_Grande_Guilherme(solveMIP, GurobiSolver(OutputFlag = 0))



runtests(solveMIP, GurobiSolver(OutputFlag = 0))
