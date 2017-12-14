include("/Users/guilhermebodin/Documents/PUC/MIPTests.jl/runmiptests.jl")
include("/Users/guilhermebodin/Documents/PUC/MIPSolver2017.jl/src/branch_and_bound.jl")
using Gurobi


runtests(solveMIP, GurobiSolver(OutputFlag = 0))
