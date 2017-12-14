include("/Users/guilhermebodin/Documents/PUC/MIPTests_Branch_Guilherme/runmiptests.jl")
include("/Users/guilhermebodin/Documents/PUC/MIPSolver2017_04.jl/src/branch_and_bound.jl")

runtests(SolveMIP, GurobiSolver())
