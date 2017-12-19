include("branch_and_bound.jl")
include("funcoes_relax.jl")

#Com Corte

H_71 = [0]
H_72 = [0]
H_20 = [0]
H_25 = [0]
H_30 = [0]
H_40 = [0]
H_50 = [0]
H_60 = [0]
H_70 = [0]
H_80 = [0]
H_90 = [0]
H_100 = [0]

for i = 1:5
    resp = test_TSPmip7_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_71 = [H_71;resp.time]
end

for i = 1:5
    resp = test_TSPbin7_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_72 = [H_72;resp.time]
end

for i = 1:5
    resp = test_TSPmip20_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_20 = [H_20;resp.time]
end

for i = 1:5
    resp = test_TSPmip25_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_25 = [H_25;resp.time]
end

for i = 1:5
    resp = test_TSPmip30_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_30 = [H_30;resp.time]
end

for i = 1:5
    resp = test_TSPmip40_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_40 = [H_40;resp.time]
end

for i = 1:5
    resp = test_TSPmip50_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_50 = [H_50;resp.time]
end

for i = 1:5
    resp = test_TSPmip60_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_60 = [H_60;resp.time]
end

for i = 1:5
    resp = test_TSPmip70_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_70 = [H_70;resp.time]
end

for i = 1:5
    resp = test_TSPmip80_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_80 = [H_80;resp.time]
end

for i = 1:5
    resp = test_TSPmip90_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_90 = [H_90;resp.time]
end

for i = 1:5
    resp = test_TSPmip100_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))
    H_100 = [H_100;resp.time]
end
