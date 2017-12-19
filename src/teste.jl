include("branch_and_bound.jl")
include("funcoes_relax.jl")


test_MIP_medio_Guilherme(solveMIP,GurobiSolver(OutputFlag=0))

function test_MIP_medio_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP médio Guilherme (TSP 15 cidades)" begin
        number_of_nodes = 15
        srand(12)
        C = 1000*rand(15,15)

        @variable(m, X[i=1:number_of_nodes,j=1:number_of_nodes], Bin)
        @variable(m, u[i=:1:number_of_nodes])#Int
        for i=1:number_of_nodes
            @constraint(m,sum(X[i,j] for j=1:number_of_nodes if j!=i) == 1 )
        end
        for j=1:number_of_nodes
            @constraint(m,sum(X[i,j] for i=1:number_of_nodes if i!=j) == 1 )
        end
        for i=2:number_of_nodes
            for j=2:number_of_nodes
              if (i!=j)
                @constraint(m, u[i] - u[j] + number_of_nodes*X[i,j] <= number_of_nodes-1)
              end
            end
        end
        @objective(m, Min, sum(C[i,j]*X[i,j] for i=1:number_of_nodes, j=1:number_of_nodes))

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 2007.2884583133053 atol = 1e-7
    end
    setoutputs!(m,solution,testresult)
    return solution
end
