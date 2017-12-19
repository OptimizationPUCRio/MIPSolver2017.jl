using JuMP, Gurobi
using MathProgBase
using Base.Test
using Combinatorics

solver = GurobiSolver()

function cutting_planes(model::JuMP.Model, VecBin::Vector{Int}, MaxIter::Int64 = 5)
    tol = 1e-5

    Astd, b, cstd, xlb, xub, solver, m, n, flag_sense = extract_data(model)

    iter = 1
    convergence = -2 #not finished
    X = 0
    z = 0

    while convergence == -2
        t_1 , t_2, t_3 = solve_LP(Astd, b, cstd, xlb, xub, solver)

        if isnan(t_2)
            break
        end
        X, z, status = t_1 , t_2, t_3


        # Finding fractional variable address
        j = find(min.(X[1:n]-floor.(X[1:n]),ceil.(X[1:n])-X[1:n]).>tol)
        # Removing index from continous variables
        tmp = min(size(VecBin)[1], size(j)[1])
        tmp_ind = []
        for i in 1:tmp
            if size(VecBin)[1] <= size(j)[1]
                if VecBin[i] in j
                    push!(tmp_ind,VecBin[i])
                end
            else
                if j[i] in VecBin
                    push!(tmp_ind,j[i])
                end
            end
        end
        j = tmp_ind

        #Test convergence
        if length(j) == 0

            convergence = 1
        elseif iter >= MaxIter
            j = j[1]
            convergence  = -1
        else
            j = j[1] #pick a fractional and binary variable address (the first one)

            NBas = find(X .<= tol)
            Bas = deleteat!(collect(1:n), NBas)

            # u = Linha do dicionario correspondente a variavel x_j
            u = find(Bas.== j)

            DB = X[Bas]
            DN = ((Astd[:,Bas]'*Astd[:,Bas])\Astd[:,Bas]')*Astd[:,NBas]

            a0 = X[Bas]
            a = zeros(size(Bas)[1],n)

            a[:,NBas] = DN
            for i = 1:size(Bas)[1]
              a[i,Bas[i]] = 1
            end
            a[find(abs.(a) .< tol)] = 0

            N1 = intersect(deleteat!(copy(VecBin),find(VecBin .== j)),NBas)
            N2 = intersect(deleteat!(collect(1:n), VecBin),NBas)

            f0 = a0[u] - floor.(a0[u])
            f0 = f0[1]

            f = a[u,union(N1,N2)] - floor.(a[u,union(N1,N2)])
            A_cut = zeros(1,n)

            A_cut[j] = 1
            if size(intersect(find( f .<= f0),N1))[1] != 0
                A_cut[intersect(find( f .<= f0),N1)] = floor.(a[u,intersect(find( f .<= f0),N1)])
            end
            if size(intersect(find( f .>  f0),N1))[1] != 0
                A_cut[intersect(find( f .>  f0),N1)] = floor.(a[u,intersect(find( f .>  f0),N1)]) + (f[intersect(find( f .>  f0),N1)]' - f0)/(1 - f0)
            end
            if size(intersect(find(a[u,:] .< 0), N2))[1] != 0
                A_cut[intersect(find(a[u,:] .< 0), N2)] = a[intersect(find(a[u,:] .< 0), N2)]*(1/(1-f0))
            end


            b_cut = floor.(a0[u])

            Astd = [[Astd zeros(m,1)]; [A_cut 1]]
            b = [b ; b_cut]

            cstd = [cstd ; 0]
            xlb = [xlb; 0]
            xub = [xub; +Inf]
            m+=1
            n+=1
            iter+=1

        end
    end
    #model_F = compose_model(Astd, b, cstd, xlb, xub, flag_sense, solver, X, z)

    solve(model)

    if isnan(z) | (z > getobjectivevalue(model))
        flag = -1
    else
        flag = 1
    end
    if flag_sense == 1
        z = -z
    end

    return z, flag
end

function extract_data(model::JuMP.Model)

    A = full(JuMP.prepConstrMatrix(model))
    c = JuMP.prepAffObjective(model)
    m, n = size(A)

    xlb = copy(model.colLower)
    xub = copy(model.colUpper)

    rowlb, rowub = JuMP.prepConstrBounds(model)

    folgas = find(rowlb .!= rowub)
    s_f = size(folgas)[1]

    Astd = [A zeros(m,s_f)]

    for i = 1:size(folgas)[1]
        if rowlb[folgas[i]] == -Inf
            Astd[folgas[i],n+i] += 1
        else
            Astd[folgas[i],n+i] -= 1
        end
    end

    cstd = [c;zeros(s_f)]

    xlb = copy(model.colLower)
    xub = copy(model.colUpper)

    nxlb = zeros(n+s_f,1)
    nxub = zeros(n+s_f,1)

    nxlb[1:size(xlb)[1]] = xlb
    nxlb[size(xlb)[1]+1 : s_f+n] = 0

    nxub[1:size(xlb)[1]] = xub
    nxub[size(xub)[1]+1 : s_f+n] = +Inf

    b = zeros(m,1)
    for i in 1:m
        if rowlb[i] == -Inf
            b[i] = rowub[i]
        else
            b[i] = rowlb[i]
        end
    end

    flag_sense = 0

    if model.objSense == :Min
    cstd = -cstd
    flag_sense = 1
    end

    m, n = size(Astd)
    solver = model.solver

    return Astd, b, cstd, nxlb, nxub, solver, m, n, flag_sense
end

function solve_LP(Astd, b, cstd, xlb, xub, solver)
    m, n = size(Astd)
    model_LP = Model(solver=GurobiSolver(OutputFlag = 0))

    @variable(model_LP, xlb[i] <= x[i=1:n] <= xub[i])
    @constraints(model_LP, begin
    constrain[i=1:m], sum(Astd[i,j]*x[j] for j=1:n) == b[i]
    end)

    @objective(model_LP, :Max, sum(cstd[j]*x[j] for j=1:n))
    status = solve(model_LP);
    return model_LP.colVal, model_LP.objVal, status
end

function compose_model(Astd, b, cstd, xlb, xub, flag_sense, solver, X, z)

    m, n = size(Astd)

    model_F = Model(solver=solver)

    @variable(model_F, xlb[i] <= x[i=1:n] <= xub[i])
    @constraints(model_F, begin
    constrain[i=1:m], sum(Astd[i,j]*x[j] for j=1:n) == b[i]
    end)

    if flag_sense == 1
        cstd = -cstd
        @objective(model_F, :Min, sum(cstd[j]*x[j] for j=1:n))
        model_F.objVal = -z
    else
        @objective(model_F, :Max, sum(cstd[j]*x[j] for j=1:n))
        model_F.objVal = z
    end

    model_F.colVal = X

    return model_F
end

"""
 Function Name   : Grasp generic mip problems
 Description     : "Greedy ramdomized adaptive search procedure" porpouse is
                    to find a primary solution to a
                    Linear Integer programim problem.
                    This inplementation was optimized to solve the
                    travelling salesman problem (TSP).

"""
function _grasp_tsp(m::JuMP.Model,maxiter::Int = 200, α::Float64 = 0.25, max_iter_improv::Int= 200,max_improvs_aux::Int=200)
    num_variables = sum(m.colCat .== :Bin)
    num_nodes = Int(sqrt(num_variables))
    num_rlc = floor(num_nodes*α)
    num_rlc = Int(min(max(1,num_rlc),num_nodes))
    c = JuMP.prepAffObjective(m)
    range_nodes = 1:num_nodes
    indVariable = collect(1:m.numCols)[m.colCat .== :Bin]

    tour = fill(-1, num_nodes)

    ## first solution ##
    # restrictions
    for i = Int.(range_nodes)
        if m.objSense == :Min
            c[i + num_nodes*(i-1)] = maximum(c[1 + num_nodes*(i-1):num_nodes*(i)])  # so it wont be selected as best cadidate
        end
    end
    for i = num_variables
        if m.colUpper[i] == 0
            row = Int(floor(i/num_nodes))
            col = Int(i - num_nodes*row)
            if m.objSense == :Min
                c[row*num_nodes+col] = maximum(c[1 + num_nodes*row:num_nodes*(row+1)]) # so it wont be selected as best cadidate
            end
        end
    end

    # rlc
    # c_aux = c[1 + num_nodes*(i-1):num_nodes*(i)]
    # sortperm(c[1 + num_nodes*(i-1):num_nodes*(i)], rev=true)
    ind_aux = collect(range_nodes)
    # ind = m.objSense == :Max ? ind_aux[c_aux .>= (1-α)*maximum(c_aux)] : ind_aux[c_aux .<= (α+1)*minimum(c_aux)]

    # first iteraction
    itr = 1
    i = ind_aux[rand(1:size(ind_aux,1))]

    tour[itr] = i

    # remaining iteractions
    proh = Int.(-1*ones(num_nodes,num_nodes))
    itr+=1
    while itr <=num_nodes-1
        rlc = m.objSense == :Max ? setdiff(setdiff(sortperm(c[1 + num_nodes*(i-1):num_nodes*(i)], rev=true),tour), proh[itr,:]):
            setdiff(setdiff(sortperm(c[1 + num_nodes*(i-1):num_nodes*(i)]),tour), proh[itr,:])
        rlc = rlc[1:min(num_rlc,size(rlc,1))]
        if size(rlc,1) == 0
            itr-=1
            proh[itr,tour[itr]] = tour[itr]
            proh[itr+1,:] .=  Int.(-1*ones(num_nodes))
            tour[itr] = -1
        else
            ind = rlc[rand(1:size(rlc,1))]
            i = ind
            tour[itr] = i
            if size(setdiff(collect(1:num_nodes),tour),1) == 0
                proh[itr,tour[itr]] = tour[itr]
                proh[itr+1,:] .=  Int.(-1*ones(num_nodes))
                tour[itr] = -1
                itr-=1
            else
                itr+=1
            end
        end
    end
    i = setdiff(collect(1:num_nodes),tour)[1]
    tour[num_nodes] = i

    # check feasibility
    sol = feasible_solution._prep_solution(m,tour)
    # check feasability

    status,obval = feasible_solution._check_feasability_solution(m,sol,indVariable)
    if status == false
        # repair
    end
    ## end first solution ##
    ## ajustments solution ##
    # Initialization
    m_aux = deepcopy(m)
    max_improvs = max(num_nodes^2,max_improvs_aux)
    best_sol = deepcopy(m.colVal)
    best_obj = deepcopy(obval)
    itr_improv = 0
    itr = 1
    neighb_param1 = 2
    neighb_param2 = 1  # shift

    neighb = collect(combinations(1:num_nodes,neighb_param1))
    while itr_improv <= max_iter_improv
        # permutation
        inds = neighb[itr] #rand(1:num_nodes,2) collect(1:5)
        tour_aux = deepcopy(tour)
        tour_aux2 = deepcopy(tour[inds])
        tour_aux[inds] = circshift(tour_aux2,neighb_param2)

        # tour_aux[inds[1]] = tour[inds[2]]
        # tour_aux[inds[2]] = tour[inds[1]]

        sol = feasible_solution._prep_solution(m_aux,tour_aux)
        status,obval = feasible_solution._check_feasability_solution(m_aux,sol,indVariable)
        if status == true
            if m.objSense == :Max && obval > best_obj
                best_sol = deepcopy(m_aux.colVal)
                best_obj = deepcopy(obval)
            elseif  m.objSense == :Min && obval < best_obj
                best_sol = deepcopy(m_aux.colVal)
                best_obj = deepcopy(obval)
            end
        end
        if itr == size(neighb,1)
            if neighb_param1 == num_nodes
                neighb_param1 = 1
                neighb_param2 +=1
            end
            itr = 0
            neighb_param1 +=1
            neighb = collect(combinations(1:num_nodes,neighb_param1))
        end
        itr+=1
        itr_improv+=1
    end

    m.colVal = best_sol
    m.objVal = best_obj
    # return status
    return :SubOptimal
end

"""
Function Name   : prep solution
Description     : constructs a solution from a tour.
"""
function _prep_solution(m::JuMP.Model,tour)
    num_variables = sum(m.colCat .== :Bin)
    num_nodes = Int(sqrt(num_variables))

    # prep solution
    sol = zeros(num_variables)
    indVariable = collect(1:m.numCols)[m.colCat .== :Bin]
    row = tour[num_nodes]
    col = tour[1]
    sol[(row-1)*num_nodes+col] = 1
    for i = 1:num_nodes-1
        row = tour[i]
        col = tour[i+1]
        sol[(row-1)*num_nodes+col] = 1
    end
    return sol
end


"""
 Function Name   : Grasp
 Description     : "Greedy ramdomized adaptive search procedure" porpouse is
                    to find a primary solution to a
                    Linear Integer programim problem.

"""
function grasp(m::JuMP.Model,flag_tsp::Bool = true, maxiter::Int = 200, α::Float64 = 0.25, max_iter_improv::Int= 200,max_improvs_aux::Int=200)
    if flag_tsp
        status = _grasp_tsp(m::JuMP.Model,maxiter, α)
    else
        return :NotImplemented
    end
end

"""
 Function Name   : Grasp generic mip problems
 Description     : "Greedy ramdomized adaptive search procedure" porpouse is
                    to find a primary solution to a
                    Linear Integer programim problem.
                    This inplementation was designed to solve any mip problem.
                    (not problem optimized consequently).

"""
function _grasp_generic_mip_problems(model::JuMP.Model,maxiter::Int = 200)
    ## First solution ##
    # Initialization #
    c = []
    variables = collect(1:model.numCols)
    indexVariables = variables[model.colCat .== :Bin]  # list of candidates
    i = 1

    # Fill solution "c" #
    while i < size(indexVariables,1) && i < maxiter
        # Initialization of restricted list of candidates (rlc)
        rlc = []
        for i in size(indexVariables,1)
            # load candidate
            candidate = deepcopy(model)
            candidate.colUpper[indVariable[i]] = 1
            candidate.colLower[indVariable[i]] = 1
            # check feasability
            if _check_feasability_variable(candidate)

            end

        end



        i += 1
    end

end

"""
Function Name   :
Description     :
"""
function _rlc_criteria(model::JuMP.Model)
  indexVariables = collect(1:model.numCols)
  selection = [model.colCat .== :Bin  model.colVal .!= 0  model.colVal .!= 1]
  for row = 1 : size(selection, 1)
    if prod(selection[row, :])
      return indexVariables[row]
    end
  end
  return 0
end

"""
Function Name   : selectInfeasableVariable
Description     : Searches for a varaible that should be binary but has a none
                  binary value in the current model solution.
"""
function _selectInfeasableVariable(model::JuMP.Model)
  indexVariables = collect(1:model.numCols)
  selection = [model.colCat .== :Bin  model.colVal .!= 0  model.colVal .!= 1]
  for row = 1 : size(selection, 1)
    if prod(selection[row, :])
      return indexVariables[row]
    end
  end
  return 0
end

"""
Function Name   : check feasability solution
Description     : Checks if the solution "c" passed as parameter is a feasible
                  solution for the given problem model when linearly relaxed.
"""
function _check_feasability_solution(m::JuMP.Model,c,indVariable)
    for i = 1:size(c,1)
        m.colUpper[indVariable[i]] = c[i]
        m.colLower[indVariable[i]] = c[i]
    end
    obj = m.obj
    obf = JuMP.prepAffObjective(m)
    @objective(m, m.objSense, 0)
    status = solve(m, relaxation=true)
    obval = obf'*m.colVal
    @objective(m, m.objSense, obj)
    if status == :Infeasible
        return false,obval
    end
    return true,obval
end

"""
Function Name   : check feasability variable
Description     : Checks if model has a linearly relaxed feasible solution
"""
function _check_feasability_variable(model::JuMP.Model)
    obj = model.obj
    @objective(m, model.objSense, 0)
    status = solve(model, relaxation=true)
    @objective(m, model.objSense, obj)
    if status == :Infeasible
        return false
    end
    return true
end

function fpump(mod::JuMP.Model, binI::Vector{Int64})
  m=deepcopy(mod)
  obf = JuMP.prepAffObjective(m)
  tam=m.numCols
  n=length(binI)

  v = Variable.(m, 1:tam)

  for i in binI
    @constraint(m,0 <= v[i])
    @constraint(m, v[i] <= 1)
  end

  # conferir se o b&b da solve no modelo antes, se n der precisa fazer aqui
  solve(m)
  xvia=m.colVal

  xint=roundbin(xvia,binI,tam)

  d=dist(xvia,xint,binI)

  cont=0
  while d>1e-7 && cont < 1e3
    cont = cont + 1

    @objective(m, Min, sum{ ifelse(xint[i] == 0  , v[i] , 0) + ifelse(xint[i] == 1  , xint[i] - v[i], 0), i in binI})

    solve(m)
    xvia = m.colVal
    d1 = m.objVal

    if d1 ≈ d
      #se o xint for igual ao anterior a distância permanecerá a mesma então criamos uma perturbação

      xint=mudaround(xvia,xint,binI,tam,n)
      d = dist(xvia,xint,binI)

    else
      xint=roundbin(xvia,binI,tam)
      d=d1
    end
  end

  objval=obf'*xvia

  if d <= 1e-7
    return xvia, objval
  end
  return false, false

end


function dist(x,xint,binI)
  dist=0
  for i  in binI
    if xint[i] == 0
      soma = x[i] - 0
    elseif xint[i] == 1
      soma = 1 - x[i]
    end
    dist=dist+soma
  end
  return dist
end

function roundbin(x::Vector, binind::Vector{Int64}, tam::Int64)
  xint=zeros(tam)
  for i  in 1:tam
    if i in binind
      xint[i] = round(x[i])
    else
      xint[i] = x[i]
    end
  end
  return xint
end

function mudaround(x::Vector, xint::Vector, binind::Vector{Int64}, tam::Int64, n::Int64)

  ind=indChange(x,xint,binind,tam,n)

  for i  in ind
    if xint[i] == 1
      xint[i] = 0
    else
      xint[i] = 1
    end
  end
  return xint
end


#essa função nos da os indices aleatorios das variaveis com maiores scores
function indChange(x::Vector, xint::Vector, binind::Vector{Int64}, tam::Int64, n::Int64)
  vet=zeros(n)

  #preenche um vetor de scores (dif entre xvia e xint das variaveis binarias)
  j=1;
  for i in binind
    vet[j] = abs(x[i] - xint[i])
    j=j+1
  end

  #p é um vetor que define a permutação dos indices de vet que o coloca em ordem decrescente
  p=sortperm(vet,rev=true)
  #como o util para mim são os índices das variaveis originais organizamos binind dessa forma
  novob=binind[p]

  #quero mudar apenas uma certa quantidade de indices, quanto menos melhor
  quant=Int.(ceil(n/400))

  if n>400
    div=0.01
  else
    div=0.1
  end

  #quero que os ind sejam das variaveis com maiores scores (uma certa proporção dos primeiros)

  percent=Int.(ceil(div*n))

  j=randperm(percent)[1:quant]
  ind=Int.(zeros(quant))
  for k in 1:quant
    ind[k]=novob[j[k]]
  end
  return ind
end

# Código original do Raphael
mutable struct node
  level::Int
  model::JuMP.Model
end

## Checks if solution is binary
function isBinary(model::JuMP.Model, binaryIndices::Vector{Int64})
  isZero = model.colVal[binaryIndices] .== 0
  isOne  = model.colVal[binaryIndices] .== 1
  if sum(isZero + isOne) == length(binaryIndices)
    return true
  end
  return false
end

## Checks if model is max: if it is, converts to min
function convertSense!(m::JuMP.Model)
  if m.objSense == :Max
    m.objSense = :Min
    m.obj = -m.obj
  else
    m.objSense = :Max
    m.obj = -m.obj
  end
end

## Strong branching -- if amount of branches == Inf --> full strong branching
function strong(currentNode::node, binaryIndices::Vector{Int64}, amountOfBranches::Int64,
                pscMatrix::Matrix{Float64})

  isZero = currentNode.model.colVal[binaryIndices] .== 0
  isOne  = currentNode.model.colVal[binaryIndices] .== 1
  indFracIndices = find(isZero .+ isOne .== 0)
  fracIndices = binaryIndices[indFracIndices]

  n = min(amountOfBranches, length(fracIndices))

  indCandidatesToBranch = randperm(length(fracIndices))[1:n]
  candidatesToBranch = fracIndices[indCandidatesToBranch]
  bounds = Array{Float64}(2,n)
  for i = 1:n
    leftModel = deepcopy(currentNode.model)
    leftModel.colUpper[candidatesToBranch[i]] = 0
    leftModel.colLower[candidatesToBranch[i]] = 0

    rightModel = deepcopy(currentNode.model)
    rightModel.colUpper[candidatesToBranch[i]] = 1
    rightModel.colLower[candidatesToBranch[i]] = 1

    solve(leftModel)
    solve(rightModel)
    bounds[1,i] = leftModel.objVal
    bounds[2,i] = rightModel.objVal

  end

  indBestFrac = ceil(Int,indmax(bounds)/2)
  indToBranch = candidatesToBranch[indBestFrac]

  leftModel = deepcopy(currentNode.model)
  leftModel.colUpper[indToBranch] = 0
  leftModel.colLower[indToBranch] = 0

  rightModel = deepcopy(currentNode.model)
  rightModel.colUpper[indToBranch] = 1
  rightModel.colLower[indToBranch] = 1

  leftDelta = abs(currentNode.model.objVal - leftModel.objVal)/currentNode.model.colVal[indToBranch]
  rightDelta = abs(currentNode.model.objVal - rightModel.objVal)/(1-currentNode.model.colVal[indToBranch])

  pscMatrix[indBestFrac, 2]+=1
  pscMatrix[indBestFrac, 1] = (pscMatrix[indBestFrac,1] + ((leftDelta + rightDelta)/2))/pscMatrix[indBestFrac,2]

  leftChild = node(currentNode.level+1, leftModel)
  rightChild = node(currentNode.level+1, rightModel)

  return leftChild, rightChild, pscMatrix
end

## Most fractional branching
function fractional(currentNode::node, binaryIndices::Vector{Int64})

  distance = abs.(currentNode.model.colVal[binaryIndices] - 0.5)
  indFrac = indmin(distance)
  indToSet = binaryIndices[indFrac]

  leftModel = deepcopy(currentNode.model)
  leftModel.colUpper[indToSet] = 0
  leftModel.colLower[indToSet] = 0

  rightModel = deepcopy(currentNode.model)
  rightModel.colUpper[indToSet] = 1
  rightModel.colLower[indToSet] = 1

  leftChild = node(currentNode.level+1, leftModel)
  rightChild = node(currentNode.level+1, rightModel)

  return leftChild, rightChild
end

## Pseudo-cost branching
function pseudocost(currentNode::node, binaryIndices::Vector{Int64}, pscMatrix::Matrix{Float64})

    isZero = currentNode.model.colVal[binaryIndices] .== 0
    isOne  = currentNode.model.colVal[binaryIndices] .== 1
    indFracIndices = find(isZero .+ isOne .== 0)
    fracIndices = binaryIndices[indFracIndices]

    indMaxPSC = indmax(pscMatrix[indFracIndices,1].*(1-currentNode.model.colVal[indFracIndices]))
    indToBranch = fracIndices[indMaxPSC]

    leftModel = deepcopy(currentNode.model)
    leftModel.colUpper[indToBranch] = 0
    leftModel.colLower[indToBranch] = 0

    rightModel = deepcopy(currentNode.model)
    rightModel.colUpper[indToBranch] = 1
    rightModel.colLower[indToBranch] = 1

    leftDelta = abs(currentNode.model.objVal - leftModel.objVal)/currentNode.model.colVal[indToBranch]
    rightDelta = abs(currentNode.model.objVal - rightModel.objVal)/(1-currentNode.model.colVal[indToBranch])

    pscMatrix[indMaxPSC, 2]+=1
    pscMatrix[indMaxPSC, 1] = (pscMatrix[indMaxPSC,1] + ((leftDelta + rightDelta)/2))/pscMatrix[indMaxPSC, 2]

    leftChild = node(currentNode.level+1, leftModel)
    rightChild = node(currentNode.level+1, rightModel)

    return leftChild, rightChild, pscMatrix
end

## Receives node and creates two children by setting a variable to 0 and 1 respectively
function branch(currentNode::node, binaryIndices::Vector{Int64}, method::Symbol;
                amountOfBranches::Int64 = 10^10, pscMatrix::Matrix{Float64} = ones(1,1))

  if method == :fractional
    leftChild, rightChild = fractional(currentNode, binaryIndices)
  elseif method == :strong
    leftChild, rightChild, pscMatrix = strong(currentNode, binaryIndices, amountOfBranches, pscMatrix)
  elseif method == :pseudocost
    leftChild, rightChild, pscMatrix = pseudocost(currentNode, binaryIndices, pscMatrix)
  elseif method == :hybrid
    if currentNode.level <= 3
      leftChild, rightChild, pscMatrix = strong(currentNode, binaryIndices, 10^10, pscMatrix)
    else
      leftChild, rightChild, pscMatrix = pseudocost(currentNode, binaryIndices, pscMatrix)
    end
  else
    println("Error on branching method defition")
  end

  return leftChild, rightChild, pscMatrix
end

function obtainBoundList(nodeList)
  boundList = Array{Float64}(length(nodeList))
  for i = 1 : length(nodeList)
    boundList[i] = nodeList[i].model.objVal
  end
  return boundList
end

## Receives a mixed binary linear JuMP model
function solveMIP(m::JuMP.Model; branchMethod = :strong, traverseMethod = :mixed, boolgrasptsp::Bool = false, cp  = false)

  tic()

  # Check if model is max; if it is, converts to min
  flagConverted = 0
  if m.objSense == :Max
    convertSense!(m)
    flagConverted = 1
  end

  # Create vector of indices of the binary variables
  binaryIndices = find(m.colCat .== :Bin)
  binarySolutions = 0

  # Best bounds: start as Inf
  if cp
      cp_lb , flag = cutting_planes(m, binaryIndices, 15)
      print("Bound cutting_planes = ", cp_lb, " e Flag = ", flag)
      if flag == -1
          bestBound = -1e200
      else
          bestBound = cp_lb
      end
  else
      bestBound = -1e200
  end

  bestVal = 1e200

  flagOpt = 0 # flag that indicates if a viable solution has been found

  # Solve using initial heuristic (GRASP)
  if boolgrasptsp
    m_aux = deepcopy(m)
    status = feasible_solution.grasp(m_aux,true,200,0.03,10000000)
    if status == :SubOptimal && m_aux.objVal < bestVal
      bestVal = m_aux.objVal
      m.colVal = m_aux.colVal
      binarySolutions += 1
      flagOpt+=1
    end
  end

  # Solve linear relaxation
  m.colCat[:] = :Cont
  status = solve(m)
  nodes = Vector{node}(0)
  if status == :Optimal && isBinary(m, binaryIndices)
    # Solution of the relaxed problem is binary: optimal solution
    if !cp
        bestBound = m.objVal
    end
    bestVal = m.objVal
    binarySolutions = 1
  else
    push!(nodes, node(0, m)) # Add root to branch and bound tree
    lastNodeLevel = 0
  end

  iter = 1 # number of visited nodes
  branched = false # flag that indicates if a branch has occurred in this iteration
  tol = 0.01 # tolerance (%)

  if traverseMethod == :depth
    traverse = -1
  else
    traverse = 1
  end

  # Initializing pseudo-cost matrix
  pscMatrix = [ones(length(binaryIndices)) zeros(length(binaryIndices))]

  time0 = time_ns()

  while !isempty(nodes) && abs((bestVal - bestBound)/bestVal) > tol && (time_ns()-time0)/1e9 < 1800

    # Change traverse method every 10 iterations for better bound discovery
    if traverseMethod == :mixed && iter%10 == 0
      traverse = (-1)*traverse
    end

    branched = false
    # Check node lower bound. If greater than current best UB, prune by limit
    if nodes[1].model.objVal <= bestVal
      status = solve(nodes[1].model)
      if status == :Optimal
        boundList = obtainBoundList(nodes)
        bestBound = minimum(boundList)
        if isBinary(nodes[1].model, binaryIndices)
          # Relaxed solution is binary: optimal solution -- don't branch
          if nodes[1].model.objVal < bestVal
            bestVal = nodes[1].model.objVal
            m.colVal = copy(nodes[1].model.colVal)
            flagOpt = 1
            binarySolutions+=1
          end
        elseif nodes[1].model.objVal <= bestVal
          # Relaxed solution is not binary and should not be pruned by limit -- branch
          leftChild, rightChild, pscMatrix = branch(nodes[1], binaryIndices, branchMethod; pscMatrix = pscMatrix)
          branched = true
        end
      end
    end

    # Check if branched -- if true, add children to the list
    if branched == true
      if traverse == 1 # breadth -- insert children at the end of the list
        push!(nodes, leftChild, rightChild)
        lastNodeLevel = nodes[1].level
        deleteat!(nodes, 1)
      else # depth -- insert children at the beginning of the list
        lastNodeLevel = nodes[1].level
        deleteat!(nodes, 1)
        unshift!(nodes, leftChild, rightChild)
      end
    else
      deleteat!(nodes, 1)
    end
    iter+=1

    if iter == 1 || iter%10 == 0
      println("UB: $bestVal")
      println("LB: $bestBound")
    end
  end

  if flagConverted == 1
    convertSense!(m)
    m.objVal = -bestVal
    m.objBound = -bestBound
  else
    m.objVal = bestVal
    m.objBound = bestBound
  end

  # Return binary variables to original state
  m.colCat[binaryIndices] = :Bin

  # Outputs
  m.ext[:status] = status
  if flagOpt == 1
    m.ext[:status] = :Optimal
  end
  m.ext[:nodes] = iter
  m.ext[:solutions] = binarySolutions
  t = toc()
  m.ext[:time] = t

  println(string("Nodes: ",m.ext[:nodes]))

  return m.ext[:status]
end

type MIPSolution
    name
    pass
    objective
    bestbound
    time
    nodes
    intsols
    status
    function MIPSolution()
        new("",false,NaN,NaN,Inf,-1,-1,:unsolved)
    end
end

function setoutputs!(m,sol::MIPSolution, test)
    sol.pass = length(test.results) == 0
    sol.objective = getobjectivevalue(m)
    sol.bestbound = m.objBound
    if haskey(m.ext,:time)
        if typeof(m.ext[:time]) <: Real
        sol.time = m.ext[:time]
        end
    end
    if haskey(m.ext,:nodes)
        if typeof(m.ext[:nodes]) <: Real
        sol.nodes = m.ext[:nodes]
        end
    end
    if haskey(m.ext,:node)
        if typeof(m.ext[:node]) <: Real
        sol.nodes = m.ext[:node]
        end
    end
    if haskey(m.ext,:intsols)
        if typeof(m.ext[:intsols]) <: Real
        sol.intsols = m.ext[:intsols]
        end
    end
    if haskey(m.ext,:solucao_inteira)
        if typeof(m.ext[:solucao_inteira]) <: Real
        sol.intsols = m.ext[:solucao_inteira]
        end
    end
    if haskey(m.ext,:solutions)
        if typeof(m.ext[:solutions]) <: Real
        sol.intsols = m.ext[:solutions]
        end
    end

    if haskey(m.ext,:status)
        if typeof(m.ext[:status]) <: Symbol
            sol.status = m.ext[:status]
        end
    end
    return nothing
end



#≈

# teste exemplo
# adicionado por Joaquim Garcia
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste examplo" begin
        @variable(m, x >=0)
        @objective(m, :Min, x)

        solveMIP(m)
        @test getobjectivevalue(m) == 0
        @test getvalue(x) == 0
    end
    setoutputs!(m,solution,testresult)
    return solution
end


#teste problema 1 da lista (mochila)
#adicionado por Guilherme Bodin
function test2(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste da Mochila" begin
        @variable(m, x[i=1:3], Bin)
        @constraint(m, 6*x[1] + 5*x[2] + 5*x[3] <= 10)
        @objective(m, Max, 6*x[1] + 4*x[2] + 3*x[3])

        solveMIP(m)
        @test getobjectivevalue(m) == 7
        @test getvalue(x) == [0, 1, 1]

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste Sudoku
# adicionado por Raphael Saavedra
function testSudoku(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste Sudoku" begin
        n = 9
        @variable(m, x[i in 1:n, j in 1:n, k in 1:n], Bin)

        fixas = [(1,3,4), (1,5,6), (1,9,2), (2,1,8), (2,3,5), (2,6,2), (2,8,3),
                (3,5,3), (3,8,6), (4,2,2), (4,3,8), (5,6,4), (6,1,7), (6,5,5),
                (6,9,9), (7,3,2), (7,6,1), (8,2,7), (8,5,4), (8,7,9), (8,8,5),
                (9,1,6), (9,8,4)]
        for idx in fixas
            @constraint(m, x[idx...] == 1)
        end
        @constraint(m, [j in 1:n, k in 1:n], sum(x[:,j,k]) == 1)
        @constraint(m, [i in 1:n, k in 1:n], sum(x[i,:,k]) == 1)
        @constraint(m, [i in 1:n, j in 1:n], sum(x[i,j,:]) == 1)
        @constraint(m, [p in [0,3,6], q in [0,3,6], k in 1:n], sum(sum(x[i+p,j+q,k] for i in 1:3) for j in 1:3) == 1)
        @objective(m, Min, 0)

        solveMIP(m)
        @test getobjectivevalue(m) == 0
        @test sum(getvalue(x)) == 81
        @test getvalue(x[1,1,1]) == 0
        @test getvalue(x[1,1,3]) == 1
        @test getvalue(x[8,1,1]) == 1
        @test getvalue(x[8,1,9]) == 0
        @test sum(getvalue(x[6,6,4:9])) == 0
        @test getvalue(x[6,6,3]) == 1
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste problema 6 da lista (Expansao da Producao)
#adicionado por Andrew Rosemberg
function test3(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    Cinv = 13.16
    M = 200
    m = Model(solver = solver)
    testresult = @testset "Teste da Expansao da Producao" begin
        @variable(m, x[i=1:2]>=0)
        @variable(m, u, Bin)
        @objective(m, Max, 4*x[1] + 3*x[2] - u*Cinv)

        @constraint(m, 2*x[1] + 1*x[2] <= 4 +u*M)
        @constraint(m, 1*x[1] + 2*x[2] <= 4 +u*M)

        @constraint(m, 1*x[1] + 0.1*x[2] <= 4 +(1-u)*M)
        @constraint(m, 0.4*x[1] + 1*x[2] <= 4 +(1-u)*M)

        sol = solveMIP(m)

        @test getobjectivevalue(m) ≈ 9.340000000000002 atol=1E-07
        @test getvalue(x) ≈ [3.75, 2.5] atol=1E-07
        @test getvalue(u) ≈ 1 atol=1E-07

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste mochila binária infeasible
# adicionado por Raphael Saavedra
function testInfeasibleKnapsack(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste Mochila infeasible" begin
        @variable(m, x[i=1:3], Bin)
        @constraint(m, 6*x[1] + 5*x[2] + 5*x[3] <= 5)
        @constraint(m, x[1] == 1)
        @objective(m, Max, 6*x[1] + 4*x[2] + 3*x[3])

        solveMIP(m)

        @test m.ext[:status] == :Infeasible
    end
    setoutputs!(m,solution,testresult)
    return solution
end


#teste problema da P1 (CVRP)
#adicionado por Eduardo Brito
function test_P1_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver=solver)
    testresult = @testset "Teste CVRP" begin
        tam = 6
        carros = 2
        c = [ 0.162329   0.934386    0.27999    0.771633  0.93058    0.757016
        0.168044   0.516489    0.561025   0.793173  0.628186   0.893351
        0.450754   0.130342    0.0550682  0.330183  0.0140276  0.5666
        0.405403   0.00683533  0.828927   0.361477  0.265888   0.437654
        0.164714   0.00988625  0.470629   0.941743  0.343984   0.906805
        0.0105587  0.825022    0.540088   0.840939  0.137927   0.637206] #Custo na rota i,j
        b = [0.3286854210610859
        0.8747455782649833
        0.2997874087044681
        0.012584530553862994
        0.5477965855045668
        0.5671227617656462] #Demandas nos vertices
        b[1] = -sum(b[2:tam])
        k = [15.2168   9.21175  46.9288   2.46427  34.4648  25.4713
        84.4671  84.7527   84.7259  30.5085   55.8006  48.8416
        59.8025  87.1999   62.5169  57.5668   47.5803  66.2564
        97.9814  53.4096   72.0769   1.71473  41.5536  80.0004
        30.9586  59.9041   65.278   88.7064   55.0077  43.9505
        46.67    64.3936   32.2418  82.8831   38.0806  68.6481] #Capacidade na rota i,j


        @variable(m,x[1:tam,1:tam] >= 0)
        @variable(m,y[1:tam,1:tam], Bin)

        @constraints(m,begin
        saida[i=2:tam],sum(y[i,j] for j = 1:tam if j!=i) == 1
        chegada[j=2:tam],sum(y[i,j] for i = 1:tam if j!=i) == 1
        s_depot, sum(y[i,1] for i = 2:tam) == carros
        c_depot, sum(y[1,j] for j = 2:tam) == carros
        end)

        @constraints(m,begin
        capacidade_rota[i=1:tam,j=1:tam], x[i,j] <= k[i,j]*y[i,j]
        capacidade_carro[i=1:tam], x[1,i] <= sum(b[2:tam])/carros + sqrt(var(b))
        demanda[i=1:tam], sum(x[j,i] for j =1:tam if j!= i) - sum(x[i,j] for j=1:tam if j!=i) == b[i]
        ciclos, sum(y[i,i] for i =1:tam) == 0
        end)

        @objective(m, Min, sum(sum(c[i,j]*y[i,j] for i = 1:tam) for j = 1:tam))

         resp_y = [0.0  -0.0   1.0   1.0  -0.0   0.0
        1.0   0.0  -0.0  -0.0  -0.0  -0.0
        -0.0  -0.0   0.0   0.0   1.0  -0.0
        -0.0   0.0  -0.0   0.0  -0.0   1.0
        0.0   1.0  -0.0  -0.0   0.0  -0.0
        1.0  -0.0  -0.0  -0.0  -0.0   0.0]

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 1.69179355 atol=exp10(-5)
        @test getvalue(y) ≈ resp_y atol=1e-3

    end
    setoutputs!(m,solution,testresult)
    return solution
end


# teste Problema da Producao (PL)
# adicionado por Eduardo Brito
function test_PL_Simples_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver=solver)
    testresult = @testset "Problema da Producao" begin
        @variable(m, x[1:2] >=0)
        @constraints(m, begin
        cons1, 2x[1] + x[2] <= 4
        cons2, x[1] + 2x[2] <= 4
        end)
        @objective(m, :Max, 4x[1] + 3x[2])

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 9.33334 atol = exp10(-5)
        @test getvalue(x) ≈ [1.3333334;1.3333334] atol = exp10(-5)

    end
    setoutputs!(m,solution,testresult)
    return solution
end


# teste Pl Infeasible
# adicionado por Eduardo Brito
function test_PL_Infeasible_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver=solver)
    testresult = @testset "Problema Infeasible" begin
        @variable(m, x[1:2] >=0)
        @constraints(m, begin
        cons1, 2x[1] + x[2] <= 4
        cons2, x[1] + 2x[2] <= 4
        cons3, x[1] + x[2] >= 5
        end)
        @objective(m, :Max, 4x[1] + 3x[2])

        solveMIP(m)

        @test m.ext[:status] == :Infeasible

    end
    setoutputs!(m,solution,testresult)
    return solution
end


# teste Pl Unbounded
# adicionado por Eduardo Brito
function test_PL_Unbounded_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver=solver)
    testresult = @testset "Problema Unbounded" begin
        @variable(m, x[1:2] >=0)
        @objective(m, :Max, 4x[1] + 3x[2])

        solveMIP(m)
        @test m.ext[:status] in [:Unbounded, :InfeasibleOrUnbounded]

    end
    setoutputs!(m,solution,testresult)
    return solution
end


# teste MIP (minimal ~5 binarias)
# adicionado por Eduardo Brito
function test_MIP_Minimal_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver=solver)
    testresult = @testset "MIP Minimal" begin
        @variable(m, x[1:5] >=0, Bin)
        @variable(m, y[1:5] >= 0)
        @constraints(m, begin
        cons1, sum(x) <= 4.5
        cons2, y[1] <= 10(x[1])
        cons3, y[2] <= 10(x[2])
        cons4, y[3] <= 10(x[3])
        cons5, y[4] <= 10(x[4])
        cons6, y[5] <= 10(x[5])
        end)
        @objective(m, :Max, 5y[1] + 4y[2] + 3y[3] + 2y[4] + 1y[5])

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 140 atol=1e-4
        @test getvalue(x) ≈ [1;1;1;1;0]
        @test getvalue(y) ≈ [10;10;10;10;0] atol = exp10(-5)

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste MIP Pequeno (~50 binarias) ~ The Assignment problem:
# adicionado por Eduardo Brito
function test_MIP_Pequeno_Brito(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Pequeno, Assignment problem" begin
        n  = 8
        if true
            c = [0.445321 0.499462 0.409753 0.471825 0.1172 0.820595 0.629809 0.333445; 0.197025 0.160481 0.00865311 0.355901 0.137367 0.199186 0.718575 0.716486;
            0.654497 0.904598 0.321483 0.171736 0.050554 0.254487 0.540093 0.724331; 0.254369 0.593379 0.205166 0.288702 0.499699 0.308233 0.869406 0.353904;
            0.854515 0.00978121 0.520072 0.985762 0.72076 0.317384 0.268573 0.315585; 0.0212753 0.754076 0.753672 0.158407 0.212617 0.403343 0.71157 0.17261;
            0.651835 0.24596 0.700141 0.989018 0.723494 0.236829 0.891181 0.568245; 0.257637 0.883802 0.0252095 0.0273074 0.450492 0.560833 0.820861 0.893546]
        end
        @variable(m, x[i=1:n,j=1:n] >=0, Bin)
        @objective(m, :Min, sum(c[i,j]*x[i,j] for i = 1:n, j = 1:n))
        @constraints(m,begin
        linhas[i=1:n], sum(x[i,j] for j = 1:n) == 1
        colunas[j=1:n], sum(x[i,j] for i = 1:n) == 1
        end)

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 1.264 atol = exp10(-5)
        @test abs.(getvalue(x)) ≈ abs.([0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0;
                                        0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0;
                                        0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0;
                                        0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0;
                                        0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0;
                                        1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                                        0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0;
                                        0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0])

        @test m.ext[:status] == :Optimal

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste n-K robusto - trabalho da P1 - caso com 10 geradores e K = 2
# adicionado por Raphael Saavedra
function testRobustCCUC(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    #------------------------------------------------------------------------------
    # Parameters
    T = collect(1:24) # periods
    N = collect(1:10) # generators
    p0 = zeros(10) # initial power output
    v0 = zeros(10) # initial on/off state
    D = 0.5*[700; 750; 850; 950; 1000; 1100; 1150; 1200; 1300; 1400; 1450; 1500; 1400;
                1300; 1200; 1050; 1000; 1100; 1200; 1400; 1300; 1100; 900; 800] # demand
    K = 2 # n-K security criterion
    Cf = [1000; 970; 700; 680; 450; 370; 480; 660; 665; 670] # fixed cost
    Cl = [16.19; 17.26; 16.6; 16.5; 19.7; 22.26; 27.74; 25.92; 27.27; 27.79] # linear cost
    Cs = 0.08*Cl # spinning reserve cost
    Cns = 0.1*Cl # non-spinning reserve cost
    Pmax = [455; 455; 130; 130; 162; 80; 85; 55; 55; 55] # generator capacity
    Pmin = [150; 150; 20; 20; 25; 20; 25; 10; 10; 10] # minimum power output
    RSmax = Pmax # maximum spinning reserve
    RNSmax = Pmax # maximum non-spinning reserve
    RD = Pmax # ramp-down limit
    RU = Pmax # ramp-up limit
    SD = RD # shutdown ramp limit
    SU = RU # startup ramp limit
    #------------------------------------------------------------------------------
    m = Model(solver = solver)
    testresult = @testset "Teste Robust n-K Unit Commitment" begin
        # Model formulation
        #----------------------------------------------------------------------------
        # Variables
        @variable(m, p[1:N[end], 1:T[end]] >= 0) # power output
        @variable(m, v[1:N[end], 1:T[end]], Bin) # 1 if generator is on, 0 otherwise
        @variable(m, 0 <= rs[i = 1:N[end], 1:T[end]] <= RSmax[i]) # spinning reserve
        @variable(m, rns[1:N[end], 1:T[end]] >= 0) # non-spinning reserve
        @variable(m, vns[1:N[end], 1:T[end]], Bin) # 1 if generator provides non-spinning reserve, 0 otherwise
        @variable(m, y[1:T[end]] >= 0) # dual variable of the n-K constraint
        @variable(m, z[1:N[end], 1:T[end]] >= 0) # dual variable of the upper bound
        #----------------------------------------------------------------------------
        # Constraints
        @constraint(m, [t in T], sum(p[i,t] for i in N) == D[t])
        @constraint(m, [t in T, i in N], p[i,t] >= Pmin[i]*v[i,t])
        @constraint(m, [t in T, i in N], p[i,t] <= Pmax[i]*v[i,t])
        @constraint(m, [t in T, i in N], p[i,t] + rs[i,t] <= Pmax[i]*v[i,t])
        @constraint(m, [t in T, i in N], rns[i,t] >= Pmin[i]*vns[i,t])
        @constraint(m, [t in T, i in N], rns[i,t] <= RNSmax[i]*vns[i,t])
        @constraint(m, [t in T, i in N], v[i,t] + vns[i,t] <= 1)
        @constraint(m, [t in 2:T[end], i in N], p[i,t-1] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v[i,t-1]-v[i,t]) + Pmax[i]*(1-v[i,t-1]))
        @constraint(m, [t in 2:T[end], i in N], p[i,t] <= p[i,t-1] + RU[i]*v[i,t-1] + SU[i]*(v[i,t]-v[i,t-1]) + Pmax[i]*(1-v[i,t]))
        @constraint(m, [t in 2:T[end], i in N], p[i,t] + rs[i,t] <= p[i,t-1] + RU[i]*v[i,t-1] + SU[i]*(v[i,t]-v[i,t-1]) + Pmax[i]*(1-v[i,t]))
        @constraint(m, [t in 2:T[end], i in N], rns[i,t] <= p[i,t-1] + RU[i]*v[i,t-1] + RNSmax[i]*(1-v[i,t-1]))
        @constraint(m, [t in 2:T[end], i in N], rns[i,t] <= SU[i]*(vns[i,t]-v[i,t-1]) + RNSmax[i]*(1-(vns[i,t]-v[i,t-1])))
        @constraint(m, [t in T], (N[end]-K)*y[t] - sum(z[i,t] for i in N) >= D[t])
        @constraint(m, [t in T, i in N], y[t] - z[i,t] <= p[i,t] + rs[i,t] + rns[i,t])

        # Constraints regarding initial state
        @constraint(m, [t in 1, i in N], p0[i] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v0[i]-v[i,t]) + Pmax[i]*(1-v0[i]))
        @constraint(m, [t in 1, i in N], p[i,t] <= p0[i] + RU[i]*v0[i] + SU[i]*(v[i,t]-v0[i]) + Pmax[i]*(1-v[i,t]))
        @constraint(m, [t in 1, i in N], p[i,t] + rs[i,t] <= p0[i] + RU[i]*v0[i] + SU[i]*(v[i,t]-v0[i]) + Pmax[i]*(1-v[i,t]))
        @constraint(m, [t in 1, i in N], rns[i,t] <= p0[i] + RU[i]*v0[i] + RNSmax[i]*(1-v0[i]))
        @constraint(m, [t in 1, i in N], rns[i,t] <= SU[i]*(vns[i,t]-v0[i]) + RNSmax[i]*(1-(vns[i,t]-v0[i])))
        #----------------------------------------------------------------------------
        # Objective function
        @objective(m, Min, sum(Cf[i]*v[i,t] for i in N, t in T) + sum(Cl[i]*p[i,t] for i in N, t in T) +
                                sum(Cs[i]*rs[i,t] for i in N, t in T) + sum(Cns[i]*rns[i,t] for i in N, t in T))
        #------------------------------------------------------------------------------
        solveMIP(m)

        @test getobjectivevalue(m) ≈ 289892.9539 rtol=1e-3
        @test getvalue(p[:,24]) ≈ [400; zeros(9)]
        @test getvalue(v[1,:]) ≈ ones(24)
        @test getvalue(v[:,1]) ≈ [1; zeros(9)]
        @test getvalue(v[:,4]) ≈ [1; zeros(4); 1; zeros(4)]
        @test sum(getvalue(v)) ≈ 43

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste Caminho mais curto
# adicionado por Carlos
function testCaminho(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste caminho mais curto" begin

        @variable(m, x[i in 1:6, j in 1:6], Bin)

        A = [0 1 1 0 0 0
             0 0 0 1 0 0
             0 1 0 1 1 0
             0 0 0 0 1 1
             0 0 0 0 0 1
             0 0 0 0 0 0]

        c = [0 2 2 0 0 0
             0 0 0 3 0 0
             0 1 0 1 3 0
             0 0 0 0 1 1
             0 0 0 0 0 2
             0 0 0 0 0 0]

        b = [1;0;0;0;0;-1]

        @constraint(m,[v=1:6], sum(A[v,j]*x[v,j] for j=1:6) - sum(A[i,v]*x[i,v] for i=1:6) == b[v])

        @objective(m, Min, sum(A[i,j]*c[i,j]*x[i,j] for i=1:6, j=1:6))

        solveMIP(m)
        @test getobjectivevalue(m) == 4
        @test getvalue(x[1,3]) == 1
        @test getvalue(x[3,4]) == 1
        @test getvalue(x[4,6]) == 1
        @test sum(getvalue(x)) == 3
    end
    setoutputs!(m,solution,testresult)
    return solution
end


#teste problema 6 da lista modificado 1 (Expansao da Producao Unbounded)
#adicionado por Andrew Rosemberg
function test3_2(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    Cinv = 13.16
    M = 200
    m = Model(solver = solver)
    testresult = @testset "Teste da Expansao da Producao Unbounded" begin
        @variable(m, x[i=1:2]>=0)
        @variable(m, u, Bin)
        @objective(m, Max, 4*x[1] + 3*x[2] - u*Cinv)

        @constraint(m, 2*x[1] + 1*x[2] >= 4 +u*M)
        @constraint(m, 1*x[1] + 2*x[2] >= 4 +u*M)

        @constraint(m, 1*x[1] + 0.1*x[2] >= 4 +(1-u)*M)
        @constraint(m, 0.4*x[1] + 1*x[2] >= 4 +(1-u)*M)

        solveMIP(m)

        @test m.ext[:status] == :Unbounded
    end
    setoutputs!(m,solution,testresult)
    return solution
end


#teste problema 6 da lista modificado 2 (Expansao da Producao Infeasible)
#adicionado por Andrew Rosemberg
function test3_3(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    Cinv = 13.16
    M = 200
    m = Model(solver = solver)
    testresult = @testset "Teste da Expansao da Infeasible" begin
        @variable(m, x[i=1:2]>=0)
        @variable(m, u, Bin)
        @objective(m, Max, 4*x[1] + 3*x[2] - u*Cinv)

        @constraint(m, 2*x[1] + 1*x[2] <= 4 -u*M)
        @constraint(m, 1*x[1] + 2*x[2] <= 4 -u*M)

        @constraint(m, 1*x[1] + 0.1*x[2] <= 4 -(1-u)*M)
        @constraint(m, 0.4*x[1] + 1*x[2] <= 4 -(1-u)*M)

        solveMIP(m)

        @test m.ext[:status] in [:Infeasible,:InfeasibleOrUnbounded]

    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste Feature selection pequeno (Viavel)
#adicionado por Andrew Rosemberg
function test_feature_selection_pequeno_viavel(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    srand(2)
    numpossiblevar = 50
    numvar = 40
    numconstraints = 50
    datamatrix = rand(numconstraints,numpossiblevar)
    weights = [rand(numvar);zeros(numpossiblevar-numvar)]
    index_vector = datamatrix*weights
    maxweight  = maximum(weights)+1
    m = Model(solver = solver)
    testresult = @testset "Teste Feature selection pequeno Viavel" begin
        @objective(m, Max, 0)

        @variable(m, w[i=1:numpossiblevar]>=0)
        @variable(m, u[i=1:numpossiblevar], Bin)

        @constraint(m, sum(u) == numvar)

        @constraints(m,begin
          dummy[i=1:numpossiblevar], w[i] <= u[i]*maxweight
        end)

        @constraints(m,begin
          linhas[i=1:numconstraints], (datamatrix*w)[i] == index_vector[i]
        end)

        solveMIP(m)

        utest = [ones(numvar);zeros(numpossiblevar-numvar)]
        @test utest == getvalue(u)
        @test getvalue(w) ≈ weights atol=1E-07

    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste Feature selection medio (Viavel)
#adicionado por Andrew Rosemberg
function test_feature_selection_medio(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    srand(1000)
    numpossiblevar = 500
    numvar = 15
    numconstraints = 100
    datamatrix = rand(numconstraints,numpossiblevar)
    weights = [rand(numvar);zeros(numpossiblevar-numvar)]
    index_vector = datamatrix*weights
    maxweight  = maximum(weights)+1
    m = Model(solver = solver)
    testresult = @testset "Teste Feature selection pequeno Viavel" begin
        @objective(m, Max, 0)

        @variable(m, w[i=1:numpossiblevar]>=0)
        @variable(m, u[i=1:numpossiblevar], Bin)

        @constraint(m, sum(u) == numvar)

        @constraints(m,begin
          dummy[i=1:numpossiblevar], w[i] <= u[i]*maxweight
        end)

        @constraints(m,begin
          linhas[i=1:numconstraints], (datamatrix*w)[i] == index_vector[i]
        end)

        solveMIP(m)

        utest = [ones(numvar);zeros(numpossiblevar-numvar)]
        @test utest == getvalue(u)
        @test getvalue(w) ≈ weights atol=1E-07

    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste Feature selection grande (Viavel)
#adicionado por Andrew Rosemberg
function test_feature_selection_grande(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    srand(1000)
    numpossiblevar = 5000
    numvar = 15
    numconstraints = 100
    datamatrix = rand(numconstraints,numpossiblevar)
    weights = [rand(numvar);zeros(numpossiblevar-numvar)]
    index_vector = datamatrix*weights
    maxweight  = maximum(weights)+1
    m = Model(solver = solver)
    testresult = @testset "Teste Feature selection pequeno Viavel" begin
        @objective(m, Max, 0)

        @variable(m, w[i=1:numpossiblevar]>=0)
        @variable(m, u[i=1:numpossiblevar], Bin)

        @constraint(m, sum(u) == numvar)

        @constraints(m,begin
          dummy[i=1:numpossiblevar], w[i] <= u[i]*maxweight
        end)

        @constraints(m,begin
          linhas[i=1:numconstraints], (datamatrix*w)[i] == index_vector[i]
        end)

        solveMIP(m)

        utest = [ones(numvar);zeros(numpossiblevar-numvar)]
        @test utest == getvalue(u)
        @test getvalue(w) ≈ weights atol=1E-07

    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste Feature selection pequeno (Inviavel)
#adicionado por Andrew Rosemberg
function test_feature_selection_pequeno_inviavel(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    srand(2)
    numpossiblevar = 50
    numvar = 40
    numconstraints = 50
    datamatrix = rand(numconstraints,numpossiblevar)
    weights = [rand(numvar);zeros(numpossiblevar-numvar)]
    index_vector = -datamatrix*weights
    maxweight  = maximum(weights)+1
    m = Model(solver = solver)
    testresult = @testset "Teste Feature selection pequeno Viavel" begin
        @objective(m, Max, 0)

        @variable(m, w[i=1:numpossiblevar]>=0)
        @variable(m, u[i=1:numpossiblevar], Bin)

        @constraint(m, sum(u) == numvar)

        @constraints(m,begin
          dummy[i=1:numpossiblevar], w[i] <= u[i]*maxweight
        end)

        @constraints(m,begin
          linhas[i=1:numconstraints], (datamatrix*w)[i] == index_vector[i]
        end)

        solveMIP(m)

        @test m.ext[:status] == :Infeasible

    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste problema ucla : https://www.math.ucla.edu/~tom/LP.pdf pg 9 (PL unbounded)
#adicionado por Andrew Rosemberg
function teste_PL_andrew_unbounded(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste PL unbounded da ucla" begin
        @variable(m, x[i=1:3]>=0)
        @variable(m, x4)
        @objective(m, Max, x[1] + 3*x[2] + 4*x[3] + 2*x4)

        @constraint(m, 4*x[1] + 2*x[2] +x[3] + 3*x4 <= 10)
        @constraint(m, x[1] - x[2] + 2*x[3] == 2)
        @constraint(m, x[1] + x[2] + x[3] + x4 >= 1)

        solveMIP(m)

        @test m.ext[:status] in [:Unbounded, :InfeasibleOrUnbounded]

    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste problema ucla modificado (PL viavel)
#adicionado por Andrew Rosemberg
function teste_PL_andrew_viavel(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste PL viavel da ucla" begin
        @variable(m, x[i=1:3]>=0)
        @variable(m, x4>=0)
        @objective(m, Max, x[1] + 3*x[2] + 4*x[3] + 2*x4)

        @constraint(m, 4*x[1] + 2*x[2] +x[3] + 3*x4 <= 10)
        @constraint(m, x[1] - x[2] + 2*x[3] == 2)
        @constraint(m, x[1] + x[2] + x[3] + x4 >= 1)

        solveMIP(m)

        @test getobjectivevalue(m) == 22
        @test getvalue(x) ≈ [0.0;3.6;2.8] atol=1E-07
        @test getvalue(x4) == 0

    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste problema ucla modificado 2 (PL inviavel)
#adicionado por Andrew Rosemberg
function teste_PL_andrew_inviavel(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste PL inviavel da ucla" begin
        @variable(m, x[i=1:3]>=0)
        @variable(m, x4>=0)
        @objective(m, Max, x[1] + 3*x[2] + 4*x[3] + 2*x4)

        @constraint(m, 4*x[1] + 2*x[2] +x[3] + 3*x4 <= -10)
        @constraint(m, x[1] - x[2] + 2*x[3] == 2)
        @constraint(m, x[1] + x[2] + x[3] + x4 >= 1)

        solveMIP(m)

        @test m.ext[:status] == :Infeasible

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste MIP unbounded
# adicionado por Raphael Saavedra
function testUnboundedKnapsack(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste Mochila unbounded" begin
        @variable(m, x[i=1:3], Bin)
        @variable(m, y >= 0)
        @constraint(m, 6*x[1] + 5*x[2] + 5*x[3] <= 5)
        @objective(m, Max, 6*x[1] + 4*x[2] + 3*x[3] + y)

        solveMIP(m)

        @test m.ext[:status] == :Unbounded

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste Infeasible Unit Commitment
# adicionado por Raphael Saavedra
function testInfeasibleUC(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    #------------------------------------------------------------------------------
    # Parameters
    T = collect(1:10) # periods
    N = collect(1:5) # generators
    p0 = zeros(N[end]) # initial power output
    v0 = zeros(N[end]) # initial on/off state
    D = [400; 390; 380; 370; 360; 350; 340; 330; 320; 500] # demand
    Cf = [100; 100; 100; 100; 100] # fixed cost
    Cl = [10; 20; 30; 40; 50] # linear cost
    Pmax = [100; 100; 100; 100; 100] # generator capacity
    Pmin = [10; 10; 10; 10; 10] # minimum power output
    RD = [10; 20; 30; 40; 50] # ramp-down limit
    RU = [10; 25; 30; 40; 50] # ramp-up limit
    SD = RD # shutdown ramp limit
    SU = RU # startup ramp limit
    #------------------------------------------------------------------------------
    m = Model(solver = solver)
    testresult = @testset "Teste Infeasible Unit Commitment" begin

        #----------------------------------------------------------------------------
        # Variables
        @variable(m, p[1:N[end], 1:T[end]] >= 0) # power output
        @variable(m, v[1:N[end], 1:T[end]], Bin) # 1 if generator is on, 0 otherwise
        #----------------------------------------------------------------------------
        # Constraints
        @constraint(m, [t in T], sum(p[i,t] for i in N) == D[t])
        @constraint(m, [t in T, i in N], p[i,t] >= Pmin[i]*v[i,t])
        @constraint(m, [t in T, i in N], p[i,t] <= Pmax[i]*v[i,t])
        @constraint(m, [t in 2:T[end], i in N], p[i,t-1] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v[i,t-1]-v[i,t]) + Pmax[i]*(1-v[i,t-1]))
        @constraint(m, [t in 2:T[end], i in N], p[i,t] <= p[i,t-1] + RU[i]*v[i,t-1] + SU[i]*(v[i,t]-v[i,t-1]) + Pmax[i]*(1-v[i,t]))
        # Constraints regarding initial state
        @constraint(m, [t in 1, i in N], p0[i] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v0[i]-v[i,t]) + Pmax[i]*(1-v0[i]))
        @constraint(m, [t in 1, i in N], p[i,t] <= p0[i] + RU[i]*v0[i] + SU[i]*(v[i,t]-v0[i]) + Pmax[i]*(1-v[i,t]))
        #----------------------------------------------------------------------------
        # Objective function
        @objective(m, Min, sum(Cf[i]*v[i,t] for i in N, t in T) + sum(Cl[i]*p[i,t] for i in N, t in T))
        #------------------------------------------------------------------------------
        solveMIP(m)

        @test m.ext[:status] == :Infeasible

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste PL simples
# adicionado por Raphael Saavedra
function test_PL_Simples_Raphael(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste PL simples" begin
        @variable(m, x[i=1:3] >= 0)
        @constraint(m, x[1] + x[2] <= 2)
        @constraint(m, x[1] + x[3] <= 2)
        @constraint(m, x[2] + x[3] <= 2)
        @objective(m, Max, x[1] + x[2] - 2*x[3])

        solveMIP(m)

        @test m.ext[:status] == :Optimal
        @test getobjectivevalue(m) == 2
    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste PL infeasible
# adicionado por Raphael Saavedra
function test_PL_Infeasible_Raphael(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste PL infeasible" begin
        @variable(m, x[i=1:2] >= 0)
        @constraint(m, x[1] + x[2] <= -1)
        @objective(m, Max, x[1] + x[2])

        solveMIP(m)

        @test m.ext[:status] == :Infeasible

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste Minimal Unit Commitment
# adicionado por Raphael Saavedra
function test_Minimal_UC(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    #------------------------------------------------------------------------------
    # Parameters
    T = collect(1:3) # periods
    N = collect(1:3) # generators
    p0 = zeros(N[end]) # initial power output
    v0 = zeros(N[end]) # initial on/off state
    D = [100; 200; 300] # demand
    Cf = [100; 100; 100] # fixed cost
    Cl = [10; 30; 50] # linear cost
    Pmax = [100; 150; 200] # generator capacity
    Pmin = [10; 10; 10] # minimum power output
    RD = [30; 50; 70] # ramp-down limit
    RU = [30; 50; 70] # ramp-up limit
    SD = RD # shutdown ramp limit
    SU = RU # startup ramp limit
    #------------------------------------------------------------------------------
    m = Model(solver = solver)
    testresult = @testset "Teste Minimal Unit Commitment" begin

        #----------------------------------------------------------------------------
        # Variables
        @variable(m, p[1:N[end], 1:T[end]] >= 0) # power output
        @variable(m, v[1:N[end], 1:T[end]], Bin) # 1 if generator is on, 0 otherwise
        #----------------------------------------------------------------------------
        # Constraints
        @constraint(m, [t in T], sum(p[i,t] for i in N) == D[t])
        @constraint(m, [t in T, i in N], p[i,t] >= Pmin[i]*v[i,t])
        @constraint(m, [t in T, i in N], p[i,t] <= Pmax[i]*v[i,t])
        @constraint(m, [t in 2:T[end], i in N], p[i,t-1] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v[i,t-1]-v[i,t]) + Pmax[i]*(1-v[i,t-1]))
        @constraint(m, [t in 2:T[end], i in N], p[i,t] <= p[i,t-1] + RU[i]*v[i,t-1] + SU[i]*(v[i,t]-v[i,t-1]) + Pmax[i]*(1-v[i,t]))
        # Constraints regarding initial state
        @constraint(m, [t in 1, i in N], p0[i] <= p[i,t] + RD[i]*v[i,t] + SD[i]*(v0[i]-v[i,t]) + Pmax[i]*(1-v0[i]))
        @constraint(m, [t in 1, i in N], p[i,t] <= p0[i] + RU[i]*v0[i] + SU[i]*(v[i,t]-v0[i]) + Pmax[i]*(1-v[i,t]))
        #----------------------------------------------------------------------------
        # Objective function
        @objective(m, Min, sum(Cf[i]*v[i,t] for i in N, t in T) + sum(Cl[i]*p[i,t] for i in N, t in T))
        #------------------------------------------------------------------------------
        solveMIP(m)

        @test m.ext[:status] == :Optimal
        @test getobjectivevalue(m) ≈ 17700 atol = 1e-5
        @test getvalue(p) ≈ [30 60 90 ; 50 100 150 ; 20 40 60] atol = 1e-5

    end
    setoutputs!(m,solution,testresult)
    return solution
end

# teste Sudoku 4x4
# adicionado por Raphael Saavedra
function testSudoku4x4(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste Sudoku 4x4" begin
        n = 4
        @variable(m, x[i in 1:n, j in 1:n, k in 1:n], Bin)

        fixas = [(1,1,1), (2,2,3), (1,3,4), (3,3,2), (4,4,4), (4,2,1)]
        for idx in fixas
            @constraint(m, x[idx...] == 1)
        end
        @constraint(m, [j in 1:n, k in 1:n], sum(x[:,j,k]) == 1)
        @constraint(m, [i in 1:n, k in 1:n], sum(x[i,:,k]) == 1)
        @constraint(m, [i in 1:n, j in 1:n], sum(x[i,j,:]) == 1)
        @constraint(m, [p in [0,2], q in [0,2], k in 1:n], sum(sum(x[i+p,j+q,k] for i in 1:2) for j in 1:2) == 1)
        @objective(m, Min, 0)

        solveMIP(m)

        M = Matrix(4,4)
        for i = 1 : 4
          for j = 1 : 4
            M[i,j] = find(getvalue(x[i,j,:]).>0)[1]
          end
        end

        @test M == [1 2 4 3; 4 3 1 2; 3 4 2 1; 2 1 3 4]
    end
    setoutputs!(m,solution,testresult)
    return solution
end


#adicionado por Rodrigo Villas
function test_rv_1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m= Model(solver = solver)
    testresult = @testset "Custo fixo" begin


        #Custo Unitário
        c = [2; 3; 2; 5; 6]

        #Custo Fixo
        f = [1; 3; 1; 5; 10]


        @variable(m, x[i=1:5]>=0)

        @variable(m,y[i=1:5], Bin)

        @constraint(m, sum(x[j] for j=1:5) >=10)

        @constraint(m,x[1]<=5*y[1])

        @constraint(m,x[2]<=4*y[2])

        @constraint(m,x[3]<=3*y[3])

        @constraint(m,x[4]<=2*y[4])

        @constraint(m,x[5]<=1*y[5])

        @objective(m, Min, sum(f[j]*y[j]+c[j]*x[j] for j=1:5))

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 27
        @test getvalue(x) ≈ [5, 2, 3, 0, 0]

    end
    setoutputs!(m,solution,testresult)
    return solution
end

#adicionado por Rodrigo Villas
function test_rv_2(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Cobertura de pontos" begin


        pontosextra= 0

        #um ponto não é coberto por nenhum subconjunto

        S1=[1 0 0 1 ones(1,pontosextra)]
        S2=[1 1 0 0 ones(1,pontosextra)]
        S3=[0 1 0 1 ones(1,pontosextra)]

        c=[4 3 2]

        A=[S1' S2' S3']

        @variable(m, x[i=1:3], Bin)
        @constraints(m, begin
          constrain[i=1:4+pontosextra], sum(A[i,j]*x[j] for j=1:3)>= 1
          end)
        @objective(m, Min, sum(c[j]*x[j] for j=1:3))
        solveMIP(m)

        @test m.ext[:status] == :Infeasible

    end
    setoutputs!(m,solution,testresult)
    return solution

end

#adicionado por Rodrigo Villas
function test_rv_3(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Cobertura de pontos" begin

        pontosextra= 50

        S1=[1 0 1 1 ones(1,pontosextra)]
        S2=[1 1 0 0 ones(1,pontosextra)]
        S3=[0 1 0 1 ones(1,pontosextra)]

        c=[4 3 2]

        A=[S1' S2' S3']

        @variable(m, x[i=1:3], Bin)
        @constraints(m, begin
          constrain[i=1:4+pontosextra], sum(A[i,j]*x[j] for j=1:3)>= 1
          end)
        @objective(m, Min, sum(c[j]*x[j] for j=1:3))

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 6
        @test getvalue(x) ≈ [1, 0, 1]

    end
    setoutputs!(m,solution,testresult)
    return solution
end

#adicionado por Rodrigo Villas
#Produção com custo fixo Inviavel
function test_rv_4(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m= Model(solver = solver)
    testresult = @testset "Custo fixo" begin

        vari=50


        @variable(m, x[i=1:vari]>=0)
        @variable(m,y[i=1:vari], Bin)
        @constraint(m, sum(x[j] for j=1:vari) >= 2*vari)
        @constraint(m, [i=1:vari], x[i] <= 1*y[i])

        @objective(m, Min, sum(j*y[j]+(vari-j)*x[j] for j=1:vari))

        solveMIP(m)
        @test m.ext[:status] == :Infeasible

    end
    setoutputs!(m,solution,testresult)
    return solution
end

#adicionado por Rodrigo Villas
function test_rv_5(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m= Model(solver = solver)
    testresult = @testset "Custo fixo" begin


        vari=500


        @variable(m, x[i=1:vari]>=0)

        @variable(m,y[i=1:vari], Bin)

        @constraint(m, sum(x[j] for j=1:vari) >= vari)

        @constraints(m, begin
          constrain[i=1:vari], x[i] <= 10*y[i]
          end)

        @objective(m, Min, sum(j*y[j]+(vari-j)*x[j] for j=1:vari))

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 36025

        #Os últimos Y's tem que estar "ligados"
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#Expansão Unbounded
 #adicionado por Rodrigo Villas

#teste P1 TSP de 7 cidades
#adicionado por Guilherme Bodin
function test_P1_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste da P1 Guilherme (TSP 7 cidades)" begin
        number_of_nodes = 7

        C = [0.0    135.484   142.801   131.0     117.154   153.473   201.022
             135.484    0.0     105.546   174.003   142.425    53.0094  105.991
             142.801  105.546     0.0      87.6641   59.0      63.8905   73.5527
             131.0    174.003    87.6641    0.0      31.9061  146.932   159.201
             117.154  142.425    59.0      31.9061    0.0     115.521   132.306
             153.473   53.0094   63.8905  146.932   115.521     0.0      55.4617
             201.022  105.991    73.5527  159.201   132.306    55.4617    0.0   ]

        ans = [0.0   0.0   0.0   1.0   0.0   0.0   0.0
               1.0   0.0   0.0   0.0   0.0   0.0   0.0
               0.0   0.0   0.0   0.0   0.0   0.0   1.0
               0.0   0.0   0.0   0.0   1.0   0.0   0.0
               0.0   0.0   1.0   0.0   0.0   0.0   0.0
               0.0   1.0   0.0   0.0   0.0   0.0   0.0
               0.0   0.0   0.0   0.0   0.0   1.0   0.0]


        @variable(m, X[i=1:number_of_nodes,j=1:number_of_nodes], Bin)
        @variable(m, u[i=:1:number_of_nodes])#int
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
        @test getobjectivevalue(m) ≈ 539.4139 atol=1e-5
        #@test getvalue(X) ≈ ans atol=1e-5
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste P1 TSP de 7 cidades formulação binária
#adicionado por Guilherme Bodin
function test_TSPbin7_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste da P1 Guilherme binário (TSP 7 cidades)" begin
        number_of_nodes = 7

        AdjacencyMatrix = [0.0    135.484   142.801   131.0     117.154   153.473   201.022
                           135.484    0.0     105.546   174.003   142.425    53.0094  105.991
                           142.801  105.546     0.0      87.6641   59.0      63.8905   73.5527
                           131.0    174.003    87.6641    0.0      31.9061  146.932   159.201
                           117.154  142.425    59.0      31.9061    0.0     115.521   132.306
                           153.473   53.0094   63.8905  146.932   115.521     0.0      55.4617
                           201.022  105.991    73.5527  159.201   132.306    55.4617    0.0   ]

             function num_edges(Matrix) # Retorna o número de arestas em uma matriz de adjacência
               n_edges = (factorial(size(Matrix,1)))/(factorial(2)*factorial(size(Matrix,1)-2))
               return round(Int64,n_edges)
             end

             function create_edges(Matrix) # Cria as arestas a partir da Matriz de Adjacencia
               number_of_edges = num_edges(Matrix)
               e = zeros(number_of_edges,2)
               k=1;
               for i=1:size(Matrix,1)
                 for j=i+1:size(Matrix,2)
                   e[k,1]=i
                   e[k,2]=j
                   k = k+1
                 end
               end
               return round(Int64,e)
             end

             function cost_of_edges(Matrix) #Cria o vetor de custo associado a cada aresta
               number_of_edges = num_edges(Matrix)
               cost = zeros(number_of_edges,1)
               k=1;
               for i=1:size(Matrix,1)
                 for j=i+1:size(Matrix,2)
                   cost[k]=Matrix[i,j]
                   k = k+1
                 end
               end
               return cost
             end

             function delta(AdjacencyMatrix,EdgesMatrix) #Devolve nas colunas as arestas ligadas à cidade i correspondente ao indice da linha
               lin=1
               delta = zeros(size(AdjacencyMatrix,1),size(AdjacencyMatrix,1)-1)
               for j=1:size(AdjacencyMatrix,1)
                 col=1
                 for i=1:size(EdgesMatrix,1)
                   if (EdgesMatrix[i,1] == j || EdgesMatrix[i,2] == j)
                     delta[lin,col]=i
                     col=col+1
                   end
                 end
                 lin=lin+1
               end
               return round(Int64,delta)
             end

             function subsets(Matrix) # Retorna todos os subconjuntos de vértices de uma Matriz de Adjacência
               number_of_nodes = size(Matrix,1)
               Subset_S = digits(1,2,number_of_nodes)
               for i = 2:2^(number_of_nodes)-2
                 aux = digits(i,2,number_of_nodes)
                 Subset_S = hcat(Subset_S,aux)
               end
               return Subset_S
             end

             function edges_no_Subconjunto(pos,Subsets,number_of_nodes) #Seleciona todas as arestas possíveis em uma dada partição do conjunto de vértices
               A = zeros(1,2)
               for i=1:number_of_nodes
                 if Subsets[i,pos]==1
                   for j=1:number_of_nodes
                     if Subsets[j,pos]==0
                       A = vcat(A,[i j])
                     end
                   end
                 end
               end
               return A
             end

             function ida_e_volta(A) # Recebe uma matriz A n x 2 e devolve uma matriz result 2*n x 2 com os elementos de A nas linhas 1 até n e o espelhamento dos elementos de A de n+1 até 2*n
               A_aux = hcat(A,A[:,1])
               A_aux = A_aux[:,2:3]
               result = vcat(A,A_aux)
               return result
             end


             function selecao_edges(e,edge_do_ciclo)
               E = []
               for i=1:size(e,1)
                 for j=1:size(edge_do_ciclo,1)
                   if edge_do_ciclo[j,1] == e[i,1]
                     if edge_do_ciclo[j,2] == e[i,2]
                       E = vcat(E,i)
                     end
                   end
                 end
               end
               return E
             end

             e=create_edges(AdjacencyMatrix)
             del = delta(AdjacencyMatrix,e)
             S = subsets(AdjacencyMatrix)
             number_of_edges = num_edges(AdjacencyMatrix)
             number_of_nodes = size(AdjacencyMatrix,1)
             C = cost_of_edges(AdjacencyMatrix)
             @variable(m, Y[i=1:number_of_edges],Bin)
             for i=1:number_of_nodes
                 @constraint(m, sum(Y[del[i,j]] for j=1:number_of_nodes-1) == 2) # Garante que todas as cidades conectadas a 2 estradas (uma de entrada e outra de saída) no caminho ótimo
             end
             for i=1:size(S,2)
                 edge_do_ciclo = ida_e_volta(edges_no_Subconjunto(i,S,number_of_nodes))
                 E = selecao_edges(e,edge_do_ciclo)
                 @constraint(m, sum(Y[E[j]] for j=1:size(E,1)) >= 2)
             end
             @objective(m,Min,sum(C[i]*Y[i] for i=1:number_of_edges))

             solveMIP(m)
             @test getobjectivevalue(m) ≈ 539.4139 atol = 1e-5
    end
    setoutputs!(m,solution,testresult)
    return solution
end


#adicionado por Guilherme Bodin
#TSP com a formulação MTZ (Muller Tucker Zemlin)
#http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.34.7256&rep=rep1&type=pdf
function test_TSPmip7_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste da P1 Guilherme mip (TSP 7 cidades)" begin
        number_of_nodes = 7

        C = [0.0    135.484   142.801   131.0     117.154   153.473   201.022
             135.484    0.0     105.546   174.003   142.425    53.0094  105.991
             142.801  105.546     0.0      87.6641   59.0      63.8905   73.5527
             131.0    174.003    87.6641    0.0      31.9061  146.932   159.201
             117.154  142.425    59.0      31.9061    0.0     115.521   132.306
             153.473   53.0094   63.8905  146.932   115.521     0.0      55.4617
             201.022  105.991    73.5527  159.201   132.306    55.4617    0.0   ]


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

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 539.4139 atol = 1e-5
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste PL Infeasible
#adicionado por Guilherme Bodin
function test_PL_Infeasible_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste PL Infeasible Guilherme" begin
        m = Model(solver = solver)
        @variable(m, x[i=1:2])
        @constraint(m, x[1] == 6)
        @constraint(m, x[2] == 6)
        @constraint(m, x[1] + x[2] <=11)
        @objective(m, Min, x[1]+x[2])

        solveMIP(m)
        @test m.ext[:status] == :Infeasible
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste MIP médio TSP de 15 cidades
#adicionado por Guilherme Bodin
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

#teste MIP grande TSP de 20 cidades
#adicionado por Guilherme Bodin
function test_TSPmip20_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Médio Guilherme (TSP 20 cidades)" begin
        number_of_nodes = 20
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

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 1.621683320972e+03 rtol = 1e-2
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste MIP grande TSP de 25 cidades
#adicionado por Guilherme Bodin
function test_TSPmip25_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Médio Guilherme (TSP 25 cidades)" begin
        number_of_nodes = 25
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

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 1862.2167124533548 rtol = 1e-2
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste MIP grande TSP de 30 cidades
#adicionado por Guilherme Bodin
function test_TSPmip30_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Médio Guilherme (TSP 30 cidades)" begin
        number_of_nodes = 30
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

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 1.810912865231e+03 rtol = 1e-2
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste MIP grande TSP de 40 cidades
#adicionado por Guilherme Bodin
function test_TSPmip40_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Médio Guilherme (TSP 40 cidades)" begin
        number_of_nodes = 40
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

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 1.847290056038e+03 rtol = 1e-2
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste MIP grande TSP de 50 cidades
#adicionado por Guilherme Bodin
function test_TSPmip50_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Médio Guilherme (TSP 50 cidades)" begin
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

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 1.868859333414e+03 rtol = 1e-2
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste MIP grande TSP de 100 cidades
#adicionado por Guilherme Bodin
function test_MIP_Grande_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Grande Guilherme (TSP 100 cidades)" begin
        number_of_nodes = 100
        srand(12)
        C = 1000*rand(100,100)

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

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 1720.190204078063 rtol = 1e-2
    end
    setoutputs!(m,solution,testresult)
    return solution
end


#teste Alocacao de portifolio P1 Andrew e Bianca (Viavel)
#adicionado por Andrew Rosemberg
function test_P1_Andrew_Bianca_viavel(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
  solution = MIPSolution()
  # params
  numD = 50 #num days
  numA = 7 #num assets
  price_rf = 5#0.0001
  W_0 = 100000
  j_robust = 3
  cost_trans = fill(8,numA+1)
  rf = (1+7.5/100).^(1/252)-1
  λ = 0.01*rf
  #init
  x_t1 = zeros(numA+1)
  x_t1[1] = W_0

  #Generate asset prices

  srand(82)
  prices = rand(1:100)+sin.(linspace(rand()*pi,numD))+randn(numD)*1

  for i= 2:numA
    srand(i)
    prices = [prices rand(1:100)+sin.(linspace(rand()*pi,numD))+randn(numD)*i]
  end
  returns= (prices[2:numD,1:numA].-prices[1:numD-1,1:numA])./prices[1:numD-1,1:numA]
  prices = [fill(price_rf, size(prices,1)) prices]

  t_1 = size(returns,1)
  ################# Predict return ######################
  r_bar_t = zeros(numA)

  r_bar_t[1] = rf
  for i = 2:numA
    r_bar_t[i] = mean(returns[:,i])
  end

  #######################################################
  myModel = Model(solver = solver)
  testresult = @testset "Alocacao de portifolio Viavel" begin

        # Decision variables
        @variable(myModel, X[1:numA]>=0)
        @variable(myModel, u_buy[1:numA]>=0 )  #Int
        @variable(myModel, u_sell[1:numA]>=0 )  #Int
        @variable(myModel, bin_buy[1:numA] , Bin)
        @variable(myModel, bin_sell[1:numA] , Bin)
        #expansao binaria
        numexp = 5
        @variable(myModel, expbin_buy[1:numA,1:numexp] , Bin)
        @variable(myModel, expbin_sell[1:numA,1:numexp] , Bin)
        @constraints(myModel, begin
          Expansaobuy[i=1:numA],  u_buy[i] == sum(expbin_buy[i,j+1]*(2^j) for j=0:numexp-1)
        end)
        @constraints(myModel, begin
          Expansaosell[i=1:numA],  u_sell[i] == sum(expbin_sell[i,j+1]*(2^j) for j=0:numexp-1)
        end)

        @constraints(myModel, begin
          link_buy[i=1:numA],  prices[i]*u_buy[i] <= bin_buy[i]*W_0
        end)
        @constraints(myModel, begin
          link_sell[i=1:numA],  prices[i]*u_sell[i] <= bin_sell[i]*W_0
        end)


        # Robust Constraints
        @constraints(myModel, begin
          robust[j=0:j_robust-1],  sum(returns[t_1-j,i-1]*X[i] for i=2:numA) >= λ*W_0
        end)

        # addapt alocation considering buy and sell costs
        @constraints(myModel, begin
          alocation1[i=1:numA],  X[i] == x_t1[i] + prices[i]*u_buy[i] - prices[i]*u_sell[i] -cost_trans[i]*(bin_buy[i]+bin_sell[i])
        end)

        # buy limit to bank reserve
        @constraint(myModel, sum(prices[i]*u_buy[i] for i=1:numA)-sum(prices[i]*u_sell[i] for i=1:numA) == 0)

        # objective function
        @objective(myModel, Max, sum(r_bar_t[i]*X[i] for i = 1:numA) - sum(cost_trans[i]*(bin_buy[i]+bin_sell[i]) for i = 1:numA)) # + sum((sum(returns[t_1-j,i-1]*X[i] for i=2:numA) - λ*W_0) for j=0:j_robust-1 )*penalizerobust )

        solveMIP(myModel)

        x =[98717.0;
           627.0;
           187.0;
             0.0;
           437.0;
             0.0;
             0.0]
        x2 = [99837.0; 67.0; 27.0; 0.0; 37.0; 0.0; 0.0]
        U_buy = [
            0.0;
            127.0;
            39.0;
            0.0;
            89.0;
            0.0;
            0.0]
        U_buy2 = [0.0; 15.0; 7.0; 0.0; 9.0; 0.0; 0.0]

        U_sell = [
            255.0;
            0.0;
            0.0;
            0.0;
            0.0;
            0.0;
            0.0]
        U_sell2 = [31.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]


        @test norm(x - getvalue(X)) < 1e-4 || norm(x2 - getvalue(X))  < 1e-4
        @test norm(U_buy - getvalue(u_buy)) < 1e-4  || norm(U_buy2 - getvalue(u_buy) )< 1e-4
        @test norm(U_sell - getvalue(u_sell))< 1e-4  || norm(U_sell2 - getvalue(u_sell) ) < 1e-4

    end
    setoutputs!(myModel,solution,testresult)
    return solution
end

#--------------------------

#adicionado por Rodrigo Villas
function test_rv_p1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)

    testresult =@testset "Bagulhão da P1 (não me pergunte pq)" begin

        QtdComp = 3
        Pcomp = [20 8 7]
        Qcomp = [5 10 10]
        Demanda = 25
        #---------#

        #Definição da sua produção#
        Pmax =100
        Pmin = 0

        custo = 0

        Qmax = 6
        Qmin = 0

        Qdiscre = 3

        Prodmáx = 10

        Bm=100

        #----------#
        deltap = (Pmax - Pmin)/(2^(Qdiscre-1))

        deltaof = (Qmax - Qmin)/(2^(Qdiscre-1))

        QtdVariaveis = 2*QtdComp + 2 + 1 +4*QtdComp

        deltasOferta = [deltaof 2*deltaof 4*deltaof]

        deltasPreço = [deltap 2*deltap 4*deltap]

        #Restrição referente a discretização do Bid#
        A = [zeros(1,QtdVariaveis-Qdiscre) deltasOferta]

        #Restrição referente a discretização do Preço#
        A = [A;zeros(1,QtdVariaveis-2*Qdiscre) deltasPreço zeros(1,Qdiscre)]

        #Total produzido = Demanda#
        A = [A;ones(1,QtdComp+1) zeros(1,QtdVariaveis-(QtdComp+1));-ones(1,QtdComp+1) zeros(1,QtdVariaveis-(QtdComp+1))]

        #Produzido <= Ofertado#
        A=[A;1 zeros(1,QtdVariaveis-1-Qdiscre) -deltasOferta]

        #Produção dos Competidores <= Oferta deles#
        A=[A;zeros(QtdComp,1) eye(QtdComp) zeros(QtdComp,QtdVariaveis-QtdComp-1)]

        #Primeira restrição dual: Spot + Dual(produção vc) - Preço de oferta <= 0#
        A=[A;zeros(1,QtdComp+1) 1 -1 zeros(1,QtdComp+2*Qdiscre) -deltasPreço zeros(1,Qdiscre)]

        #Dual da produção deles#
        A=[A;zeros(QtdComp,QtdComp+1) ones(QtdComp,1) zeros(QtdComp,1) -eye(QtdComp) zeros(QtdComp,4*Qdiscre)]

        #Primal = Dual#
        DP=[Pmin Pcomp -Demanda Qmin Qcomp deltasPreço deltasOferta zeros(1,2*Qdiscre)]
        DP=[DP;-DP]
        A=[A;DP]


        #Ordem: Gvc G1 G2 G3 Spot PIvc PI1 PI2 PI3 z1 z2 z3 w1 w2 w3 x1 x2 x3 y1 y2 y3 #

        #Restrições referentes a linearização binária#
        #Se Xk = 1, Zk = produção#
        A=[A;ones(Qdiscre,1) zeros(Qdiscre,1+2*QtdComp+1) -eye(Qdiscre) zeros(Qdiscre,Qdiscre) Bm*eye(Qdiscre) zeros(Qdiscre,Qdiscre)]

        #Se Xk = 0, Zk = 0#
        A=[A;zeros(Qdiscre,1) zeros(Qdiscre,1+2*QtdComp+1) eye(Qdiscre) zeros(Qdiscre,Qdiscre) -Bm*eye(Qdiscre) zeros(Qdiscre,Qdiscre)]

        #Se yk = 1, wk = dual da produção#
        A=[A;zeros(Qdiscre,1) zeros(Qdiscre,1+QtdComp) ones(Qdiscre,1) zeros(Qdiscre,QtdComp) zeros(Qdiscre,Qdiscre) -eye(Qdiscre) zeros(Qdiscre,Qdiscre) Bm*eye(Qdiscre)]

        #Se yk = 0, wk = 0#
        A=[A;zeros(Qdiscre,1) zeros(Qdiscre,1+2*QtdComp+1) zeros(Qdiscre,Qdiscre) -eye(Qdiscre) zeros(Qdiscre,Qdiscre) -Bm*eye(Qdiscre)]

        #Bin menor que 1#
        A=[A;zeros(2*QtdComp,QtdVariaveis-2*Qdiscre) eye(2*QtdComp)]

        #Ordem: Gvc G1 G2 G3 Spot PIvc PI1 PI2 PI3 z1 z2 z3 w1 w2 w3 x1 x2 x3 y1 y2 y3 #

        #Primal#
        b = [Qmax-Qmin;Pmax-Pmin;Demanda;-Demanda;Qmin;Qcomp']

        #Dual#
        b=[b;Pmin;Pcomp';0;0]

        #Binários#
        b=[b;Bm*ones(Qdiscre,1);zeros(Qdiscre,1);Bm*ones(Qdiscre,1);zeros(Qdiscre,1)]

        #Bin <1#
        b=[b;ones(2*QtdComp,1)]


        #Ordem: Gvc G1 G2 G3 Spot PIvc PI1 PI2 PI3 z1 z2 z3 w1 w2 w3 x1 x2 x3 y1 y2 y3 #
        c = [Pmin-custo; zeros(QtdComp+1,1); Qmin; zeros(QtdComp,1); deltasPreço' ;deltasOferta' ;zeros(2*Qdiscre,1)]

        @variable(m, y[i=QtdVariaveis-2*Qdiscre+1:QtdVariaveis]>=0, Bin)

        @variable(m, dual[i=QtdComp+3:2*QtdComp+3]>=0)

        p, n = size(A)
        Astd = zeros(p,p+n)
        cstd = zeros(p+n)
        Astd = [A eye(p)]
        cstd = [c ; zeros(p)]

        #---------------------------

        p,k = size(Astd)
        @variable(m, x[i=1:QtdComp+2]>=0)
        @variable(m, xc[i=2*QtdComp+4:QtdVariaveis-2*Qdiscre]>=0)
        @variable(m, z[i=QtdVariaveis+1:k]>=0)
        @constraints(m, begin
        constrain[i=1:p], sum(Astd[i,j]*x[j] for j=1:QtdComp+2)+sum(Astd[i,k]*y[k] for k=QtdVariaveis-2*Qdiscre+1:QtdVariaveis)+sum(Astd[i,l]*z[l] for l=QtdVariaveis+1:k)+sum(Astd[i,h]*dual[h] for h=QtdComp+3:2*QtdComp+3)+sum(Astd[i,u]*xc[u] for u=2*QtdComp+4:QtdVariaveis-2*Qdiscre)<= b[i]
        end)
        @objective(m, Max, sum(cstd[j]*x[j] for j=1:QtdComp+2)+sum(cstd[k]*y[k] for k=QtdVariaveis-2*Qdiscre+1:QtdVariaveis)+sum(cstd[h]*dual[h] for h=QtdComp+3:2*QtdComp+3)+sum(cstd[u]*xc[u] for u=2*QtdComp+4:QtdVariaveis-2*Qdiscre))


        sol = solveMIP(m)
        @test getobjectivevalue(m) ≈ 90  atol = exp10(-5)
        # vc tem que produzir 4.5, confiram
    end
    setoutputs!(m,solution,testresult)
    return solution
end


 #adicionado por Rodrigo Villas
function test_rv_7(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "PL da minha cabeça" begin

        @variable(m, x[i=1:4]>=0)

        @constraint(m, x[1]+x[2]+x[3]<=3)
        @constraint(m, x[4]+2*x[1]+6*x[3]<=10)
        @constraint(m, 4*x[3]+x[1]+3*x[2]<=5)

        @objective(m, Max, 4*x[1]+5*x[2]+2*x[3]-3*x[4])

        solveMIP(m)
        @test getobjectivevalue(m) ≈ 13
        @test getvalue(x) ≈ [2, 1, 0, 0]

    end
    setoutputs!(m,solution,testresult)
    return solution
end

 #adicionado por Rodrigo Villas
function test_rv_8(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "PL da minha cabeça" begin

        @variable(m, x[i=1:4]>=0)

        @constraint(m, x[1]+x[2]+x[3]<=3)
        @constraint(m, x[4]+2*x[1]+6*x[3]<=10)
        @constraint(m, 4*x[3]+x[1]+3*x[2]<=5)
        @constraint(m, x[3]==2)
        @objective(m, Max, 4*x[1]+5*x[2]+2*x[3]-3*x[4])

        solveMIP(m)
        @test m.ext[:status] == :Infeasible

    end
    setoutputs!(m,solution,testresult)
    return solution
end

function test_optimal_dispatch(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Despacho ótimo" begin
        # Conjuntos

        T = collect(1:2) # Periodos
        L =collect(1:2) # conjunto de linhas
        N = collect(1:2) # conjunto de barras
        I = collect(1:5) # conjunto de usinas
        C = collect(1:7) # conjunto de contingencias

        # Constantes

        x = 1 # reatancia
        c = [10;20;50;200;300] # custo de geração
        u = [10;15;5;20;25] # custo da reserva de subida
        d = [10;15;5;1;1] # custo da reserva de descida
        CU = [100;100;50;10;10] # custo fixo de acoplamento
        CD = CU # custo fixo de desacoplamento
        Pmin = [50;40;10;5;1] # Produção minima de cada usinas
        Pmax = [80;70;60;60;60] # produção maxima de cada usina
        RU = [5;15;25;50;60] # Rampa de subida
        RD = [5;15;50;50;60] # Rampa de descida
        F = [50;50] # capacidade de cada linha
        K = 1 # criterio de segurança

        # Lei de Kirshoff

        A = [1  1;
            -1 -1] # Matriz de incidencia

        S = (A*(1/x)*diagm(ones(length(L))))'

        # Balanço de potencia

        B = [1 1 0 0 1;
            0 0 1 1 0]

        # Demanda

        D = [50;50;50;60;70;80;81;82;83;84;85;80;75;70;80;90;100;100;100;90;80;70;60;50]
        Db=0.5*(D.*ones(24,2))'
        # Contingencias

        a = ones(7,7)-diagm(ones(7))

        # Estado inicial do sistema

        p0 = zeros(5)
        v0 = zeros(5)

        #------------------------------------------------------------------------------
        # Formulação do modelo
        #------------------------------------------------------------------------------

        # Variáveis

        @variable(m, p[1:length(I),1:length(T)] >= 0)
        @variable(m, v[1:length(I),1:length(T)], Bin)
        @variable(m, f[1:length(L),1:length(T)])
        @variable(m, θ[1:length(N),1:length(T)])
        @variable(m, ru[1:length(I),1:length(T)] >= 0)
        @variable(m, rd[1:length(I),1:length(T)] >= 0)
        @variable(m, pk[1:length(I),1:length(T),1:7] >= 0)
        @variable(m, fk[1:length(L),1:length(T),1:7])
        @variable(m, θk[1:length(N),1:length(T),1:7])
        @variable(m, 0 <= vu[1:length(I),1:length(T)] <= 1)
        @variable(m, 0 <= vd[1:length(I),1:length(T)] <= 1)

        # Função objetivo

        @objective(m, Min,sum(c[i]*p[i,t]+CU[i]*vu[i,t]+CD[i]*vd[i,t]+u[i]*ru[i,t]+d[i]*rd[i,t] for i in I, t in T))

        # Restições

        @constraint(m, [i in I,t in 2:T[end]], vu[i,t] - vd[i,t] == v[i,t] - v[i,t-1])
        @constraint(m, [i in I,t in T], vu[i,t] <= v[i,t])
        @constraint(m, [i in I,t in T], vd[i,t] <= 1- v[i,t])
        @constraint(m, [i in I,t in T], Pmin[i]*v[i,t] <= p[i,t])
        @constraint(m, [i in I,t in T], p[i,t] <= Pmax[i]*v[i,t])
        @constraint(m, [i in I, t in T[2:end]], p[i,t] - p[i,t-1] <= RU[i]*v[i,t-1] + Pmax[i]*vu[i,t])
        @constraint(m, [i in I, t in T[2:end]], p[i,t-1] - p[i,t] <= RD[i]*v[i,t] + Pmax[i]*vd[i,t])
        @constraint(m, pre[b in N,t in T], sum(A[b,l]*f[l,t] for l in L) + sum(B[b,i]*p[i,t] for i in I) == Db[b,t])
        @constraint(m, [l in L,t in T], f[l,t] == sum(S[l,b]*θ[b,t] for b in N))
        @constraint(m, [l in L,t in T], -F[l] <= f[l,t])
        @constraint(m, [l in L,t in T], f[l,t] <= F[l])
        @constraint(m, [i in I,t in T], p[i,t] + ru[i,t] <= Pmax[i]*v[i,t])
        @constraint(m, [i in I,t in T], p[i,t] - rd[i,t] >= Pmin[i]*v[i,t])
        @constraint(m, [i in I,t in T], ru[i,t] <= 0.9*RU[i])
        @constraint(m, [i in I,t in T], rd[i,t] <= 0.9*RD[i])


        if K!=0
            @constraint(m, pos[b in N,t in T,k in C], sum(A[b,l]*fk[l,t,k] for l in L) + sum(B[b,i]*pk[i,t,k] for i in I) == Db[b,t])
            @constraint(m, [l in L,t in T,k in C], fk[l,t,k] == a[k,l+5]*sum(S[l,b]*θk[b,t,k] for b in N))
            @constraint(m, [l in L,t in T,k in C], -F[l]*a[k,l+5] <= fk[l,t,k])
            @constraint(m, [l in L,t in T,k in C], fk[l,t,k] <= F[l]*a[k,l+5])
            @constraint(m, [i in I,t in T,k in C], a[k,i]*(p[i,t]-rd[i,t]) <= pk[i,t,k])
            @constraint(m, [i in I,t in T,k in C], pk[i,t,k] <= a[k,i]*(p[i,t]+ru[i,t]))
        end


        @constraint(m, [i in I,t in 1], vu[i,t] - vd[i,t] == v[i,t] - v0[i])
        @constraint(m, [i in I, t in 1], p[i,t] - p0[i] <= RU[i]*v0[i] + Pmax[i]*vu[i,t])
        @constraint(m, [i in I, t in 1], p0[i] - p[i,t] <= RD[i]*v[i,t] + Pmax[i]*vd[i,t])
        sol = solveMIP(m)

        @test getobjectivevalue(m) == 5150.0
        @test sum(getvalue(p[2,:])) == 98.0
        @test sum(getvalue(v[1,:])) == 0.00
        @test sum(getvalue(v)) == 4.00
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#PL Unbounded
#Adicionado por Guilherme Bodin
function test_PL_Unbounded_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste PL Unbounded Guilherme" begin
        m = Model(solver = solver)
        @variable(m, x[i=1:3])
        @objective(m, Min, x[1]+x[2] + x[3])

        solveMIP(m)
        @test m.ext[:status] == :Unbounded
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste MIP Unbounded
#adicionado por Guilherme Bodin
function test_MIP_Unbounded_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Unbounded Guilherme " begin
        number_of_nodes = 7

        C = [0.0    135.484   142.801   131.0     117.154   153.473   201.022
             135.484    0.0     105.546   174.003   142.425    53.0094  105.991
             142.801  105.546     0.0      87.6641   59.0      63.8905   73.5527
             131.0    174.003    87.6641    0.0      31.9061  146.932   159.201
             117.154  142.425    59.0      31.9061    0.0     115.521   132.306
             153.473   53.0094   63.8905  146.932   115.521     0.0      55.4617
             201.022  105.991    73.5527  159.201   132.306    55.4617    0.0   ]

        ans = [0.0   0.0   0.0   1.0   0.0   0.0   0.0
               1.0   0.0   0.0   0.0   0.0   0.0   0.0
               0.0   0.0   0.0   0.0   0.0   0.0   1.0
               0.0   0.0   0.0   0.0   1.0   0.0   0.0
               0.0   0.0   1.0   0.0   0.0   0.0   0.0
               0.0   1.0   0.0   0.0   0.0   0.0   0.0
               0.0   0.0   0.0   0.0   0.0   1.0   0.0]


        @variable(m, X[i=1:number_of_nodes,j=1:number_of_nodes], Bin)
        @variable(m, u[i=:1:number_of_nodes])#Int
        @objective(m, Min, sum(C[i,j]*u[i] for i=1:number_of_nodes, j=1:number_of_nodes))
        solveMIP(m)
        @test m.ext[:status] == :Unbounded
    end
    setoutputs!(m,solution,testresult)
    return solution
end


#teste MIP Infeasible Minimal
#adicionado por Guilherme Bodin
function test_MIP_Infeasible_Minimal_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Infeasible Minimal Guilherme" begin
        m = Model(solver = solver)
        @variable(m, x[i=1:5], Bin)
        @constraint(m, x[1] == 1)
        @constraint(m, x[2] == 1)
        @constraint(m, x[3] == 1)
        @constraint(m, x[4] == 1)
        @constraint(m, x[5] == 1)
        @constraint(m, x[1] + x[2] + x[3] + x[4] + x[5] <= 4)
        @objective(m, Min, x[1]+x[2])

        solveMIP(m)
        @test m.ext[:status] == :Infeasible
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#teste MIP Infeasible Pequeno
#adicionado por Guilherme Bodin
function test_MIP_Infeasible_Pequeno_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Infeasible Guilherme " begin
        number_of_nodes = 7

        C = [0.0    135.484   142.801   131.0     117.154   153.473   201.022
             135.484    0.0     105.546   174.003   142.425    53.0094  105.991
             142.801  105.546     0.0      87.6641   59.0      63.8905   73.5527
             131.0    174.003    87.6641    0.0      31.9061  146.932   159.201
             117.154  142.425    59.0      31.9061    0.0     115.521   132.306
             153.473   53.0094   63.8905  146.932   115.521     0.0      55.4617
             201.022  105.991    73.5527  159.201   132.306    55.4617    0.0   ]


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
        @constraint(m, sum(X[i,j] for i=1:number_of_nodes, j=1:number_of_nodes) == 50)
        @objective(m, Min, sum(C[i,j]*X[i,j] for i=1:number_of_nodes, j=1:number_of_nodes))

        solveMIP(m)
        @test m.ext[:status] == :Infeasible
    end
    setoutputs!(m,solution,testresult)
    return solution
end

#PL
#Adicionado por Guilherme Bodin
function test_PL_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste PL Unbounded Guilherme" begin
        m = Model(solver = solver)
        @variable(m, x[i=1:50])
        @constraint(m, sum(x[i] for i=1:50) >= 30)
        @objective(m, Min, sum(x[i] for i=1:50))
        solveMIP(m)
        @test m.objVal == 30
    end
    setoutputs!(m,solution,testresult)
    return solution
end




function test_TSPmip60_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Médio Guilherme (TSP 50 cidades)" begin
        number_of_nodes = 60
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

        solveMIP(m)
    end
    setoutputs!(m,solution,testresult)
    return solution
end

function test_TSPmip70_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Médio Guilherme (TSP 50 cidades)" begin
        number_of_nodes = 70
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

        solveMIP(m)
    end
    setoutputs!(m,solution,testresult)
    return solution
end

function test_TSPmip80_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Médio Guilherme (TSP 50 cidades)" begin
        number_of_nodes = 80
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

        solveMIP(m)
    end
    setoutputs!(m,solution,testresult)
    return solution
end

function test_TSPmip90_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Médio Guilherme (TSP 50 cidades)" begin
        number_of_nodes = 90
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

        solveMIP(m)
    end
    setoutputs!(m,solution,testresult)
    return solution
end

function test_TSPmip100_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    solution = MIPSolution()
    m = Model(solver = solver)
    testresult = @testset "Teste MIP Médio Guilherme (TSP 50 cidades)" begin
        number_of_nodes = 100
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

        solveMIP(m)
    end
    setoutputs!(m,solution,testresult)
    return solution
end
