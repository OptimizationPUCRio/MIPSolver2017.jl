"""
 Module Name    : Primary solutions by heuristics functions
 Description    : The module focus on delivering a primary solution upon its
                  calling.
 Authors        : Andrew Rosemberg and Bianca Lacê

 Initial date   : 04/12/2017
"""
module feasible_solution
using JuMP


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
    rlc = m.objSense == :Max ? setdiff(sortperm(c[1 + num_nodes*(i-1):num_nodes*(i)], rev=true),tour):
        setdiff(sortperm(c[1 + num_nodes*(i-1):num_nodes*(i)]),tour)
    rlc = rlc[1:min(num_rlc,size(rlc,1))]
    ind = rlc[rand(1:size(rlc,1))]
    i = ind
    # remaining iteractions
    for  itr= 2:num_nodes-2
        tour[itr] = i
        rlc = m.objSense == :Max ? setdiff(sortperm(c[1 + num_nodes*(i-1):num_nodes*(i)], rev=true),tour) :
            setdiff(sortperm(c[1 + num_nodes*(i-1):num_nodes*(i)]),tour)
        rlc = rlc[1:min(num_rlc,size(rlc,1))]
        ind = rlc[rand(1:size(rlc,1))]
        i = ind
    end
    tour[num_nodes-1] = i
    i = setdiff(indVariable,tour)[1]
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

    max_improvs = max(num_nodes^2,max_improvs_aux)
    best_sol = m.colVal
    best_obj = obval
    improv_count = 0
    itr = 0
    # neighbourhood search
    while itr <= min(max_iter_improv,max_improvs)
        # permutation
        itr+=1
        inds =rand(1:num_nodes,2)
        tour_aux = deepcopy(tour)
        tour_aux[inds[1]] = tour[inds[2]]
        tour_aux[inds[2]] = tour[inds[1]]

        sol = feasible_solution._prep_solution(m,tour)
        status,obval = feasible_solution._check_feasability_solution(m,sol,indVariable)
        if status == true
            if m.objSense == :Max && obval > best_obj
                best_sol = m.colVal
                best_obj = obval
            elseif  m.objSense == :Min && obval < best_obj
                best_sol = m.colVal
                best_obj = obval
            end
        end
    end
    m.colVal = best_sol
    m.objVal = obval
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


"""
Function Name   : Feasability Pump
"""
function feasability_pump()

end





end
