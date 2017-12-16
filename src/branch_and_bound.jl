# ------------------------------------------------------------------
# Includes para os outros arquivos do resto da turma
#include("funcoes_relax.jl")
#include("feasible_solution.jl")



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

## Strong branching
function strong(currentNode::node, binaryIndices::Vector{Int64}, amountOfBranches::Int64)

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
    bounds[1,i] = leftModel.objVal
    solve(rightModel)
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

  leftChild = node(currentNode.level+1, leftModel)
  rightChild = node(currentNode.level+1, rightModel)

  return leftChild, rightChild
end

## Most fractioned branching
function fractioned(currentNode::node, binaryIndices::Vector{Int64})

  distance = abs(currentNode.model.colVal[binaryIndices] - 0.5)
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

## Receives node and creates two children by setting a variable to 0 and 1 respectively
function branch(currentNode::node, binaryIndices::Vector{Int64}, method::Symbol)

  if method == :fractioned
    leftChild, rightChild = fractioned(currentNode, binaryIndices)
  elseif method == :strong
    amountOfBranches = 10
    leftChild, rightChild = strong(currentNode, binaryIndices, amountOfBranches)
  else
    println("Error on branching method defition")
  end

  return leftChild, rightChild
end

function obtainBoundList(nodeList)
  boundList = Array{Float64}(length(nodeList))
  for i = 1 : length(nodeList)
    boundList[i] = nodeList[i].model.objVal
  end

  return boundList
end

## Receives a mixed binary linear JuMP model
function solveMIP(m::JuMP.Model)

  tic()

  # Check if model is max; if it is, converts to min
  flagConverted = 0
  if m.objSense == :Max
    convertSense!(m)
    flagConverted = 1
  end

  # Best bounds: start as Inf
  bestBound = -1e200
  bestVal = 1e200

  # Create vector of indices of the binary variables
  binaryIndices = find(m.colCat .== :Bin)
  binarySolutions = 0

  # Solve linear relaxation
  m.colCat[:] = :Cont
  status = solve(m)
  nodes = Vector{node}(0)
  if status == :Optimal && isBinary(m, binaryIndices)
    # Solution of the relaxed problem is binary: optimal solution
    bestBound = m.objVal
    bestVal = m.objVal
    binarySolutions = 1
  else
    push!(nodes, node(0, m)) # Add root to branch and bound tree
    lastNodeLevel = 0
  end

  iter = 1 # number of visited nodes
  flagOpt = 0 # flag that indicates if a viable solution has been found
  branched = false # flag that indicates if a branch has occurred in this iteration
  traverse = 1 # 1 for breadth, -1 for depth
  tol = 0.01 # tolerance (%)
  time0 = time_ns()

  while !isempty(nodes) && abs((bestVal - bestBound)/bestVal) > tol && (time_ns()-time0)/1e9 < 600

    # Change traverse method every 10 iterations for better bound discovery
    if iter%10 == 0
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
          leftChild, rightChild = branch(nodes[1], binaryIndices, :strong)
          branched = true
        end
      end
    end

    # Check if branched -- if true, add children to the list
    if branched == true
      if traverse == 1 # breadth -- insert children at the end of the list
        push!(nodes, leftChild, rightChild)
      else # depth -- insert children at the beginning of the list
        unshift!(nodes, leftChild, rightChild)
      end
    end
    lastNodeLevel = nodes[1].level
    deleteat!(nodes, 1)
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

  return m.ext[:status]
end
