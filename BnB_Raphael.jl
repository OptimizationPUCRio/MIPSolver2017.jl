# ------------------------------------------------------------------
# Includes para os outros arquivos do resto da turma
include("funcoes_relax.jl")


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

## Receives node and creates two children by setting a variable to 0 and 1 respectively
function branch(currentNode::node, binaryIndices::Vector{Int64})

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

function obtainBoundList(nodeList::Vector{node})
  boundList = Array{Float64}(length(nodes))
  for i = 1 : length(nodes)
    boundList[i] = nodes[i].model.objVal
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
  status = solve(m, relaxation=true)
  if status == :Optimal && isBinary(m, binaryIndices)
    # Solution of the relaxed problem is binary: optimal solution
    nodes = Vector{node}(0)
    bestBound = m.objVal
    bestVal = m.objVal
    binarySolutions = 1
  else
    # Create branch and bound tree
    nodes = Vector{node}(1)
    nodes[1] = node(0, m) # root
    lastNodeLevel = 0
  end

  iter = 0
  flagOpt = 0
  tol = 0.01 # Tolerance (%)
  time0 = time_ns()
  while !isempty(nodes) && abs((bestVal - bestBound)/bestVal) > tol && (time_ns()-time0)/1e9 < 600
    if iter == 0
      iter = 1
      continue
    end
    # Check node lower bound. If greater than current best UB, prune by limit
    if nodes[1].model.objVal <= bestVal
      status = solve(nodes[1].model)
      if status == :Optimal
        bestBound = minimum(obtainBoundList(nodes))
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
          (leftChild, rightChild) = branch(nodes[1], binaryIndices)
          push!(nodes, leftChild)
          push!(nodes, rightChild)
        end
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
