using JuMP

#essa função recebe um modelo relaxado e os indices das variaveis que deveriam ser binarias
#apartir disso retorna uma solução viável para o problema, no caso de min isso nos da um upperbound
function fpump(mod::JuMP.Model, binI::Vector{Int64})
  m=deepcopy(mod)
  tam=m.numCols
  n=length(binI)

  # conferir se o b&b da solve no modelo antes, se n der precisa fazer aqui
  solve(m)
  xotim=m.colVal

  xint=roundbin(xotim,binI,tam)

  v = Variable.(m, 1:tam)

  d=dist(xotim,xint,binI)

  cont=0
  while d>1e-7 && cont < 1e3
    cont = cont + 1

    @objective(m, Min, sum{ ifelse(xint[i] == 0  , v[i] , 0) + ifelse(xint[i] == 1  , xint[i] - v[i], 0), i in binI})

    solve(m)
    xotim=m.colVal
    d1 = dist(xotim,xint,binI)
    if d1 == d
      xint=mudaround(xint,binI,tam,n)
      d = dist(xotim,xint,binI)
    else
      xint=roundbin(xotim,binI,tam)
      d=d1
    end
  end

  if d <= 1e-7
    return true, xotim
  end
  return false, 0

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

function mudaround(xint::Vector, binind::Vector{Int64}, tam::Int64, n::Int64)
  quant=Int.(ceil(n/400))
  ind = rand(binind,quant)
  for i  in ind
    if xint[i] == 1
      xint[i] = 0
    else
      xint[i] = 1
    end
  end
  return xint
end
