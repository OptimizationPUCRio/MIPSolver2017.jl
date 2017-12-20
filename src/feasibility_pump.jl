module feasibility_pump

using JuMP

#essa função recebe um problema relaxado e os indices das variaveis que deveriam ser binarias
#apartir disso retorna uma solução viável para o problema, no caso de min isso nos da um upperbound
function fpump(mod::JuMP.Model, binI::Vector{Int64})
  m=deepcopy(mod)
  obf = JuMP.prepAffObjective(m)
  tam=m.numCols
  n=length(binI)

  v = Variable.(m, 1:tam)

  # conferir se o b&b da solve no modelo antes, se n der precisa fazer aqui
  solve(m)
  xrel=m.colVal

  xint=roundbin(xrel,binI,tam)

  if xint == false
    for i in binI
      @constraint(m,0 <= v[i])
      @constraint(m, v[i] <= 1)
    end
    xint=roundbin(xrel,binI,tam)
  end


  d=dist(xrel,xint,binI)

  cont=0
  while d>1e-7 && cont < 1e3
    cont = cont + 1

    @objective(m, Min, sum{ ifelse(xint[i] == 0  , v[i] , 0) + ifelse(xint[i] == 1  , xint[i] - v[i], 0), i in binI})

    solve(m)
    xrel = m.colVal
    d1 = m.objVal

    if d1 ≈ d
      #se o xint for igual ao anterior a distância permanecerá a mesma então criamos uma perturbação

      xint=mudaround(xrel,xint,binI,tam,n)
      d = dist(xrel,xint,binI)

    else
      xint=roundbin(xrel,binI,tam)
      d=d1
    end
  end

  objval=obf'*xrel

  if d <= 1e-7
    return xrel, objval
  end
  return false, false

end


function dist(x,xint,binI)
  dist=0
  for i  in binI
    if xint[i] ≈ 0
      soma = x[i] - 0
    elseif xint[i] ≈ 1
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
      if xint[i] != 1 && xint[i] != 0
        return false
      end
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

  #preenche um vetor de scores (dif entre xrel e xint das variaveis binarias)
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


end
