using JuMP


function cutting_planes(model::JuMP.Model, VecBin::Vector{Int}, MaxIter::Int64 = 5)
    tol = 1e-5

    Astd, b, cstd, xlb, xub, solver, m, n, flag_sense = extract_data(model)

    iter = 1
    convergence = -2 #not finished
    H = zeros(1,4+n) #Iterations history
    X = 0
    z = 0
    status = ""

    while convergence == -2
        X, z, status = solve_LP(Astd, b, cstd, xlb, xub, solver)

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

            #return X[1:n], z
            convergence = 1
        elseif iter >= MaxIter
            j = j[1]

            convergence  = -1
        else
            j = j[1] #pick a fractional variable address (the first one)

            #Calculating Chvatal-Gomory coefficients
            NBas = find(X.<=tol) #finding nonbasic variable addresses
            Bas = deleteat!(collect(1:n), NBas) #finding basic variable addresses
            i = find(Bas.==j) # finding dictionary line related to the fractional variable

            D=((Astd[:,Bas]'*Astd[:,Bas])\Astd[:,Bas]')*Astd[:,NBas]
            q0 =  floor.(X[j])
            qi = zeros(1,n)
            qi[j] = 1
            qi[NBas] = floor.(D[i,:])

            #Adding the cut as new a constraint into the problem stand. form
            Astd = [[Astd zeros(m,1)]; [qi 1]] #new constrait and new slack
            b = [b;q0]
            cstd = [cstd;0]
            xlb = [xlb; 0]
            xub = [xub; +Inf]
            m+=1
            n+=1
            iter+=1
        end
    end

    println(X)
    println(z)
    println(status)
    print(iter)
    model_F = compose_model(Astd, b, cstd, xlb, xub, flag_sense, solver, X, z, status)

    return model_F
end

function extract_data(model::JuMP.Model)
    A = full(JuMP.prepConstrMatrix(model))
    c = JuMP.prepAffObjective(model)
    m, n = size(A)

    Astd = zeros(m,m+n)
    cstd = zeros(m+n)
    b = zeros(n,1)

    Astd = [A eye(m)]
    cstd = [c ; zeros(m)]

    xlb = copy(model.colLower)
    xub = copy(model.colUpper)

    nxlb = zeros(n+m,1)
    nxub = zeros(n+m,1)
    nxlb[1:size(xlb)[1]] = xlb
    nxlb[size(xlb)[1]+1 : n+m] = 0
    nxub[1:size(xlb)[1]] = xub
    nxub[size(xub)[1]+1 : n+m] = +Inf


    rowlb, rowub = JuMP.prepConstrBounds(model)
    solver = model.solver

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

function compose_model(Astd, b, cstd, xlb, xub, flag_sense, solver, X, z, status)
    m, n = size(Astd)

    model_F = Model(solver=solver)

    @variable(model_F, xlb[i] <= x[i=1:n] <= xub[i])
    @constraints(model_F, begin
    constrain[i=1:m], sum(Astd[i,j]*x[j] for j=1:n) == b[i]
    end)

    if flag_sense == 1
        cstd = -cstd
        @objective(model_F, :Min, sum(cstd[j]*x[j] for j=1:n))
    else
        @objective(model_F, :Max, sum(cstd[j]*x[j] for j=1:n))
    end

    model_F.colVal = X
    model_F.objVal = z
    model_F.ext[:status] = status


    return model_F
end
