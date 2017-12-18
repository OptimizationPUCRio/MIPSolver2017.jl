function cutting_planes(model::JuMP.Model, VecBin::Vector{Int}, MaxIter::Int64 = 5)
    tol = 1e-5

    Astd, b, cstd, xlb, xub, solver, m, n, flag_sense = extract_data(model)

    iter = 1
    convergence = -2 #not finished
    X = 0
    z = 0
    status = ""
    println(" VecBin ", VecBin)
    while convergence == -2
        X, z, status = solve_LP(Astd, b, cstd, xlb, xub, solver)
        println("Astd = ", Astd)
        println("b = ", b)
        println("cstd = ", cstd)
        println("xlb = ", xlb)
        println("xub = ", xub)
        println("x = ", X)
        println("m = " , m)
        println("n = " , n)

        # Finding fractional variable address
        j = find(min.(X[1:n]-floor.(X[1:n]),ceil.(X[1:n])-X[1:n]).>tol)
        println("j fracionario = ", j)
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
        println("j binario e fracionario = ", j)

        #Test convergence
        if length(j) == 0

            #return X[1:n], z
            convergence = 1
        elseif iter >= MaxIter
            j = j[1]

            convergence  = -1
        else
            j = j[1] #pick a fractional and binary variable address (the first one)
            print("j escolhido = ", j)
            #Calculating Chvatal-Gomory coefficients
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

            N1 = intersect(deleteat!(copy(VecBin),j),NBas)
            N2 = intersect(deleteat!(collect(1:n), VecBin),NBas)

            f0 = a0[u] - floor.(a0[u])
            f0 = f0[1]


            f = a[u,union(N1,N2)] - floor.(a[u,union(N1,N2)])

            A_cut = zeros(1,n)
            A_cut[intersect(find(f .<= f0), N1)] = f[intersect(find(f .<= f0), N1)]
            A_cut[intersect(find(f .> f0), N1)] = f0*(1-f[intersect(find(f .> f0), N1)])/(1-f0)
            A_cut[intersect(find(a[u,:] .> 0), N2)] = a[u,intersect(find(a[u,:] .> 0), N2)]
            A_cut[intersect(find(a[u,:] .< 0), N2)] = f0/(1-f0) * a[u,intersect(find(a[u,:] .< 0), N2)]

            b_cut = f0

            Astd = [[Astd zeros(m,1)]; [A_cut -1]]
            b = [b ; b_cut]

            cstd = [cstd ; 0]
            xlb = [xlb; 0]
            xub = [xub; +Inf]
            m+=1
            n+=1
            iter+=1
            println(" ")
            println(" ")
        end
    end
    println(" ")
    println(" ")

    println("X = ", X)
    println("Z = ", z)
    println("Status = ", status)
    print("Iter = " , iter)
    model_F = compose_model(Astd, b, cstd, xlb, xub, flag_sense, solver, X, z, status)

    return model_F
end

function extract_data(model::JuMP.Model)
    A = full(JuMP.prepConstrMatrix(model))
    c = JuMP.prepAffObjective(model)
    m, n = size(A)

    Astd = zeros(m,m+n)
    cstd = zeros(m+n)
    b = zeros(m,1)

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
