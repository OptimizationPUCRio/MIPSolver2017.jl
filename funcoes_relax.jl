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
        X, z, status = solve_LP(Astd, b, cstd, xlb, xub, solver, m, n)

        # Finding fractional variable address
        j = find(min(X[1:n]-floor(X[1:n]),ceil(X[1:n])-X[1:n]).>tol)

        # Removing index from continous variables
        for i in VecBin
            deleteat!(j,find(j .== i))
        end

        #Test convergence
        if length(j) == 0
            H = [H; iter z -999 status X[1:n]']
            #return X[1:n], z
            convergence = 1
        elseif iter >= MaxIter
            j = j[1]
            H = [H; iter z j status X[1:n]']
            convergence  = -1
        else

            j = j[1] #pick a fractional variable address (the first one)

            H = [H; iter z j status X[1:n]']
            #Calculating Chvatal-Gomory coefficients
            NBas = find(X.<tol)[1:n] #finding nonbasic variable addresses
            Bas = deleteat!(collect(1:(n+m)), NBas) #finding basic variable addresses
            i = find(Bas.==j) # finding dictionary line related to the fractional variable
            #D = \(Astd[:,Bas],Astd[:,NBas])
            D=((Astd[:,Bas]'*Astd[:,Bas])\Astd[:,Bas]')*Astd[:,NBas]
            q0 =  floor(X[j])
            qi = zeros(1,n+m)
            qi[j] = 1
            qi[NBas] = floor(D[i,:])

            #Adding the cut as new a constraint into the problem stand. form
            Astd = [[Astd zeros(m,1)];[qi 1]] #new constrait and new slack
            b = [b;q0]
            cstd = [cstd;0]
            m+=1
            iter+=1
        end
    end


    model = compose_model(model, m_origin, n_origin, Astd, b, cstd, xlb, xub, m, n)

    return model
end

function extract_data(model::JuMP.Model)
    A = full(JuMP.prepConstrMatrix(model))
    m, n = size(A)

    Astd = zeros(m,m+n)
    cstd = zeros(m+n)
    b = zeros(n,1)

    Astd = [A eye(m)]
    cstd = [c ; zeros(m)]

    xlb = copy(model.colLower)
    xub = copy(model.colUpper)
    rowlb, rowub = JuMP.prepConstrBounds(model)
    solver = model.solver

    for i in 1:n
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


    return Astd, b, cstd, xlb, xub, solver, m, n, flag_sense
end

function solve_LP(Astd, b, cstd, xlb, xub, solver, m, n)
    model_LP = Model(solver=solver)

    @variable(model_LP, xlb <= x[1:n] <= xub)
    @constraints(model_LP, begin
    constrain[i=1:m], sum(Astd[i,j]*x[j] for j=1:n) == b[i]
    end)

    @objective(model, :Max, sum(cstd[j]*x[j] for j=1:k))
    status = solve(model_LP);
    return model_LP.colVal, model_LP.objVal, status
end

function compose_model(model, m_origin, n_origin, Astd, b, cstd, xlb, xub, m, n)

    #@addConstraint(model, begin
    #cut[i = 1:(m-m_origin)], A)

end
