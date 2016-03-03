# wrapper to convert Nonlinear solver into Conic solver
# The translation is lossy...
# Authors: Emre Yamangil and Miles Lubin

type NonlinearToConicBridge <: MathProgBase.AbstractConicModel
    # SOLUTION DATA
    solution::Vector{Float64}                       # Vector containing solution
    status                                          # Termination status of the nlp_solver
    objval::Float64                                 # Objective value of optimal solution

    # SOLVER DATA
    nlp_solver::MathProgBase.AbstractMathProgSolver # Choice of nonlinear solver
    remove_single_rows                              # Preprocessing singleton rows flag

    # PROBLEM DATA
    x                                               # Variables in nonlinear model
    numVar                                          # Number of variables in nonlinear model
    numConstr                                       # Number of constraints in nonlinear model
    nlp_model                                       # Reference to nonlinear model
    A_ini                                           # Initial constraint matrix
    b                                               # Right hand side vector
    constr_cones_ini                                # Initial constraint cones
    var_cones_ini                                   # Initial variable cones

    # CONSTRUCTOR
    function NonlinearToConicBridge(nlp_solver,remove_single_rows)
        m = new()
        m.nlp_solver = nlp_solver
        m.remove_single_rows = remove_single_rows
        return m
    end
end 

export ConicNLPWrapper
immutable ConicNLPWrapper <: MathProgBase.AbstractMathProgSolver
    nlp_solver::MathProgBase.AbstractMathProgSolver
    remove_single_rows
end
ConicNLPWrapper(;nlp_solver=nothing,remove_single_rows=false) = ConicNLPWrapper(nlp_solver,remove_single_rows)
MathProgBase.ConicModel(s::ConicNLPWrapper) = NonlinearToConicBridge(s.nlp_solver,s.remove_single_rows)

function MathProgBase.loadproblem!(
    m::NonlinearToConicBridge, c, A, b, constr_cones, var_cones)

    if m.nlp_solver == nothing
        error("NLP solver is not specified.")
    end

    nlp_model = Model(solver=m.nlp_solver)
    numVar = length(c) # number of variables
    numConstr = length(b) # number of constraints
    m.A_ini = A
    m.b = b
    m.constr_cones_ini = constr_cones
    m.var_cones_ini = var_cones

    # b - Ax \in K => b - Ax = s, s \in K
    new_var_cones = Any[x for x in var_cones]
    new_constr_cones = Any[]
    copy_constr_cones = copy(constr_cones)
    lengthSpecCones = 0
    # ADD SLACKS FOR ONLY SOC AND EXP
    A_I, A_J, A_V = findnz(A)
    slack_count = numVar+1
    for (cone, ind) in copy_constr_cones
        if cone == :SOC || cone == :ExpPrimal
            lengthSpecCones += length(ind)
            slack_vars = slack_count:(slack_count+length(ind)-1)
            append!(A_I, ind)
            append!(A_J, slack_vars)
            append!(A_V, ones(length(ind)))
            
            push!(new_var_cones, (cone, slack_vars))
            push!(new_constr_cones, (:Zero, ind))
            slack_count += length(ind)
        else
            push!(new_constr_cones, (cone, ind))
        end
    end
    A = sparse(A_I,A_J,A_V, numConstr, numVar + lengthSpecCones)

    m.numVar = size(A,2)
    m.numConstr = numConstr 
    c = [c;zeros(m.numVar-numVar)]

    # LOAD NLP MODEL
    @defVar(nlp_model, x[i=1:m.numVar], start = 1)
    
    @setObjective(nlp_model, Min, dot(c,x))

    for (cone, ind) in new_var_cones
        if cone == :Zero
            for i in ind
                setLower(x[i], 0.0)
                setUpper(x[i], 0.0)
            end
        elseif cone == :Free
            # do nothing
        elseif cone == :NonNeg
            for i in ind
                setLower(x[i], 0.0)
            end
        elseif cone == :NonPos
            for i in ind
                setUpper(x[i], 0.0)
            end
        elseif cone == :SOC
            @addNLConstraint(nlp_model, sqrt(sum{x[i]^2, i in ind[2:length(ind)]}) <= x[ind[1]])
            setLower(x[ind[1]], 0.0)
        elseif cone == :ExpPrimal
            @addNLConstraint(nlp_model, x[ind[2]] * exp(x[ind[1]]/x[ind[2]]) <= x[ind[3]])
            setLower(x[ind[2]], 0.0)
            setLower(x[ind[3]], 0.0)
        end
    end

    # *************** PREPROCESS *******************
    constr_cones_map = [:NoCone for i in 1:numConstr]
    for (cone, ind) in new_constr_cones
        constr_cones_map[ind] = cone
    end

    nonZeroElements = [Any[] for i in 1:numConstr] # by row
    for i in 1:length(A_I)
        push!(nonZeroElements[A_I[i]], (A_J[i], A_V[i]))
    end
    rowIndicator = [true for i in 1:numConstr]
    if m.remove_single_rows
        for i in 1:numConstr
            if length(nonZeroElements[i]) == 1
                (ind, val) = nonZeroElements[i][1]
                if constr_cones_map[i] == :Zero
                    setLower(x[ind], b[i]/val)
                    setUpper(x[ind], b[i]/val)
                elseif constr_cones_map[i] == :NonNeg
                    if val < 0.0
                        setLower(x[ind], b[i]/val)
                    else
                        setUpper(x[ind], b[i]/val)
                    end
                elseif constr_cones_map[i] == :NonPos
                    if val < 0.0
                        setUpper(x[ind], b[i]/val)
                    else
                        setLower(x[ind], b[i]/val)
                    end
                else
                    error("special cone $(constr_cones_map[i]) in constraint cones after preprocess.")
                end
                rowIndicator[i] = false
            end
        end
    end

    rowIndicator = [true for i in 1:numConstr]
    for (cone,ind) in new_constr_cones
        for i in 1:length(ind)
            if rowIndicator[ind[i]]
                if cone == :Zero
                    @addConstraint(nlp_model, A[ind[i]:ind[i],:]*x .== b[ind[i]])
                elseif cone == :NonNeg
                    @addConstraint(nlp_model, A[ind[i]:ind[i],:]*x .<= b[ind[i]])
                elseif cone == :NonPos
                    @addConstraint(nlp_model, A[ind[i]:ind[i],:]*x .>= b[ind[i]])
                else
                    error("unrecognized cone $cone")
                end
            end
        end
    end

    m.x = x
    m.numVar = numVar
    m.nlp_model = nlp_model

end

function MathProgBase.optimize!(m::NonlinearToConicBridge)
 
    m.status = solve(m.nlp_model)
    m.objval = getObjectiveValue(m.nlp_model)
    if (m.status != :Infeasible)
        m.solution = getValue(m.x)
    end   

end

MathProgBase.supportedcones(s::ConicNLPWrapper) = [:Free,:Zero,:NonNeg,:NonPos,:SOC,:ExpPrimal]

function MathProgBase.setwarmstart!(m::NonlinearToConicBridge, x) 

    x_expanded = copy(x)
    val = m.b - m.A_ini*x
    nonlinear_cones = 0
    for (cone, ind) in m.constr_cones_ini
        if cone == :SOC || cone == :ExpPrimal
            append!(x_expanded, val[ind])
            nonlinear_cones += 1
        end
    end
    m.solution = x_expanded
    setValue(m.x, m.solution)
end

function MathProgBase.freemodel!(m::NonlinearToConicBridge)
    if applicable(MathProgBase.freemodel!,m.nlp_model.internalModel)
        MathProgBase.freemodel!(m.nlp_model.internalModel)
    end
end

MathProgBase.setvartype!(m::NonlinearToConicBridge, v::Vector{Symbol}) = (m.vartype = v)

MathProgBase.status(m::NonlinearToConicBridge) = m.status
MathProgBase.getobjval(m::NonlinearToConicBridge) = m.objval
MathProgBase.getsolution(m::NonlinearToConicBridge) = m.solution
