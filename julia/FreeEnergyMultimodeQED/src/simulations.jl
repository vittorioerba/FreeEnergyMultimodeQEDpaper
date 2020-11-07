
function simulateGS(p,q;file=dataPath*"/raw/local_minima.csv")
    open(file, "a") do io
        writedlm(io, [findMinimum(p,q,optimizer=:cplex)], '|')
    end
end

# simulate minima using gradient descent, and appends them in minima.csv
function simulateMinima(p,q,iter;file=dataPath*"/raw/local_minima.csv")
    result = []
    prog = Progress(iter,1)
    @threads for i in 1:iter
        push!(result, findMinimum(p,q))
        next!(prog)
    end
    open(file, "a") do io
        writedlm(io, result, '|')
    end
end


# simulate heating and cooling
function heatAndCool(p,q,ts;eta::Float64=1e-2)

    # record energies
    result = []

    # compute ground state at T=0
    GSconfiguration = findMinimum(p,q; optimizer=:cplex)[5]
    push!(result, [0., GSconfiguration])

    # change the temperature of the system and minimize it locally
    @showprogress for t in ts
        state = findMinimum(p,q; optimizer=:gradient,
                            temperature=t,
                            startingCondition=addNoise(result[end][2];eta=eta)
                           )[5]
        push!(result,[t,state])
    end

    return result
end


