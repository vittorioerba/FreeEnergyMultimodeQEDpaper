"""
    findMinimum(p::Int, q::Int)
    findMinimum(p::Int, q::Int; optimizer::Symbol) where T <: Real

Finds a minimum of modelFreeEnergy in the box ```[-1,1]^modelDimension(q)``` at parameter (p,q).

Returns a Tuple with (p, q, temperature, optimizer, minimum configuration, starting configuration, timestamp , timespan of optimization, status).

Optional arguments:
-   ```optimizer```: ```:gradient``` for gradient descent, ```:ipopt``` for interior point method, ```:cplex``` for cplex optimizer. Default value: ```:gradient```.

```:gradient``` finds a local minimum.
Optional arguments for ```:gradient```.
-   ```startingCondition```: initialization for the Ipopt optimizer. Default value: random point in the box ```[-1,1]^modelDimension(q)```;
-   ```temperature```: temperature of the model. Default value: ```0```;
-   ```step```: step-size of gradient descent update. Default value: ```1e-3```;
-   ```maxIter```: maximum number of gradient descent steps allowed. Default value: ```10^6```;
-   ```threshold```: threshold for the norm of the gradient. If ```norm(grad) < threshold```, the gradient is considered null. Default value: 1e-3. 

**CURRENTLY BROKEN!!**
```:ipopt``` finds a local minimum.
Optional arguments for ```:ipopt```:
-   ```verbose```: if ```true``` activates JuMP verbose mode. Default value: ```false```;
-   ```startingCondition```: initialization for the Ipopt optimizer. Default value: random point in the box ```[-1,1]^modelDimension(q)```;
-   ```temperature```: temperature of the model. Default value: ```0```.

```:cplex``` finds the global minimum at zero temperature.
Optional arguments for ```:cplex```:
-   ```verbose```: if ```true``` activates JuMP verbose mode. Default value: ```false```;
-   ```size```: to solve the problem, magnetizations are quantized in 1/size steps. Default value: ```100```.


"""
function modelFindMinimum(p::Int64,q::Int64;
                     optimizer::Symbol =:gradient,
                     kwargs...
                     ) where T <: Real

    id = Dates.format(now(), "yymmddHHMMSS")
    result = Tuple{Vector{Float64},Vector{Float64}}

    if optimizer == :gradient
        threshold = 1e-3
        if haskey(kwargs, :threshold)
            threshold = kwargs[:threshold]
        end
        maxIter = 1000000
        if haskey(kwargs, :maxIter)
            maxIter = kwargs[:maxIter]
        end
        step = 1e-4
        if haskey(kwargs, :step)
            step = kwargs[:step]
        end
        temperature = 0
        if haskey(kwargs, :temperature)
            temperature = kwargs[:temperature]
        end
        startingCondition = rand(Float64,modelDimension(q)) .|> u -> 2*u-1
        if haskey(kwargs, :startingCondition)
           startingCondition = kwargs[:startingCondition]
        end

        res = @timed localMinimum_gradientDescent(p,q,
                                                  startingCondition=startingCondition,
                                                  temperature=temperature,
                                                  threshold=threshold,
                                                  maxIter=maxIter,
                                                  step=step
                                                 )
        minimum,starting,status = res[:value]
        time = res[:time]

    elseif optimizer == :ipopt
        verbose = false
        if haskey(kwargs, :verbose)
            verbose = kwargs[:verbose]
        end
        temperature = 0
        if haskey(kwargs, :temperature)
            temperature = kwargs[:temperature]
        end
        startingCondition = rand(Float64,modelDimension(q)) .|> u -> 2*u-1
        if haskey(kwargs, :startingCondition)
           startingCondition = kwargs[:startingCondition]
        end

        minimum, starting, time, status = localMinimum_interiorPoint(p,q,
                                                                     startingCondition=startingCondition,
                                                                     temperature=temperature,
                                                                     verbose=verbose
                                                                    )
    elseif optimizer == :cplex
        size = 100
        if haskey(kwargs, :size)
            size = kwargs[:size]
        end
        verbose = false
        if haskey(kwargs, :verbose)
            verbose = kwargs[:verbose]
        end
        temperature = 0

        minimum, time, status = globalMinimum(p,q,verbose=verbose, size=size)
        starting = []

    else
        println("Wrong optimizer!")
        minimum = []
        starting = []
        time = -1
    end

    return ( p, q, temperature, string(optimizer), minimum, starting, id, time, status )
end

"""
    globalMinimum(p::Int,q::Int) 
    globalMinimum(p::Int,q::Int; verbose::Bool, size::Int)

Finds the global minimum of modelFreeEnergy in the box ```[-1,1]^modelDimension(q)``` at zero temperature.

Returns a Tuple with (optimal magnetizations, time elapsed, JuMP status).

Optional arguments:
-   ```verbose```: if ```true``` activates JuMP verbose mode. Default value: ```false```;
-   ```size```: to solve the problem, magnetizations are quantized in 1/size steps. Default value: ```100```.
"""
function globalMinimum(p::Int,q::Int;
                       verbose::Bool=false, 
                       size::Int=100)

    q2 = modelDimension(q)

    # initialize quadratic interaction matrix
    J = modelInteraction(p,q)

    # set optimizer to find global minima 
    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPXPARAM_OptimalityTarget", 3)

    # set verbosity of optimizer
    if !verbose
        set_silent(model)
    end

    # create quantized variables 
    @variable(model, -size <= M[i=1:q2] <= size, Int)
    vars= JuMP.all_variables(model)

    # set quadratic objective function
    @objective(model, Min, dot(vars,J,vars))

    # optimize
    optimize!(model)

    # a priori there may be more optimal solutions
    # here we ignore this possibility

    GSconfiguration = (value.(vars))./size        
    solveTime = solve_time(model)                
    status = primal_status(model)                 |> string

    return (GSconfiguration,solveTime,status)
end


"""
    localMinimum_interiorPoint(p::Int,q::Int) 
    localMinimum_interiorPoint(p::Int,q::Int; verbose::Bool, startingCondition::Vector{T}, temperature::Real) where T <: Real

CURRENTLY BROKEN AS IT EVALUATES THE ENTROPY OUTSIDE OF THE DOMAIN BY VERY SMALL AMOUNTS

Finds a local minimum of modelFreeEnergy in the box ```[-1,1]^modelDimension(q)``` at temperature t using the interior point method provided by Ipopt.

Returns a Tuple with (optimal magnetizations, starting condition, time elapsed, JuMP status).

Optional arguments:
-   ```verbose```: if ```true``` activates JuMP verbose mode. Default value: ```false```;
-   ```startingCondition```: initialization for the Ipopt optimizer. Default value: random point in the box ```[-1,1]^modelDimension(q)```;
-   ```temperature```: temperature of the model. Default value: ```0```.
"""
function localMinimum_interiorPoint(p::Int64, q::Int64;
                                    startingCondition::Vector{T} = rand(Float64,modelDimension(q)) .|> u -> 2*u-1,
                                    temperature::Real=0,
                                    verbose::Bool=false) where T <: Real

    q2 = modelDimension(q)
    if length(startingCondition) != q2
        throw(ArgumentError("Length of startingCondition incompatible with q!"))
    end

    # initialize quadratic interaction matrix
    J = modelInteraction(p,q)

    # set optimizer
    model = Model(Ipopt.Optimizer)

    register(model, :magnetizationEntropy, 1, magnetizationEntropy, autodiff=true)

    # set verbosity of optimizer
    if !verbose
        set_silent(model)
    end

    # create continuous variables and set starting condition
    @variable(model, -1 <= M[i=1:q2] <= 1)
    for i in 1:q2
        set_start_value(M[i],startingCondition[i])
    end
    vars= JuMP.all_variables(model)

    # initialize quadratic interaction matrix
    # and compute energy and entropy
    energy = @NLexpression(model, 
                           sum(vars[i]*J[i,j]*vars[j] for i in 1:q2, j in 1:q2)
                          )
    entropy = @NLexpression(model, 
                            2*magnetizationEntropy(vars[1])
                            + 4*sum( 
                                    magnetizationEntropy(vars[i])
                                    for i in 2:q2-1)
                            + ( q%2==0 ? 2 : 4 )*magnetizationEntropy(vars[q2])
                           )

    # set objective function
    @NLobjective(model, Min, energy-temperature*entropy)

    # optimize
    optimize!(model)

    GSconfiguration = (value.(vars))
    solveTime = solve_time(model)
    status = primal_status(model) |> string

    return (GSconfiguration,startingCondition,solveTime,status)
end



"""
    localMinimum_gradientDescent(p::Int,q::Int) 
    localMinimum_gradientDescent(p::Int,q::Int; startingCondition::Vector{T}, temperature::Real, step::Real, maxIter::Int, threshold::Real) where T <: Real

Finds a local minimum of modelFreeEnergy in the box ```[-1,1]^modelDimension(q)``` at temperature t using a custom gradient descent method that respects the box boundary conditions.

Returns a Tuple with (optimal magnetizations, starting condition, status).
Status can be:
-   ```"good"```: optimization ended due to null gradient;
-   ```"max iter reached": no convergence was reached before the maximum number of allower iterations;
-   ```"stuck in position"```: the descent is stuck in a position while the gradient is non null.
Only ```"good"``` should be considered as a successful optimization.

Optional arguments:
-   ```startingCondition```: initialization for the Ipopt optimizer. Default value: random point in the box ```[-1,1]^modelDimension(q)```;
-   ```temperature```: temperature of the model. Default value: ```0```;
-   ```step```: step-size of gradient descent update. Default value: ```1e-3```;
-   ```maxIter```: maximum number of gradient descent steps allowed. Default value: ```10^6```;
-   ```threshold```: threshold for the norm of the gradient. If ```norm(grad) < threshold```, the gradient is considered null. Default value: 1e-3. 
"""
function localMinimum_gradientDescent(p::Int,q::Int;
                                      startingCondition::Vector{T} = rand(Float64,modelDimension(q)) .|> u -> 2*u-1,
                                      temperature::Real =0.,
                                      threshold::Real = 1e-3,
                                      maxIter::Int = 1000000,
                                      step::Real = 1e-3
                                     ) where T <: Real

    q2 = modelDimension(q)
    if length(startingCondition) != q2
        throw(ArgumentError("Length of startingCondition incompatible with q!"))
    end
    # set starting point and starting gradient
    old_position = startingCondition.+0.1 # add 0.1 to make it different from actual_postition before starting the while cycle
    actual_position = startingCondition
    grad = gradModelFreeEnergy(actual_position,p,q,temperature)
    iter = 0

    # while gradient not too small, keep following the steepest descent
    while norm(grad) > threshold && iter < maxIter && actual_position != old_position
        iter += 1
        # update position with gradient
        old_position = actual_position
        actual_position -= step * grad
        # make sure that the new position is in the box
        # otherwise push it back
        for i in 1:q2
            if actual_position[i] > 1
                actual_position[i] = 1
            elseif actual_position[i] < -1
                actual_position[i] = -1        
            end
        end
        # compute new gradient
        grad = gradModelFreeEnergy(actual_position,p,q,temperature)
        # if gradient descent direction (i.e. -gradient) points out of the box, put that component to zero
        for i in 1:q2
            if ( actual_position[i] == 1 && grad[i] < 0 ) || ( actual_position[i] == -1 && grad[i] > 0 )
                grad[i] = 0
            end
        end
    end

    if iter == maxIter
        status = "max iter reached"
    elseif actual_position == old_position
        status = "stuck in position"
    else
        status = "good"
    end

    return (actual_position, startingCondition, status)
end
