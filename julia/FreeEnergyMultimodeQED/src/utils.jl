# added absolute value to avoid numerical problems with small negative x
function xlog(x)
    x == 0 ? 0 : x * log(abs(x))
end

# generate a list of fractions [numerator, denominator] between [pMin,qMin] and [pMax,qMax]
function coprimeFractionsGenerator(pMin,qMin,pMax,qMax)
    result = []
    for q in qMin:qMax
        pmin = ( q==qMin ? pMin : 1   )
        pmax = ( q==qMax ? pMax : q-1 )
        for p in pmin:pmax
            if gcd(p,q)==1
                push!(result,[p,q])
            end
        end
    end
    result
end

# round vector components
function roundVector(vec;digits::Int64=4)
    map( u->round(u,digits=digits) , vec)
end

# noise compatible with boundary conditions
function addBoxNoise(vars;
                  eta::Float64=1e-2)
    q = length(vars)
    result = vars .+ eta * ( 2*rand(Float64,q).-1 ) 
    for i in 1:q
        r = result[i]
        if r <= -1
            result[i] = -1.
        elseif r >= 1
            result[i] = 1.
        end
    end
    result
end

# flips or not the sign of the vector so that the first non null component is positive
function flipZ2(vars)
    if length(vars)==0
        return vars
    elseif vars[1] == 0
        return vcat(0,flipZ2(vars[2:end]))
    elseif vars[1] < 0
        return -vars
    else
        return vars
    end
end
