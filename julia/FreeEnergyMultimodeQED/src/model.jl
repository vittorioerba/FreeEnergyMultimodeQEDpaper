"""
    modelDimension(q::Int)

Returns the minimal dimension of the model at parameter q.
"""
function modelDimension(q::Int)
    q/2 + 1 |> floor |> Int64
end

"""
    modelInteraction(p::Int,q::Int)

Returns the interaction matrix J of the model at parameter (p,q).
Notice:
- some terms appear with a factor 1/2 with respect to Equation (S8). This is because Equation (S8) is not completely symmetric in the magnatizations;
- ```dot(ms,J,ms)``` is the energy of Equation (S8) divided by ell.
"""
function modelInteraction(p::Int,q::Int)
    q2 = modelDimension(q)
    J = zeros(q2,q2)
    J[1,1] = 1
    for i in 2:q2
        J[1,i] = 2
        J[i,1] = 2

        for j in i:q2
            J[i,j] = 4 * cos( 2 * pi * p / q * (i-1) * (j-1) )
            J[j,i] = 4 * cos( 2 * pi * p / q * (i-1) * (j-1) )
        end
    end

    if q%2==0
        J[q2,q2] = cos(pi*p*q/2)
        J[1,q2] = 1
        J[q2,1] = 1

        for j in 2:q2-1
            J[j,q2] = 2 * cos( pi * p * (j-1) )
            J[q2,j] = 2 * cos( pi * p * (j-1) )
        end
    end

    J /= q
end

"""
    magnetizationEntropy(m::Real)

Returns the entropy per degree of freedom associated to a magnetization m; see Equation (S9)
"""
function magnetizationEntropy(m::Real)
    if m < -1 || m > 1
        throw(ArgumentError(string("m = ", m, " outside of [-1,1] domain!")))
    else
        log(2) - 1/2*xlog(1-m) - 1/2*xlog(1+m)
    end
end

"""
    gradMagnetizationEntropy(m::Real)
    gradMagnetizationEntropy(m::Real; regularizer::Float64)

Returns the derivative of the entropy per degree of freedom associated to a magnetization m; see Equation (S9).
Optional arguments:
-   ```regularizer```: small constant to deal with the logarithmic divergences at m=1 and m=-1. Default value ```1e-3```.
"""
function gradMagnetizationEntropy(m::Real;
                                  regularizer::Real=1e-3)
    if m < -1 || m > 1
        throw(ArgumentError("m outside of [-1,1] domain!"))
    else
        1/2*log(1-m+regularizer) - 1/2*log(1+m+regularizer)
    end
end

"""
    modelEntropy(ms::Vector{T}, q::Int) where T <: Real

Returns the entropy of the model (divided by ell) at magnetizations ms, and at parameter q; see Equation (S11).
"""
function modelEntropy(ms::Vector{T},q::Int) where T <: Real
    q2 = length(ms)
    if q2 == modelDimension(q)
        entropy = 2*magnetizationEntropy(ms[1]) 
        if q2 - 1 >= 2
            entropy += 4*sum( magnetizationEntropy(ms[i]) for i in 2:q2-1)
        end
        return entropy += ( q%2==0 ? 2 : 4 )*magnetizationEntropy(ms[q2])
    else
        throw(ArgumentError("Length of ms incompatible with q!"))
    end
end

"""
    modelFreeEnergy(ms::Vector{T}, p::Int, q::Int, t::Real) where T <: Real

Returns the free energy of the model (divided by ell) at magnetizations ms, at parameters (p,q) and at temperature t; see Equation (S8) and (S11).
"""
function modelFreeEnergy(ms::Vector{T},p::Int,q::Int,t::Real) where T <: Real
    q2 = length(ms)
    if q2 == modelDimension(q)
        J = modelInteraction(p,q)
        return dot(ms,J,ms) - t * modelEntropy(ms,q)
    else
        throw(ArgumentError("Length of ms incompatible with q!"))
    end
end


"""
    gradModelFreeEnergy(ms::Vector{T}, p::Int, q::Int, t::Real) where T <: Real

Returns the gradient of the free energy of the model (divided by ell) at magnetizations ms, at parameters (p,q) and at temperature t; see Equation (S8) and Equation (S11).
"""
function gradModelFreeEnergy(ms::Vector{T},p::Int,q::Int,t::Real) where T <: Real
    q2 = length(ms)
    if q2 == modelDimension(q)
        J = modelInteraction(p,q)

        gradF = 2*(J*ms)
        gradF[1] += -t * 2 * gradMagnetizationEntropy(ms[1])
        for i in 2:q2-1
            gradF[i] += -t * 4 * gradMagnetizationEntropy(ms[i])
        end
        gradF[q2] += -t * ( q%2==0 ? 2 : 4 ) * gradMagnetizationEntropy(ms[q2])
        return gradF
    else
        throw(ArgumentError("Length of ms incompatible with q!"))
    end
end
