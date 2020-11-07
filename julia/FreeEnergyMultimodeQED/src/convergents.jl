################################################# CONVERGENTS
function continuedFraction(r::Float64,n)
    if n == 0
        return Vector{Int64}()
    end
    i = floor(r) |> Int64
    f = r - i
    if f == 0
        return r |> Int64
    else
        return vcat(i,continuedFraction(1/f,n-1))
    end
end

function continuedFraction(r::Rational{Int64})
    i = floor(r) |> Int64
    f = r - i
    if f == 0
        return r |> Int64
    else
        return vcat(i,continuedFraction(1//f))
    end
end

function convergents(cf::Vector{Int64})
    [ foldl( (u,v) -> 1//u + Rational(v) , cf[1:i] |> reverse )
     for i in 1:length(cf) ]
end

function convergents(r::Float64; n=20)
    convergents(continuedFraction(r,n))
end

function convergents(r::Rational{Int64})
    convergents(continuedFraction(r))
end


