
# read which fractions have data
function whichQP()
    res = []
    for file in readdir("/run/media/vittorio/Extreme SSD/ConfocalCavity/data/parsed2")[1:end-1] 
        push!(res, Tuple(split(file[2:end-4],"_p") .|> u->parse(Int64,u) )) 
    end
    return res
end

function computeTotal(p,q)
    file = string("/run/media/vittorio/Extreme SSD/ConfocalCavity/data/parsed2/q",q,"_p",p,".csv")
    read(pipeline(`wc -l $file`,`cut -d ' ' -f 1 `), String) |> u -> parse(Int,u)
end

# reads csv file of minima into dataframe
function readMinima(p,q)
    file = string("/run/media/vittorio/Extreme SSD/ConfocalCavity/data/parsed2/q",q,"_p",p,".csv")
    raw = readdlm(file,'|') 
    names = map( Symbol , readdlm("/run/media/vittorio/Extreme SSD/ConfocalCavity/data/qq_pp.csv", '|') ) |> vec
    raw[:,4] = map( Symbol , raw[:,4])
    raw[:,5] = map( u -> parse.(Float64,split(u[2:end-1],", ")) , raw[:,5])
    raw[:,6] = map( u -> u=="Any[]" ? Vector{Float64}() : parse.(Float64,split(u[2:end-1],", ")) , raw[:,6])
    raw[:,7] = map( string , raw[:,7])
    raw[:,9] = map( string |> u->u=="FEASIBLE_POINT" ? "good" : u , raw[:,9])
    tup = Vector{Tuple}()
    for row in eachrow(raw)
        if row[9] == "FEASIBLE_POINT"
            row[9] = "good"
        end
        push!(tup,Tuple(row))
    end
    data = DataFrame(tup)
    names!(data,names)
    data
end

# preprocess dataframe
function preProcess(db)
    # select gradient
    db = db[ db.optimizer .|> u->occursin("gradient",u |> string ) , :]
    # select good simulations
    db = db[db.status .== "good", :]
    # select zero temperature
    db = db[db.temperature .== 0., :]
    # totalTrials
    tot = nrow(db)
    db.totalTrials = ones(tot)*tot
    # add info on ground state configuration (rounded to 1 decimal, I'm interested in 0,+-1 basically) normalized to factor out the Z2 symmetry
    db.cleanMinimum = map( u -> roundVec(u[:minimum],digits=1) |> normalizeZ2 , db |> eachrow)
    # and info on free energy of the exact minimum found, rounded to the second decimal
    db.energy = map( u -> roundVec(freeEnergy(u[:minimum],u[:p],u[:q],u[:temperature]),digits=2) , db |> eachrow)
    # and info on entropy of the clean minimum found, rounded to the second decimal
    db.entropy = map( u -> roundVec(entropy(u[:minimum],u[:p],u[:q]),digits=2) , db |> eachrow)

    return db
end

function loadProcessed(q,p,type)
    if type == :GS
        return load(string("/run/media/vittorio/Extreme SSD/ConfocalCavity/data/processed/GS_q",q,"_p",p,".jld2")) 
    elseif type == :META
        return load(string("/run/media/vittorio/Extreme SSD/ConfocalCavity/data/processed/META_q",q,"_p",p,".jld2")) 
    else
        println("Type ", type, " unknown!")
        return 0
    end
end


# get preprocesssed dataframe and perform statistics on the GS energy
function analyzeGS(db)
    
    # get useful parameters for dataset
    p = db[1,:p]
    q = db[1,:q]
    
    # compute GS quantities
    groundStateEnergy = min(db.energy...)
    groundStates = db[db.energy .== groundStateEnergy,:].cleanMinimum |> unique
    # groundStateDegeneracy = groundStates |> length
    # groundStateEntropies = groundStates .|> s -> entropy(s,p,q)

    save(
         string("/run/media/vittorio/Extreme SSD/ConfocalCavity/data/processed/GS_q",q,"_p",p,".jld2"),
         Dict(
            # "p"          => p,
            # "q"          => q,
            "energy"     => groundStateEnergy,
            # "degeneracy" => groundStateDegeneracy,
            # "entropies"  => groundStateEntropies,
            "states"     => groundStates
           )
        )
end

# TODO: move saturation in another function
function analyzeMetastable(db)
    
    # get useful parameters for dataset
    p = db[1,:p]
    q = db[1,:q]
    
    # add useful columns to db
    db.distance = map( u -> norm(u[1]-u[2]), zip(db.cleanMinimum, db.starting))

    # states 
    tmp = db.cleanMinimum |> countmap
    states = tmp |> keys |> collect

    l = length(states)
    basinNumber = Vector{Int64}(undef,l)
    basinDistances = Vector{Vector{Float64}}(undef,l)
    for (i,s) in enumerate(states)
        basinNumber[i] = tmp[s]
        basinDistances[i] = db.distance[db.cleanMinimum .== [s]]
    end
    
    # basin volume proxy
    # [basinNumber = states .|> s -> [tmp[s],db.distance[db.cleanMinimum .== [s]]]
    # basinNumber = states .|> s -> length( db.cleanMinimum[db.cleanMinimum .== [s]] )
    
    # basin linear distances
    # basinDistance = states .|> s -> db.distance[db.cleanMinimum .== [s]]

    # number of distinct metastable states
    # number = states |> length

    # entropies
    # entropies = states .|> u->entropy(u,p,q) 

    # energies
    # energies = states .|> u->freeEnergy(u,p,q,0)

    # basin volume proxy
    # basinNumber = states .|> s -> tmp[s]
    # basinNumber = states .|> s -> length( db.cleanMinimum[db.cleanMinimum .== [s]] )

    # basin linear distances
    # basinDistance = states .|> s -> db.distance[db.cleanMinimum .== [s]]

    # saturation of metastable states vs total trials
    # vec = db.cleanMinimum
    # distinctMin = Set{Vector{Float64}}()
    # history = Matrix{Float64}(undef,length(vec),2)
    # perm = randperm(length(vec))
    # for i in 1:length(vec)
    #     if !(vec[perm[i]] in distinctMin)
    #         push!(distinctMin, vec[perm[i]])
    #     end
    #     history[i,:] = [ i,length(distinctMin)]
    # end 
    save(
         string("/run/media/vittorio/Extreme SSD/ConfocalCavity/data/processed/META_q",q,"_p",p,".jld2"),
         Dict(
            # number        => number,
            "states"        => states,
            # entropies     => entropies,
            # energies      => energies,
            # history       => history,
            "basinDistances" => basinDistances,
            "basinNumber"   => basinNumber
           )
        )
end

# history
function history(db)
    # saturation of metastable states vs total trials
    vec = db.cleanMinimum
    distinctMin = Set{Vector{Float64}}()
    history = Matrix{Float64}(undef,length(vec),2)
    perm = randperm(length(vec))
    for i in 1:length(vec)
        if !(vec[perm[i]] in distinctMin)
            push!(distinctMin, vec[perm[i]])
        end
        history[i,:] = [ i,length(distinctMin)]
    end 
    history
end

