# load information from raypaths.p and p_tstar.dat, including stations.lst
# output raypaths.jld and traces.jld for MCMC input
using DelimitedFiles,Interpolations,JLD, YAML

include("./scripts/model.jl")
include("./scripts/stats.jl")
include("./scripts/load.jl")
include("./scripts/utils.jl")


file_name = "par_sph_test.yml"
@time par = load_par_from_yml(file_name)
@time itp = load_UUP07(par)
# @time itp = load_lauvel(par)
@time Find_Qp = load_PREM(par)
make_dir(par, file_name)

function intersect_trace_raypath()
    """
    Load raypaths.dat and p_tstar.dat, and store the data with intersected 
    event-station pairs as the new matched raypaths.dat/p_tstar.dat.
    """
    loadatten = readlines(par["DataDir"] * par["tstar_file"])
    oridset1 = Set{String}()
    sta1 = Dict{String, Vector{String}}()

    for line in loadatten
        tokens = split(line)
        push!(oridset1, tokens[1])
        haskey(sta1, tokens[1]) ? push!(sta1[tokens[1]],tokens[2]) : sta1[tokens[1]] = [tokens[2]]
    end
    # oridlst = collect(oridset1)

    
    loadrayp = readlines(par["DataDir"] * par["rayp_file"])
    oridset2 = Set{String}()
    sta2 = Dict{String, Vector{String}}()

    for line in loadrayp
        tokens = split(line)
        if length(tokens) > 4
            push!(oridset2, tokens[1])
            haskey(sta2, tokens[1]) ? push!(sta2[tokens[1]],tokens[2]) : sta2[tokens[1]] = [tokens[2]]
        end
    end

    # Find the intersection of the two sets of origin IDs.
    oridset = intersect(oridset1, oridset2)

    # For each origin ID in the intersected set,
    # find and store the intersection of the station vectors.
    sta = Dict{String, Vector{String}}()
    for orid in oridset
        if haskey(sta1, orid) && haskey(sta2, orid)
            sta[orid] = intersect(sta1[orid], sta2[orid])
        end
    end

    #Rewrite raypaths_matched.dat/p_tstar_matched.dat based on the intersected dataset
    f1 = open(par["base_dir"] * "data/p_tstar.dat", "w")
    f2 = open(par["base_dir"] * "data/raypaths.dat", "w")

    # Filter and write matched records from loadatten to p_tstar_matched.dat
    for line in loadatten
        tokens = split(line)
        if tokens[1] in oridset && tokens[2] in sta[tokens[1]]
            println(f1, line)
        end
    end

    # Filter and write matched records from loadrayp to raypaths_matched.dat
    flag = false  # This is used to ensure subsequent lines of a matched record are also written
    i = 1
    while i <= length(loadrayp)
        line = loadrayp[i]
        tokens = split(line)
        if length(tokens) > 4  # This appears to identify header lines in the dataset
            if tokens[1] in oridset && tokens[2] in sta[tokens[1]]
                flag = true  # Match found, so set flag to true
                println(f2, line)
            else
                flag = false  # No match, so set flag to false
            end
        else
            # This writes non-header lines only if the previous header was a match
            if flag
                println(f2, line)
            end
        end 
        i += 1
    end          

    # Close the output files
    close(f1)
    close(f2)



end

function load_sph_raypath()
    """
    Load 3D raytracing results in spherical coordinate system.
    The current raypaths are from Fan Wang, utilizing 3D fast marching method.
    """
    loadraypath = readlines(par["base_dir"] * "data/raypaths.dat")

    X, Y, Z, U = [], [], [], []
    ix, iy, iz, iu = [], [], [], []
    n=0

    ENV["GKSwstype"] = "nul"

    # p = plot(legend=false)
    for i in 1:length(loadraypath)
        points = split(loadraypath[i])
        if length(points) > 4
            n+=1
            if ix == []
                continue
            end

            # Plot the difference between velocity obtained from
            # 1) trilinear interpolation 2) fm3d ray tracing
            # iu1 = itp.(ix,iy,iz)
            # if n < 10000
            #     plot!(p, iz,1 ./ iu1 .- 1 ./ iu)
            # elseif n == 10000
            #     hline!(p, [0])
            #     title!(p,"vel differece (interp-fm3d)")
            #     xlabel!(p,"depth")
            #     ylabel!(p,"vel diff (km/s)")
            #     savefig(p,"vel diff")
            # end

            push!(X,ix); push!(Y,iy); push!(Z,iz); push!(U,iu)

            ix, iy, iz, iu = [], [], [], []
        else
            iix = parse(Float64,points[3])
            iiy = parse(Float64,points[2])
            iiz = parse(Float64,points[1])
            if length(points) > 3
                iiu = 1/parse(Float64,points[4])
            else
                iiu = itp(iix,iiy,iiz)
            end
            append!(ix,iix); append!(iy,iiy); append!(iz,iiz); append!(iu, iiu)
        end
    end
    if !isempty(ix)
        push!(X,ix); push!(Y,iy); push!(Z,iz); push!(U,iu)
    end
        
    #interpolate the ray if individual segment is shorter than seg_len [km]
    # 09/29/23 yurong: potential modifications
    # if Spherical: calc the seg lengths
    if par["ray_interpolate"]
        seg_len = par["seg_len"]
        for i in 1:length(U)
            for j in 2:length(U[i])
                ilength = sph_dist(X[i][j], Y[i][j], Z[i][j], X[i][j-1], Y[i][j-1], Z[i][j-1])
                if ilength > seg_len && Z[i][j] < par["interp_depth"]
                    #note: use round() instead of ceil() since there are many 10.5 km segments when 10 km is set to the threshold
                    seg_num     = Int(ceil(ilength/seg_len))
                    insert_x    = range(X[i][j-1],stop=X[i][j],length=seg_num+1)
                    insert_y    = range(Y[i][j-1],stop=Y[i][j],length=seg_num+1)
                    insert_z    = range(Z[i][j-1],stop=Z[i][j],length=seg_num+1)
                    insert_u    = itp.(insert_x,insert_y,insert_z)
                    for k in 2:seg_num
                        insert!(X[i],j+k-2,insert_x[k])
                        insert!(Y[i],j+k-2,insert_y[k])
                        insert!(Z[i],j+k-2,insert_z[k])
                        insert!(U[i],j+k-2,insert_u[k])
                    end
                end
            end
        end
    end

    #rearrange the raypaths and according slowness to be a matrix
    #fill the empty cell with NaN
    x = fill(NaN,(maximum(length.(U)),length(U)))
    y = fill(NaN,(maximum(length.(U)),length(U)))
    z = fill(NaN,(maximum(length.(U)),length(U)))
    u = fill(NaN,(maximum(length.(U)),length(U)))

    for i in 1:length(U)
        x[1:length(U[i]),i]=X[i]
        y[1:length(U[i]),i]=Y[i]
        z[1:length(U[i]),i]=Z[i]
        u[1:length(U[i]),i]=U[i]
    end

    save(par["base_dir"] * "data/raypaths.jld","x",x,"y",y,"z",z,"u",u)  
 
    #release memory
    X, Y, Z, U, loadraypath = nothing, nothing, nothing, nothing, nothing
    x, y, z, u = nothing, nothing, nothing, nothing
    lon, lat = nothing,nothing

end

function load_cart_raypath()
    """
    Load 3D raytracing results in cartesian coordinate system.
    """
    loadraypath = readlines(par["DataDir"] * par["rayp_file"])

    #load raypaths from raypaths
    #interpolate 3D slowness for each point in a ray
    X, Y, Z, U = [], [], [], []
    ix, iy, iz, iu = [], [], [], []
    n=0
    for i in 1:length(loadraypath)
        sign = split(loadraypath[i])[1]
        if sign == "1234567"
            n+=1
            if ix == []
                continue
            end
            iu = itp.(ix,iy,iz)
            push!(X,ix); push!(Y,iy); push!(Z,iz); push!(U,iu)
            ix, iy, iz, iu = [], [], [], []
        else
            iix = parse(Float64,split(loadraypath[i])[1])
            iiy = parse(Float64,split(loadraypath[i])[2])
            iiz = parse(Float64,split(loadraypath[i])[3])
            append!(ix,iix); append!(iy,iiy); append!(iz,iiz)
        end
        
        if i == length(loadraypath)
            iu = itp.(ix,iy,iz)
            push!(X,ix); push!(Y,iy); push!(Z,iz); push!(U,iu)
        end
    end

    #interpolate the ray if individual segment is shorter than seg_len [km]
    # 09/29/23 yurong: potential modifications
    # if Spherical: calc the seg lengths
    if par["ray_interpolate"]
        seg_len = par["seg_len"]
        for i in 1:length(U)
            for j in 2:length(U[i])
                ilength = cart_dist(X[i][j], Y[i][j], Z[i][j], X[i][j-1], Y[i][j-1], Z[i][j-1])
                if ilength > seg_len && Z[i][j] < par["interp_depth"]
                    #note: use round() instead of ceil() since there are many 10.5 km segments when 10 km is set to the threshold
                    seg_num     = Int(ceil(ilength/seg_len))
                    insert_x    = range(X[i][j-1],stop=X[i][j],length=seg_num+1)
                    insert_y    = range(Y[i][j-1],stop=Y[i][j],length=seg_num+1)
                    insert_z    = range(Z[i][j-1],stop=Z[i][j],length=seg_num+1)
                    insert_u    = itp.(insert_x,insert_y,insert_z)
                    for k in 2:seg_num
                        insert!(X[i],j+k-2,insert_x[k])
                        insert!(Y[i],j+k-2,insert_y[k])
                        insert!(Z[i],j+k-2,insert_z[k])
                        insert!(U[i],j+k-2,insert_u[k])
                    end
                end
            end
        end
    end

    #rearrange the raypaths and according slowness to be a matrix
    #fill the empty cell with NaN
    x = fill(NaN,(maximum(length.(U)),length(U)))
    y = fill(NaN,(maximum(length.(U)),length(U)))
    z = fill(NaN,(maximum(length.(U)),length(U)))
    u = fill(NaN,(maximum(length.(U)),length(U)))

    for i in 1:length(U)
        x[1:length(U[i]),i]=X[i]
        y[1:length(U[i]),i]=Y[i]
        z[1:length(U[i]),i]=Z[i]
        u[1:length(U[i]),i]=U[i]
    end

    save(par["base_dir"] * "data/raypaths.jld","x",x,"y",y,"z",z,"u",u) 

        
    #release memory
    X, Y, Z, U, loadraypath = nothing, nothing, nothing, nothing, nothing
    x, y, z, u = nothing, nothing, nothing, nothing
    lon, lat = nothing,nothing
end

function load_sph_traceinfo()
    """
    Load 3D raytracing results in spherical coordinate system.
    """
    loadatten = readlines(par["base_dir"] * "data/p_tstar.dat")
    loadsta = readlines(par["DataDir"] * par["sta_info"])

    #load station latitude and longitude
    stalat, stalon, staele = Dict(), Dict(), Dict()
    for line in loadsta
        stalat[split(line)[1]] = parse(Float64,split(line)[2])
        stalon[split(line)[1]] = parse(Float64,split(line)[3])
        staele[split(line)[1]] = parse(Float64,split(line)[4])
    end

    #load all the information from p_tstar.dat (from t* inversion or synthetic dataset)
    stations, EventLatitudes, EventLongitudes, EventDepths = String[], Float64[], Float64[], Float64[]
    latitudes, longitudes, elevations, tStars, errors, aveattens = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    # Determine the systematic error addition 
    systematic_error = par["add_systematic_tsterr"] ? par["systematic_tsterr"] : 0.0

    for line in loadatten
        tokens = split(line)
        # orid
        push!(stations, tokens[2])
        push!(EventLatitudes, parse(Float64, tokens[3]))
        push!(EventLongitudes, mod(parse(Float64, tokens[4]), 360))
        push!(EventDepths, parse(Float64, tokens[5]))
        push!(latitudes, stalat[tokens[2]])
        push!(longitudes, mod(stalon[tokens[2]], 360))
        push!(elevations, staele[tokens[2]])
        push!(tStars, parse(Float64, tokens[6]))
        push!(errors, parse(Float64, tokens[7]) + systematic_error)
        push!(aveattens, parse(Float64, tokens[9]))
    end
    
    save(par["base_dir"] * "data/traces.jld","station",reshape(stations,(length(stations),1)),
    "EventLatitude",reshape(EventLatitudes,(length(EventLatitudes),1)),
    "EventLongitude",reshape(EventLongitudes,(length(EventLongitudes),1)),
    "EventDepth",reshape(EventDepths,(length(EventDepths),1)),
    "latitude",reshape(latitudes,(length(latitudes),1)),
    "longitude",reshape(longitudes,(length(longitudes),1)),
    "elevation",reshape(elevations,(length(elevations),1)),
    "tStar",reshape(tStars,(length(tStars),1)),
    "error",reshape(errors,(length(errors),1)),
    "aveatten",reshape(aveattens,(length(aveattens),1)))
    
end

function load_cart_traceinfo()
    """
    Load 3D raytracing results in cartesian coordinate system.
    """
    loadatten = readlines(par["DataDir"] * par["tstar_file"])
    loadsta = readlines(par["DataDir"] * par["sta_info"])

    #load station latitude and longitude
    stalat, stalon, staele = Dict(), Dict(), Dict()
    for line in loadsta
        stalat[split(line)[1]] = parse(Float64,split(line)[2])
        stalon[split(line)[1]] = parse(Float64,split(line)[3])
        staele[split(line)[1]] = parse(Float64,split(line)[4])
    end

    #load all the information from p_tstar.dat (from t* inversion or synthetic dataset)
    stations, EventLatitudes, EventLongitudes, EventDepths = String[], Float64[], Float64[], Float64[]
    latitudes, longitudes, elevations, tStars, errors, aveattens = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    # Determine the systematic error addition 
    systematic_error = par["add_systematic_tsterr"] ? par["systematic_tsterr"] : 0.0

    for line in loadatten
        tokens = split(line)
        push!(stations, tokens[1])
        push!(EventLatitudes, parse(Float64, tokens[2]))
        push!(EventLongitudes, mod(parse(Float64, tokens[3]), 360))
        push!(EventDepths, parse(Float64, tokens[4]))
        push!(latitudes, stalat[tokens[1]])
        push!(longitudes, mod(stalon[tokens[1]], 360))
        push!(elevations, staele[tokens[1]])
        push!(tStars, parse(Float64, tokens[5]))
        push!(errors, parse(Float64, tokens[6]) + systematic_error)
        push!(aveattens, parse(Float64, tokens[8]))
    end
    
    save(par["base_dir"] * "data/traces.jld","station",reshape(stations,(length(stations),1)),
    "EventLatitude",reshape(EventLatitudes,(length(EventLatitudes),1)),
    "EventLongitude",reshape(EventLongitudes,(length(EventLongitudes),1)),
    "EventDepth",reshape(EventDepths,(length(EventDepths),1)),
    "latitude",reshape(latitudes,(length(latitudes),1)),
    "longitude",reshape(longitudes,(length(longitudes),1)),
    "elevation",reshape(elevations,(length(elevations),1)),
    "tStar",reshape(tStars,(length(tStars),1)),
    "error",reshape(errors,(length(errors),1)),
    "aveatten",reshape(aveattens,(length(aveattens),1)))
    
end

if par["coordinates"] == 1
    @time intersect_trace_raypath()
    @time load_sph_raypath()
    @time load_sph_traceinfo()
elseif par["coordinates"] == 2
    @time load_cart_raypath()
    @time load_cart_traceinfo()
end

