using CSV, DataFrames, Interpolations

function load_UUP07(par::Dict{String,Any})
    """
    Return interp function based on U-UP07 velocity model (absolute).
        https://www.atlas-of-the-underworld.org/  

    Parameters
    ----------
    - `par` : Dict{String,Any}

    Returns
    -------
    - `itp` : interpolate((lon ::Vector{Float64},lat ::Vector{Float64},z ::Vector{Float64}), sn ::Array{Float64, 3})
        Function to interpolate slowness (s/m) based on a given (lon, lat, z).
    """
    vel = readlines(par["DataDir"] * par["vel_UUP07"])
    bld, nlon, nlat, nz = parse(Float64, split(vel[1])[1]), parse(Int, split(vel[1])[2]),
                            parse(Int, split(vel[1])[3]), parse(Int, split(vel[1])[4])

    lons    = [parse(Float64,i) for i in split(vel[2])]
    lats    = [parse(Float64,i) for i in split(vel[3])]
    z       = [parse(Float64,i) for i in split(vel[4])]
    
    lons[lons.<0] .+= 360
    
    sn, vp = Array{Float64,3}(undef,nlon,nlat,nz), Array{Float64,3}(undef,nlon,nlat,nz)
    
    for iz = 1:nz
        for ilat = 1:nlat
            for ilon = 1:nlon
                line_index = (iz-1)*nlat + (ilat-1) + 5
                vp[ilon, ilat, iz] = parse(Float64, split(vel[line_index])[ilon])
                sn[ilon, ilat, iz] = 1/vp[ilon, ilat, iz]
            end
        end
    end

    if par["coordinates"] == 1
        # itp = interpolate((round.(lon;digits=2)[:,1],round.(lat;digits=2)[1,:],z),sn[1,:,:,:],BSpline(Linear()), OnGrid()) # for uniformly gridded data
        itp = interpolate((lons, lats, z), sn[:,:,:],Gridded(Linear()))
        #08/23/23 yurong: dangerous! extrapolation is included!
        itp = extrapolate(itp, Flat())

        lon_bounds = (minimum(lons), maximum(lons))
        lat_bounds = (minimum(lats), maximum(lats))
        z_bounds = (first(z), last(z))

        println("Lon bounds: ", lon_bounds)
        println("Lat bounds: ", lat_bounds)
        println("Z bounds: ", z_bounds)
        return itp

    elseif par["coordinates"] == 2
        par["lat0"] = -23.1000
        dataLat = repeat(lats', length(lons), 1)
        dataLon = repeat(lons, 1, length(lats))
        (x, y) = lonlat2xy(par["lon0"], par["lat0"], par["beta"], dataLon, dataLat)

        itp = interpolate((round.(x;digits=2)[:,1],round.(y;digits=2)[1,:],z),sn[:,:,:],Gridded(Linear()))
        #08/23/23 yurong: dangerous! extrapolation is included!
        itp = extrapolate(itp, Flat())

        x_bounds = (minimum(x), maximum(x))
        y_bounds = (minimum(y), maximum(y))
        z_bounds = (first(z), last(z))
        println("X bounds: ", x_bounds)
        println("Y bounds: ", y_bounds)
        println("Z bounds: ", z_bounds)
        return itp
    else
        error("\"coordinates\" options are 1(Spherical coordinate system) and 2(Cartesian coordinate system).")
    end
end

function load_lauvel(par::Dict{String,Any})
    """
    Return interp function based on lau.vel velocity model.
        
    Parameters
    ----------
    - `par` : Dict{String,Any}

    Returns
    -------
    - `itp` : interpolate((lon ::Vector{Float64},lat ::Vector{Float64},z ::Vector{Float64}), sn ::Array{Float64, 3})
        Function to interpolate slowness (s/m) based on a given (lon, lat, z).
    """
    vel = readlines(par["DataDir"] * par["vel_lau"])
    nnx, nny, nnz = parse(Int,split(vel[1])[1]), parse(Int,split(vel[1])[2]), parse(Int,split(vel[1])[3])
    lat0, lon0, beta = parse(Float64,split(vel[2])[1]), parse(Float64,split(vel[2])[2]), parse(Float64,split(vel[2])[3])
    lat, lon = Array{Float64,2}(undef,nnx,nny), Array{Float64,2}(undef,nnx,nny)
    for i = 1:nnx
        for j = 1:nny
            lat[i,j] = parse(Float64,split(vel[(i-1)*nny+j+2])[1])
            lon[i,j] = parse(Float64,split(vel[(i-1)*nny+j+2])[2]) 
        end
    end

    z = [parse(Float64,i) for i in split(vel[nnx*nny+3])]
    sn, vps = Array{Float64,4}(undef,2,nnx,nny,nnz), Array{Float64,4}(undef,2,nnx,nny,nnz)
    for p = 1:2
        for i = 1:nnx
            for j = 1:nny
                for k = 1:nnz
                    vps[p,i,j,k] = parse(Float64,split(vel[(i-1+p*nnx)*nny+j+3])[k])
                    sn[p,i,j,k] = 1/vps[p,i,j,k]
                end
            end
        end
    end

    if par["coordinates"] == 1
        # itp = interpolate((round.(lon;digits=2)[:,1],round.(lat;digits=2)[1,:],z),sn[1,:,:,:],BSpline(Linear()), OnGrid()) # for uniformly gridded data
        itp = interpolate((round.(lon;digits=2)[:,1],round.(lat;digits=2)[1,:],z),sn[1,:,:,:],Gridded(Linear()))
        #08/23/23 yurong: dangerous! extrapolation is included!
        itp = extrapolate(itp, Flat())

        lon_bounds = (minimum(lon), maximum(lon))
        lat_bounds = (minimum(lat), maximum(lat))
        z_bounds = (first(z), last(z))

        println("Lon bounds: ", lon_bounds)
        println("Lat bounds: ", lat_bounds)
        println("Z bounds: ", z_bounds)
        return itp

    elseif par["coordinates"] == 2
        (dataX, dataY) = lonlat2xy(lon0, lat0, beta, lon, lat)

        itp = interpolate((round.(dataX;digits=2)[:,1],round.(dataY;digits=2)[1,:],z),sn[1,:,:,:],Gridded(Linear()))
        #08/23/23 yurong: dangerous! extrapolation is included!
        itp = extrapolate(itp, Flat())

        x_bounds = (first(round.(dataX;digits=2)[:,1]), last(round.(dataX;digits=2)[:,1]))
        y_bounds = (first(round.(dataY;digits=2)[1,:]), last(round.(dataY;digits=2)[1,:]))
        z_bounds = (first(z), last(z))

        println("X bounds: ", x_bounds)
        println("Y bounds: ", y_bounds)
        println("Z bounds: ", z_bounds)
        return itp
    else
        error("\"coordinates\" options are 1(Spherical coordinate system) and 2(Cartesian coordinate system).")
    end
end

function load_par_from_yml(file_name::String)
    """
    Return a Dict{String,Any} from the YAML file.
    """
    par = YAML.load_file(file_name;dicttype=Dict{String,Any})
    par["base_dir"] = "./" * file_name[1:end-4] * "/"
    par["min_depth"] = par["add_litho"] ? par["litho_thickness"] : par["min_depth"]
    return par
end

function load_data_Tonga(par::Dict{String,Any})
    """
    Return dataStruct and RayTraces.

    Parameters
    ----------
    par : Dict{String,Any}

    Returns
    -------
    dataStruct : Dict{String,Any}
        "tS" : Array{Float64,1}
            tStar for each event-station pair.
        "allaveatten" : Array{Float64,1}
            Path-average attenuation for each event-station pair.
        "allLats" : Array{Float64,1}
            Latitude for the stations in each event-station pair.
        "allLons" : Array{Float64,1}
            Longitude for the stations in each event-station pair.
        "allEles" : Array{Float64,1}
            Elevation for the stations in each event-station pair.
        "allSig" : Array{Float64,1}
            Error for each t*.
        "dataX" : Array{Float64,1}
            x (cartesion coordinate)/ lon (spherical coordinate) for stations.
        "dataY" : Array{Float64,1}
            y (cartesion coordinate)/ lat (spherical coordinate) for stations.
        "xVec" : StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
            x (cartesion coordinate)/ lon (spherical coordinate) for regular grids for plotting.
        "yVec" : StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
            y (cartesion coordinate)/ lat (spherical coordinate) for regular grids for plotting.
        "zVec" : StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
            z (depth) for regular grids for plotting.
        "elonsX" : Array{Float64,2}
            x (cartesion coordinate) for events.
        "elatsY" : Array{Float64,2}
            y (cartesion coordinate) for events.
        "elons" : Array{Float64,2}
            Longitude for events.
        "elats" : Array{Float64,2}
            Latitude for events.
        "edep" : Array{Float64,2}
            Depth for events.

    RayTraces : Dict{String,Any}
        "rayX" : Array{Float64,2}
            x (cartesian coordinates)/ lon (spherical coordinates) for each ray data point.
        "rayY" : Array{Float64,2}
            y (cartesian coordinates)/ lat (spherical coordinates) for each ray data point.
        "rayZ" : Array{Float64,2}
            depth for each ray data point.
        "rayL" : Array{Float64,2}
            raylength for each ray segment.
        "rayU"::Array{Float64,2}
            slowness for each ray segment.
        "U" : Array{Float64,2}
            slowness for each ray data point.   
    """
    allLats, allLons    = [], [] #station latitudes and longitudes
    allTS, allSig, allaveatten  = [], [], []
    elats, elons, edep          = [], [], [] #event latitudes and longitudes
    
    traces      = load(par["base_dir"] * "data/traces.jld")
    elats       = Array{Float64}(traces["EventLatitude"])
    elons       = Array{Float64}(traces["EventLongitude"])
    edep        = Array{Float64}(traces["EventDepth"])
    allLats     = Array{Float64}(traces["latitude"])
    allLons     = Array{Float64}(traces["longitude"])
    allEles     = Array{Float64}(traces["elevation"])
    allTS       = Array{Float64}(traces["tStar"])
    allaveatten = Array{Float64}(traces["aveatten"])
    allSig      = Array{Float64}(traces["error"])

    #drop to singleton dimensions
    allTS       = Vector{Float64}(vec(allTS))
    allLats     = Vector{Float64}(vec(allLats))
    allLons     = Vector{Float64}(vec(allLons))
    allEles     = Vector{Float64}(vec(allEles))
    allaveatten = Vector{Float64}(vec(allaveatten))
    allSig      = Vector{Float64}(vec(allSig))
    
    if par["coordinates"] == 1
        #stations and events coordicates in spherical coordinate system
        dataX, dataY    = allLons, allLats
        elonsX, elatsY  = elons, elats

        #study area
        minX = minimum(dataX) - par["buffer_sph"]
        maxX = maximum(dataX) + par["buffer_sph"]
        minY = minimum(dataY) - par["buffer_sph"]
        maxY = maximum(dataY) + par["buffer_sph"]

        xVec = minX:par["latlonnodeSpacing"]:maxX
        yVec = minY:par["latlonnodeSpacing"]:maxY
        zVec = par["min_depth"]:par["ZnodeSpacing"]:par["max_depth"]

    elseif par["coordinates"] == 2
        #stations and events coordicates in cartesian system
        lat0 = par["lat0"]
        lon0 = par["lon0"]
        beta = par["beta"]
        (dataX, dataY)      = lonlat2xy(lon0, lat0, beta, allLons, allLats)
        (elonsX, elatsY)    = lonlat2xy(lon0, lat0, beta, elons, elats)


        #study area
        minX = minimum(dataX) - par["buffer_cart"]
        maxX = maximum(dataX) + par["buffer_cart"]
        minY = minimum(dataY) - par["buffer_cart"]
        maxY = maximum(dataY) + par["buffer_cart"]

        xVec = minX:par["XYnodeSpacing"]:maxX
        yVec = minY:par["XYnodeSpacing"]:maxY
        zVec = par["min_depth"]:par["ZnodeSpacing"]:par["max_depth"]

    end


    #load raypaths
    raypath=load(par["base_dir"] * "data/raypaths.jld")
    x = raypath["x"]
    y = raypath["y"]
    z = raypath["z"]
    U = raypath["u"]

    #raylength and slowness for each segment
    if par["coordinates"] == 1
        rayl = sph_dist.(x[1:end - 1,:], y[1:end - 1,:], z[1:end - 1,:],
                x[2:end,:], y[2:end,:], z[2:end,:])
        
    elseif par["coordinates"] == 2
        rayl = cart_dist.(x[1:end - 1,:], y[1:end - 1,:], z[1:end - 1,:],
                x[2:end,:], y[2:end,:], z[2:end,:])

    end
    rayu = 0.5 .* (U[1:end - 1,:] + U[2:end,:])

    dataStruct = Dict(
        "tS" => allTS   ::Array{Float64,1},
        "allaveatten" => allaveatten   ::Array{Float64,1},
        "allLats" => allLats        ::Array{Float64,1},    #Latitude for the stations in each event-station pair
        "allLons" => allLons        ::Array{Float64,1},    #Longitude for the stations in each event-station pair
        "allEles" => allEles        ::Array{Float64,1},    #Elevation for the stations in each event-station pair
        "allSig" => allSig      ::Array{Float64,1},      #ATTENTION!!! remain unknown!
        "dataX" => dataX        ::Array{Float64,1},        #Station position in each event-station pair
        "dataY" => dataY        ::Array{Float64,1},
        "xVec" => xVec      ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
        "yVec" => yVec      ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
        "zVec" => zVec      ::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
        "elonsX" => elonsX  ::Array{Float64,2},
        "elatsY" => elatsY  ::Array{Float64,2},
        "elons" => elons    ::Array{Float64,2},
        "elats" => elats    ::Array{Float64,2},
        "edep" => edep      ::Array{Float64,2},
    )

    RayTraces = Dict(
        "rayX" => x         ::Array{Float64,2},
        "rayY" => y         ::Array{Float64,2},
        "rayZ" => z         ::Array{Float64,2},
        "rayL" => rayl      ::Array{Float64,2},
        "rayU" => rayu      ::Array{Float64,2},
        "U" => U            ::Array{Float64,2}
    )

    return dataStruct, RayTraces
end

function load_PREM(par::Dict{String,Any})
    """
    Return Qp interp function based on PREM model.
        PREM header: radius,depth,density,Vpv,Vph,Vsv,Vsh,eta,Qµ,Qk
        http://ds.iris.edu/ds/products/emc-prem/

        Assuming isotropic velocity model:
        Vs = (Vsh+Vsv)/2, Vp = (Vph+Vpv)/2.

        1/Qp = (1-L)/Qk + L/Qµ
        1/Qs = 1/Qµ
        L = 4Vs^2/3Vp^2
        (Wei and Wiens, 2020)

    Parameters
    ----------
    par : Dict{String,Any}

    Returns
    -------
    Find_Qp : interpolate(z ::Vector{Float64})
        Function to interpolate Qp based on a given depth.
    """

    prem = CSV.read(par["DataDir"] * par["PREM"], DataFrame)
    z = prem[1:end, 2]
    Vp = (prem[1:end, 4] + prem[1:end, 5]) ./ 2
    Vs = (prem[1:end, 6] + prem[1:end, 7]) ./ 2
    Qk = prem[1:end, 10]
    Qµ = prem[1:end, 9]
    L = 4 .* Vs.^2 ./ (3 .* Vp.^2)
    Qp = 1 ./ ((1 .- L) ./ Qk .+ L ./ Qµ)
    Find_Qp = interpolate((z,), Qp, Gridded(Linear()))

    return Find_Qp
end

