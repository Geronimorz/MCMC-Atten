using Base
function lonlat2xy( 
    lon0::Float64, 
    lat0::Float64, 
    beta::Float64, 
    lon1, 
    lat1)
    # function lonlat2xy(lon0::Float64, lat0::Float64, beta::Float64, lon1::Array{Float64,2}, lat1::Array{Float64,2})
    # Convert lon, lat to x, y
    # lon0, lat0: reference location
    # beta: angle of x axis from east (counterclockwise)
    re = 6371
    # pi = 4.*atan(1.)
    r2d = 180.0 / pi
    # yurong 05/12/23 high latitude correction
    # xx = (lon1 .- lon0) .* re ./ r2d
    xx = cos.(lat1 ./ r2d) .* (lon1 .- lon0) .* re ./ r2d
    yy = (lat1 .- lat0) .* re ./ r2d
    x1 = (xx .- yy .* tan.(beta)) .* cos.(beta)
    y1 = x1 .* tan.(beta) .+ yy ./ cos.(beta)

    ###################CANNOT WORK!!!! HOW TO DEAL WITH IT?!!##############
    # if abs.(x1) < 0.01
    #     x1 = 0
    # end
    # if abs.(y1) < 0.01
    #     y1 = 0
    # end
    return x1, y1
end

function xy2lonlat(
    lon0::Float64, 
    lat0::Float64, 
    beta::Float64, 
    x2, 
    y2
    )
    # Convert x, y to lon, lat 
    # lon0, lat0: reference location
    # beta: angle of x axis from east (counterclockwise)
    re = 6371
    # pi = 4.*atan(1.)
    r2d = 180.0 / pi
    yy = (y2 .- x2 .* tan.(beta)) .* cos.(beta)
    xx = yy .* tan.(beta) .+ x2 ./ cos.(beta)
    # yurong 05/12/23 high latitude correction
    # lon2 = xx .* r2d ./ re .+ lon0
    lat2 = yy .* r2d ./ re .+ lat0
    lon2 = xx .* r2d ./ re ./ cos.(lat2 ./ r2d) .+ lon0
    
    return lon2, lat2
end

function sph_dist(
    lon1, lat1, z1, lon2, lat2, z2, 
    Re::Float64=6371.0
    )
    """
    Calculate the distance between (lon1, lat1, z1) and (lon2, lat2, z2)

    Parameters
    ----------
    Float64/list

    Returns
    -------
    Distance : Float64/list
        Distance between the 2 points [km]
    """
    #=
    Equation is from
    https://math.stackexchange.com/questions/833002/distance-between-two-points-in-spherical-coordinates
    =#
    
    # Convert degrees to radians
    θ1 = lat1 .+ 90
    θ2 = lat2 .+ 90
    ϕ1, θ1, ϕ2, θ2 = map(x -> deg2rad.(x), (lon1, θ1, lon2, θ2))

    r1 = Re .- z1
    r2 = Re .- z2
    return sqrt.(r1.^2 .+ r2.^2 .- 2 .* r1 .* r2 .* (sin.(θ1).*sin.(θ2).*cos.(ϕ2 - ϕ1) + cos.(θ1).*cos.(θ2)))
end

function cart_dist(
    x1, y1, z1, x2, y2, z2, 
    )
    """
    Calculate the distance between (x1, y1, z1) and (x2, y2, z2)

    Parameters
    ----------
    Float64/list

    Returns
    -------
    Distance : Float64/list
        Distance between the 2 points [km]
    """

    return sqrt.((x1 - x2).^2 + (y1 - y2).^2 + (z1 - z2).^2)
end

function deg2rad(deg)
    return deg * π / 180
end

function hypotenuse(a, b)
    return sqrt(a^2 + b^2)
end

function normalize(arr::Array, min_val::Number, max_val::Number)
    min_arr, max_arr = extrema(arr)
    return min_val .+ (max_val - min_val) .* ((max_arr .- arr) ./ (max_arr - min_arr))
end

function generate_new_zeta(zeta::Float64, sig_zeta::Float64)
    """
    Generate a new zeta value based on Gaussian distribution with the given standard deviation
    """
    return rand(Normal(zeta, sig_zeta), 1)[1]
end

function calculate_alpha_for_birth(
    par::Dict{String,Any},
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    model::Model,
    modeln::Model, 
    czeta::Float64, 
    zetanew::Float64, 
    sig_zeta::Float64)
    """
    Calculate the acceptance function (α) for the birth action.

    Parameters
    ----------
    - `par`
    - `dataStruct`
    - `RayTraces`
    - `model`:      The current state of the model before the birth action.
    - `modeln`:     The proposed state of the model after the birth action.
    - `czeta`:      The current value of zeta at the new cell nuclei.
    - `zetanew`:    The proposed new value of zeta.
    - `sig_zeta`:   The standard deviation for zeta.

    Returns
    -------
    - `α`:          The calculated acceptance ratio for the birth action.
    - `valid`:      Whether the proposed state is valid.

    """
    α, valid = 0, 1
    if par["prior"] == 1 # Uniform
        if zetanew > 0 && zetanew < par["zeta_scale"]
            # α(birth), Byrnes and Bezada, 2020, eq. 16
            α = ((model.nCells) / (model.nCells + 1)) * ((sig_zeta * sqrt(2 * pi)) / (par["zeta_scale"])) * 
            exp(((czeta - zetanew)^2) / (2 * sig_zeta^2) - (modeln.phi - model.phi) / 2)
            α = min([1 α]...)
        else
            valid = 0
        end
    elseif par["prior"] == 2 # Normal
        α = ((model.nCells) / (model.nCells + 1)) * (sig_zeta / par["zeta_scale"]) * 
        exp(-(zetanew)^2 / (par["zeta_scale"]^2) + (czeta - zetanew)^2 / (2 * sig_zeta^2) - (modeln.phi - model.phi) / 2)
        α = min([1 α]...)
    elseif par["prior"] == 3 # Exponential
        if zetanew > 0
            α = ((model.nCells) / (model.nCells + 1)) * (sqrt(2 * pi) * sig_zeta / par["zeta_scale"]) * 
            exp(-zetanew / par["zeta_scale"] + (czeta - zetanew)^2 / (2 * sig_zeta^2) - (modeln.phi - model.phi) / 2)
            α = min([1 α]...)  
        else
            valid = 0
        end
    end
    return α, valid
end

function calculate_alpha_for_death(
    par::Dict{String,Any},
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    model::Model,
    modeln::Model, 
    zetanew::Float64, 
    sig_zeta::Float64,
    kill::Int64
    )
    """
    Calculate the acceptance function (α) for the death action.

    Parameters
    ----------
    - `par`
    - `dataStruct`
    - `RayTraces`
    - `model`:      The current state of the model before the death action.
    - `modeln`:     The proposed state of the model after the death action.
    - `zetanew`:    The proposed new value of zeta.
    - `sig_zeta`:   The standard deviation for zeta.
    - `kill`:       The index of the cell to be deleted.

    Returns
    -------
    - `α`:          The calculated acceptance ratio for the death action.
    - `valid`:      Whether the proposed state is valid.

    """
    α, valid = 0, 1
    if par["prior"] == 1 # Uniform
        α = ((model.nCells) / (model.nCells - 1)) * ((par["zeta_scale"]) / (sig_zeta * sqrt(2 * pi))) * 
        exp(-((model.zeta[kill] - zetanew)^2) / (2 * sig_zeta^2) - (modeln.phi - model.phi) / 2)
        α = min([1 α]...)           
    elseif par["prior"] == 2 # Normal
        (modeln, dataStructn, valid) = sph_evaluate(modeln, dataStruct, RayTraces, par)
        α = ((model.nCells) / (model.nCells - 1)) * (par["zeta_scale"] / sig_zeta) * 
        exp((model.zeta[kill]^2) / (2 * par["zeta_scale"]^2) - (model.zeta[kill] - zetanew)^2 / (2 * sig_zeta^2) - 
        (modeln.phi - model.phi) / 2)
        α = min([1 α]...)
    elseif par["prior"] == 3 # Exponential
        if zetanew > 0
            α = ((model.nCells) / (model.nCells - 1)) * (par["zeta_scale"] / (sqrt(2 * pi) * sig_zeta)) * 
            exp(model.zeta[kill] / par["zeta_scale"] - (model.zeta[kill] - zetanew)^2 / (2 * sig_zeta^2) - 
            (modeln.phi - model.phi) / 2)
            α = min([1 α]...)  
        else
            valid = 0
        end
    end       
    return α, valid
end

function calculate_alpha_for_change(
    par::Dict{String,Any},
    model::Model,
    modeln::Model, 
    change::Int64,
    )
    """
    Calculate the acceptance function (α) for the change action.

    Parameters
    ----------
    - `par`
    - `model`:  The current state of the model before the change action.
    - `modeln`: The proposed state of the model after the change action.
    - `change`: The index of the cell to be changed.

    Returns
    -------
    - `α`:      The calculated acceptance ratio for the change action.
    - `valid`:  Whether the proposed state is valid.

    """
    α, valid = 0, 1
    if par["prior"] == 1 # Uniform
        if modeln.zeta[change] > 0 && modeln.zeta[change] < par["zeta_scale"]
            α = exp(-(modeln.phi - model.phi) / 2)
            α = min([1 α]...)
        else
            valid = 0;
        end      
    elseif par["prior"] == 2 # Normal
        α = exp( (model.zeta[change]^2 - modeln.zeta[change]^2) / (2 * par["zeta_scale"]^2) -
        (modeln.phi - model.phi) / 2)
        α = min([1 α]...)
    elseif par["prior"] == 3 # Exponential
        if modeln.zeta[change] > 0
            α = exp( (model.zeta[change] - modeln.zeta[change]) / par["zeta_scale"] -
            (modeln.phi - model.phi) / 2)
            α = min([1 α]...)  
        else
            α = 0
        end
    end   
    return α, valid
end

function calculate_alpha_for_move(
    modeln::Model,
    model::Model
    )
    """
    Calculate the acceptance function (α) for the move action.

    Parameters
    ----------
    - `model`:  The current state of the model before the move action.
    - `modeln`: The proposed state of the model after the move action.

    Returns
    -------
    - `α`:      The calculated acceptance ratio for the move action.

    """
    # α(move), Byrnes and Bezada, 2020, eq. 14
    α = exp(-(modeln.phi - model.phi) / 2)
    α = min([1 α]...) 
    return α
end