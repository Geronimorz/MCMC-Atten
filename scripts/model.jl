mutable struct Model
    nCells      ::Float64
    xCell       ::Vector{Float64}
    yCell       ::Vector{Float64}
    zCell       ::Vector{Float64}
    zeta        ::Vector{Float64}
    phi         ::Float64
    # ptS         ::Vector{Float64}
    # tS          ::Vector{Float64}
    # predicted_traveltime
    # observed_traveltime
    likelihood  ::Float64
    action      ::Int64
    accept      ::Int64
    zeta_xz     ::Float64
    zeta_xy     ::Float64
end

function build_starting(
    par::Dict{String,Any}, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}
    )

    # log uniform distribution, in Byrnes and Bezada, 2020, eq. 11

    xVec = dataStruct["xVec"]
    yVec = dataStruct["yVec"]
    zVec = dataStruct["zVec"]
    # 11/14/23 yurong: change max_cells to prior_max_cells
    # nCells = floor.(exp.(rand(1) * log(par["max_cells"] / 
    #         par["min_cells"]) .+ log(par["min_cells"])))
    nCells = floor.(exp.(rand(1) * log(par["prior_max_cells"] / 
            par["min_cells"]) .+ log(par["min_cells"])))
    # nCells = floor.(rand(1)*(par.max_cells- par.min_cells)) .+ par.min_cells

    xCell = min(xVec...) .+ (max(xVec...) - min(xVec...)) .* rand(Int(nCells[1]))
    yCell = min(yVec...) .+ (max(yVec...) - min(yVec...)) .* rand(Int(nCells[1]))
    zCell = min(zVec...) .+ (max(zVec...) - min(zVec...)) .* rand(Int(nCells[1]))


    if par["prior"] == 1 
        # Uniform 
            # for absolute t*
        zeta = rand(Int(nCells[1])) .* par["zeta_scale"]
            # for relative t*
        # zeta = rand(Int(nCells[1])) .* par["zeta_scale"] * 2 .- par["zeta_scale"]  
    elseif par["prior"] == 2
        # Normal 
        zeta = rand(Normal(0, par["zeta_scale"]), Int(nCells[1]))
    elseif par["prior"] == 3
        # Exponential 
        zeta = -log.(rand(Int(nCells[1]))) .* par["zeta_scale"]
    end

    model = Model(nCells[1],
        xCell, yCell, zCell,
        zeta,
        -1, -1, -1, -1, -1, -1
    )
    # setindex!(dataStruct, model["allSig"] .* ones(length(dataStruct["tS"])), "allSig")
    if par["coordinates"] == 1
        (model, dataStruct, valid) = sph_evaluate(model, dataStruct, RayTraces, par)
    elseif par["coordinates"] == 2
        (model, dataStruct, valid) = cart_evaluate(model, dataStruct, RayTraces, par)
    end
    
    valid = 1 # Valid is designed to deal with discontinuity and other features. Will work on it afterwards.
    return model, dataStruct, valid
end

function cart_evaluate(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    par::Dict{String,Any}
    )
    valid       = 1
    likelihood  = 1
    phi         = 1
    model.phi = phi
    model.likelihood = likelihood

    if par["debug_prior"] == 1
        return model, dataStruct, valid
    end

    (m, n) = size(RayTraces["rayX"])
    ptS = zeros(n)
    observed_traveltime = zeros(n,1)
    predicted_traveltime = zeros(n,1)
    Threads.@threads for i in 1:n
        zeta0 = cart_interpolation(par, model, RayTraces["rayX"][:,i], RayTraces["rayY"][:,i], RayTraces["rayZ"][:,i]) 
        rayzeta = []
 
        endpoint = length(zeta0)
        rayzeta = 0.5 .* (zeta0[1:endpoint - 1] + zeta0[2:endpoint])
    
        rayl = RayTraces["rayL"][:,i]
        index = findall(x -> isnan(x), rayl)
        
        if isempty(index)
            ptS[i] = sum(rayl .* RayTraces["rayU"][:,i] .* (rayzeta ./ 1000))
            predicted_traveltime[i] = sum(rayl .* RayTraces["rayU"][:,i])
        else
            rayl = RayTraces["rayL"][1:index[1]-1,i]  
            rayu = RayTraces["rayU"][1:index[1]-1,i]
        # compare travel time
            ptS[i] = sum(rayl .* rayu .* (rayzeta ./ 1000))
            predicted_traveltime[i] = sum(rayl .* rayu)
        end
        observed_traveltime[i] = 1000*dataStruct["tS"][i]/dataStruct["allaveatten"][i]
    end

    (m, n) = size(RayTraces["rayZ"])

    tS = dataStruct["tS"]

    C = 0
    for k in 1:length(dataStruct["allSig"])
        C += (ptS - tS)[k].^2 .* 1.0 / dataStruct["allSig"][k][1]^2
    end
    model.phi = C
    # model.ptS = ptS
    # model.tS = tS
    # model.predicted_traveltime = predicted_traveltime
    # model.observed_traveltime = observed_traveltime

    # yurong 03/02/23 change llh fuction based on central limit theorem
    likelihood = sum(-log.(dataStruct["allSig"] * sqrt(2 * pi))) - 
    sum(0.5 * ((vec(ptS) .- vec(tS)) ./ vec(dataStruct["allSig"])).^2)
    # likelihood = sum(-log.(dataStruct["allSig"] * sqrt(2 * pi)) * length(tS)) - 
    # sum(0.5 * ((vec(ptS) .- vec(tS)) ./ vec(dataStruct["allSig"])).^2)


    model.likelihood = likelihood
    
    return model, dataStruct, valid
end

function sph_evaluate(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    par::Dict{String,Any}
    )
    valid       = 1
    likelihood  = 1
    phi         = 1
    model.phi = phi
    model.likelihood = likelihood


    if par["debug_prior"] == 1
        return model, dataStruct, valid
    end

    (m, n) = size(RayTraces["rayX"])
    ptS = zeros(n)
    observed_traveltime = zeros(n,1)
    predicted_traveltime = zeros(n,1)
    Threads.@threads for i in 1:n
        zeta0 = sph_interpolation(par, model, RayTraces["rayX"][:,i], RayTraces["rayY"][:,i], RayTraces["rayZ"][:,i]) 
        rayzeta = []
 
        endpoint = length(zeta0)
        rayzeta = 0.5 .* (zeta0[1:endpoint - 1] + zeta0[2:endpoint])
    
        rayl = RayTraces["rayL"][:,i]
        index = findall(x -> isnan(x), rayl)
        
        if isempty(index)
            ptS[i] = sum(rayl .* RayTraces["rayU"][:,i] .* (rayzeta ./ 1000))
            predicted_traveltime[i] = sum(rayl .* RayTraces["rayU"][:,i])
        else
            rayl = RayTraces["rayL"][1:index[1]-1,i]  
            rayu = RayTraces["rayU"][1:index[1]-1,i]
        # compare travel time
            ptS[i] = sum(rayl .* rayu .* (rayzeta ./ 1000))
            predicted_traveltime[i] = sum(rayl .* rayu)
        end
        observed_traveltime[i] = 1000*dataStruct["tS"][i]/dataStruct["allaveatten"][i]
    end

    (m, n) = size(RayTraces["rayZ"])

    tS = dataStruct["tS"]

    C = 0
    for k in 1:length(dataStruct["allSig"])
        C += (ptS - tS)[k].^2 .* 1.0 / dataStruct["allSig"][k][1]^2
    end
    model.phi = C
    # model.ptS = ptS
    # model.tS = tS
    # model.predicted_traveltime = predicted_traveltime
    # model.observed_traveltime = observed_traveltime

    # yurong 03/02/23 change llh fuction based on central limit theorem
    likelihood = sum(-log.(dataStruct["allSig"] * sqrt(2 * pi))) - 
    sum(0.5 * ((vec(ptS) .- vec(tS)) ./ vec(dataStruct["allSig"])).^2)
    # likelihood = sum(-log.(dataStruct["allSig"] * sqrt(2 * pi)) * length(tS)) - 
    # sum(0.5 * ((vec(ptS) .- vec(tS)) ./ vec(dataStruct["allSig"])).^2)


    model.likelihood = likelihood
    
    return model, dataStruct, valid
end

function update_model_birth!(
    modeln::Model, 
    xNew::Float64, 
    yNew::Float64, 
    zNew::Float64, 
    zetanew::Float64)
    """
    Update the `modeln` with the newly generated cell coordinates (xNew, yNew, zNew) and the associated zeta value (zetanew).

    Parameters
    ----------
    - `modeln`:             The model to be updated.
    - `xNew`,`yNew`,`zNew`: The new coordinates to add to the model.
    - `zetanew`:            The new zeta value associated with the new coordinates.
    """
    append!(modeln.xCell, xNew)
    append!(modeln.yCell, yNew)
    append!(modeln.zCell, zNew)
    append!(modeln.zeta, zetanew)
    modeln.nCells += 1
end

# function update_model_birth(
#     modeln::Model, 
#     xNew::Float64, 
#     yNew::Float64, 
#     zNew::Float64, 
#     zetanew::Float64)
#     """
#     Update the `modeln` with the newly generated cell coordinates (xNew, yNew, zNew) and the associated zeta value (zetanew).

#     Parameters
#     ----------
#     - `modeln`:             The model to be updated.
#     - `xNew`,`yNew`,`zNew`: The new coordinates to add to the model.
#     - `zetanew`:            The new zeta value associated with the new coordinates.
#     """
#     append!(modeln.xCell, xNew)
#     append!(modeln.yCell, yNew)
#     append!(modeln.zCell, zNew)
#     append!(modeln.zeta, zetanew)
#     modeln.nCells += 1

#     return modeln
# end

function update_model_death!(
    modeln::Model, 
    kill::Int64
    )
    """
    Update the `modeln` by deleting the randomly selected cell.

    Parameters
    ----------
    - `modeln`: The model to be updated.
    - `kill`:   The index of the cell to be deleted.
    """
    deleteat!(modeln.xCell, kill)
    deleteat!(modeln.yCell, kill)
    deleteat!(modeln.zCell, kill)
    deleteat!(modeln.zeta, kill)
    modeln.nCells -= 1
end

function update_model_change!(
    modeln::Model, 
    zetanew::Float64,
    change::Int64
    )
    """
    Update the `modeln` by changing the zeta value of the randomly selected cell.

    Parameters
    ----------
    - `modeln`:     The model to be updated.
    - `zetanew`:    The new zeta value of the cell.
    - `change`:     The index of the cell to be changed.
    """
    modeln.zeta[change] = zetanew
end

function update_model_move!(
    modeln::Model, 
    move::Int64,
    xCell_n::Float64, 
    yCell_n::Float64, 
    zCell_n::Float64
    )
    """
    Update the `modeln` by moving the randomly selected cell.

    Parameters
    ----------
    - `modeln`: The model to be updated.
    - `move`:   The index of the cell to be moved.
    - `xCell_n`,`yCell_n`,`zCell_n`:
                The coordinate of the cell after moving.
    """
    modeln.xCell[move] = xCell_n
    modeln.yCell[move] = yCell_n
    modeln.zCell[move] = zCell_n
end

function sph_perform_birth_action(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    par::Dict{String,Any},
    sig_zeta::Float64, 
    xVec::AbstractVector{Float64},
    yVec::AbstractVector{Float64},
    zVec::AbstractVector{Float64}
    )
    if model.nCells < par["max_cells"]
        xNew, yNew, zNew = generate_new_coordinates(xVec, yVec, zVec)
        czeta = sph_interpolation(par, model, xNew, yNew, zNew)[1]
        zetanew = generate_new_zeta(czeta, sig_zeta)

        modeln = deepcopy(model)
        update_model_birth!(modeln, xNew, yNew, zNew, zetanew)
        (modeln, dataStructn, valid) = sph_evaluate(modeln, dataStruct, RayTraces, par)
        α, valid = calculate_alpha_for_birth(par, dataStruct, RayTraces, model, modeln, czeta, zetanew, sig_zeta)
        a1 = rand()
        if a1 < α && valid == 1
            modeln.accept = 1
            return modeln, dataStructn
        end
    end
    return model, dataStruct
end

function sph_perform_death_action(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    par::Dict{String,Any},
    sig_zeta::Float64
    )
    if model.nCells > par["min_cells"]
        kill = Int(rand(1:model.nCells))

        modeln =  deepcopy(model)
        update_model_death!(modeln, kill)
        (modeln, dataStructn, valid) = sph_evaluate(modeln, dataStruct, RayTraces, par)

        zetanew = sph_interpolation(par, modeln, [model.xCell[kill]], [model.yCell[kill]], [model.zCell[kill]])[1]

        α, valid = calculate_alpha_for_death(par, dataStruct, RayTraces, model, modeln, zetanew, sig_zeta, kill)
        a1 = rand()

        if a1 < α && valid == 1
            modeln.accept = 1
            return modeln, dataStructn
        end
    end
    return model, dataStruct 
end

function sph_perform_change_action(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    par::Dict{String,Any},
    sig_zeta::Float64
    )

    change = Int(rand(1:model.nCells))
    zetanew = generate_new_zeta(model.zeta[change], sig_zeta)
    modeln =  deepcopy(model)
    update_model_change!(modeln, zetanew, change)

    (modeln, dataStructn, valid) = sph_evaluate(modeln, dataStruct, RayTraces,par)

    α, valid = calculate_alpha_for_change(par, model, modeln, change)

    if rand() < α && valid == 1
        modeln.accept = 1
        return modeln, dataStructn
    end
    
    return model, dataStruct 
end

function sph_perform_move_action(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    par::Dict{String,Any},
    xr::Float64,
    yr::Float64,
    zr::Float64,
    xVec::AbstractVector{Float64},
    yVec::AbstractVector{Float64},
    zVec::AbstractVector{Float64}
    )
    if model.nCells > 0
        move = Int(rand(1:model.nCells))
        modeln = deepcopy(model)
        xCell_n, yCell_n, zCell_n = move_new_coordinates(model,move,xr,yr,zr)
        # define α in advance
        α = NaN
        if xCell_n >= min(xVec...) && xCell_n <= max(xVec...) &&
            yCell_n >= min(yVec...) && yCell_n <= max(yVec...) &&
            zCell_n >= min(zVec...) && zCell_n <= max(zVec...)
            update_model_move!(modeln,move,xCell_n, yCell_n, zCell_n)
            (modeln, dataStructn, valid) = sph_evaluate(modeln, dataStruct, RayTraces, par)

            α = calculate_alpha_for_move(modeln, model) 
        else
            valid = 0 
        end
        if rand() < α && valid == 1
            modeln.accept = 1
            return modeln, dataStructn
        end
    end
    return model, dataStruct 
end

function cart_perform_birth_action(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    par::Dict{String,Any},
    sig_zeta::Float64, 
    xVec::AbstractVector{Float64},
    yVec::AbstractVector{Float64},
    zVec::AbstractVector{Float64}
    )
    if model.nCells < par["max_cells"]
        xNew, yNew, zNew = generate_new_coordinates(xVec, yVec, zVec)
        czeta = cart_interpolation(par, model, xNew, yNew, zNew)[1]
        zetanew = generate_new_zeta(czeta, sig_zeta)

        modeln = deepcopy(model)
        update_model_birth!(modeln, xNew, yNew, zNew, zetanew)

        (modeln, dataStructn, valid) = cart_evaluate(modeln, dataStruct, RayTraces, par)
           
        α, valid = calculate_alpha_for_birth(par, dataStruct, RayTraces, model, modeln, czeta, zetanew, sig_zeta)

        if rand() < α && valid == 1
            modeln.accept = 1
            return modeln, dataStructn
        end
    end
    return model, dataStruct
end

function cart_perform_death_action(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    par::Dict{String,Any},
    sig_zeta::Float64
    )
    if model.nCells > par["min_cells"]
        kill = Int(rand(1:model.nCells))

        modeln =  deepcopy(model)
        update_model_death!(modeln, kill)

        (modeln, dataStructn, valid) = cart_evaluate(modeln, dataStruct, RayTraces, par)

        zetanew = cart_interpolation(par, modeln, [model.xCell[kill]], [model.yCell[kill]], [model.zCell[kill]])[1]

        α, valid = calculate_alpha_for_death(par, dataStruct, RayTraces, model, modeln, zetanew, sig_zeta, kill)

        if rand() < α && valid == 1
            modeln.accept = 1
            return modeln, dataStructn
        end
    end
    return model, dataStruct 
end

function cart_perform_change_action(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    par::Dict{String,Any},
    sig_zeta::Float64
    )

    change = Int(rand(1:model.nCells))
    zetanew = generate_new_zeta(model.zeta[change], sig_zeta)
    modeln =  deepcopy(model)
    update_model_change!(modeln, zetanew, change)

    (modeln, dataStructn, valid) = cart_evaluate(modeln, dataStruct, RayTraces,par)

    α, valid = calculate_alpha_for_change(par, model, modeln, change)

    if rand() < α && valid == 1
        modeln.accept = 1
        return modeln, dataStructn
    end
    
    return model, dataStruct 
end

function cart_perform_move_action(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N},
    RayTraces::Dict{String,Array{Float64,2}}, 
    par::Dict{String,Any},
    xr::Float64,
    yr::Float64,
    zr::Float64,
    xVec::AbstractVector{Float64},
    yVec::AbstractVector{Float64},
    zVec::AbstractVector{Float64}
    )
    if model.nCells > 0
        move = Int(rand(1:model.nCells))
        modeln = deepcopy(model)
        xCell_n, yCell_n, zCell_n = move_new_coordinates(model,move,xr,yr,zr)
        α = 0
        if xCell_n >= min(xVec...) && xCell_n <= max(xVec...) &&
            yCell_n >= min(yVec...) && yCell_n <= max(yVec...) &&
            zCell_n >= min(zVec...) && zCell_n <= max(zVec...)
            update_model_move!(modeln,move,xCell_n, yCell_n, zCell_n)
            (modeln, dataStructn, valid) = cart_evaluate(modeln, dataStruct, RayTraces, par)

            α = calculate_alpha_for_move(modeln, model)
        else
            valid = 0 
        end
        if rand() < α && valid == 1
            modeln.accept = 1
            return modeln, dataStructn
        end
    end
    return model, dataStruct 
end

function perform_sigma_action(dataStruct,RayTraces,sig_sig,par,model)
    """
    Change the standard deviation of sigma.
    ! To be updated.
    """
    sig_n = rand(Normal(dataStruct["allSig"][1], sig_sig), 1)[1] 
    allSig_n = sig_n * ones(size(dataStruct["tS"]))
    if sig_n > 0 && sig_n < par["max_sig"]
        dataStructn        = deepcopy(dataStruct)
        dataStructn["allSig"] = allSig_n
        (modeln, dataStructn, valid) = sph_evaluate(model, dataStructn, RayTraces, par)
        # α(error), Byrnes and Bezada, 2020, eq. 18
        α = log(dataStruct["allSig"][1] / sig_n) * n - (modeln.phi - model.phi) / 2
        α = min([log(1) α]...)

        if rand() < α && valid == 1
            modeln.accept = 1
            return modeln, dataStructn
        end
    end
    return model,dataStruct
end