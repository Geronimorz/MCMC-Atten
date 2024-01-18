function make_dir(par::Dict{String,Any}, file_name::String)
    """
    Make directories for figures and models.
    """
    base_dir = file_name[1:end-4]
    if !isdir(base_dir)
        mkdir(base_dir)
    end

    FigDirLst = ["models", "figures", "data", "figures/nCells", "figures/llh",
                 "figures/xzUncertainty", "figures/xzMasked", "figures/xzMean",
                 "figures/xyUncertainty", "figures/xyMasked", "figures/xyMean",
                 "figures/xyContour", "figures/xzContour",
                 "figures/phi", "figures/nCells"]

    # Append directories based on parameters
    for i in par["z0"]
        push!(FigDirLst, joinpath("figures", "xyVoronoi_"*string(i)))
    end

    for i in par["y0"]
        push!(FigDirLst, joinpath("figures", "xzVoronoi_"*string(i)))
    end    

    # Create directories
    for iDir in FigDirLst
        full_path = joinpath(base_dir, iDir)
        if !isdir(full_path)
            mkdir(full_path)
        end
    end

    # copy the par file to the base_dir
    if !isfile(joinpath(base_dir, file_name))
        cp(file_name, joinpath(base_dir, file_name))
    end
end


function initialize_par(
    par::Dict{String,Any}, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}
    )
    """
    Return initialized parameters.

    Parameters
    ----------
    - `par` 
    - `dataStruct`

    Returns
    -------
    - `modelCount`:         Current number of models in a chain after burn_in # of iterations.
    - `savedModelCount`:    Current number of models in a chain saved in model_hist.
    - `makeplot`:           Whether to plot the voronoi diagram of each individual model. 
                            Takes a lot of time and memory, recommend to use during simple tests.
    - `sig_zeta`:           Standard deviation of zeta.
    - `sig_sig`:            Stardard deviation of sigma in perturbation #5. Currently unenabled.
    - `totalModelCount` :   In principle, the number of saved number of a chain.
    - `xVec`,`yVec`,`zVec`: The regular grids for interpolation in the model space.
    - `n`:                  Number of traces.
    - `xr`, `yr`, `zr`:     Standard deviation of the nuclei coordinate in perturbation #4.
    """
    modelCount, savedModelCount = 0, 0
    makeplot = par["plot_voronoi"]

    sig_zeta = par["zeta_scale"] * par["sig"] / 100
    sig_sig = par["max_sig"] * par["sig"] / 100

    totalModelCount = Int64((par["n_iter"] - par["burn_in"]) / par["keep_each"])

    xVec, yVec, zVec = dataStruct["xVec"], dataStruct["yVec"], dataStruct["zVec"]

    n = length(dataStruct["dataX"])
    xr = (par["sig"] / 100) * (max(xVec...) - min(xVec...))
    yr = (par["sig"] / 100) * (max(yVec...) - min(yVec...))
    zr = (par["sig"] / 100) * (max(zVec...) - min(zVec...))

    return modelCount, savedModelCount, makeplot, sig_zeta, sig_sig, totalModelCount, 
            xVec, yVec, zVec, n, xr, yr, zr
end

function load_latest_checkpoint(model_checkpoint_lists, chain::Int64)
    """
    Extract iteration numbers from the checkpoint filenames
    """
    iter_nums = [parse(Int, split(split(filename, "_")[end], "%")[1]) for filename in model_checkpoint_lists]
    
    # Find the checkpoint with the highest iteration number
    latest_checkpoint_index = argmax(iter_nums)
    latest_checkpoint = model_checkpoint_lists[latest_checkpoint_index]
    
    println("----Loading Chain $(chain)------")
    @time load_model = load(latest_checkpoint)
    
    return load_model, latest_checkpoint_index
end

function initialize_model_data(
    par::Dict{String,Any}, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}}
    )
    """
    Return the starting model and dataStruct.
    """
    valid = 0
    model = nothing
    @time while valid == 0
        (model, dataStruct, valid) = build_starting(par, dataStruct, RayTraces)
    end
    return model, dataStruct, valid, 1, [], [], [], 0, 0
end

function cleanup_checkpoints(model_checkpoint_lists, latest_checkpoint_index)
    """
    Delete all the checkpoints except the latest one.
    """
    for (index, checkpoint) in enumerate(model_checkpoint_lists)
        if index != latest_checkpoint_index
            rm(checkpoint, force=true)
        end
    end
end

function handle_checkpoints(
    par::Dict{String,Any}, 
    chain::Int64,
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}}
    )
    """
    Load the latest checkpoint if it exists, otherwise build a new starting model.
    """
    model_checkpoint_lists = glob(par["base_dir"] * "models/chain$(chain)_*")
    
    if isempty(model_checkpoint_lists)
        println("-------------build starting model-------------")
        return initialize_model_data(par, dataStruct, RayTraces)
    else
        load_model, latest_checkpoint_index = load_latest_checkpoint(model_checkpoint_lists, chain)
        cleanup_checkpoints(model_checkpoint_lists, latest_checkpoint_index)

        model = load_model["model"]
        dataStruct = load_model["dataStruct"]
        cellnumber_list = load_model["nCells"]
        phi_list = load_model["phi"]
        valid = 1
        iter_ind = Int(load_model["iter"])
        model_hist = load_model["burnin"] ? load_model["model_hist"] : []
        savedModelCount = load_model["burnin"] ? load_model["saved_#"] : 0
        modelCount = load_model["burnin"] ? load_model["modelCount"] : 0

        return model, dataStruct, valid, iter_ind, model_hist, cellnumber_list, phi_list, savedModelCount, modelCount
    end
end

function generate_new_coordinates(
    xVec::AbstractVector{Float64},
    yVec::AbstractVector{Float64},
    zVec::AbstractVector{Float64}
    )
    """
    Return the coordinate for a new cell.

    Parameters
    ----------
    - `xVec`,`yVec`,`zVec`: The regular grids in the model space.

    Returns
    -------
    - `xNew`,`yNew`,`zNew`: The coordinate of the new cell.
    """
    xNew   = rand() .* (max(xVec...) - min(xVec...)) .+ min(xVec...)
    yNew   = rand() .* (max(yVec...) - min(yVec...)) .+ min(yVec...)
    zNew   = rand() .* (max(zVec...) - min(zVec...)) .+ min(zVec...)
    return xNew, yNew, zNew
end
    
function move_new_coordinates(
    model::Model,
    move::Int64,
    xr::Float64,
    yr::Float64,
    zr::Float64
    )
    """
    Return the coordinate after randomly moving the cell.

    Parameters
    ----------
    - `model`:          The current state of the model.
    - `move`:           The index of the cell to be moved.
    - `xr`,`yr`,`zr`:   The standard deviation of the cell nuclei coordinate.

    Returns
    -------
    - `xCell_n`,`yCell_n`,`zCell_n`: The coordinate of the cell after moving.
    """
    xCell_n = rand(Normal(model.xCell[move], xr), 1)[1]
    yCell_n = rand(Normal(model.yCell[move], yr), 1)[1]
    zCell_n = rand(Normal(model.zCell[move], zr), 1)[1]
    return xCell_n, yCell_n, zCell_n
end

function delete_old_checkpoint(par, chain_id, iteration, percentage, percent_interval)
    """
    Delete the old checkpoint file.
    """
    old_checkpoint_pattern = joinpath(par["base_dir"], "models/chain$(chain_id)_iter*_$(percentage - percent_interval)%.jld")
    old_checkpoints = glob(old_checkpoint_pattern)
    for old_checkpoint in old_checkpoints
        rm(old_checkpoint)
    end
end

