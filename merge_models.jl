using JLD, HDF5
using Distributed 
using Glob, Plots, YAML

@everywhere include("./scripts/model.jl")
@everywhere include("./scripts/utils.jl")
@everywhere include("./scripts/stats.jl")
@everywhere include("./scripts/plot_sub.jl")
@everywhere include("./scripts/interp.jl")
@everywhere include("./scripts/load.jl")
@everywhere include("./scripts/inversion_function.jl")

function merge_models(file_name)
    """
    This function merges the saved models from different chains and interpolates them to a grid
    1: Arithmetic average the models from different chains
    2: Bayesian average the models from different chains
    3: Inverse variance weighted average the models from different chains

    Parameters
    ----------
    - `file_name` : string
        The name of the yaml file that contains the parameters for the inversion
    """
    println("--------Loading Data-------")
    println("Loading parameters from file: "*file_name)
    par = load_par_from_yml(file_name)
    (dataStruct, RayTraces) = load_data_Tonga(par)
    make_dir(par, file_name)
    println(string(length(dataStruct["tS"]))*" t* loaded")
    println("--------Data Loaded-------\n")

    std_threshold = par["std_threshold"]
    transparency_threshold = par["transparency_threshold"]
    lat0 = par["lat0"]
    lon0 = par["lon0"]
    beta = par["beta"]

    if isdir(par["base_dir"] * "TongaAttenData") == false
        mkdir(par["base_dir"] * "TongaAttenData")
    end

    println("--------Loading Models-------")
    model_hist = []
    for chain in 1:par["n_chains"]
        model_checkpoint_lists = glob(par["base_dir"] * "models/chain" * string(chain) * "_*")
        if length(model_checkpoint_lists) == 0
            println("ERROR: Couldn't Find Models For Chain" * string(chain))
        else
            # load the newest model
            split1      = split.(model_checkpoint_lists,"%")
            split2      = [split1[i][1] for i in 1:length(split1)]
            split3      = split.(split2,"_")
            split4      = [split3[i][end] for i in 1:length(split3)]
            model_ind   =  findmax(parse.(Float64,split4))[2]
            load_model  = load(model_checkpoint_lists[model_ind])
            for irm in 1:length(model_checkpoint_lists)
                if irm == model_ind
                    continue
                end
                rm(sort(model_checkpoint_lists)[irm])
            end
            push_models = load_model["model_hist"]
            push!(model_hist,push_models)
            println("Chain " * string(chain) * " Loaded with " * string(length(push_models)) * " Models")
        end
    end
    println("--------Models Loaded-------\n")

    println("--------Interpolating Models-------")
    if par["merge_method"] == 1
        println("Calculating Arithmetic Average")
    elseif par["merge_method"] == 2
        println("Calculating Bayesian Model Average")
    elseif par["merge_method"] == 3
        println("Calculating Inverse Variance Weighted Average")
    end
    println("Models with std > " * string(std_threshold) * " will be masked out")


    if par["xyMap"] == true
        println("Interpolating to a grid in Map View at depths of " * string(par["z0"]) * " km")
        if par["coordinates"] == 1
            for l0 in par["z0"]
                m_xy = []
                m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
                mcount = 0

                for i = 1:length(model_hist)

                    for j = 1:length(model_hist[i])

                        m  = [ sph_v_nearest(xs, ys, l0, model_hist[i][j])
                            for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]
                        append!(m_xy,[m])

                    end
                end

                if par["merge_method"] == 1
                    (model_mean_xy, poststd_xy, mask_model_xy) = merge_arithmetic(m_xy, std_threshold)
                end
                output_merged_models(par, dataStruct, model_mean_xy, poststd_xy, mask_model_xy, std_threshold, transparency_threshold, l0)
            end
        elseif par["coordinates"] == 2
            for l0 in par["z0"]
                m_xy = []
                m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
                mcount = 0

                for i = 1:length(model_hist)

                    for j = 1:length(model_hist[i])

                        m  = [ cart_v_nearest(xs, ys, l0, model_hist[i][j])
                            for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]
                        append!(m_xy,[m])

                    end
                end

                if par["merge_method"] == 1
                    (model_mean_xy, poststd_xy, mask_model_xy) = merge_arithmetic(m_xy, std_threshold)
                end
                output_merged_models(par, dataStruct, model_mean_xy, poststd_xy, mask_model_xy, std_threshold, transparency_threshold, l0)
            end
        end
    end

    if par["xzMap"] == true
        open(par["base_dir"] * "TongaAttenData/Interpolated_Tonga_Atten_Model.txt", "w") do io
            write(io, "Longitude\t Latitude\t 1000/Q\t Uncertainty\t 1000/Q with Mask\n")
        end
        if par["coordinates"] == 1
            println("Interpolating to a grid of the model with a resolution of " * string(par["latlonnodeSpacing"]) * " degree in Longitude and Latitude 
            and " * string(par["ZnodeSpacing"]) * " km in Depth")
            for z0 in dataStruct["zVec"]
                m_iz0 = []
                m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
        
                for i = 1:length(model_hist)
        
                    for j = 1:length(model_hist[i])
        
                        m  = [ sph_v_nearest(xs, ys, z0, model_hist[i][j])
                            for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]
                        append!(m_iz0,[m])
        
                    end
                end
                (model_mean, poststd, mask_model) = merge_arithmetic(m_iz0, std_threshold)
                open(par["base_dir"] * "TongaAttenData/Interpolated_Tonga_Atten_Model.txt", "a") do io
                    for i in 1:length(dataStruct["xVec"])
                        for j in 1:length(dataStruct["yVec"])
                            write(io, string(dataStruct["xVec"][i]) *"\t "* string(dataStruct["yVec"][j]) *
                                "\t " * string(z0) *"\t " * string(model_mean[i,j]) *"\t "* 
                                string(poststd[i,j]) *"\t " * string(mask_model[i,j]) *"\n")
                        end
                    end
                end
            end
        elseif par["coordinates"] == 2
            println("Interpolating to a grid of the model with a resolution of " * string(par["XYnodeSpacing"]) * " km in Longitude and Latitude 
            and " * string(par["ZnodeSpacing"]) * " km in Depth")
            for z0 in dataStruct["zVec"]
                m_iz0 = []
                m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
        
                for i = 1:length(model_hist)
        
                    for j = 1:length(model_hist[i])
        
                        m  = [ cart_v_nearest(xs, ys, z0, model_hist[i][j])
                            for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]
                        append!(m_iz0,[m])
        
                    end
                end
                (model_mean, poststd, mask_model) = merge_arithmetic(m_iz0, std_threshold)
                open(par["base_dir"] * "TongaAttenData/Interpolated_Tonga_Atten_Model.txt", "a") do io
                    for i in 1:length(dataStruct["xVec"])
                        for j in 1:length(dataStruct["yVec"])
                            (lon,lat) = xy2lonlat(lon0,lat0,beta,dataStruct["xVec"][i],dataStruct["yVec"][j])
                            write(io, string(lon) *"\t "* string(lat) *"\t " * string(z0) *"\t "
                                * string(model_mean[i,j]) *"\t "* string(poststd[i,j]) *"\t " 
                                * string(mask_model[i,j]) *"\n")
                        end
                    end
                end
            end
        end
    end
    println("--------Models Interpolated-------\n")


end


function output_merged_models(par, dataStruct, 
    model_mean_xy, poststd_xy, mask_model_xy, 
    std_threshold, transparency_threshold, l0)
    """ 
    This function outputs the merged models to files

    Parameters
    ----------
    - `par` : Dict
    - `dataStruct` : Dict
    - `model_mean_xy` : Array{Float64,2}
        The mean of the models
    - `poststd_xy` : Array{Float64,2}
        The standard deviation of the models
    - `mask_model_xy` : Array{Float64,2}
        The mask of the models
    - `l0` : Float64
        The depth of the models
    """
    open(par["base_dir"] * "TongaAttenData/Tonga_Map_Mask_"*string(l0)*".txt", "w") do io
        for i in 1:length(dataStruct["xVec"])
            for j in 1:length(dataStruct["yVec"])
                write(io, string(dataStruct["xVec"][i]) *"  "* string(dataStruct["yVec"][j]) *
                        "  " * string(mask_model_xy[i,j]) *"\n")
            end
        end
    end

    open(par["base_dir"] * "TongaAttenData/Tonga_Map_Model_"*string(l0)*".txt", "w") do io
        for i in 1:length(dataStruct["xVec"])
            for j in 1:length(dataStruct["yVec"])
                write(io, string(dataStruct["xVec"][i]) *"  "* string(dataStruct["yVec"][j]) *
                        "  " * string(model_mean_xy[i,j]) *"\n")
            end
        end
    end

    open(par["base_dir"] * "TongaAttenData/Tonga_Map_Uncertainty_"*string(l0)*".txt", "w") do io
        for i in 1:length(dataStruct["xVec"])
            for j in 1:length(dataStruct["yVec"])
                write(io, string(dataStruct["xVec"][i]) *"  "* string(dataStruct["yVec"][j]) *
                        "  " * string(poststd_xy[i,j]) *"\n")
            end
        end
    end

    open(par["base_dir"] * "TongaAttenData/Tonga_Map_Transparency_"*string(l0)*".txt", "w") do io
        for i in 1:length(dataStruct["xVec"])
            for j in 1:length(dataStruct["yVec"])
                if poststd_xy[i,j]<std_threshold
                    transparency = 0
                elseif poststd_xy[i,j]<transparency_threshold
                    transparency = (poststd_xy[i,j]-std_threshold)/(transparency_threshold-std_threshold)
                else
                    transparency = 1
                end
                write(io, string(dataStruct["xVec"][i]) *"  "* string(dataStruct["yVec"][j]) *
                        "  " * string(transparency) *"\n")
            end
        end
    end
    println("- Map view at depth " * string(l0) * " km saved.")

end
