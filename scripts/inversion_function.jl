using JLD, HDF5
using Dates, MAT
using Glob
using Random, Distributions

function perform_inversion(
    par::Dict{String,Any}, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}}
    )
    """
    Function to perform inversion based on coordinate system
    """
    coord_system = par["coordinates"]
    inv_func = coord_system == 1 ? sph_inversion_function : cart_inversion_function
    println("Inverting data in $(coord_system == 1 ? "Spherical" : "Cartesian") Coordinate System")

    return pmap(x -> inv_func(par, dataStruct, RayTraces, x), 1:par["n_chains"])
end

function delete_old_checkpoint(chain_id, iteration, percentage, percent_interval)
    # Find and remove the old checkpoint file
    old_checkpoint_pattern = "./models/chain$(chain_id)_iter*_$(percentage - percent_interval)%.jld"
    
    old_checkpoints = glob(old_checkpoint_pattern)
    for old_checkpoint in old_checkpoints
        rm(old_checkpoint)
    end
end

function sph_inversion_function(
    par::Dict{String,Any}, 
    dataStruct1::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}}, 
    chain::Int64
    )
    
    # This line is needed to make sure that the chains are different from one another
    Random.seed!(chain * Int(round(mod(Dates.millisecond(now()), 1e3)))) 
    dataStruct = deepcopy(dataStruct1)

    # Initialze parameters needed for the iterations ####
    modelCount, savedModelCount, makeplot, sig_zeta, sig_sig, totalModelCount, 
        xVec, yVec, zVec, n, xr, yr, zr = initialize_par(par, dataStruct)

    valid = 0 
    model, dataStruct, valid, iter_ind, model_hist, cellnumber_list, phi_list, savedModelCount, modelCount = 
        handle_checkpoints(par,chain,dataStruct,RayTraces)
        
    α = 0
    for iter in iter_ind:par["n_iter"]
        # By using t* error from t* inversion, we no longer use ACTION5 to perturb "allSig"
        action = Int(rand(1:4))
        # 11/14/23 yurong: lower the possibility of action 1
        # actions = [1,2,2,3,3,4,4]
        # action = rand(actions)

        model.action = action
        model.accept = 0
        if action == 1 # Birth
            model, dataStruct = sph_perform_birth_action(model, dataStruct, RayTraces, par, sig_zeta, xVec, yVec, zVec)
        elseif action == 2 # Death
            model, dataStruct = sph_perform_death_action(model, dataStruct, RayTraces, par, sig_zeta)
        elseif action == 3 # change ζ in a cell
            model, dataStruct = sph_perform_change_action(model, dataStruct, RayTraces, par, sig_zeta)
        elseif action == 4 # move
            model, dataStruct = sph_perform_move_action(model, dataStruct, RayTraces, par, xr, yr, zr, xVec, yVec, zVec)
        elseif action == 5 # change sigma
            model, dataStruct = perform_sigma_action(dataStruct,RayTraces,sig_sig,par,model)
        end

        # plot voronoi diagram over iterations
        if par["plot_voronoi"] == true
            ENV["GKSwstype"] = "nul"
            plot_voronoi(model, dataStruct, par, chain, iter)
        end

        push!(cellnumber_list, model.nCells)
        push!(phi_list, model.phi)
        
        (m, n) = size(RayTraces["rayZ"])

        if iter >= par["burn_in"]
            modelCount += 1
            if mod(modelCount, par["keep_each"]) == 0
                savedModelCount += 1 
                push!(model_hist, model)
            end

            if -1e-9 < mod(iter*100/par["n_iter"],par["save_percent"]) < 1e-9
                CurrentModel = model
                percentage = Int(100 * iter / par["n_iter"])
                modelname = par["base_dir"] * "models/chain" * string(chain) * "_iter" * string(Int(iter)) * "_" * string(percentage)* "%.jld"
                println("----Saving Chain" * string(chain) * " "* string(Int(100 * iter / par["n_iter"]))* "% models------")
                @time save(modelname,"model",CurrentModel,"dataStruct",dataStruct,"iter",iter,"saved_#",savedModelCount,"modelCount",modelCount,"model_hist",model_hist,"burnin",true,"nCells",cellnumber_list,"phi",phi_list)
                delete_old_checkpoint(par, chain, iter, percentage, par["save_percent"])
                CurrentModel = nothing
            end
            
        elseif -1e-9 < mod(iter*100/par["n_iter"],par["save_percent"]) < 1e-9
            CurrentModel = model
            percentage = Int(100 * iter / par["n_iter"])
            modelname = par["base_dir"] * "models/chain" * string(chain) * "_iter" * string(Int(iter)) * "_" * string(percentage)* "%.jld"
            println("----Saving Chain" * string(chain) * " "* string(Int(100 * iter / par["n_iter"]))* "% models------")
            @time save(modelname,"model",CurrentModel,"dataStruct",dataStruct,"iter",Int64(iter),"burnin",false,"nCells",cellnumber_list,"phi",phi_list)
            delete_old_checkpoint(par, chain, iter, percentage, par["save_percent"])
            CurrentModel = nothing
        end

        if mod(iter, par["print_each"]) == 0
            println("Chain #", chain, " at ", 100 * iter / par["n_iter"], "% with a phi of ", model.phi)
        end

        

    end
    
    return model_hist
end

function cart_inversion_function(
    par::Dict{String,Any}, 
    dataStruct1::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}}, 
    chain::Int64
    )
    
    # This line is needed to make sure that the chains are different from one another
    Random.seed!(chain * Int(round(mod(Dates.millisecond(now()), 1e3)))) 
    dataStruct = deepcopy(dataStruct1)

    # Initialze parameters needed for the iterations ####
    modelCount, savedModelCount, makeplot, sig_zeta, sig_sig, totalModelCount, 
        xVec, yVec, zVec, n, xr, yr, zr = initialize_par(par, dataStruct)

    valid = 0 
    model, dataStruct, valid, iter_ind, model_hist, cellnumber_list, phi_list, savedModelCount, modelCount = 
        handle_checkpoints(par,chain,dataStruct,RayTraces)
        
    α = 0
    for iter in iter_ind:par["n_iter"]
        # By using t* error from t* inversion, we no longer use ACTION5 to perturb "allSig"
        action = Int(rand(1:4))
        # 11/14/23 yurong: lower the possibility of action 1
        # actions = [1,2,2,3,3,4,4]
        # action = rand(actions)

        model.action = action
        model.accept = 0
        if action == 1 # Birth
            model, dataStruct = cart_perform_birth_action(model, dataStruct, RayTraces, par, sig_zeta, xVec, yVec, zVec)
        elseif action == 2 # Death
            model, dataStruct = cart_perform_death_action(model, dataStruct, RayTraces, par, sig_zeta)
        elseif action == 3 # change ζ in a cell
            model, dataStruct = cart_perform_change_action(model, dataStruct, RayTraces, par, sig_zeta)
        elseif action == 4 # move
            model, dataStruct = cart_perform_move_action(model, dataStruct, RayTraces, par, xr, yr, zr, xVec, yVec, zVec)
        elseif action == 5 # change sigma
            model, dataStruct = perform_sigma_action(dataStruct,RayTraces,sig_sig,par,model)

        end

        # plot voronoi diagram over iterations
        if par["plot_voronoi"] == true
            ENV["GKSwstype"] = "nul"
            plot_voronoi(model, dataStruct, par, chain, iter)
        end

        push!(cellnumber_list, model.nCells)
        push!(phi_list, model.phi)
        
        (m, n) = size(RayTraces["rayZ"])
        if iter >= par["burn_in"]
            modelCount += 1
            if mod(modelCount, par["keep_each"]) == 0
                savedModelCount += 1 
                push!(model_hist, model)
            end

            if -1e-9 < mod(iter*100/par["n_iter"],par["save_percent"]) < 1e-9
                CurrentModel = model
                percentage = Int(100 * iter / par["n_iter"])
                modelname = par["base_dir"] * "models/chain" * string(chain) * "_iter" * string(Int(iter)) * "_" * string(percentage)* "%.jld"
                println("----Saving Chain" * string(chain) * " "* string(Int(100 * iter / par["n_iter"]))* "% models------")
                @time save(modelname,"model",CurrentModel,"dataStruct",dataStruct,"iter",iter,"saved_#",savedModelCount,"modelCount",modelCount,"model_hist",model_hist,"burnin",true,"nCells",cellnumber_list,"phi",phi_list)
                delete_old_checkpoint(par, chain, iter, percentage, par["save_percent"])
                CurrentModel = nothing
            end
            
        elseif -1e-9 < mod(iter*100/par["n_iter"],par["save_percent"]) < 1e-9
            CurrentModel = model
            percentage = Int(100 * iter / par["n_iter"])
            modelname = par["base_dir"] * "models/chain" * string(chain) * "_iter" * string(Int(iter)) * "_" * string(percentage)* "%.jld"
            println("----Saving Chain" * string(chain) * " "* string(Int(100 * iter / par["n_iter"]))* "% models------")
            @time save(modelname,"model",CurrentModel,"dataStruct",dataStruct,"iter",Int64(iter),"burnin",false,"nCells",cellnumber_list,"phi",phi_list)
            delete_old_checkpoint(par, chain, iter, percentage, par["save_percent"])
            CurrentModel = nothing
        end

        if mod(iter, par["print_each"]) == 0
            println("Chain #", chain, " at ", 100 * iter / par["n_iter"], "% with a phi of ", model.phi)
        end

        

    end
    
    return model_hist
end