function Plot_model(
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}}, 
    model_mean::Array{Float64,2},  #used to be Array{Float64} since model_mean is a 1 dimensional vector, but it is now a 2D matrix
    TD_parameters::Dict{String,Any},
    cmax::Float64,
    l0::Int64,
    ModelType::String,
    CrossSection::String
    )

    ENV["GKSwstype"] = "nul"
    gr()
    closeenough = 2.0
    if ModelType == "Masked" 
        cmap = :jet
        cbtitle = "1000/Qp"
    elseif ModelType == "Uncertainty"
        cmap = :bone
        cbtitle = "σ"
    elseif ModelType == "Mean"
        cmap = :jet
        cbtitle = "1000/Qp"
    end

    dirname = "./figures/"*CrossSection*ModelType

    if TD_parameters["add_yVec"] == 0 # Tonga 2D
        println("****2D Plotting****")
        tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
        tmpy = vcat(zeros(size(dataStruct["dataX"]))', vec(dataStruct["edep"])') 
        tmpy = Array{Float64,2}(tmpy)
        p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model_mean), xlabel="distance(km)", ylabel="depth(km)", yflip=true, c=cmap)

        title!("model_mean_2d_Tonga")
        savefig(p, "model_mean_2d_Tonga")

        title!("model_mean_2d_Tonga_0-300km")
        ylims!(p, (0, 300))
        savefig(p, "model_mean_2d_Tonga_0-300km") 

        ylims!(p, (0, 660))
        plot!(p, tmpx, tmpy, color=:white, legend=false)

        scatter!(p, dataStruct["dataX"], zeros(size(dataStruct["dataX"])), marker=:utri, color=:pink, label="station", markersize=6)
        scatter!(p, dataStruct["elonsX"], vec(dataStruct["edep"]), marker=:o, color=:lightblue, label="events", markersize=4)
        title!("model_mean_2d_Tonga_ray")       
        savefig(p, dirname*"/model_mean_2d_Tonga_ray")   

        title!("model_mean_2d_Tonga__0-300km_ray")
        ylims!(p, (0, 300))
        savefig(p, dirname*"/model_mean_2d_Tonga__0-300km_ray") 
    else # 3d
        if CrossSection == "xz"
            (m, n) = size(RayTraces["rayY"])
            nearrays = zeros(n)
            for i in 1:n
                yvec = RayTraces["rayY"][:,i] .- l0    
                if sum(abs.(yvec)) - abs(sum(yvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(yvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(zeros(size(dataStruct["dataX"]))', vec(dataStruct["edep"])') 
            tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]

            p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model_mean), linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", yflip=true, clims=(0,cmax), c=cmap, colorbar_title = cbtitle)
           
            if ModelType == "Uncertainty"
                title!("Model uncertainty on cross-section")
            else
                title!(ModelType*" model on cross-section")
            end
            savefig(p, dirname*"/Model_"*ModelType*"_xzMap" * string(l0) * "km")

            title!(ModelType*" model on cross-section")
            # ylims!(p, (0, 300))
            # savefig(p, "model_mean_xzMap_0-300km" * string(TD_parameters["y0"][1]) * "km") 

            ylims!(p, (0, 660))
            scatter!(p, dataStruct["elonsX"], vec(dataStruct["edep"]), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
            savefig(p, dirname*"/Model_"*ModelType*"_xzMap_events" * string(l0) * "km") 

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct["dataX"], zeros(size(dataStruct["dataX"])), marker=:utri, color=:pink, label="station", markersize=6)
            
            title!("Model_"*ModelType*"_xzMap_ray" * string(l0) * "km")       
            savefig(p, dirname*"/Model_"*ModelType*"_xzMap_ray" * string(l0) * "km")   
  
            # title!("model_mean_xzMap_0-300km_ray" * string(TD_parameters["y0"][1]) * "km")
            # ylims!(p, (0, 300))
            # savefig(p, "model_mean_xzMap_0-300km_ray" * string(TD_parameters["y0"][1]) * "km") 
        end

        if CrossSection == "xy"
            (m, n) = size(RayTraces["rayZ"])
            nearrays = zeros(n)
            for i in 1:n
                zvec = RayTraces["rayZ"][:,i] .- l0    
                if sum(abs.(zvec)) - abs(sum(zvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(zvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(vec(dataStruct["dataY"])', vec(dataStruct["elatsY"])') 
            # tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]
        
            p = contourf(dataStruct["xVec"], dataStruct["yVec"], vec(model_mean), xlabel="X(km)", ylabel="Y(km)",linewidth=0.001, clims=(0,cmax), c=cmap, colorbar_title = cbtitle)
            title!("Model_"*ModelType*"_xyMap" * string(l0) * "km")  
            if ModelType == "Uncertainty"
                title!("Model uncertainty in map view (" * string(l0) * " km)")
            else
                title!(ModelType*" model in map view (" * string(l0) * " km)")
            end      
            savefig(p, dirname*"/Model_"*ModelType*"_xyMap" * string(l0) * "km")

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct["dataX"], dataStruct["dataY"], shape=:utri, color=:pink, label="station", markersize=6)
            scatter!(p, dataStruct["elonsX"], dataStruct["elatsY"], shape=:o, color=:lightblue, label="events", markersize=4)
            title!("Model_"*ModelType*"_xyMap_ray" * string(l0) * "km")
        
            savefig(p, dirname*"/Model_"*ModelType*"_xyMap_ray" * string(l0) * "km")       
        end

    end


end

function Plot_model_with_uncertainty(
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}}, 
    model_mean::Array{Float64,2},  #used to be Array{Float64} since model_mean is a 1 dimensional vector, but it is now a 2D matrix
    model_poststd::Array{Float64,2},
    TD_parameters::Dict{String,Any},
    cmax::Float64,
    l0::Int64,
    CrossSection::String
    )

    ENV["GKSwstype"] = "nul"
    gr()
    closeenough = 20.0
    ModelType = "Contour" 
    cmap = :jet
    cbtitle = "1000/Qp"
    threshold = 5.0
    dirname = "./figures/"*CrossSection*ModelType

    if TD_parameters["add_yVec"] == 0 # Tonga 2D
        println("****2D Plotting****")
        tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
        tmpy = vcat(zeros(size(dataStruct["dataX"]))', vec(dataStruct["edep"])') 
        tmpy = Array{Float64,2}(tmpy)

        p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model_mean), xlabel="distance(km)", ylabel="depth(km)", yflip=true, c=cmap)

        title!("model_mean_2d_Tonga")
        savefig(p, "model_mean_2d_Tonga")

        title!("model_mean_2d_Tonga_0-300km")
        ylims!(p, (0, 300))
        savefig(p, "model_mean_2d_Tonga_0-300km") 

        ylims!(p, (0, 660))
        plot!(p, tmpx, tmpy, color=:white, legend=false)

        scatter!(p, dataStruct["dataX"], zeros(size(dataStruct["dataX"])), marker=:utri, color=:pink, label="station", markersize=6)
        scatter!(p, dataStruct["elonsX"], vec(dataStruct["edep"]), marker=:o, color=:lightblue, label="events", markersize=4)
        title!("model_mean_2d_Tonga_ray")       
        savefig(p, dirname*"/model_mean_2d_Tonga_ray")   

        title!("model_mean_2d_Tonga__0-300km_ray")
        ylims!(p, (0, 300))
        savefig(p, dirname*"/model_mean_2d_Tonga__0-300km_ray") 
    else # 3d
        if CrossSection == "xz"
            (m, n) = size(RayTraces["rayY"])
            nearrays = zeros(n)
            for i in 1:n
                yvec = RayTraces["rayY"][:,i] .- l0    
                if sum(abs.(yvec)) - abs(sum(yvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(yvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(zeros(size(dataStruct["dataX"]))', vec(dataStruct["edep"])') 
            tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]
            
            alpha_mask = [if zi >= threshold zi else 1 end for zi in model_poststd]
            alpha_mask = normalize(alpha_mask, 0, 1)
            # alpha_mask = zeros(size(alpha_mask))

            p = contourf(
                vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model_mean), 
                linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", 
                yflip=true, clims=(0,cmax), c=cmap, alpha = alpha_mask, colorbar_title = cbtitle)
            
            contour!(
                p, vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(model_poststd), 
                levels = 20, linewidth=1, cbar = false, contour_labels = true)

            title!(ModelType*" model on cross-section")

            savefig(p, dirname*"/Model_"*ModelType*"_xzMap" * string(l0) * "km")

            title!(ModelType*" model on cross-section")
            # ylims!(p, (0, 300))
            # savefig(p, "model_mean_xzMap_0-300km" * string(TD_parameters["y0"][1]) * "km") 

            ylims!(p, (0, 660))
            scatter!(p, dataStruct["elonsX"], vec(dataStruct["edep"]), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
            savefig(p, dirname*"/Model_"*ModelType*"_xzMap_events" * string(l0) * "km") 

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct["dataX"], zeros(size(dataStruct["dataX"])), marker=:utri, color=:pink, label="station", markersize=6)
            
            title!("Model_"*ModelType*"_xzMap_ray" * string(l0) * "km")       
            savefig(p, dirname*"/Model_"*ModelType*"_xzMap_ray" * string(l0) * "km")   
  
            # title!("model_mean_xzMap_0-300km_ray" * string(TD_parameters["y0"][1]) * "km")
            # ylims!(p, (0, 300))
            # savefig(p, "model_mean_xzMap_0-300km_ray" * string(TD_parameters["y0"][1]) * "km") 
        end

        if CrossSection == "xy"
            (m, n) = size(RayTraces["rayZ"])
            nearrays = zeros(n)
            for i in 1:n
                zvec = RayTraces["rayZ"][:,i] .- l0    
                if sum(abs.(zvec)) - abs(sum(zvec)) > 1e-7
                    nearrays[i] = 1
                end  
                if min((abs.(zvec) .- closeenough)...) < 1e-7
                    nearrays[i] = 1
                end
            end
            nearrays = Array{Bool,1}(nearrays)
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(vec(dataStruct["dataY"])', vec(dataStruct["elatsY"])') 
            # tmpy = Array{Float64,2}(tmpy)
            nearraysX = tmpx[:,nearrays]
            nearraysY = tmpy[:,nearrays]
        
            alpha_mask = [if zi >= threshold zi else 1 end for zi in model_poststd]
            alpha_mask = normalize(alpha_mask, 0, 1)

            p = contourf(
                dataStruct["xVec"], dataStruct["yVec"], vec(model_mean), 
                xlabel="X(km)", ylabel="Y(km)",linewidth=0.001, 
                clims=(0,cmax), c=cmap, alpha = alpha_mask, colorbar_title = cbtitle)
            contour!(
                p, vec(dataStruct["xVec"]), vec(dataStruct["yVec"]), vec(model_poststd), 
                levels = 10, linewidth=1, cbar = false, contour_labels = true)

            title!("Model_"*ModelType*"_xyMap" * string(l0) * "km")  
            if ModelType == "Uncertainty"
                title!("Model uncertainty in map view (" * string(l0) * " km)")
            else
                title!(ModelType*" model in map view (" * string(l0) * " km)")
            end      
            savefig(p, dirname*"/Model_"*ModelType*"_xyMap" * string(l0) * "km")

            plot!(p, tmpx, tmpy, color=:white, legend=false)
            plot!(p, nearraysX, nearraysY, color=:forestgreen)
            scatter!(p, dataStruct["dataX"], dataStruct["dataY"], shape=:utri, color=:pink, label="station", markersize=6)
            scatter!(p, dataStruct["elonsX"], dataStruct["elatsY"], shape=:o, color=:lightblue, label="events", markersize=4)
            title!("Model_"*ModelType*"_xyMap_ray" * string(l0) * "km")
        
            savefig(p, dirname*"/Model_"*ModelType*"_xyMap_ray" * string(l0) * "km")       
        end

    end


end

function plot_model_hist(
    model_hist, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}},
    TD_parameters::Dict{String,Any}, 
    cmax::Float64
    )

    if TD_parameters["xzMap"] == true
        for l0 in TD_parameters["y0"]
            m_xz = []
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["zVec"])))

            mcount = 0

            for i = 1:length(model_hist)

                for j = 1:length(model_hist[i])

                    m  = [ cart_v_nearest(xs, l0, zs,
                        model_hist[i][j].xCell, model_hist[i][j].yCell, model_hist[i][j].zCell, model_hist[i][j].zeta)
                        for xs in vec(dataStruct["xVec"]), zs in vec(dataStruct["zVec"]) ]
                    append!(m_xz,[m])

                end
            end

            model_mean_xz   = mean(m_xz)
            poststd_xz      = std(m_xz)
            mask_xz         = ones(size(poststd_xz))
            for i = 1:length(poststd_xz)
                if poststd_xz[i] > 5
                    mask_xz[i] = NaN
                end
            end
            mask_model_xz = mask_xz .* model_mean_xz

            # @time Plot_model(dataStruct, RayTraces, model_mean_xz, TD_parameters, cmax, l0, "Mean", "xz")
            # @time Plot_model(dataStruct, RayTraces, poststd_xz, TD_parameters, maximum(poststd_xz), l0, "Uncertainty", "xz")
            # @time Plot_model(dataStruct, RayTraces, mask_model_xz, TD_parameters, cmax, l0, "Masked", "xz")
            @time Plot_model_with_uncertainty(
                dataStruct, RayTraces, model_mean_xz, poststd_xz, 
                TD_parameters, cmax, l0,  "xz" 
                )
        end
    end

    if TD_parameters["xyMap"] == true
        for l0 in TD_parameters["z0"]
            m_xy = []
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
            mcount = 0

            for i = 1:length(model_hist)

                for j = 1:length(model_hist[i])

                    m  = [ cart_v_nearest(xs, ys, l0,
                        model_hist[i][j].xCell, model_hist[i][j].yCell, model_hist[i][j].zCell, model_hist[i][j].zeta)
                        for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]
                    append!(m_xy,[m])

                end
            end

            model_mean_xy   = mean(m_xy)
            poststd_xy      = std(m_xy)
            mask_xy         = ones(size(poststd_xy))
            for i = 1:length(poststd_xy)
                if poststd_xy[i] > 5
                    mask_xy[i] = NaN
                end
            end
            mask_model_xy = mask_xy .* model_mean_xy

            # @time Plot_model(dataStruct, RayTraces, model_mean_xy, TD_parameters, cmax, l0, "Mean", "xy")
            # @time Plot_model(dataStruct, RayTraces, poststd_xy, TD_parameters, maximum(poststd_xy), l0, "Uncertainty", "xy")
            # @time Plot_model(dataStruct, RayTraces, mask_model_xy, TD_parameters, cmax, l0, "Masked", "xy")
            @time Plot_model_with_uncertainty(
                dataStruct, RayTraces, model_mean_xy, poststd_xy, 
                TD_parameters, cmax, l0, "xy" 
                )
        end
    end

end


function plot_model_hist_weighted(
    model_hist, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    RayTraces::Dict{String,Array{Float64,2}},
    TD_parameters::Dict{String,Any}, 
    cmax::Float64
    )

    if TD_parameters["xzMap"] == true
        for l0 in TD_parameters["y0"]
            m_xz = []
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["zVec"])))

            mcount = 0

            llh = []
            for i = 1:length(model_hist)
                
                for j = 1:length(model_hist[i])

                    m  = [ cart_v_nearest(xs, l0, zs,
                        model_hist[i][j].xCell, model_hist[i][j].yCell, model_hist[i][j].zCell, model_hist[i][j].zeta)
                        for xs in vec(dataStruct["xVec"]), zs in vec(dataStruct["zVec"]) ]
                    append!(m_xz,[m])

                    push!(llh, model_hist[i][j].likelihood)
                
                end

            end

            weight_llh      = llh / sum(llh)
            # yurong 03/13/23 could still try wsum in package StatsBase to calculate the weighted average
            # and try weighted std
            # weight_llh      = [ones(size(m_xz[i]))*weight_llh[i] for i in 1:length(weight_llh)]
            # model_mean_xz   = wsum(m_xz, weight_llh)
            model_mean_xz   = zeros(size(m_xz[1]))
            for i in 1:length(m_xz)
                model_mean_xz += m_xz[i]*weight_llh[i]
            end
            poststd_xz      = std(m_xz)
            mask_xz         = ones(size(poststd_xz))
            for i = 1:length(poststd_xz)
                if poststd_xz[i] > 5
                    mask_xz[i] = NaN
                end
            end
            mask_model_xz = mask_xz .* model_mean_xz

            @time Plot_model(dataStruct, RayTraces, model_mean_xz, TD_parameters, cmax, l0, "Mean", "xz")
            @time Plot_model(dataStruct, RayTraces, poststd_xz, TD_parameters, maximum(poststd_xz), l0, "Uncertainty", "xz")
            @time Plot_model(dataStruct, RayTraces, mask_model_xz, TD_parameters, cmax, l0, "Masked", "xz")
        end
    end

    if TD_parameters["xyMap"] == true
        
        for l0 in TD_parameters["z0"]
            m_xy = []
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
            mcount = 0

            llh = []
            for i = 1:length(model_hist)

                for j = 1:length(model_hist[i])

                    m  = [ cart_v_nearest(xs, ys, l0,
                        model_hist[i][j].xCell, model_hist[i][j].yCell, model_hist[i][j].zCell, model_hist[i][j].zeta)
                        for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]
                    append!(m_xy,[m])

                    push!(llh, model_hist[i][j].likelihood)

                end

            end

            weight_llh      = llh / sum(llh)
            model_mean_xy   = zeros(size(m_xy[1]))
            for i in 1:length(m_xy)
                model_mean_xy += m_xy[i]*weight_llh[i]
            end
            poststd_xy      = std(m_xy)
            mask_xy         = ones(size(poststd_xy))
            for i = 1:length(poststd_xy)
                if poststd_xy[i] > 5
                    mask_xy[i] = NaN
                end
            end
            mask_model_xy = mask_xy .* model_mean_xy

            @time Plot_model(dataStruct, RayTraces, model_mean_xy, TD_parameters, cmax, l0, "Mean", "xy")
            @time Plot_model(dataStruct, RayTraces, poststd_xy, TD_parameters, maximum(poststd_xy), l0, "Uncertainty", "xy")
            @time Plot_model(dataStruct, RayTraces, mask_model_xy, TD_parameters, cmax, l0, "Masked", "xy")
        end
    end

end

function plot_voronoi(
    model::Model, 
    dataStruct::Dict{String,AbstractArray{Float64,N} where N}, 
    TD_parameters::Dict{String,Any}, 
    chain::Int64, 
    iter::Float64
    )
    ENV["GKSwstype"] = "nul"

    cmap    = :jet
    cbtitle = "1000/Qp"
    cmax    = TD_parameters["cmax"]

    if TD_parameters["xzMap"] == true
        for l0 in TD_parameters["y0"]
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["zVec"])))
            m  = [ cart_v_nearest(xs, l0, zs,
                model.xCell, model.yCell, model.zCell, model.zeta)
                for xs in vec(dataStruct["xVec"]), zs in vec(dataStruct["zVec"]) ]

            p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["zVec"]), vec(m), 
                 linewidth=0.001, xlabel="distance(km)", ylabel="depth(km)", yflip=true, 
                 clims=(0,cmax), c=cmap, colorbar_title = cbtitle)
            title!("Cross section Y="*string(l0)*";Chain"*string(chain)*"_"*string(Int(iter))*
            ";llh="*string(round(model.likelihood,digits=2)))
            savefig(p, "./figures/xzVoronoi_"*string(l0)*"/Chain"*string(chain)*"_"*string(Int(iter)))

            # add rays, EQs, and stations
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(zeros(size(dataStruct["dataX"]))', vec(dataStruct["edep"])') 
            tmpy = Array{Float64,2}(tmpy)
            plot!(p, tmpx, tmpy, color=:white, legend=false)
            scatter!(p, dataStruct["elonsX"], vec(dataStruct["edep"]), marker=:o, color=:lightblue, label="events", markersize=4,legend=false)
            scatter!(p, dataStruct["dataX"], zeros(size(dataStruct["dataX"])), marker=:utri, color=:pink, label="station", markersize=6)
            title!("Cross section Y="*string(l0)*";Chain"*string(chain)*"_"*string(Int(iter))*
            ";llh="*string(round(model.likelihood,digits=2)))
            savefig(p, "./figures/xzVoronoi_"*string(l0)*"/Chain"*string(chain)*"_"*string(Int(iter))*"ray")

        end
    end

    if TD_parameters["xyMap"] == true
        for l0 in TD_parameters["z0"]
            m = zeros(length(vec(dataStruct["xVec"])), length(vec(dataStruct["yVec"])))
            m  = [ cart_v_nearest(xs, ys, l0,
                model.xCell, model.yCell, model.zCell, model.zeta)
                for xs in vec(dataStruct["xVec"]), ys in vec(dataStruct["yVec"]) ]

            p = contourf(vec(dataStruct["xVec"]), vec(dataStruct["yVec"]), vec(m),
                 xlabel="X(km)", ylabel="Y(km)",linewidth=0.001, clims=(0,cmax), 
                 c=cmap, colorbar_title = cbtitle)
            title!("Map View Z="*string(l0)*";Chain"*string(chain)*"_"*string(Int(iter))*
            ";llh="*string(round(model.likelihood,digits=2)))
            savefig(p, "./figures/xyVoronoi_"*string(l0)*"/Chain"*string(chain)*"_"*string(Int(iter)))

            # add rays, EQs, and stations
            tmpx = vcat(vec(dataStruct["dataX"])', vec(dataStruct["elonsX"])')
            tmpy = vcat(vec(dataStruct["dataY"])', vec(dataStruct["elatsY"])')
            plot!(p, tmpx, tmpy, color=:white, legend=false)
            scatter!(p, dataStruct["dataX"], dataStruct["dataY"], shape=:utri, color=:pink, label="station", markersize=6)
            scatter!(p, dataStruct["elonsX"], dataStruct["elatsY"], shape=:o, color=:lightblue, label="events", markersize=4)
            title!("Cross section Y="*string(l0)*";Chain"*string(chain)*"_"*string(Int(iter))*
            ";llh="*string(round(model.likelihood,digits=2)))
            savefig(p, "./figures/xyVoronoi_"*string(l0)*"/Chain"*string(chain)*"_"*string(Int(iter))*"ray")
        end
    end

end

function plot_convergence(TD_parameters::Dict{String,Any})
# plot the convergence of nCells and phi value over iterations for each chain
# saved in ./figures/nCells and ./figures/phi
    ENV["GKSwstype"] = "nul"

    for chain in 1:TD_parameters["n_chains"]
        model_checkpoint_lists = glob("./models/chain" * string(chain) * "_*")
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
            # plot nCells and phi over iterations
            cellnumber_list = load_model["nCells"]
            phi_list        = load_model["phi"]

            p1 = plot(1:length(cellnumber_list),cellnumber_list)
            title!(p1,"nCells Convergence in Chain" *string(chain))
            xlabel!(p1,"Iterations")
            ylabel!(p1,"Number of Cells")
            p2 = plot(1:length(phi_list),phi_list)
            title!(p2,"Phi Convergence in Chain" *string(chain))
            xlabel!(p2,"Iterations")
            ylabel!(p2,"Phi")
            savefig(p1,"./figures/nCells/nCells_chain" * string(chain))
            savefig(p2,"./figures/phi/phi_chain" * string(chain))
 
        end
    end
end

function plot_llh(models)
    ENV["GKSwstype"] = "nul"
    
    llh_all = []
    for i in 1:length(models)
        llh = []
        for j in 1:length(models[i])
            push!(llh, models[i][j].likelihood)
        end
        push!(llh_all, llh)
        p = histogram(llh,bins=10,xlabel="likelihood", ylabel="normalized ratio",normalize=true)
        title!(p,"llh_chain"*string(i))
        savefig(p,"./figures/llh/chain"*string(i)*".png")
    end
    p = histogram(llh_all,xlabel="likelihood", ylabel="normalized ratio",normalize=true)

    title!(p,"llh_all")
    savefig(p,"./figures/llh/llh_all.png")

end
