using JLD, HDF5
using Distributed 
using Glob,YAML
@everywhere include("./scripts/model.jl")
@everywhere include("./scripts/utils.jl")
@everywhere include("./scripts/stats.jl")
@everywhere include("./scripts/plot_sub.jl")
@everywhere include("./scripts/interp.jl")
@everywhere include("./scripts/load.jl")
@everywhere include("./scripts/inversion_function.jl")

file_name = "par_cart_test.yml"
@time par = load_par_from_yml(file_name)
@time (dataStruct, RayTraces) = load_data_Tonga(par)
make_dir(par, file_name)
println("--------Data Loaded-------")

@time models = perform_inversion(par, dataStruct, RayTraces)
println("--------Inversion End-------")

# println("--------Plotting Models-------")
# @time plot_model_hist(models, dataStruct, par, 30.0)
# @time plot_convergence(par)
# @time save("model.jld","model",models)
# model_checkpoint_lists = glob("chain*jld")
# [rm(model_checkpoint_lists[i]) for i in 1:length(model_checkpoint_lists)]

