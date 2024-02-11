include("prepare_data.jl")
include("main_inversion.jl")
include("merge_models.jl")

file_name = "plot_sph_test.yml" 
@timev prepare_data(file_name)
@timev main_inversion(file_name)
@timev merge_models(file_name)
