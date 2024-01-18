include("prepare_data.jl")
include("main_inversion.jl")
# include("Plot_MapView.jl")

file_name = "par_cart_1e4.yml" 
@timev prepare_data(file_name)
@timev main_inversion(file_name)
# @timev plot_mapview(file_name)