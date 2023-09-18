"""

################################################################################                                     #
################################################################################
# Based on StochasticPowerModels.jl                                            #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# This example is for Numerical Illustration 
################################################################################

# variables
"""

using Pkg
Pkg.activate(".")
Pkg.instantiate()
using JuMP
using Ipopt
using PowerModels
using StochasticPowerModels
using PowerModelsDistribution
using JSON
using DataFrames
using CSV
using Statistics
using Plots
# constants 
const PM = PowerModels
const SPM = StochasticPowerModels
const PMD = PowerModelsDistribution

# solvers
ipopt_solver = Ipopt.Optimizer
# input
deg  = 2
aux  = true
red  = false

feeder ="Pola/1076069_1274129_mod_configuration.json" #feeder with 7 consumer for test
# feeder= "All_feeder/65019_74478_configuration.json"
# data
file  = joinpath(BASE_DIR, "test/data/Spanish/")
load_file= "beta_lm_2016_8_6_60min.csv"
pv_file = "Normal_pm_2021.csv"
# plot!(xticks=([j for j=4:4:40], [7, 8, 9, 10,11,12,13,14,15,16]))
"""
# Code to plot the result of Fig. 5 without the scenario based method #
"""
k=Dict()
r=Dict()
for i=7:17
    data  = SPM.build_mathematical_model_mc(file, feeder,load_file, pv_file, t_s=i)
    [data["bus"]["$i"]["vmax"]= [1.05, 1.05, 1.05] for i=1:length(data["bus"])]
    [data["bus"]["$i"]["vmin"]= [0.95, 0.95, 0.95] for i=1:length(data["bus"])]
    [data["pv"]["$i"]["p_size"]= 3 for i in [6]]   
    [data["pv"]["$i"]["p_size"]= 6 for i in [1, 3]]
    [data["pv"]["$i"]["p_size"]= 10 for i in [2, 4, 7]]
    [data["pv"]["$i"]["p_size"]= 12 for i in [5]]    
    sdata = SPM.build_stochastic_data_mc(data,deg)
    r["$i"]= SPM.solve_sopf_iv_mc_doe(data, PMD.IVRUPowerModel, ipopt_solver)
    #, aux=aux, deg=deg, red=red, stochastic=false)
end
df=[r["$j"]["solution"]["nw"]["1"]["pv"]["$k"]["p_c"] for j=7:17, k=1:7]
df1=DataFrame(df, ["D$(i)" for i=1:7])
fig=heatmap(df,xlabel="Device", ylabel="hrs",xticks=([i for i=1:7], ["D1", "D2", "D3", "D4", "D5", "D6", "D7"]), yticks=([j for j=1:11],[j+6 for j=1:11]), clims=(0,1))
# savefig(fig, "min curtailment.pdf")
# CSV.write("min curtailment.csv", df1)
# savefig(fig, "Qos relative_0o0_new.pdf")
# CSV.write("Qos relative_0o0_new.csv", df1)
# savefig(fig, "Equality relative_0o6_v02.pdf")
# CSV.write("Equality relative_0o6_v02.csv", df1)
savefig(fig, "alpha relative_0o0_v04.pdf")
CSV.write("alpha relative_0o0_v04.csv", df1)