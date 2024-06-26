"""

################################################################################                                     #
################################################################################
# Based on StochasticPowerModels.jl                                            #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# This example is for Numerical Illustration  in 
#Koirala, Arpan; Geth, Frederik, Van Aacker, Tom (2023): 
#Determining Dynamic Operating Envelopes Using Stochastic Unbalanced Optimal Power Flow
#Submitted to PSCC 2024
################################################################################

# variables
"""

using Pkg
Pkg.activate(".")
# Pkg.instantiate()
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
feeder ="Pola/1076069_1274129_bal_configuration.json" #feeder with 7 consumer for test
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
for i=6:19
    data  = SPM.build_mathematical_model_mc(file, feeder,load_file, pv_file, t_s=i)
    [data["bus"]["$i"]["vmax"]= [1.05, 1.05, 1.05] for i=1:length(data["bus"])]
    [data["bus"]["$i"]["vmin"]= [0.95, 0.95, 0.95] for i=1:length(data["bus"])]
    [data["bus"]["$i"]["λvmax"]= 1.03643 for i=1:length(data["bus"])]
    [data["pv"]["$i"]["p_size"]= 2 for i in [6]]   
    [data["pv"]["$i"]["p_size"]= 6 for i in [1, 3]]
    [data["pv"]["$i"]["p_size"]= 8 for i in [2, 4, 7]]
    [data["pv"]["$i"]["p_size"]= 15 for i in [5]]  
    # [data["pv"]["$i"]["p_size"]= 15 for i in [1,2,3,4,5,6,7]]  
    sdata = SPM.build_stochastic_data_mc(data,deg)
    r["$i"]= SPM.solve_sopf_iv_mc_doe(data, PMD.IVRUPowerModel, ipopt_solver)
    #, aux=aux, deg=deg, red=red, stochastic=false)
end
df=[r["$j"]["solution"]["nw"]["1"]["pv"]["$k"]["p_c"] for j=6:19, k=1:7]
df1=DataFrame(df, ["D$(i)" for i=1:7])

inj1=[r["$j"]["solution"]["nw"]["1"]["pv"]["$k"]["n_inj"] for j=6:19, k=1:7]
inj2=DataFrame(inj1, ["D$(i)" for i=1:7])
fig_i=heatmap(inj1,xlabel="Device", ylabel="hrs",xticks=([i for i=1:7], ["D1", "D2", "D3", "D4", "D5", "D6", "D7"]), guidefontsize=11,tickfontsize=11, yticks=([j for j=1:134],[j+5 for j=1:14]), clims=(0,1), size=(300,350))

result_dir= "c:\\Users\\aka\\OneDrive - Ampner Oy\\Documents\\collaboration\\DOE\\voltage_bound\\"
fig=heatmap(df,xlabel="Device", ylabel="hrs",xticks=([i for i=1:7], ["D1", "D2", "D3", "D4", "D5", "D6", "D7"]), guidefontsize=11,tickfontsize=11,yticks=([j for j=1:14],[j+5 for j=1:14]), clims=(0,1), size=(300,350))

"""Codes for saving figures"""
savefig(fig, result_dir*"maxinj_curt_10.pdf")
CSV.write(result_dir*"maxinj_curt_bal_10.csv", df1)
savefig(fig_i,result_dir*"max inj_bal_10.pdf")
CSV.write(result_dir*"max inj_bal_10.csv", inj2)
# savefig(fig, result_dir*"Qos_curt.pdf")
# CSV.write(result_dir*"Qos_curt.csv", df1)
# savefig(fig_i, result_dir*"QoS_inj.pdf")
# CSV.write(result_dir*"QoS_inj.csv", inj2)
# savefig(fig, result_dir*"Equality_curt_3.pdf")
# CSV.write(result_dir*"Equality_curt_3.csv", df1)
# savefig(fig_i, result_dir*"Equality_inj_3.pdf")
# CSV.write(result_dir*"Equality_inj_3.csv", inj2)
# savefig(fig, result_dir*"alpha relative_10o0_curt.pdf")
# CSV.write(result_dir*"alpha relative_10o0_curt.csv", df1)
# savefig(fig_i, result_dir*"alpha relative_10o0_inj.pdf")
# CSV.write(result_dir*"alpha relative_10o0_inj.csv", inj2)