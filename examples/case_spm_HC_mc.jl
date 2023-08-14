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
feeder= "All_feeder/65019_74478_configuration.json"
# data
file  = joinpath(BASE_DIR, "test/data/Spanish/")
load_file= "beta_lm_2016_8_6.csv"
pv_file = "beta_pm_2016_8_6.csv"

"""
# Code to plot the result of Fig. 5 without the scenario based method #
"""
r=Dict()
for i=61:70
    data  = SPM.build_mathematical_model_mc(file, feeder,load_file, pv_file, t_s=i)
    [data["bus"]["$i"]["vmax"]= [1.05, 1.05, 1.05] for i=1:length(data["bus"])]
    [data["bus"]["$i"]["vmin"]= [0.95, 0.95, 0.95] for i=1:length(data["bus"])]    
    sdata = SPM.build_stochastic_data_mc(data,deg)
    r["$i"]= SPM.solve_sopf_iv_mc_hc(data, PMD.IVRUPowerModel, ipopt_solver)#, aux=aux, deg=deg, red=red, stochastic=false)
    s2 = Dict("output" => Dict("duals" => true))
end
# result_hc_2= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)
# scatter([1,2,3,4,5,6,7],[result_hc_2["solution"]["nw"]["1"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="gPC-CC-OPF HC", figsize=(28,8))
# scatter!([1,2,3,4,5,6,7],[result_hc["solution"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="OPF HC", figsize=(28,8))
# plot!(xlabel="Device Id")
# plot!(ylabel="PV size [kWp]")
# display(plot!(title="Fig. 5: deterministic vs stochastic OPF"))

SPM.print_summary_mc(result_hc["solution"])

# print variables for a specific index k
k=1
SPM.print_summary_mc(result_hc["solution"]["nw"]["$k"])

# get polynomial chaos coefficients for specific component
pg_coeff = pce_coeff(result_hc, "gen", 1, "pg") 

# obtain 10 samples of the generator active power output variable
pg_sample = sample(result_hc, "gen", 1, "pg"; sample_size=10) 

v_sample_1 = sqrt.(sample(r["57"], "bus", 7, "vms", phase=1; sample_size=1000))
v_sample_2 = sqrt.(sample(r["57"], "bus", 7, "vms", phase=2; sample_size=1000)) 
v_sample_3 = sqrt.(sample(r["57"], "bus", 7, "vms", phase=3; sample_size=1000))

histogram((v_sample_1))
histogram!((v_sample_2))
histogram!((v_sample_3))

# p_sample_1 = sample(result_hc, "branch", 6, "cmss", phase=1; sample_size=1000)
# p_sample_2 = sample(result_hc, "branch", 6, "cmss", phase=2; sample_size=1000) 
# p_sample_3 = sample(result_hc, "branch", 6, "cmss", phase=3; sample_size=1000)

# histogram(p_sample_3)
# histogram!(p_sample_1)
# histogram!(p_sample_2)

p_sample_1 = sample(result_hc, "gen", 1, "pg", phase=1; sample_size=1000)
p_sample_2 = sample(result_hc, "gen", 1, "pg", phase=2; sample_size=1000) 
p_sample_3 = sample(result_hc, "gen", 1, "pg", phase=3; sample_size=1000)

histogram(p_sample_3)
histogram!(p_sample_1)
histogram!(p_sample_2)



