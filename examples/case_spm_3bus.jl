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

feeder =  "case3_unbalanced.dss"
file = joinpath(BASE_DIR, "test/data/opendss/"*feeder)
data_eng = PMD.parse_file(file)
data_eng["voltage_source"]["source"]["rs"]*=0
data_eng["voltage_source"]["source"]["xs"]*=0
data_eng["load"]["l1"]["pd_nom"]=0.5*data_eng["load"]["l1"]["pd_nom"]
org_base=data_eng["settings"]["sbase_default"]
data_eng["settings"]["sbase_default"]=1
ratio=data_eng["settings"]["sbase_default"]/org_base

data=PMD.transform_data_model(data_eng)

[data["bus"]["$i"]["vmax"]= [1.1, 1.1, 1.1] for i=1:length(data["bus"])]
[data["bus"]["$i"]["vmin"]= [0.9, 0.9, 0.9] for i=1:length(data["bus"])]
data["bus"]["2"]["vmax"]= [1.0, 1.0, 1.0] 
data["bus"]["2"]["vmin"]= [1.0, 1.0, 1.0] 
[data["bus"]["$i"]["λvmin"]= 1.65 for i=1:length(data["bus"])]
[data["bus"]["$i"]["λvmax"]= 1.65 for i=1:length(data["bus"])]
[data["branch"]["$i"]["λcmax"]= 1.65 for i=1:length(data["branch"])]
[data["gen"]["$i"]["pgmin"]= [-10, -10, -10]./ratio for i=1:length(data["gen"])]
[data["gen"]["$i"]["pgmax"]= [20, 20, 200]./ratio for i=1:length(data["gen"])]
[data["gen"]["$i"]["qgmin"]= [-100, -100, -100]./ratio for i=1:length(data["gen"])]
[data["gen"]["$i"]["qgmax"]= [100, 100, 100]./ratio for i=1:length(data["gen"])]
[data["gen"]["$i"]["λpmax"]= 1.65 for i=1:length(data["gen"])]
[data["gen"]["$i"]["λpmin"]= 1.65 for i=1:length(data["gen"])]
[data["gen"]["$i"]["λqmax"]= 1.65 for i=1:length(data["gen"])]
[data["gen"]["$i"]["λqmin"]= 1.65 for i=1:length(data["gen"])]



data["gen"]["1"]["cost"]= data["gen"]["1"]["cost"]
sdata = SPM.build_stochastic_data_mc_dss(data,deg)

result_hc= PMD.solve_mc_model(sdata, IVRUPowerModel, ipopt_solver, SPM.build_sopf_iv_mc; multinetwork=true, solution_processors=[SPM.sol_data_model!])
result_hc["mop"] = sdata["mop"]

v_sample_1 = sample(result_hc, "bus", 3, "vms", phase=1; sample_size=1000)
v_sample_2 = sample(result_hc, "bus", 3, "vms", phase=2; sample_size=1000) 
v_sample_3 = sample(result_hc, "bus", 3, "vms", phase=3; sample_size=1000)

histogram((v_sample_1))
histogram!((v_sample_2))
histogram!((v_sample_3))


# result = PMD.solve_mc_opf(data, IVRUPowerModel, ipopt_solver)


# feeder ="Pola/1076069_1274129_mod_configuration.json" #feeder with 7 consumer for test
# # data
# file  = joinpath(BASE_DIR, "test/data/Spanish/")
# load_file= "beta_lm_2016_8_6.csv"
# pv_file = "beta_pm_2016_8_6.csv"

# """
# # Code to plot the result of Fig. 5 without the scenario based method #
# """
# data  = SPM.build_mathematical_model_mc(file, feeder,load_file, pv_file)
# sdata = SPM.build_stochastic_data_mc(data,deg)
# result_hc= solve_sopf_iv_mc(data, PMD.IVRUPowerModel, ipopt_solver)#, aux=aux, deg=deg, red=red, stochastic=false)
# s2 = Dict("output" => Dict("duals" => true))
# result_hc_2= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)
# scatter([1,2,3,4,5,6,7],[result_hc_2["solution"]["nw"]["1"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="gPC-CC-OPF HC", figsize=(28,8))
# scatter!([1,2,3,4,5,6,7],[result_hc["solution"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="OPF HC", figsize=(28,8))
# plot!(xlabel="Device Id")
# plot!(ylabel="PV size [kWp]")
# display(plot!(title="Fig. 5: deterministic vs stochastic OPF"))

# SPM.print_summary_mc(result_hc["solution"])

# # print variables for a specific index k
# k=1
# SPM.print_summary_mc(result_hc["solution"]["nw"]["$k"])

# # get polynomial chaos coefficients for specific component
# pg_coeff = pce_coeff(result_hc, "gen", 1, "pg") 

# # obtain 10 samples of the generator active power output variable
# pg_sample = sample(result_hc, "gen", 1, "pg"; sample_size=10) 



p_sample_1 = sample(result_hc, "branch", 2, "cmss", phase=1; sample_size=1000)
p_sample_2 = sample(result_hc, "branch", 2, "cmss", phase=2; sample_size=1000) 
p_sample_3 = sample(result_hc, "branch", 2, "cmss", phase=3; sample_size=1000)

histogram(p_sample_3)
histogram!(p_sample_1)
histogram!(p_sample_2)

p_sample_1 = sample(result_hc, "gen", 1, "pg", phase=1; sample_size=1000)
p_sample_2 = sample(result_hc, "gen", 1, "pg", phase=2; sample_size=1000) 
p_sample_3 = sample(result_hc, "gen", 1, "pg", phase=3; sample_size=1000)

histogram(p_sample_3)
histogram!(p_sample_1)
histogram!(p_sample_2)

# # obtain an kernel density estimate of the generator active power output variable
# # pg_density = _PC.density(result_ivr, "gen", 1, "pg"; sample_size=10) 

# # it = Iterators.product(ntuple(_ -> [1.6,2], 11)...)
# # p=collect(it);
# # d=[]
# # for x in p
# #     if sum(x)==20
# #         push!(d, [l for l in x])
# #     end
# # end

# # # e=DataFrame(d)
# # writedlm("C:\\Users\\karpan\\Documents\\GitHub_ak\\HostingCapacityLVDS\\combination2.csv", d, ',')
