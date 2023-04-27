"""
    StochasticPowerModels.extend_matlab_file(path::String)
"""
function extend_matlab_file(path::String)
    # data
    data = _PM.parse_file(path)

    # general data
    baseMVA = data["baseMVA"]

    # bus data 
    Nb   = length(data["bus"])
    μ, σ = zeros(Nb), zeros(Nb)
    for (l,load) in data["load"]
        bus = load["load_bus"]
        μ[bus] = load["pd"] * baseMVA
        σ[bus] = load["pd"] * baseMVA * 0.10
    end
    for (b,bus) in data["bus"]
        bus["dst_id"]   = 0
        bus["μ"]        = μ[parse(Int,b)]
        bus["σ"]        = σ[parse(Int,b)]
        bus["λvmin"]    = 1.03643
        bus["λvmax"]    = 1.03643
    end

    # generator data
    for (g,gen) in data["gen"]
        gen["λpmin"] = 1.03643
        gen["λpmax"] = 1.03643
        gen["λqmin"] = 1.03643
        gen["λqmax"] = 1.03643
    end

    # branch data
    for (l,branch) in data["branch"]
        branch["λcmax"] = 1.03643
    end

    # export file
    _PM.export_file(path[1:end-2] * "_spm.m", data)
end


"Converts JSON file of three phase DN to single phase equivalent"
function build_mathematical_model_single_phase(dir, config_file_name, load_dist_csv, pv_dist_csv; t_s=52, pd = 0.0, qd = 0.0, scale_factor = 1.0, curt=0.0, cross_area_fact=1.0)
#configuration = "star"

"""
Specify voltage and power base, voltage base should be the phase-to-ground voltage
of the feeder studied (kV), the power base can be arbitrairly chosen (MW)
"""
voltage_base = 0.230  # (kV)
power_base = 0.5  # (MW)
Z_base = voltage_base^2/power_base # (Ohm)
current_base = power_base/(voltage_base*1e-3) # (A)

mwpu = 1/power_base
kwpu = (1e-3)/power_base

network_model = Dict{String,Any}()
configuration_json_dict = Dict{Any,Any}()
device_df=CSV.read(dir*config_file_name[1:length(config_file_name)-19]*".csv", DataFrame)

dist_lv=CSV.read(dir*load_dist_csv, DataFrame)
dist_pv=CSV.read(dir*pv_dist_csv, DataFrame)
# dist_pv=CSV.read(dir*"beta_pm_2022_181"*".csv", DataFrame)
dist_pv_ts= dist_pv[in([t_s]).(dist_pv.timeslot),:]
dist_lv_ts=dist_lv[in([t_s]).(dist_lv.timeslot),:]

dist_lv_ts_feeder = dist_lv_ts[in(unique(device_df.category)).(dist_lv_ts.cluster),:]
s_dict=Dict()
i=1
for dist in eachrow(dist_lv_ts_feeder)
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = dist["cluster"]
    s["pa"]= dist["alpha"]
    s["pb"]= dist["beta"]
    s["pc"]= dist["lower"]
    s["pd"]= dist["lower"]+dist["upper"]
    s_dict[string(i)] = s
    i=i+1
end


##add Irradiance if day time or there is some Irradiance
if dist_pv_ts.upper[1]>0
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = 55
    s["pa"]= dist_pv_ts[!,"alpha"][1]
    s["pb"]= dist_pv_ts[!,"beta"][1]
    s["pc"]= dist_pv_ts[!,"lower"][1]
    s["pd"]= dist_pv_ts[!,"lower"][1]+dist_pv_ts[!,"upper"][1]
    s_dict[string(i)] = s
end


#network_model["is_kron_reduced"] = true
network_model["p_factor"] = 0.95
network_model["q_factor"] = sqrt(1-0.95^2)
network_model["dcline"] = Dict{String,Any}()
network_model["switch"] = Dict{String,Any}()
#network_model["is_projected"] = true
network_model["per_unit"] = true
#network_model["data_model"] = MATHEMATICAL
network_model["shunt"] = Dict{String,Any}()
network_model["transformer"] = Dict{String,Any}()
network_model["bus"] = Dict{String,Any}()
network_model["map"] = Dict{String,Any}()
#network_model["conductors"] = 1
network_model["baseMVA"] =  power_base
network_model["basekv"] =  voltage_base
network_model["bus_lookup"] = Dict{Any,Int64}()
network_model["run_type"] = 1
network_model["load"] = Dict{String,Any}()
network_model["gen"] = Dict{String,Any}("1" => Dict{String,Any}(
"pg"            => 0.2,
"model"         => 2,
#"connections"   => [1, 2, 3],
"shutdown"      => 0.0,
"startup"       => 0.0,
#"configuration" => WYE,
"name"          => "virtual_generator",
"qg"            => 0.0,
"gen_bus"       => 1,
"vbase"         =>  voltage_base,
"source_id"     => Any["gen",1],
"index"         => 1,
"cost"          => [20000.0, 1400.0, 0.0],
"gen_status"    => 1,
"qmax"          => 1.275,
"qmin"          => -1.275,
"pmax"          => 1.5,
"pmin"          => -1.5,
"ncost"         => 3,
"λpmin"         => 1.65, #1.03643 ,
"λpmax"         => 1.65, #1.03643 ,
"λqmin"         => 1.65, #1.03643 ,
"λqmax"         => 1.65 #1.03643
))
network_model["settings"] = Dict{String,Any}(
"sbase_default"        => power_base,
"vbases_default"       => Dict{String,Any}(), #No default is specified for now, since default is never used
"voltage_scale_factor" => 1E3, #Voltages are thus expressed in kV
"sbase"                => power_base,
"power_scale_factor"   => 1E6, #Power is expressed in MW
"base_frequency"       => 50.0 #Hertz
)
network_model["branch"] = Dict{String,Any}()
network_model["storage"] = Dict{String,Any}()
open(dir * config_file_name,"r") do io
configuration_json_dict = JSON.parse(io)
end;
#voltage_base = configuration_json_dict["gridConfig"]["basekV"]
#power_base = configuration_json_dict["gridConfig"]["baseMVA"]
sub_dir=splitpath(config_file_name)[1]
configuration = configuration_json_dict["gridConfig"]["connection_configuration"]
branches_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["branches_file"])[2]
buses_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["buses_file"])[2]
devices_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["devices_file"])[2]


open(dir * buses_file_name,"r") do io
buses_json_dict = JSON.parse(io)
    for bus in buses_json_dict
        id = bus["busId"] + 1 #Indexing starts at one in Julia
        id_s = string(id)
        network_model["bus_lookup"][id_s] = id
        network_model["settings"]["vbases_default"][id_s] =  voltage_base

        if id == 1 #Settings for slack bus
            network_model["bus"][id_s] = Dict{String,Any}(
                "name"      => "slack",
                "bus_type"  => 3,
                ##"grounded"  => Bool[0, 0, 0],
                #"terminals" => [1, 2, 3],
                "vbase"     =>  voltage_base,
                "index"     => id,
                "bus_i"     => id,
                "λvmin"     => 1.65, #1.03643,
                "λvmax"     => 1.65, #1.03643,
                "vmin"      => 1.0,
                "vmax"      => 1,
                "va"        => 0.0,
                "vm"        => 1, 
                #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"run_type"  => 1
                )
        else
            network_model["bus"][id_s] = Dict{String,Any}(
                "name"      => id_s,
                "bus_type"  => 1,
                #"grounded"  => Bool[0, 0, 0],
                #"terminals" => [1, 2, 3],
                "vbase"     =>  voltage_base,
                "index"     => id,
                "bus_i"     => id,
                "λvmin"     => 1.65,#1.03643,
                "λvmax"     => 1.65, #1.03643,
                "vmin"      =>  0.95,
                "vmax"      => 1.05, 
                #"LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                #"LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                "run_type"  => 1
                )
        end;
    end;
end;

#print(device_df)
open(dir * devices_file_name,"r") do io
devices_json_dict = JSON.parse(io)
  for device in devices_json_dict["LVcustomers"]
    id = device["deviceId"] + 1 #Indexing starts at one in Julia
    d=device_df[in(id-1).(device_df.dev_id),:]
    id_s = string(id)
    μ = dist_lv_ts_feeder[in(d[!,"category"][1]).(dist_lv_ts_feeder.cluster),:][!,"lower"][1]
    σ  = dist_lv_ts_feeder[in(d[!,"category"][1]).(dist_lv_ts_feeder.cluster),:][!,"upper"][1] 
    cons = convert(Float64,device["yearlyNetConsumption"])
    network_model["load"][id_s] = Dict{String,Any}(
        #"connections"   => vec(Int.(device["phases"])),
        "name"          => id_s*"-"*device["coded_ean"],
        "status"        => 1,
        "vbase"         =>  voltage_base,
        "vnom_kv"       => 1.0,
        "source_id"     => device["coded_ean"],
        "load_bus"      => device["busId"] + 1,
        "dispatchable"  => 0,
        "index"         => id,
        "yearlyNetConsumption" => cons,
        #"phases"        => device["phases"],
        "pd"            => max(0.1, μ)/1e3/power_base/ 3,
        "qd"            =>  max(0.01,μ)/1e3/ power_base/ 3/10,
        "p_inj"         => 0.0,
        "q_inj"         => 0.0,
        "conn_cap_kW"   => device["connectionCapacity"],
        "dst_id" => d[!,"category"][1],
        "cluster_id"  => findall(x->x==1,[s_dict["$i"]["dst_id"]==d[!,"category"][1] for i=1:length(s_dict)])[1],
        "μ"  => μ,
        "σ"  => σ 

    )
    end;
end;

open(dir * branches_file_name,"r") do io
branches_json_dict = JSON.parse(io)
impedance_dict = Dict{String,Any}(
"BT - Desconocido BT" => [0.21, 0.075],
"BT - MANGUERA" => [0.3586, 0.089], #[1.23, 0.08],
"BT - RV 0,6/1 KV 2*16 KAL" => [2.14, 0.09], #16 = 20 A
"BT - RV 0,6/1 KV 2*25 KAL" => [1.34, 0.097],
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => [0.2309, 0.085],
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => [0.1602, 0.079],
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => [0.1602, 0.079],
"BT - RV 0,6/1 KV 4*25 KAL" => [1.34, 0.097],
"BT - RV 0,6/1 KV 4*50 KAL" => [0.71849, 0.093],
"BT - RV 0,6/1 KV 4*95 KAL" => [0.3586, 0.089],
"BT - RX 0,6/1 KV 2*16 Cu" => [1.23, 0.08],
"BT - RX 0,6/1 KV 2*2 Cu" => [9.9, 0.075],
"BT - RX 0,6/1 KV 2*4 Cu" => [4.95, 0.075],
"BT - RX 0,6/1 KV 2*6 Cu" => [3.3, 0.075],
"BT - RZ 0,6/1 KV 2*16 AL" => [2.14, 0.09],
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => [0.2309, 0.85],
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => [0.2309, 0.85],
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => [1.34, 0.097],
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => [0.9073, 0.095],
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => [0.718497, 0.093],
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => [0.4539, 0.091],
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => [0.3586, 0.089],
"BT - RZ 0,6/1 KV 4*16 AL" => [2.14, 0.09],
"aansluitkabel" => [1.15, 0.150]
)

open(dir * branches_file_name,"r") do io
branches_json_dict = JSON.parse(io)
currentmax_dict = Dict{String,Any}(
"BT - Desconocido BT" => 200,
"BT - MANGUERA" => 150, #200 certain  40.18#150
"BT - RV 0,6/1 KV 2*16 KAL" => 75,
"BT - RV 0,6/1 KV 2*25 KAL" => 100,
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => 305,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => 344,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => 344,
"BT - RV 0,6/1 KV 4*25 KAL" => 100,
"BT - RV 0,6/1 KV 4*50 KAL" => 150,
"BT - RV 0,6/1 KV 4*95 KAL" => 230,
"BT - RX 0,6/1 KV 2*16 Cu" => 75,#95,
"BT - RX 0,6/1 KV 2*2 Cu" => 40, #30,
"BT - RX 0,6/1 KV 2*4 Cu" => 60,#40,
"BT - RX 0,6/1 KV 2*6 Cu" => 80, #50,
"BT - RZ 0,6/1 KV 2*16 AL" => 20, #75,
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => 305, #264,
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => 305,#264,
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => 100, #78.98,
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => 120,
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => 150, #118.47,
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => 160,
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => 230, # 182.21,
"BT - RZ 0,6/1 KV 4*16 AL" => 75,
"aansluitkabel" => 120 #200
)


    for branch in branches_json_dict
        id = branch["branchId"] +1
        id_s = string(id)
        network_model["branch"][id_s] = Dict{String,Any}(
            "shift"         => 0.0,
            #"f_connections" => [1, 2, 3],
            "name"          => id_s,
            "switch"        => false,
            "g_to"          => 0.0,
            "c_rating_a"    => currentmax_dict[branch["cableType"]],
            "vbase"         =>  voltage_base,
            "g_fr"          => 0.0,
            #"t_connections" => [1, 2, 3],
            "f_bus"         => branch["upBusId"]+1,
            "b_fr"          => 0.0,
            "c_rating_b"    => currentmax_dict[branch["cableType"]],
            "br_status"     => 1,
            "t_bus"         => branch["downBusId"]+1,
            "b_to"          => 0.0,
            "index"         => id,
            "angmin"        => -1.0472,
            "angmax"        => 1.0472,
            "transformer"   => false,
            "tap"           => 1.0,
            "c_rating_c"    => currentmax_dict[branch["cableType"]], 
            "λcmax"         => 1.65 #1.65 #2.5 #1.03643    
            )

        if haskey(impedance_dict,branch["cableType"])
            # network_model["branch"][id_s]["br_x"] = (cross_area_fact^(1/4)).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
            # network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base
            network_model["branch"][id_s]["br_x"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
            network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base  
        end;
        

        if haskey(currentmax_dict,branch["cableType"])
            network_model["branch"][id_s]["rate_a"] = cross_area_fact.*((currentmax_dict[branch["cableType"]]*voltage_base)/1e3)/power_base
            
            #network_model["branch"][id_s]["I_rating"] = currentmax_dict[branch["cableType"]]/current_base

        end;
    end;
end;
end;


network_model["sdata"]= s_dict 
network_model["curt"]= curt

network_model["PV"]=deepcopy(network_model["load"]);
[network_model["PV"][d]["μ"]=s_dict[string(length(s_dict))]["pc"] for d in   keys(network_model["PV"])]
[network_model["PV"][d]["σ"]=s_dict[string(length(s_dict))]["pd"] for d in   keys(network_model["PV"])]
[network_model["PV"][d]["pd"]=s_dict[string(length(s_dict))]["pd"]/1e6/ power_base / 3 for d in   keys(network_model["PV"])]
return network_model
end;

"""
Function that builds a network model in the PowerModelsDistribution format
(mathematical) from a JSON file name & location. Three phase active & reactive power
for all devices in the network can be set using pd & qd respectively (in pu).
"""

function build_mathematical_model_mc(dir, config_file_name, load_dist_csv, pv_dist_csv; t_s=52, pd = [0.0, 0.0, 0.0], qd = [0.0, 0.0, 0.0], scale_factor = 1.0, curt=0.0, cross_area_fact=1.0)


"""
Specify voltage and power base, voltage base should be the phase-to-ground voltage
of the feeder studied (kV), the power base can be arbitrairly chosen (MW)
"""
voltage_base = 0.230  # (kV)
power_base = 0.5  # (MW)
Z_base = voltage_base^2/power_base # (Ohm)
current_base = power_base/(voltage_base*1e-3) # (A)
configuration = "star"
mwpu = 1/power_base
kwpu = (1e-3)/power_base

network_model = Dict{String,Any}()
configuration_json_dict = Dict{Any,Any}()
device_df=CSV.read(dir*config_file_name[1:length(config_file_name)-19]*".csv", DataFrame)

dist_lv=CSV.read(dir*load_dist_csv, DataFrame)
dist_pv=CSV.read(dir*pv_dist_csv, DataFrame)
# dist_pv=CSV.read(dir*"beta_pm_2022_181"*".csv", DataFrame)
dist_pv_ts= dist_pv[in([t_s]).(dist_pv.timeslot),:]
dist_lv_ts=dist_lv[in([t_s]).(dist_lv.timeslot),:]

dist_lv_ts_feeder = dist_lv_ts[in(unique(device_df.category)).(dist_lv_ts.cluster),:]
s_dict=Dict()
i=1
for dist in eachrow(dist_lv_ts_feeder)
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = dist["cluster"]
    s["pa"]= dist["alpha"]
    s["pb"]= dist["beta"]
    s["pc"]= dist["lower"]
    s["pd"]= dist["lower"]+dist["upper"]
    s_dict[string(i)] = s
    i=i+1
end


##add Irradiance if day time or there is some Irradiance
if dist_pv_ts.upper[1]>0
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = 55
    s["pa"]= dist_pv_ts[!,"alpha"][1]
    s["pb"]= dist_pv_ts[!,"beta"][1]
    s["pc"]= dist_pv_ts[!,"lower"][1]
    s["pd"]= dist_pv_ts[!,"lower"][1]+dist_pv_ts[!,"upper"][1]
    s_dict[string(i)] = s
end

network_model["multinetwork"]   = false
network_model["is_kron_reduced"] = true
network_model["conductor_ids"]   = [1, 2, 3]
network_model["p_factor"] = 0.95
network_model["q_factor"] = sqrt(1-0.95^2)
network_model["dcline"] = Dict{String,Any}()
network_model["switch"] = Dict{String,Any}()
network_model["is_projected"] = true
network_model["per_unit"] = true
network_model["data_model"] = _PMD.MATHEMATICAL
network_model["shunt"] = Dict{String,Any}()
network_model["transformer"] = Dict{String,Any}()
network_model["bus"] = Dict{String,Any}()
network_model["map"] = Dict{String,Any}()
network_model["conductors"] = 3
network_model["baseMVA"] =  power_base
network_model["basekv"] =  voltage_base
network_model["bus_lookup"] = Dict{Any,Int64}()
# network_model["run_type"] = 1
network_model["load"] = Dict{String,Any}()
network_model["gen"] = Dict{String,Any}("1" => Dict{String,Any}(
"pg"            => [1.0, 1.0, 1.0],
"model"         => 2,
"connections"   => [1, 2, 3],
"shutdown"      => 0.0,
"startup"       => 0.0,
"configuration" => WYE,
"name"          => "virtual_generator",
"qg"            => [1.0, 1.0, 1.0],
"gen_bus"       => 1,
"vbase"         =>  voltage_base,
"source_id"     => Any["gen",1],
"index"         => 1,
"cost"          => [20000.0, 1400.0],
"gen_status"    => 1,
"qgmax"          => [1.275, 1.275, 1.275],
"qgmin"          => [-1.275, -1.275, -1.275],
"pgmax"          => [5, 5, 5],
"pgmin"          => [-5, -5, -5],
"ncost"         => 2,
"λpmin"         => 1.65, #1.03643 ,
"λpmax"         => 1.65, #1.03643 ,
"λqmin"         => 1.65, #1.03643 ,
"λqmax"         => 1.65 #1.03643
))
network_model["settings"] = Dict{String,Any}(
"sbase_default"        => power_base,
"vbases_default"       => Dict{String,Any}(), #No default is specified for now, since default is never used
"voltage_scale_factor" => 1E3, #Voltages are thus expressed in kV
"sbase"                => power_base,
"power_scale_factor"   => 1E6, #Power is expressed in MW
"base_frequency"       => 50.0 #Hertz
)
network_model["branch"] = Dict{String,Any}()
network_model["storage"] = Dict{String,Any}()
open(dir * config_file_name,"r") do io
configuration_json_dict = JSON.parse(io)
end;
#voltage_base = configuration_json_dict["gridConfig"]["basekV"]
#power_base = configuration_json_dict["gridConfig"]["baseMVA"]
sub_dir=splitpath(config_file_name)[1]
configuration = configuration_json_dict["gridConfig"]["connection_configuration"]
branches_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["branches_file"])[2]
buses_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["buses_file"])[2]
devices_file_name = sub_dir*"/"*splitpath(configuration_json_dict["gridConfig"]["devices_file"])[2]


open(dir * buses_file_name,"r") do io
buses_json_dict = JSON.parse(io)
    for bus in buses_json_dict
        id = bus["busId"] + 1 #Indexing starts at one in Julia
        id_s = string(id)
        network_model["bus_lookup"][id_s] = id
        network_model["settings"]["vbases_default"][id_s] =  voltage_base

        if id == 1 #Settings for slack bus
            network_model["bus"][id_s] = Dict{String,Any}(
                "name"      => "slack",
                "bus_type"  => 3,
                "grounded"  => Bool[0, 0, 0],
                "rg"        => [5, 5, 5],
                "xg"        => [0.1, 0.1, 0.1],
                "terminals" => [1, 2, 3],
                "vbase"     =>  voltage_base,
                "index"     => id,
                "bus_i"     => id,
                "λvmin"     => 1.65, #1.03643,
                "λvmax"     => 1.65, #1.03643,
                "vmin"      => [0.9, 0.9, 0.9],
                "vmax"      => [1.1, 1.1, 1.1],
                "va"        => [0.0, -2.0944, 2.0944],
                # "vm"        => [1.0, 1.0, 1.0], 
                "LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                "LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                "LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                "LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                "run_type"  => 1
                )
        else
            network_model["bus"][id_s] = Dict{String,Any}(
                "name"      => id_s,
                "bus_type"  => 1,
                "grounded"  => Bool[0, 0, 0],
                "terminals" => [1, 2, 3],
                "vbase"     =>  voltage_base,
                "index"     => id,
                "bus_i"     => id,
                "λvmin"     => 1.65,#1.03643,
                "λvmax"     => 1.65, #1.03643,
                "vmin"      => [0.9, 0.9, 0.9],
                "vmax"      => [1.1, 1.1, 1.1], 
                # "LPp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                # "LPq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                # "LQp"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                # "LQq"       => Dict(c => Dict(b => 0.0 for b in 1:length(keys(buses_json_dict))) for c in 1:3),
                "run_type"  => 1
                )
        end;
    end;
end;

#print(device_df)
open(dir * devices_file_name,"r") do io
devices_json_dict = JSON.parse(io)
  for device in devices_json_dict["LVcustomers"]
    id = device["deviceId"] + 1 #Indexing starts at one in Julia
    d=device_df[in(id-1).(device_df.dev_id),:]
    id_s = string(id)
    μ = dist_lv_ts_feeder[in(d[!,"category"][1]).(dist_lv_ts_feeder.cluster),:][!,"lower"][1]
    σ  = dist_lv_ts_feeder[in(d[!,"category"][1]).(dist_lv_ts_feeder.cluster),:][!,"upper"][1] 
    cons = convert(Float64,device["yearlyNetConsumption"])
    pd_temp=5e-3
    print(pd_temp.*500)
    network_model["load"][id_s] = Dict{String,Any}(
        "model"         => POWER,
        "configuration" => configuration=="star" ? WYE : DELTA,
        "connections"   => vec(Int.(device["phases"])),
        "name"          => id_s*"-"*device["coded_ean"],
        "status"        => 1,
        "vbase"         =>  voltage_base,
        "vnom_kv"       => 1.0,
        "source_id"     => device["coded_ean"],
        "load_bus"      => device["busId"] + 1,
        "dispatchable"  => 0,
        "index"         => id,
        "yearlyNetConsumption" => cons,
        "phases"        => device["phases"],
        "pd"            => length(device["phases"])==1 ? [pd_temp] : fill(pd_temp./3, 3),
        "qd"            => length(device["phases"])==1 ? [pd_temp./10] : fill(pd_temp./30, 3),
        # "p_inj"         => fill(0.0, 3),
        # "q_inj"         => fill(0.0, 3),
        "conn_cap_kW"   => device["connectionCapacity"],
        "dst_id" => d[!,"category"][1],
        "cluster_id"  => findall(x->x==1,[s_dict["$i"]["dst_id"]==d[!,"category"][1] for i=1:length(s_dict)])[1],
        "μ"  => μ,
        "σ"  => σ 

    )
    end;
end;

open(dir * branches_file_name,"r") do io
branches_json_dict = JSON.parse(io)
impedance_dict = Dict{String,Any}(
"BT - Desconocido BT" => [0.21, 0.075],
"BT - MANGUERA" => [0.3586, 0.089], #[1.23, 0.08],
"BT - RV 0,6/1 KV 2*16 KAL" => [2.14, 0.09], #16 = 20 A
"BT - RV 0,6/1 KV 2*25 KAL" => [1.34, 0.097],
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => [0.2309, 0.085],
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => [0.1602, 0.079],
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => [0.1602, 0.079],
"BT - RV 0,6/1 KV 4*25 KAL" => [1.34, 0.097],
"BT - RV 0,6/1 KV 4*50 KAL" => [0.71849, 0.093],
"BT - RV 0,6/1 KV 4*95 KAL" => [0.3586, 0.089],
"BT - RX 0,6/1 KV 2*16 Cu" => [1.23, 0.08],
"BT - RX 0,6/1 KV 2*2 Cu" => [9.9, 0.075],
"BT - RX 0,6/1 KV 2*4 Cu" => [4.95, 0.075],
"BT - RX 0,6/1 KV 2*6 Cu" => [3.3, 0.075],
"BT - RZ 0,6/1 KV 2*16 AL" => [2.14, 0.09],
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => [0.2309, 0.85],
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => [0.2309, 0.85],
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => [1.34, 0.097],
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => [0.9073, 0.095],
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => [0.718497, 0.093],
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => [0.4539, 0.091],
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => [0.3586, 0.089],
"BT - RZ 0,6/1 KV 4*16 AL" => [2.14, 0.09],
"aansluitkabel" => [1.15, 0.150]
)

open(dir * branches_file_name,"r") do io
branches_json_dict = JSON.parse(io)
currentmax_dict = Dict{String,Any}(
"BT - Desconocido BT" => 200,
"BT - MANGUERA" => 150, #200 certain  40.18#150
"BT - RV 0,6/1 KV 2*16 KAL" => 75,
"BT - RV 0,6/1 KV 2*25 KAL" => 100,
"BT - RV 0,6/1 KV 3(1*150 KAL) + 1*95 KAL" => 305,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*150 KAL" => 344,
"BT - RV 0,6/1 KV 3(1*240 KAL) + 1*95 KAL" => 344,
"BT - RV 0,6/1 KV 4*25 KAL" => 100,
"BT - RV 0,6/1 KV 4*50 KAL" => 150,
"BT - RV 0,6/1 KV 4*95 KAL" => 230,
"BT - RX 0,6/1 KV 2*16 Cu" => 75,#95,
"BT - RX 0,6/1 KV 2*2 Cu" => 40, #30,
"BT - RX 0,6/1 KV 2*4 Cu" => 60,#40,
"BT - RX 0,6/1 KV 2*6 Cu" => 80, #50,
"BT - RZ 0,6/1 KV 2*16 AL" => 20, #75,
"BT - RZ 0,6/1 KV 3*150 AL/80 ALM" => 305, #264,
"BT - RZ 0,6/1 KV 3*150 AL/95 ALM" => 305,#264,
"BT - RZ 0,6/1 KV 3*25 AL/54,6 ALM" => 100, #78.98,
"BT - RZ 0,6/1 KV 3*35 AL/54,6 ALM" => 120,
"BT - RZ 0,6/1 KV 3*50 AL/54,6 ALM" => 150, #118.47,
"BT - RZ 0,6/1 KV 3*70 ALM/54,6 AL" => 160,
"BT - RZ 0,6/1 KV 3*95 AL/54,6 ALM" => 230, # 182.21,
"BT - RZ 0,6/1 KV 4*16 AL" => 75,
"aansluitkabel" => 120 #200
)


    for branch in branches_json_dict
        id = branch["branchId"] +1
        id_s = string(id)
        network_model["branch"][id_s] = Dict{String,Any}(
            "shift"         => [0.0, 0.0, 0.0],
            "f_connections" => [1, 2, 3],
            "name"          => id_s,
            "switch"        => false,
            "g_to"          => [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
            "c_rating_a"    => fill(currentmax_dict[branch["cableType"]], 3),
            "vbase"         =>  voltage_base,
            "g_fr"          => [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
            "t_connections" => [1, 2, 3],
            "f_bus"         => branch["upBusId"]+1,
            "b_fr"          => [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
            "c_rating_b"    => currentmax_dict[branch["cableType"]],
            "br_status"     => 1,
            "t_bus"         => branch["downBusId"]+1,
            "b_to"          => [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
            "index"         => id,
            "angmin"        => [-1.0472, -1.0472, -1.0472],
            "angmax"        => [1.0472, 1.0472, 1.0472],
            "transformer"   => false,
            "tap"           => [1.0, 1.0, 1.0],
            "c_rating_c"    => fill(currentmax_dict[branch["cableType"]], 3), 
            "λcmax"         => 1.65 #1.65 #2.5 #1.03643    
            )

        if haskey(impedance_dict,branch["cableType"])
            # network_model["branch"][id_s]["br_x"] = (cross_area_fact^(1/4)).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base 
            # network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base
            network_model["branch"][id_s]["br_x"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][2] .* (branch["cableLength"]/1000+1E-6)./  Z_base.* ([2.0 1.0 1.0; 1.0 2.0 1.0; 1.0 1.0 2.0]) 
            network_model["branch"][id_s]["br_r"] = (1/cross_area_fact).* impedance_dict[branch["cableType"]][1] .* (branch["cableLength"]/1000+1E-6)./  Z_base  .* ([2.0 1.0 1.0; 1.0 2.0 1.0; 1.0 1.0 2.0])
        end;
        

        if haskey(currentmax_dict,branch["cableType"])
            network_model["branch"][id_s]["rate_a"] = cross_area_fact.*fill(((currentmax_dict[branch["cableType"]]*voltage_base)/1e3)/power_base, 3)
            network_model["branch"][id_s]["I_rating"] = fill(currentmax_dict[branch["cableType"]]/current_base, 3)

            
            #network_model["branch"][id_s]["I_rating"] = currentmax_dict[branch["cableType"]]/current_base

        end;
    end;
end;
end;


network_model["sdata"]= s_dict 
network_model["curt"]= curt

network_model["PV"]=deepcopy(network_model["load"]);
[network_model["PV"][d]["μ"]=s_dict[string(length(s_dict))]["pc"] for d in   keys(network_model["PV"])]
[network_model["PV"][d]["σ"]=s_dict[string(length(s_dict))]["pd"] for d in   keys(network_model["PV"])]
[network_model["PV"][d]["pd"]=fill(s_dict[string(length(s_dict))]["pd"]/1e6/ power_base,length(network_model["PV"][d]["connections"])) for d in   keys(network_model["PV"])]
return network_model
end;


function build_stochastic_data_mc(data::Dict{String,Any}, deg::Int)
    # add maximum current
    for (nb, branch) in data["branch"]
        f_bus = branch["f_bus"]
        branch["cmax"] = branch["rate_a"] / data["bus"]["$f_bus"]["vmin"]
    end

    # build mop
    opq = [parse_dst(ns[2]["dst"], ns[2]["pa"], ns[2]["pb"], deg) for ns in data["sdata"]]
    mop = _PCE.MultiOrthoPoly(opq, deg)

    # build load matrix
    Nd, Npce = length(data["load"]), mop.dim
    nc= [length(data["load"]["$j"]["connections"]) for j=1:Nd]
    pd, qd = zeros(Nd, Npce), zeros(Nd, Npce)
    pd_g, qd_g = zeros(Nd, Npce), zeros(Nd, Npce)
    for nd in 1:Nd 
        # reactive power
        qd[nd,1] = data["load"]["$nd"]["qd"][1]
        # active power
        nb = data["load"]["$nd"]["load_bus"]
        ni = data["load"]["$nd"]["cluster_id"]
        if ni == 55
            pd[nd,1] = data["load"]["$nd"]["pd"][1]
        else
            base = data["baseMVA"]
            μ, σ = data["load"]["$nd"]["μ"] /1e3/ base, data["load"]["$nd"]["σ"] /1e3/ base
            if mop.uni[ni] isa _PCE.GaussOrthoPoly
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            else
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            end
        end

        np = length(opq)
        base = data["baseMVA"]
        μ, σ = data["PV"]["1"]["μ"]/1e6 / base , data["PV"]["1"]["σ"] /1e6/ base 
        
            if mop.uni[np] isa _PCE.GaussOrthoPoly
                pd_g[nd,[1,np+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np])
            else
                pd_g[nd,[1,np+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np])
            end
    end

    # replicate the data
    sdata = _PM.replicate(data, Npce)

    # add the stochastic data 
    sdata["T2"] = _PCE.Tensor(2,mop)
    sdata["T3"] = _PCE.Tensor(3,mop)
    sdata["T4"] = _PCE.Tensor(4,mop)
    sdata["mop"] = mop
    for nw in 1:Npce, nd in 1:Nd
        if nc[nd] == 1
            sdata["nw"]["$nw"]["load"]["$nd"]["pd"][1] = pd[nd,nw]
            sdata["nw"]["$nw"]["load"]["$nd"]["qd"][1] = qd[nd,nw]
        else
            [sdata["nw"]["$nw"]["load"]["$nd"]["pd"][j] = pd[nd,nw]/3 for j=1:nc[nd]]
            [sdata["nw"]["$nw"]["load"]["$nd"]["qd"][j] = qd[nd,nw]/3 for j=1:nc[nd]]
        end

    end
    for nw in 1:Npce, nd in 1:Nd
        if nc[nd] == 1
            sdata["nw"]["$nw"]["PV"]["$nd"]["pd"][1] = pd_g[nd,nw]
            sdata["nw"]["$nw"]["PV"]["$nd"]["qd"][1] = qd_g[nd,nw]
        else
            [sdata["nw"]["$nw"]["PV"]["$nd"]["pd"][j] = pd[nd,nw]/3 for j=1:nc[nd]]
            [sdata["nw"]["$nw"]["PV"]["$nd"]["qd"][j] = qd[nd,nw]/3 for j=1:nc[nd]]
        end

    end
    sdata["data_model"] =data["data_model"]
    [sdata["nw"]["$k"]["per_unit"] = data["per_unit"] for k=1:length(sdata["nw"])]
    return sdata
end

function build_stochastic_data_mc_dss(data, deg::Int)

    s_dict=Dict()
    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = 1
    s["pa"]= 0.5
    s["pb"]= 2
    s["pc"]= 0
    s["pd"]= data["load"]["1"]["pd"][1]
    s_dict["1"] = s

    s=Dict()
    s["dst"]= "Beta"
    s["dst_id"] = 2
    s["pa"]= 1
    s["pb"]= 2
    s["pc"]= 0
    s["pd"]= data["load"]["3"]["pd"][1]
    s_dict["2"] = s
    
    data["sdata"]=s_dict
    data["load"]["1"]["cluster_id"]=1
    data["load"]["2"]["cluster_id"]=1
    data["load"]["3"]["cluster_id"]=2



    opq = [parse_dst(ns[2]["dst"], ns[2]["pa"], ns[2]["pb"], deg) for ns in data["sdata"]]
    mop = _PCE.MultiOrthoPoly(opq, deg)

    # add maximum current
    for (nb, branch) in data["branch"]
        f_bus = branch["f_bus"]
        branch["cmax"] = branch["c_rating_a"] ./ data["bus"]["$f_bus"]["vmin"]
    end

    # build mop

    # build load matrix
    Nd, Npce = length(data["load"]), mop.dim
    nc= [length(data["load"]["$j"]["connections"]) for j=1:Nd]
    pd, qd = zeros(Nd, Npce), zeros(Nd, Npce)
    pd_g, qd_g = zeros(Nd, Npce), zeros(Nd, Npce)
    for nd in 1:Nd 
        # reactive power
        qd[nd,1] = data["load"]["$nd"]["qd"][1]
        # active power
        nb = data["load"]["$nd"]["load_bus"]
        ni = data["load"]["$nd"]["cluster_id"]
        if ni == 55
            pd[nd,1] = data["load"]["$nd"]["pd"][1]
        else
            # base = data["baseMVA"]
            μ, σ = data["sdata"]["$ni"]["pc"], data["sdata"]["$ni"]["pd"]
            if mop.uni[ni] isa _PCE.GaussOrthoPoly
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            else
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            end
        end
#For adding PV
        # np = length(opq)
        # base = data["baseMVA"]
        # μ, σ = data["PV"]["1"]["μ"]/1e6 / base , data["PV"]["1"]["σ"] /1e6/ base 
        
        #     if mop.uni[np] isa _PCE.GaussOrthoPoly
        #         pd_g[nd,[1,np+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np])
        #     else
        #         pd_g[nd,[1,np+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[np])
        #     end
    end

    # replicate the data
    sdata = _PM.replicate(data, Npce)

    # add the stochastic data 
    sdata["T2"] = _PCE.Tensor(2,mop)
    sdata["T3"] = _PCE.Tensor(3,mop)
    sdata["T4"] = _PCE.Tensor(4,mop)
    sdata["mop"] = mop
    for nw in 1:Npce, nd in 1:Nd
        if nc[nd] == 1
            sdata["nw"]["$nw"]["load"]["$nd"]["pd"][1] = pd[nd,nw]
            sdata["nw"]["$nw"]["load"]["$nd"]["qd"][1] = qd[nd,nw]
        else
            [sdata["nw"]["$nw"]["load"]["$nd"]["pd"][j] = pd[nd,nw]/3 for j=1:nc[nd]]
            [sdata["nw"]["$nw"]["load"]["$nd"]["qd"][j] = qd[nd,nw]/3 for j=1:nc[nd]]
        end

    end
    # for nw in 1:Npce, nd in 1:Nd
    #     if nc[nd] == 1
    #         sdata["nw"]["$nw"]["PV"]["$nd"]["pd"][1] = pd_g[nd,nw]
    #         sdata["nw"]["$nw"]["PV"]["$nd"]["qd"][1] = qd_g[nd,nw]
    #     else
    #         [sdata["nw"]["$nw"]["PV"]["$nd"]["pd"][j] = pd[nd,nw]/3 for j=1:nc[nd]]
    #         [sdata["nw"]["$nw"]["PV"]["$nd"]["qd"][j] = qd[nd,nw]/3 for j=1:nc[nd]]
    #     end

    # end
    sdata["data_model"] =data["data_model"]
    [sdata["nw"]["$k"]["per_unit"] = data["per_unit"] for k=1:length(sdata["nw"])]
    return sdata
end