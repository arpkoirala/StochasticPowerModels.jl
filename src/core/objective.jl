################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

"expected cost of active power generation"
function objective_min_expected_generation_cost(pm::AbstractPowerModel; kwargs...)
    gen_cost = Dict()

    T2 = pm.data["T2"]

    for (g, gen) in _PM.ref(pm, :gen, nw=1)
        pg = Dict(nw => _PM.var(pm, nw, :pg, g) for nw in _PM.nw_ids(pm))

        if length(gen["cost"]) == 1
            gen_cost[g] = gen["cost"][1]
        elseif length(gen["cost"]) == 2
            gen_cost[g] = gen["cost"][1]*pg[1] + 
                          gen["cost"][2]
        elseif length(gen["cost"]) == 3
            gen_cost[g] = gen["cost"][1]*sum(T2.get([n-1,n-1]) * pg[n]^2 for n in _PM.nw_ids(pm)) + 
                          gen["cost"][2]*pg[1] + 
                          gen["cost"][3]
        else
            gen_cost[g] = 0.0
        end
    end

    return JuMP.@objective(pm.model, Min,
            sum(gen_cost[g] for g in _PM.ids(pm, :gen, nw=1))
    )
end

"expected cost of active power generation"
function objective_mc_min_expected_generation_cost(pm::AbstractUnbalancedPowerModel; kwargs...)
    gen_cost = Dict()

    T2 = pm.data["T2"]
    for nw=1:1
        for (g, gen) in _PMD.ref(pm, :gen, nw=1)
            pg =  Dict(nw => sum(_PMD.var(pm, nw, :pg, g)[c] for c in gen["connections"]) for nw in _PMD.nw_ids(pm))
                if length(gen["cost"]) == 1
                    gen_cost[g] = gen["cost"][1]
                elseif length(gen["cost"]) == 2
                    gen_cost[g] = gen["cost"][1]*pg[1] + 
                                gen["cost"][2]
                elseif length(gen["cost"]) == 3
                    gen_cost[g][c] = gen["cost"][1]*sum(T2.get([n-1,n-1]) * pg[n]^2 for n in _PMD.nw_ids(pm)) + 
                                gen["cost"][2]*pg[1] + 
                                gen["cost"][3]
                else
                    gen_cost[g] = 0.0
                end
        end
    end
    return JuMP.@NLobjective(pm.model, Min,
            sum(gen_cost[g] for g in _PMD.ids(pm, :gen, nw=1))
    )
end

"expected max PV generation"
function objective_max_PV_mc(pm::AbstractUnbalancedPowerModel; kwargs...)
    p_size = Dict()


    for (p, pv) in _PMD.ref(pm, :pv, nw=1)
        p_size[p] = Dict(nw => _PMD.var(pm, nw, :p_size, p) for nw in [1])
    end

    return JuMP.@objective(pm.model, Max,
            sum(p_size[p][1] for p in _PMD.ids(pm, :pv, nw=1))
    )
end

# "expected absolute curtailment"
# function objective_min_PV_curtail_absolute(pm::AbstractUnbalancedPowerModel; kwargs...)
#     curtailment = Dict()

#     T2 = pm.data["T2"]
#     for nw=1:1
#         for (p, pv) in _PMD.ref(pm, :pv, nw=1)
#             p_curtail=  _PMD.var(pm, nw, :p_c, p)
#             p_size= _PMD.ref(pm, 1, :pv, p ,"p_size")
#             curtailment[p] = p_curtail*p_size
#         end
#     end
#    display(sum(curtailment[p] for p in _PMD.ids(pm, :pv, nw=1)))
#     return JuMP.@NLobjective(pm.model, Max,
#             sum(curtailment[p] for p in _PMD.ids(pm, :pv, nw=1))
#     )
# end

"expected relative curtailment"
function objective_min_PV_curtail(pm::AbstractUnbalancedPowerModel; kwargs...)
    inj = Dict()
    active=[]
    T2 = pm.data["T2"]
    for nw=1:1
        for (p, pv) in _PMD.ref(pm, :pv, nw=1)
            bus = ref(pm, nw, :bus, pv["load_bus"])
            load = ref(pm, nw, :load, p) #load of the consumer with pv p
            bus = ref(pm, nw,:bus, load["load_bus"])
            p_size= _PMD.ref(pm, 1, :pv,p,"p_size")
            a_load, alpha_load, b_load, beta_load = _PMD._load_expmodel_params(load, bus)
            a_pv, alpha_pv, b_pv, beta_pv = _PMD._load_expmodel_params(pv, bus)
            n_inj=  _PMD.var(pm, nw, :n_inj, p)
            inj_max= sum((a_pv[c]*p_size-a_load[c]) for c in 1:length(pv["connections"]))
            if inj_max>0
                inj[p] = n_inj*inj_max
                push!(active,p)
            else
                inj[p]= 0
            end
        end
    end
   display(sum(inj[p] for p in active))
    return JuMP.@NLobjective(pm.model, Max,
            sum(inj[p] for p in active)
    )
end

"expected relative curtailment"
function objective_max_PV_curtail(pm::AbstractUnbalancedPowerModel; kwargs...)
    curt = Dict()

    T2 = pm.data["T2"]
    for nw=1:1
        for (p, pv) in _PMD.ref(pm, :pv, nw=1)
           
            bus = ref(pm, nw, :bus, pv["load_bus"])
            load = ref(pm, nw, :load, p) #load of the consumer with pv p
            bus = ref(pm, nw,:bus, load["load_bus"])
            a_load, alpha_load, b_load, beta_load = _PMD._load_expmodel_params(load, bus)
            a_pv, alpha_pv, b_pv, beta_pv = _PMD._load_expmodel_params(pv, bus)
            n_curt=  _PMD.var(pm, nw, :p_c, p)*_PMD.ref(pm, 1, :pv, p ,"p_size")*sum(a_pv)
            curt[p] = n_curt
        end
    end
   display(sum(curt[p] for p in _PMD.ids(pm, :pv, nw=1)))
    return JuMP.@NLobjective(pm.model, Max,
            sum(curt[p] for p in _PMD.ids(pm, :pv, nw=1))
    )
end

"expected relative curtailment"
function objective_min_PV_curtail_absolute(pm::AbstractUnbalancedPowerModel; kwargs...)
    curtailment = Dict()

    T2 = pm.data["T2"]
    for nw=1:1
        for (p, pv) in _PMD.ref(pm, :pv, nw=1)
            p_curtail =  _PMD.var(pm, nw, :p_c, p)
            curtailment[p] = p_curtail*_PMD.ref(pm, 1, :pv,p,"p_size")
        end
    end
   display(sum(curtailment[p] for p in _PMD.ids(pm, :pv, nw=1)))
    return JuMP.@NLobjective(pm.model, Max,
            sum(curtailment[p] for p in _PMD.ids(pm, :pv, nw=1))
    )
end

"expected cost of active power generation"
function objective_max_PV_injection(pm::AbstractUnbalancedPowerModel; kwargs...)
    injection=Dict()
    T2 = pm.data["T2"]
    for nw=1:1
        for (p, pv) in _PMD.ref(pm, :pv, nw=1)
            bus = ref(pm, nw, :bus, pv["load_bus"])
            load = ref(pm, nw, :load, p) #load of the consumer with pv p
            bus = ref(pm, nw,:bus, load["load_bus"])
            a_load, alpha_load, b_load, beta_load = _PMD._load_expmodel_params(load, bus)
            a_pv, alpha_pv, b_pv, beta_pv = _PMD._load_expmodel_params(pv, bus)
            p_curtail =  _PMD.var(pm, nw, :p_c, p)
            p_size= _PMD.ref(pm, 1, :pv,p,"p_size")
            # if (a_pv*p_size)>a_load
            inj_max= sum((a_pv[c]*p_size-a_load[c]) for c in 1:length(pv["connections"]))
            injection[p]= sum((p_curtail*a_pv[c]*p_size-a_load[c]) for c in 1:length(pv["connections"]))*inj_max
            # else
                # injection[p] = 0
            # end
        end
    end
   display(sum(injection[p] for p in _PMD.ids(pm, :pv, nw=1)))
    return JuMP.@NLobjective(pm.model, Max,
    sum(injection[p] for p in _PMD.ids(pm, :pv, nw=1))
    )
end

"expected cost of active power generation"
function objective_equality_PV_injection(pm::AbstractUnbalancedPowerModel; kwargs...)
    injection=Dict()
    active=[]
    T2 = pm.data["T2"]
    for nw=1:1
        for (p, pv) in _PMD.ref(pm, :pv, nw=1)
            bus_pv = _PMD.ref(pm, nw, :bus, pv["load_bus"])
            load = _PMD.ref(pm, nw, :load, p) #load of the consumer with pv p
            bus_load = _PMD.ref(pm, nw,:bus, load["load_bus"])
            a_load, alpha_load, b_load, beta_load = _PMD._load_expmodel_params(load, bus_load)
            a_pv, alpha_pv, b_pv, beta_pv = _PMD._load_expmodel_params(pv, bus_pv)
            n_inj =  _PMD.var(pm, nw, :n_inj, p)
            p_size= _PMD.ref(pm, 1, :pv,p,"p_size")
            # p_size= _PMD.ref(pm, 1, :pv,id,"p_size")
            inj_max= sum((a_pv[c]*p_size-a_load[c]) for c in 1:length(pv["connections"]))
            if inj_max>0
                injection[p]= n_inj
                push!(active,p)
            else
                injection[p]= 0    
            end        
        end
    end
    n=length(active)
    ϵ=1e-3
    # u_mean=1.0
    u_mean=sum(injection[p] for p in active)/length(active)
    # JuMP.@expression(pm.model,  u_mean = sum(injection)/length(injection))
    return JuMP.@NLobjective(pm.model, Max,
    -(sum((injection[p]-u_mean)^2 for p in active)+ϵ)^(1/2)/(u_mean+ϵ)/n
    )
end


# "expected cost of active power generation"
# function objective_Qualityofservice_PV_injection(pm::AbstractUnbalancedPowerModel; kwargs...)
#     injection=Dict()
#     T2 = pm.data["T2"]
#     for nw=1:1
#         for (p, pv) in _PMD.ref(pm, :pv, nw=1)
#             bus_pv = _PMD.ref(pm, nw, :bus, pv["load_bus"])
#             load = _PMD.ref(pm, nw, :load, p) #load of the consumer with pv p
#             bus_load = _PMD.ref(pm, nw,:bus, load["load_bus"])
#             a_load, alpha_load, b_load, beta_load = _PMD._load_expmodel_params(load, bus_load)
#             a_pv, alpha_pv, b_pv, beta_pv = _PMD._load_expmodel_params(pv, bus_pv)
#             p_curtail =  _PMD.var(pm, nw, :p_c, p)
#             p_size= _PMD.ref(pm, 1, :pv,p,"p_size")
#             # p_size= _PMD.ref(pm, 1, :pv,id,"p_size")
            
#             injection[p]= sum((p_curtail*a_pv[c]*p_size-a_load[c]) for c in 1:length(pv["connections"]))/(sum((a_pv[c]*p_size)-a_load[c] for c in 1:length(pv["connections"])))
            
#         end
#     end

#    u_mean=0.8
#    σ = 0.5
#    display(u_mean)
#    display(injection)
#     return JuMP.@NLobjective(pm.model, Min,
#     (((sum((((injection[p]-u_mean))^2)/7 for p in _PMD.ids(pm, :pv, nw=1)))/σ))
#     )
#     # return JuMP.@NLobjective(pm.model, Min,
#     # (sum((injection[p]) for p in _PMD.ids(pm, :pv, nw=1))^2)/(7* sum((injection[p])^2 for p in _PMD.ids(pm, :pv, nw=1)))
#     # )
# end

"QoS"
function objective_Qualityofservice_PV_injection(pm::AbstractUnbalancedPowerModel; kwargs...)
    injection=Dict()
    active=[]
    for nw=1:1
        for (p, pv) in _PMD.ref(pm, :pv, nw=1)
            bus_pv = _PMD.ref(pm, nw, :bus, pv["load_bus"])
            load = _PMD.ref(pm, nw, :load, p) #load of the consumer with pv p
            bus_load = _PMD.ref(pm, nw,:bus, load["load_bus"])
            a_load, alpha_load, b_load, beta_load = _PMD._load_expmodel_params(load, bus_load)
            a_pv, alpha_pv, b_pv, beta_pv = _PMD._load_expmodel_params(pv, bus_pv)
            n_inj =  _PMD.var(pm, nw, :n_inj, p)
            p_size= _PMD.ref(pm, 1, :pv,p,"p_size")
            # p_size= _PMD.ref(pm, 1, :pv,id,"p_size")
            inj_max= sum((a_pv[c]*p_size-a_load[c]) for c in 1:length(pv["connections"]))
            if inj_max>=0
                injection[p]= n_inj
                push!(active,p)
            else
                injection[p]= 0    
            end        
        end
    end
    display(injection)
    u_mean=sum(injection[p] for p in active)/length(active)
   σ = 0.5
#    u_mean=1
   ϵ=1e-2
   n=length(active)
   display(u_mean)
#    return JuMP.@NLobjective(pm.model, Max,
#    sum(((injection[p]+ϵ)^(1-alpha))/(1-alpha) for p in _PMD.ids(pm, :pv, nw=1))
#    )
    return JuMP.@NLobjective(pm.model, Max,
    u_mean/100+1-((sum((((injection[p]-u_mean)^2)/n)/(σ^2) for p in active))+ϵ)^(1/2)
    )
    # return JuMP.@NLobjective(pm.model, Min,
    # (sum((injection[p]) for p in _PMD.ids(pm, :pv, nw=1))^2)/(7* sum((injection[p])^2 for p in _PMD.ids(pm, :pv, nw=1)))
    # )
end

# "expected cost of active power generation"
# function objective_alpha_PV_injection(pm::AbstractUnbalancedPowerModel; kwargs...)
#     injection=Dict()
#     T2 = pm.data["T2"]
#     for nw=1:1
#         for (p, pv) in _PMD.ref(pm, :pv, nw=1)
#             bus_pv = _PMD.ref(pm, nw, :bus, pv["load_bus"])
#             load = _PMD.ref(pm, nw, :load, p) #load of the consumer with pv p
#             bus_load = _PMD.ref(pm, nw,:bus, load["load_bus"])
#             a_load, alpha_load, b_load, beta_load = _PMD._load_expmodel_params(load, bus_load)
#             a_pv, alpha_pv, b_pv, beta_pv = _PMD._load_expmodel_params(pv, bus_pv)
#             p_curtail =  _PMD.var(pm, nw, :p_c, p)
#             p_size= _PMD.ref(pm, 1, :pv,p,"p_size")
#             # p_size= _PMD.ref(pm, 1, :pv,id,"p_size")
            
#             injection[p]= sum((p_curtail*a_pv[c]*p_size-a_load[c]) for c in 1:length(pv["connections"]))/(sum((a_pv[c]*p_size-a_load[c]) for c in 1:length(pv["connections"])))
            
#         end
#     end

#    alpha = 0
#    display(injection)
#     return JuMP.@NLobjective(pm.model, Max,
#     sum(((injection[p])^(1-alpha))/(1-alpha) for p in _PMD.ids(pm, :pv, nw=1))
#     )
# end

"Alpha fairness"
function objective_alpha_PV_injection(pm::AbstractUnbalancedPowerModel; kwargs...)
    injection=Dict()
    active=[]
    T2 = pm.data["T2"]
    for nw=1:1
        for (p, pv) in _PMD.ref(pm, :pv, nw=1)
            bus_pv = _PMD.ref(pm, nw, :bus, pv["load_bus"])
            load = _PMD.ref(pm, nw, :load, p) #load of the consumer with pv p
            bus_load = _PMD.ref(pm, nw,:bus, load["load_bus"])
            a_load, alpha_load, b_load, beta_load = _PMD._load_expmodel_params(load, bus_load)
            a_pv, alpha_pv, b_pv, beta_pv = _PMD._load_expmodel_params(pv, bus_pv)
            n_inj =  _PMD.var(pm, nw, :n_inj, p)
            p_size= _PMD.ref(pm, 1, :pv,p,"p_size")
            # p_size= _PMD.ref(pm, 1, :pv,id,"p_size")
            inj_max= sum((a_pv[c]*p_size-a_load[c]) for c in 1:length(pv["connections"]))
            if inj_max>0
                injection[p]= n_inj
                push!(active,p)
            else
                injection[p]= 0    
            end        
        end
    end
    display(injection)
   
   alpha = 10
    ϵ=1e-2
    return JuMP.@NLobjective(pm.model, Max,
    sum(((injection[p]+ϵ)^(1-alpha))/(1-alpha) for p in active)
    )
end