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
    display(sum(gen_cost[g] for g in _PMD.ids(pm, :gen, nw=1)))
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