################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# bus
""
function variable_bus_voltage(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_bus_voltage_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
"variable: `vms[i]` for `i` in `bus`"
function variable_bus_voltage_magnitude_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vms = _PM.var(pm, nw)[:vms] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_vms",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "vms_start", 1.0)
    )

    if bounded
        for (i, bus) in _PM.ref(pm, nw, :bus)
            JuMP.set_lower_bound(vms[i], -2.0 * bus["vmax"]^2)
            JuMP.set_upper_bound(vms[i],  2.0 * bus["vmax"]^2)
        end
    end

    report && _PM.sol_component_value(pm, nw, :bus, :vms, _PM.ids(pm, nw, :bus), vms)
end
""
function variable_mc_bus_voltage(pm::AbstractUnbalancedACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMD.variable_mc_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMD.variable_mc_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    variable_mc_bus_voltage_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end


"variable: `vms[i]` for `i` in `bus`"
function variable_mc_bus_voltage_magnitude_squared(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    terminals = Dict(i => bus["terminals"] for (i,bus) in ref(pm, nw, :bus))
    vms = var(pm, nw)[:vms] = Dict(i => JuMP.@variable(pm.model,
            [t in terminals[i]], base_name="$(nw)_vms_$(i)",
            start = comp_start_value(ref(pm, nw, :bus, i), "vms_start", t, 1.0)
        ) for i in ids(pm, nw, :bus)
    )

    if bounded
        for (i, bus) in _PMD.ref(pm, nw, :bus)
            if haskey(bus, "vmax")
                for (idx, t) in enumerate(terminals[i])
                    JuMP.set_lower_bound(vms[i][t], 0 * bus["vmax"][idx]^2)
                    JuMP.set_upper_bound(vms[i][t],  2.0 * bus["vmax"][idx]^2)
        
                end
            end
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :bus, :vms, _PMD.ids(pm, nw, :bus), vms)
end

# branch
"variable: `cmss[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_branch_series_current_magnitude_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cmss = _PM.var(pm, nw)[:cmss] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :branch)], base_name="$(nw)_cmss",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :branch, l), "cmss_start", 0.0)
    )

    if bounded
        bus = _PM.ref(pm, nw, :bus)
        branch = _PM.ref(pm, nw, :branch)

        for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
            b = branch[l]
            ub = Inf
            if haskey(b, "rate_a")
                rate = b["rate_a"] * b["tap"]
                y_fr = abs(b["g_fr"] + im * b["b_fr"])
                y_to = abs(b["g_to"] + im * b["b_to"])
                shunt_current = max(y_fr * bus[i]["vmax"]^2, y_to * bus[j]["vmax"]^2)
                series_current = max(rate / bus[i]["vmin"], rate / bus[j]["vmin"])
                ub = series_current + shunt_current
            end
            if haskey(b, "c_rating_a")
                total_current = b["c_rating_a"]
                y_fr = abs(b["g_fr"] + im * b["b_fr"])
                y_to = abs(b["g_to"] + im * b["b_to"])
                shunt_current = max(y_fr * bus[i]["vmax"]^2, y_to * bus[j]["vmax"]^2)
                ub = total_current + shunt_current
            end

            if !isinf(ub)
                JuMP.set_lower_bound(cmss[l], -2.0 * ub^2)
                JuMP.set_upper_bound(cmss[l],  2.0 * ub^2)
            end
        end
    end

    report && _PM.sol_component_value(pm, nw, :branch, :cmss, _PM.ids(pm, nw, :branch), cmss)
end


"variable: `cmss[l,i,j]` for `(l,i,j)` in `arcs_from`"
function variable_mc_branch_series_current_magnitude_squared(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict((l,i,j) => connections for (bus,entry) in ref(pm, nw, :bus_arcs_conns_branch) for ((l,i,j), connections) in entry)
    cmss = _PMD.var(pm, nw)[:cmss] = Dict(l => JuMP.@variable(pm.model,
            [c in connections[(l,i,j)]], base_name="$(nw)_cmss_$(l)",
            start = comp_start_value(_PMD.ref(pm, nw, :branch, l), "cmss_start", c, 0.0)
        ) for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch)
    )

    if bounded
        bus = _PMD.ref(pm, nw, :bus)
        branch = _PMD.ref(pm, nw, :branch)

        for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch_from)
            # b = branch[l]
            # ub = Inf
            cmax = _PMD._calc_branch_series_current_max(ref(pm, nw, :branch, l), ref(pm, nw, :bus, i), ref(pm, nw, :bus, j))
            for (idx,c) in enumerate(connections[(l,i,j)])
                set_upper_bound(cmss[l][c],  cmax[idx]^2)
                set_lower_bound(cmss[l][c], 0*cmax[idx]^2)
            end
            # if haskey(b, "rate_a")
            #     rate = b["rate_a"] * b["tap"]
            #     y_fr = abs(b["g_fr"] + im * b["b_fr"])
            #     y_to = abs(b["g_to"] + im * b["b_to"])
            #     shunt_current = max(y_fr * bus[i]["vmax"]^2, y_to * bus[j]["vmax"]^2)
            #     series_current = max(rate / bus[i]["vmin"], rate / bus[j]["vmin"])
            #     ub = series_current + shunt_current
            # end
            # if haskey(b, "c_rating_a")
            #     total_current = b["c_rating_a"]
            #     y_fr = abs(b["g_fr"] + im * b["b_fr"])
            #     y_to = abs(b["g_to"] + im * b["b_to"])
            #     shunt_current = max(y_fr * bus[i]["vmax"]^2, y_to * bus[j]["vmax"]^2)
            #     ub = total_current + shunt_current
            # end

            # if !isinf(ub)
            #     JuMP.set_lower_bound(cmss[l], -2.0 * ub^2)
            #     JuMP.set_upper_bound(cmss[l],  2.0 * ub^2)
            # end
        end
    end

    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :branch, :cmss, _PMD.ids(pm, nw, :branch), cmss)
end

# generator
""
function variable_gen_power(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_gen_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_gen_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# PV size
"variable: `p_size` for `j` in `load`"
function variable_mc_pv_size(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    # connections = Dict(i => pv["connections"] for (i,pv) in _PMD.ref(pm, nw, :pv))
    p_size = _PMD.var(pm, nw)[:p_size] = Dict(i => JuMP.@variable(pm.model, base_name="$(nw)_p_size$(i)",
            start = comp_start_value(_PMD.ref(pm, nw, :pv, i), "crd_pv_start", i, 0.0)
        ) for i in ids(pm, nw, :pv)
        )
        if bounded
            for (i, pv) in _PMD.ref(pm, nw, :pv)
                if haskey(pv, "p_max") & haskey(pv,"p_min")
                    JuMP.set_lower_bound(p_size[i], pv["p_min"])
                    JuMP.set_upper_bound(p_size[i], pv["p_max"]) #2*PV["conn_cap_kW"])
                else
                    JuMP.set_lower_bound(p_size[i], 0)
                    JuMP.set_upper_bound(p_size[i], 15) #2*PV["conn_cap_kW"])
                end
            end
        end     
        _IM.sol_component_value(pm, pmd_it_sym, nw, :pv, :p_size, ids(pm, nw, :pv), p_size)
    end