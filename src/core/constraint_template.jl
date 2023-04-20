################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# general constraints
## reference
""
function constraint_bus_voltage_ref(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    constraint_bus_voltage_ref(pm, nw, i)
end

""
function constraint_mc_bus_voltage_ref(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    va_ref = ref(pm, nw, :bus, i, "va")
    constraint_mc_bus_voltage_ref(pm, nw, i, va_ref)
end

## bus
""
function constraint_current_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    if !haskey(_PM.con(pm, nw), :kcl_cr)
        _PM.con(pm, nw)[:kcl_cr] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(_PM.con(pm, nw), :kcl_ci)
        _PM.con(pm, nw)[:kcl_ci] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_current_balance(pm, nw, i, bus_arcs, bus_gens, bus_loads, bus_gs, bus_bs)
end

""
function constraint_mc_gp_current_balance(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    # if !haskey(_PM.con(pm, nw), :kcl_cr)
    #     _PM.con(pm, nw)[:kcl_cr] = Dict{Int,JuMP.ConstraintRef}()
    # end
    # if !haskey(_PM.con(pm, nw), :kcl_ci)
    #     _PM.con(pm, nw)[:kcl_ci] = Dict{Int,JuMP.ConstraintRef}()
    # end

    bus = _PMD.ref(pm, nw, :bus, i)
    bus_arcs = _PMD.ref(pm, nw, :bus_arcs_conns_branch, i)
    bus_gens = _PMD.ref(pm, nw, :bus_conns_gen, i)
    bus_loads = _PMD.ref(pm, nw, :bus_conns_load, i)
    bus_shunts = _PMD.ref(pm, nw, :bus_conns_shunt, i)

    # bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    # bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_mc_gp_current_balance(pm, nw, i, bus["terminals"], bus["grounded"], bus_arcs, bus_gens, bus_loads, bus_shunts)
end
""
function constraint_power_balance(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus_arcs   = _PM.ref(pm, nw, :bus_arcs, i)
    bus_gens   = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads  = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)

    bus_pd = Dict(k => _PM.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PM.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance(pm, nw, i, bus_arcs, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end

# galerkin projection constraints
## bus
""
function constraint_gp_bus_voltage_magnitude_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_bus_voltage_magnitude_squared(pm, nw, i, T2, T3)
end

function constraint_mc_gp_bus_voltage_magnitude_squared(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    bus = _PMD.ref(pm, nw, :bus, i)
    constraint_mc_gp_bus_voltage_magnitude_squared(pm, nw, i,bus["terminals"], T2, T3)
end

## branch
""
function constraint_gp_branch_series_current_magnitude_squared(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_branch_series_current_magnitude_squared(pm, nw, b, T2, T3)
end

""
function constraint_mc_gp_branch_series_current_magnitude_squared(pm::AbstractUnbalancedPowerModel, b::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    branch = _PMD.ref(pm, nw, :branch, b)

    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (b, f_bus, t_bus)
    t_idx = (b, t_bus, f_bus)

    constraint_mc_gp_branch_series_current_magnitude_squared(pm, nw, b, f_bus, t_bus, f_idx, t_idx, branch["f_connections"], branch["t_connections"],T2, T3)
end
""
function constraint_gp_power_branch_to(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_to = branch["g_to"]
    b_to = branch["b_to"]
    tm = branch["tap"]

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_power_branch_to(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, T2, T3)
end
""
function constraint_gp_power_branch_from(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
    branch = _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = _PM.calc_branch_y(branch)
    tr, ti = _PM.calc_branch_t(branch)
    g_fr = branch["g_fr"]
    b_fr = branch["b_fr"]
    tm = branch["tap"]
    
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_power_branch_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
end

## generator
""
function constraint_gp_gen_power(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    i   = _PM.ref(pm, nw, :gen, g, "gen_bus")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_gen_power_real(pm, nw, i, g, T2, T3)
    constraint_gp_gen_power_imaginary(pm, nw, i, g, T2, T3)
end

## generator
""
function constraint_mc_gp_gen_power(pm::AbstractUnbalancedPowerModel, id::Int; nw::Int=nw_id_default,bounded::Bool=true, report::Bool=true)
    # i   = _PMD.ref(pm, nw, :gen, g, "gen_bus")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    generator = _PMD.ref(pm, nw, :gen, id)
    bus = _PMD.ref(pm, nw, :bus, generator["gen_bus"])

    N = length(generator["connections"])
    pmin = get(generator, "pmin", fill(-Inf, N))
    pmax = get(generator, "pmax", fill( Inf, N))
    qmin = get(generator, "qmin", fill(-Inf, N))
    qmax = get(generator, "qmax", fill( Inf, N))
    display(get(generator, "configuration", WYE))
    if get(generator, "configuration", WYE) == WYE
        constraint_mc_gp_generator_power_wye(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax, T2, T3; report=report, bounded=bounded)
    else
        constraint_mc_generator_power_delta(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax, T2, T3; report=report, bounded=bounded)
    end
    # constraint_mc_gp_gen_power_real(pm, nw, i, g, T2, T3)
    # constraint_mc_gp_gen_power_imaginary(pm, nw, i, g, T2, T3)
end

## load
""
function constraint_gp_load_power(pm::AbstractPowerModel, l::Int; nw::Int=nw_id_default)
    i   = _PM.ref(pm, nw, :load, l, "load_bus") 

    pd  = _PM.ref(pm, nw, :load, l, "pd")
    qd  = _PM.ref(pm, nw, :load, l, "qd")

    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    

    constraint_gp_load_power_real(pm, nw, i, l, pd, T2, T3)
    constraint_gp_load_power_imaginary(pm, nw, i, l, qd, T2, T3)
end


function constraint_mc_gp_load_power(pm::AbstractUnbalancedPowerModel, id::Int; nw::Int=nw_id_default, report::Bool=true)::Nothing
    load = ref(pm, nw, :load, id)
    bus = ref(pm, nw,:bus, load["load_bus"])
    
    configuration = load["configuration"]
    
    a, alpha, b, beta = _PMD._load_expmodel_params(load, bus)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]
    T4  = pm.data["T4"]
    
    if configuration==WYE
        constraint_mc_gp_load_power_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta, T2,T3, T4; report=report)
    else
        # constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    end
    nothing
end

# chance constraint limit
## bus
""
function constraint_cc_bus_voltage_magnitude_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PM.ref(pm, nw, :bus, i, "vmin")
    vmax = _PM.ref(pm, nw, :bus, i, "vmax")
    
    λmin = _PM.ref(pm, nw, :bus, i, "λvmin")
    λmax = _PM.ref(pm, nw, :bus, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_bus_voltage_magnitude_squared(pm, i, vmin, vmax, λmin, λmax, T2, mop)
end

## branch
""
function constraint_cc_branch_series_current_magnitude_squared(pm::AbstractPowerModel, b::Int; nw::Int=nw_id_default)
    cmax = _PM.ref(pm, nw, :branch, b, "cmax")
    λmax = _PM.ref(pm, nw, :branch, b, "λcmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_branch_series_current_magnitude_squared(pm, b,  cmax, λmax, T2, mop)
end

## generator
""
function constraint_cc_gen_power(pm::AbstractPowerModel, g::Int; nw::Int=nw_id_default)
    pmin = _PM.ref(pm, nw, :gen, g, "pmin")
    pmax = _PM.ref(pm, nw, :gen, g, "pmax")
    qmin = _PM.ref(pm, nw, :gen, g, "qmin")
    qmax = _PM.ref(pm, nw, :gen, g, "qmax")

    λpmin = _PM.ref(pm, nw, :gen, g, "λpmin")
    λpmax = _PM.ref(pm, nw, :gen, g, "λpmax")
    λqmin = _PM.ref(pm, nw, :gen, g, "λqmin")
    λqmax = _PM.ref(pm, nw, :gen, g, "λqmax")

    T2  = pm.data["T2"]
    mop = pm.data["mop"]
    
    constraint_cc_gen_power_real(pm, g, pmin, pmax, λpmin, λpmax, T2, mop)
    constraint_cc_gen_power_imaginary(pm, g, qmin, qmax, λqmin, λqmax, T2, mop)
end


# mc chance constraint limit
## bus
""
function constraint_mc_cc_bus_voltage_magnitude_squared(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)
    vmin = _PMD.ref(pm, nw, :bus, i, "vmin")
    vmax = _PMD.ref(pm, nw, :bus, i, "vmax")
    
    λmin = _PMD.ref(pm, nw, :bus, i, "λvmin")
    λmax = _PMD.ref(pm, nw, :bus, i, "λvmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_mc_cc_bus_voltage_magnitude_squared(pm, i, _PMD.ref(pm, nw, :bus, i)["terminals"], vmin, vmax, λmin, λmax, T2, mop)
end

## branch
""
function constraint_mc_cc_branch_series_current_magnitude_squared(pm::AbstractUnbalancedPowerModel, b::Int; nw::Int=nw_id_default)
    cmax = _PMD.ref(pm, nw, :branch, b, "cmax")
    λmax = _PMD.ref(pm, nw, :branch, b, "λcmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_mc_cc_branch_series_current_magnitude_squared(pm, b, _PMD.ref(pm, nw, :branch, b)["t_connections"], cmax, λmax, T2, mop)
end

## generator
""
function constraint_mc_cc_gen_power(pm::AbstractUnbalancedPowerModel, g::Int; nw::Int=nw_id_default)
    pmin = _PMD.ref(pm, nw, :gen, g, "pgmin")
    pmax = _PMD.ref(pm, nw, :gen, g, "pgmax")
    qmin = _PMD.ref(pm, nw, :gen, g, "qgmin")
    qmax = _PMD.ref(pm, nw, :gen, g, "qgmax")

    λpmin = _PMD.ref(pm, nw, :gen, g, "λpmin")
    λpmax = _PMD.ref(pm, nw, :gen, g, "λpmax")
    λqmin = _PMD.ref(pm, nw, :gen, g, "λqmin")
    λqmax = _PMD.ref(pm, nw, :gen, g, "λqmax")

    T2  = pm.data["T2"]
    mop = pm.data["mop"]
    
    constraint_mc_cc_gen_power_real(pm, g,_PMD.ref(pm, nw, :gen, g)["connections"], pmin, pmax, λpmin, λpmax, T2, mop)
    constraint_mc_cc_gen_power_imaginary(pm, g, _PMD.ref(pm, nw, :gen, g)["connections"],qmin, qmax, λqmin, λqmax, T2, mop)
end