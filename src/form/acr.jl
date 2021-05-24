################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

#sorted_nw_ids(pm) = sort(collect(_PMs.nw_ids(pm)))

# variables
""
function variable_bus_voltage(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_bus_voltage_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_bus_voltage_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    #variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_bus_voltage_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

""
function variable_branch_power(pm::AbstractACRModel;nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_branch_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_branch_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    #variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_branch_current(pm::AbstractACRModel;nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_series_current_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end



function variable_gen_power(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMs.variable_gen_power_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMs.variable_gen_power_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# constraints

""
function constraint_voltage_ref(pm::AbstractACRModel, i::Int , nw::Int, vm)
    vr  = _PMs.var(pm, nw, :vr, i) 
    vi  = _PMs.var(pm, nw, :vi, i)
   
        JuMP.@constraint(pm.model, _PMs.var(pm, nw, :vi)[i] == 0)
    if nw ==1
        JuMP.@constraint(pm.model, _PMs.var(pm, nw, :vr)[i] == vm)
    else
        JuMP.@constraint(pm.model, _PMs.var(pm, nw, :vr)[i] == 0)

    end
end
""
function constraint_gp_bus_voltage_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    vs  = _PMs.var(pm, n, :vs, i)
    vr  = Dict(nw => _PMs.var(pm, nw, :vr, i) for nw in _PMs.nw_ids(pm))
    vi  = Dict(nw => _PMs.var(pm, nw, :vi, i) for nw in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vs 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * vr[n2] + vi[n1] * vi[n2]) 
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end

""
function constraint_theta_ref(pm::AbstractACRModel, n::Int, i::Int)
    vr = var(pm, n, :vr, i)
    vi = var(pm, n, :vi, i)

    vn = ifelse(n == 1, 1.0, 0.0)

    JuMP.@constraint(pm.model, vr == vn)
    JuMP.@constraint(pm.model, vi == 0.0)
end
""
function expression_branch_power_ohms_yt_from(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    vr_fr = _PMs.var(pm, n, :vr, f_bus)
    vr_to = _PMs.var(pm, n, :vr, t_bus)
    vi_fr = _PMs.var(pm, n, :vi, f_bus)
    vi_to = _PMs.var(pm, n, :vi, t_bus)

    _PMs.var(pm, n, :p)[f_idx] =  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
    _PMs.var(pm, n, :q)[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
end


""
function expression_branch_power_ohms_yt_to(pm::AbstractACRModel, n::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    vr_fr = _PMs.var(pm, n, :vr, f_bus)
    vr_to = _PMs.var(pm, n, :vr, t_bus)
    vi_fr = _PMs.var(pm, n, :vi, f_bus)
    vi_to = _PMs.var(pm, n, :vi, t_bus)

    _PMs.var(pm, n, :p)[t_idx] =  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
    _PMs.var(pm, n, :q)[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
end


""
function constraint_gp_load_power_real(pm::AbstractACRModel, n::Int, i, l, pd, T2, T3)
    vr  = Dict(nw => _PMs.var(pm, nw, :vr, i) for nw in _PMs.nw_ids(pm))
    vi  = Dict(nw => _PMs.var(pm, nw, :vi, i) for nw in _PMs.nw_ids(pm))

    crd = Dict(nw => _PMs.var(pm, nw, :crd, l) for nw in _PMs.nw_ids(pm))
    cid = Dict(nw => _PMs.var(pm, nw, :cid, l) for nw in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pd
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vr[n1] * crd[n2] + vi[n1] * cid[n2])
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end


""
function constraint_gp_load_power_imaginary(pm::AbstractACRModel, n::Int, i, l, qd, T2, T3)
    vr  = Dict(n => _PMs.var(pm, n, :vr, i) for n in _PMs.nw_ids(pm))
    vi  = Dict(n => _PMs.var(pm, n, :vi, i) for n in _PMs.nw_ids(pm))

    crd = Dict(n => _PMs.var(pm, n, :crd, l) for n in _PMs.nw_ids(pm))
    cid = Dict(n => _PMs.var(pm, n, :cid, l) for n in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qd
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crd[n2] - vr[n1] * cid[n2])
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end

""
function constraint_gp_power_branch_from(pm::AbstractACRModel, n::Int,f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
    #vs  = _PMs.var(pm, n, :vs, i)
    #bus = _PMs.ref(pm, n, :bus, i)
    #bus_arcs = _PMs.ref(pm, n, :bus_arcs, i)
    vr_fr = Dict(nw => _PMs.var(pm, nw, :vr, f_bus) for nw in _PMs.nw_ids(pm))
    vr_to = Dict(nw => _PMs.var(pm, nw, :vr, t_bus) for nw in _PMs.nw_ids(pm))
    vi_fr = Dict(nw => _PMs.var(pm, nw, :vi, f_bus) for nw in _PMs.nw_ids(pm))
    vi_to = Dict(nw => _PMs.var(pm, nw, :vi, t_bus) for nw in _PMs.nw_ids(pm))

    #_PMs.var(pm, n, :p)[f_idx] =  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
    #_PMs.var(pm, n, :q)[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)

    cstr_p = JuMP.@constraint(pm.model,
        _PMs.var(pm, n, :p)[f_idx]* T2.get([n-1,n-1])
        ==
        sum(T3.get([n1-1,n2-1,n-1]) *
        (g+g_fr)/tm^2*(vr_fr[n1]*vr_fr[n2] + vi_fr[n1]*vi_fr[n1]) + (-g*tr+b*ti)/tm^2*(vr_fr[n1]*vr_to[n2] + vi_fr[n1]*vi_to[n2]) + (-b*tr-g*ti)/tm^2*(vi_fr[n1]*vr_to[n2] - vr_fr[n1]*vi_to[n2])
        for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
    )

    cstr_q = JuMP.@constraint(pm.model,
    _PMs.var(pm, n, :q)[f_idx]*T2.get([n-1,n-1])
                    ==
                    sum(T3.get([n1-1,n2-1,n-1])*
                    -(b+b_fr)/tm^2*(vr_fr[n1]*vr_fr[n2] + vi_fr[n1]*vi_fr[n2]) - (-b*tr-g*ti)/tm^2*(vr_fr[n1]*vr_to[n2] + vi_fr[n1]*vi_to[n2]) + (-g*tr+b*ti)/tm^2*(vi_fr[n1]*vr_to[n2] - vr_fr[n1]*vi_to[n2])
                        for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                            )
end


""
function constraint_gp_power_branch_to(pm::AbstractACRModel, n::Int,f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, T2, T3)
    #vs  = _PMs.var(pm, n, :vs, i)
    #bus = _PMs.ref(pm, n, :bus, i)
    #bus_arcs = _PMs.ref(pm, n, :bus_arcs, i)
    vr_fr = Dict(nw => _PMs.var(pm, nw, :vr, f_bus) for nw in _PMs.nw_ids(pm))
    vr_to = Dict(nw => _PMs.var(pm, nw, :vr, t_bus) for nw in _PMs.nw_ids(pm))
    vi_fr = Dict(nw => _PMs.var(pm, nw, :vi, f_bus) for nw in _PMs.nw_ids(pm))
    vi_to = Dict(nw => _PMs.var(pm, nw, :vi, t_bus) for nw in _PMs.nw_ids(pm))

    #_PMs.var(pm, n, :p)[f_idx] =  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
    #_PMs.var(pm, n, :q)[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)

    cstr_p = JuMP.@constraint(pm.model,
        _PMs.var(pm, n, :p)[t_idx]* T2.get([n-1,n-1])
        ==
        sum(T3.get([n1-1,n2-1,n-1]) *
        (g+g_to)*(vr_to[n1]*vr_to[n2] + vi_to[n1]*vi_to[n2]) + (-g*tr-b*ti)/tm^2*(vr_fr[n1]*vr_to[n2] + vi_fr[n1]*vi_to[n2]) + (-b*tr+g*ti)/tm^2*(-(vi_fr[n1]*vr_to[n2] - vr_fr[n1]*vi_to[n2]))
            for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
    )

    cstr_q = JuMP.@constraint(pm.model,
    _PMs.var(pm, n, :q)[t_idx]*T2.get([n-1,n-1])
                    ==
                    sum(T3.get([n1-1,n2-1,n-1])*
                    -(b+b_to)*(vr_to[n1]*vr_to[n2] + vi_to[n1]*vi_to[n2]) - (-b*tr+g*ti)/tm^2*(vr_fr[n1]*vr_to[n2] + vi_fr[n1]*vi_to[n2]) + (-g*tr-b*ti)/tm^2*(-(vi_fr[n1]*vr_to[n2] - vr_fr[n1]*vi_to[n2]))
                        for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                            )
end

"""
function constraint_gp_branch_series_current_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    vs  = _PMs.var(pm, n, :vs, i)
    vr  = Dict(nw => _PMs.var(pm, nw, :vr, i) for nw in _PMs.nw_ids(pm))
    vi  = Dict(nw => _PMs.var(pm, nw, :vi, i) for nw in _PMs.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vs 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * vr[n2] + vi[n1] * vi[n2]) 
                                    for n1 in _PMs.nw_ids(pm), n2 in _PMs.nw_ids(pm))
                    )
end
"""

# chance constraints
""
function constraint_bus_voltage_squared_cc_limit(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop)
    vs  = [_PMs.var(pm, n, :vs, i) for n in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model,  _PCE.var(vs,T2)
                                <=
                                ((_PCE.mean(vs,mop) - vmin^2) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(vs,T2)
                                <=
                                ((vmax^2 - _PCE.mean(vs,mop)) / λmax)^2
                    )
end

""
function constraint_branch_series_current_squared_cc_limit(pm::AbstractACRModel, b, imax, λmax, T2, mop)
    css  = [_PMs.var(pm, nw, :css, b) for nw in sorted_nw_ids(pm)]
    display(css)
    JuMP.@constraint(pm.model,  _PCE.var(css,T2)
                                <=
                                ((imax^2 - _PCE.mean(css,mop)) / λmax)^2
                    )
end

""
function constraint_gen_power_real_cc_limit(pm::AbstractACRModel, g, pmin, pmax, λmin, λmax, T2, mop)
    pg  = [_PMs.var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model,  _PCE.var(pg,T2)
                                <=
                                ((_PCE.mean(pg,mop) - pmin) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(pg,T2)
                                <=
                                ((pmax - _PCE.mean(pg,mop)) / λmax)^2
                    )
end

""
function constraint_gen_power_imaginary_cc_limit(pm::AbstractACRModel, g, qmin, qmax, λmin, λmax, T2, mop)
    qg  = [_PMs.var(pm, nw, :qg, g) for nw in sorted_nw_ids(pm)]

    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                                <=
                                ((_PCE.mean(qg,mop) - qmin) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                                <=
                                ((qmax - _PCE.mean(qg,mop)) / λmax)^2
                    )
end

function constraint_power_balance(pm::AbstractACRModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vr = _PMs.var(pm, n, :vr, i)
    vi = _PMs.var(pm, n, :vi, i)
    p    = _PMs.get(_PMs.var(pm, n),    :p, Dict()); _PMs._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = _PMs.get(_PMs.var(pm, n),    :q, Dict()); _PMs._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = _PMs.get(_PMs.var(pm, n),   :pg, Dict()); _PMs._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = _PMs.get(_PMs.var(pm, n),   :qg, Dict()); _PMs._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = _PMs.get(_PMs.var(pm, n),   :ps, Dict()); _PMs._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = _PMs.get(_PMs.var(pm, n),   :qs, Dict()); _PMs._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = _PMs.get(_PMs.var(pm, n),  :psw, Dict()); _PMs._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = _PMs.get(_PMs.var(pm, n),  :qsw, Dict()); _PMs._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = _PMs.get(_PMs.var(pm, n), :p_dc, Dict()); _PMs._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = _PMs.get(_PMs.var(pm, n), :q_dc, Dict()); _PMs._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")


    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*(vr^2 + vi^2)
    )
    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*(vr^2 + vi^2)
    )

end
