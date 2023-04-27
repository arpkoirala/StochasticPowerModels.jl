################################################################################
#  Copyright 2023, Arpan Koirala                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# variable
## branch 
""
function variable_mc_branch_current(pm::AbstractUnbalancedIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PMD.variable_mc_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMD.variable_mc_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMD.variable_mc_branch_current_series_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMD.variable_mc_branch_current_series_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_branch_series_current_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_mc_variable_branch_current_real(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    p = Dict()
    q = Dict()
    bus = _PMD.ref(pm, nw, :bus)
    branch = _PMD.ref(pm, nw, :branch)
    for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch_from)
        b = branch[l]
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]
        f_connections = _PMD.ref(pm, nw, :branch, l, "f_connections")
        t_connections = _PMD.ref(pm, nw, :branch, l, "t_connections")

        vr_fr = [_PMD.var(pm, nw, :vr, i)[c] for c in f_connections]
        vi_fr = [_PMD.var(pm, nw, :vi, i)[c] for c in f_connections]
        cr_fr = [_PMD.var(pm, nw, :cr, (l,i,j))[c] for c in f_connections]
        ci_fr = [_PMD.var(pm, nw, :ci, (l,i,j))[c] for c in f_connections]

        vr_to = [_PMD.var(pm, nw, :vr, j)[c] for c in t_connections]
        vi_to = [_PMD.var(pm, nw, :vi, j)[c] for c in t_connections]
        cr_to = [_PMD.var(pm, nw, :cr, (l,j,i))[c] for c in t_connections]
        ci_to = [_PMD.var(pm, nw, :ci, (l,j,i))[c] for c in t_connections]

        cr[(l,i,j)] = (g_sh_fr .* vr_fr - b_sh_fr .* vi_fr)
        cr[(l,j,i)] = -csr_fr + g_sh_to .* vr_to - b_sh_to .* vi_to
    end

    report && _IM.sol_component_value_edge(pm, _PMD.pm_it_sym, nw, :branch, :cr_fr, :cr_to, _PMD.ref(pm, nw, :arcs_from), _PMD.ref(pm, nw, :arcs_to), cr)
end
"variable: `ci[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_mc_variable_branch_current_imaginary(pm::AbstractUnbalancedPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    bus = _PMD.ref(pm, nw, :bus)
    branch = _PMD.ref(pm, nw, :branch)

    for (l,i,j) in _PMD.ref(pm, nw, :arcs_branch_from)
        b = branch[l]
       
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]
        vr_fr = [_PMD.var(pm, nw, :vr, i)[c] for c in f_connections]
        vi_fr = [_PMD.var(pm, nw, :vi, i)[c] for c in f_connections]
        cr_fr = [_PMD.var(pm, nw, :crt, (l,i,j))[c] for c in f_connections]
        ci_fr = [_PMD.var(pm, nw, :cit, (l,i,j))[c] for c in f_connections]

        vr_to = [_PMD.var(pm, nw, :vr, j)[c] for c in t_connections]
        vi_to = [_PMD.var(pm, nw, :vi, j)[c] for c in t_connections]
        cr_to = [_PMD.var(pm, nw, :crt, (l,j,i))[c] for c in t_connections]
        ci_to = [_PMD.var(pm, nw, :cit, (l,j,i))[c] for c in t_connections]


        ci[(l,i,j)] .= (g_sh_fr * vi_fr + b_sh_fr * vr_fr)
        ci[(l,j,i)] .= -csi_fr + g_sh_to * vi_to + b_sh_to * vr_to
    end

    report && _IM.sol_component_value_edge(pm, _PMD.pm_it_sym, nw, :branch, :ci_fr, :ci_to, _PMD.ref(pm, nw, :arcs_from), _PMD.ref(pm, nw, :arcs_to), ci)
end

## load
""
function variable_mc_load_current(pm::AbstractUnbalancedIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_load_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_load_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end


"""
	function variable_mc_load_current_real(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates load real current variables `:crd` for models with explicit neutrals
"""
function variable_mc_load_current_real(pm::AbstractUnbalancedIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => load["connections"] for (i,load) in _PMD.ref(pm, nw, :load))
    crd = _PMD.var(pm, nw)[:crd] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_crd_$(i)",
            start = comp_start_value(_PMD.ref(pm, nw, :load, i), "crd_start", c, 0.0)
        ) for i in ids(pm, nw, :load)
    )
    _PMD.var(pm, nw)[:crd_bus] = Dict{Int, Any}()
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :crd, ids(pm, nw, :load), crd)
end


"""
	function variable_mc_load_current_imaginary(
		pm::ExplicitNeutralModels;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true
	)

Creates load imaginary current variables `:cid` for models with explicit neutrals
"""
function variable_mc_load_current_imaginary(pm::AbstractUnbalancedIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    connections = Dict(i => load["connections"] for (i,load) in _PMD.ref(pm, nw, :load))
    cid = _PMD.var(pm, nw)[:cid] = Dict(i => JuMP.@variable(pm.model,
            [c in connections[i]], base_name="$(nw)_cid_$(i)",
            start = comp_start_value(_PMD.ref(pm, nw, :load, i), "cid_start", c, 0.0)
        ) for i in ids(pm, nw, :load)
    )
    _PMD.var(pm, nw)[:cid_bus] = Dict{Int, Any}()
    report && _IM.sol_component_value(pm, pmd_it_sym, nw, :load, :cid, _PMD.ids(pm, nw, :load), cid)
end


# "variable: `crd[j]` for `j` in `load`"
# function variable_mc_load_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     crd = _PM.var(pm, nw)[:crd] = JuMP.@variable(pm.model,
#         [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_crd",
#         start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "crd_start")
#     )

#     report && _PM.sol_component_value(pm, nw, :load, :crd, _PM.ids(pm, nw, :load), crd)
# end
# "variable: `cid[j]` for `j` in `load`"
# function variable_mc_load_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
#     cid = _PM.var(pm, nw)[:cid] = JuMP.@variable(pm.model,
#         [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_cid",
#         start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "cid_start")
#     )

#     report && _PM.sol_component_value(pm, nw, :load, :cid, _PM.ids(pm, nw, :load), cid)
# end


## generator
""
function variable_mc_generator_current(pm::AbstractUnbalancedIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_mc_gen_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_mc_gen_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PMD.var(pm, nw)[:crg_bus] = Dict{Int, Any}()
    _PMD.var(pm, nw)[:cig_bus] = Dict{Int, Any}()
end

"""
Creates gen real current variables `:crg` for models with explicit neutrals
    """
function variable_mc_gen_current_real(pm::AbstractUnbalancedIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
        connections = Dict(i => gen["connections"] for (i,gen) in _PMD.ref(pm, nw, :gen))
        crg = _PMD.var(pm, nw)[:crg] = Dict(i => JuMP.@variable(pm.model,
                [c in connections[i]], base_name="$(nw)_crg_$(i)",
                start = comp_start_value(_PMD.ref(pm, nw, :gen, i), "crg_start", c, 0.0)
            ) for i in ids(pm, nw, :gen)
        )
        # _PMD.var(pm, nw)[:crg_bus] = Dict{Int, Any}()
        report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :crg, _PMD.ids(pm, nw, :gen), crg)
end

    """
Creates gen imaginary current variables `:cig` for models with explicit neutrals
    """
function variable_mc_gen_current_imaginary(pm::AbstractUnbalancedIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
        connections = Dict(i => gen["connections"] for (i,gen) in _PMD.ref(pm, nw, :gen))
        cig = _PMD.var(pm, nw)[:cig] = Dict(i => JuMP.@variable(pm.model,
                [c in connections[i]], base_name="$(nw)_cig_$(i)",
                start = comp_start_value(_PMD.ref(pm, nw, :gen, i), "cig_start", c, 0.0)
            ) for i in ids(pm, nw, :gen)
        )
        # _PMD.var(pm, nw)[:cig_bus] = Dict{Int, Any}()
        report && _IM.sol_component_value(pm, pmd_it_sym, nw, :gen, :cig, _PMD.ids(pm, nw, :gen), cig)
end


    # general constraints
## bus
""
function constraint_mc_gp_current_balance(pm::AbstractUnbalancedIVRModel, n::Int, i, terminals::Vector{Int}, grounded::Vector{Bool}, bus_arcs::Vector{Tuple{Tuple{Int,Int,Int},Vector{Int}}}, bus_gens::Vector{Tuple{Int,Vector{Int}}}, bus_loads::Vector{Tuple{Int,Vector{Int}}}, bus_shunts::Vector{Tuple{Int,Vector{Int}}})
    vr = _PMD.var(pm, n, :vr, i)
    vi = _PMD.var(pm, n, :vi, i)
    
    cr    = get(_PMD.var(pm, n),    :cr, Dict()); _PMD._check_var_keys(cr, bus_arcs, "real current", "branch")
    ci    = get(_PMD.var(pm, n),    :ci, Dict()); _PMD._check_var_keys(ci, bus_arcs, "imaginary current", "branch")
    crd   = get(_PMD.var(pm, n),   :crd_bus, Dict()); _PMD._check_var_keys(crd, bus_loads, "real current", "load")
    cid   = get(_PMD.var(pm, n),   :cid_bus, Dict()); _PMD._check_var_keys(cid, bus_loads, "imaginary current", "load")
    crg   = get(_PMD.var(pm, n),   :crg_bus, Dict()); _PMD._check_var_keys(crg, bus_gens, "real current", "generator")
    cig   = get(_PMD.var(pm, n),   :cig_bus, Dict()); _PMD._check_var_keys(cig, bus_gens, "imaginary current", "generator")
    # crs   = get(_PMD.var(pm, n),   :crs, Dict()); _PMD._check_var_keys(crs, bus_storage, "real currentr", "storage")
    # cis   = get(_PMD.var(pm, n),   :cis, Dict()); _PMD._check_var_keys(cis, bus_storage, "imaginary current", "storage")
    # crsw  = get(_PMD.var(pm, n),  :crsw, Dict()); _PMD._check_var_keys(crsw, bus_arcs_sw, "real current", "switch")
    # cisw  = get(_PMD.var(pm, n),  :cisw, Dict()); _PMD._check_var_keys(cisw, bus_arcs_sw, "imaginary current", "switch")
    # crt   = get(_PMD.var(pm, n),   :crt, Dict()); _PMD._check_var_keys(crt, bus_arcs_trans, "real current", "transformer")
    # cit   = get(_PMD.var(pm, n),   :cit, Dict()); _PMD._check_var_keys(cit, bus_arcs_trans, "imaginary current", "transformer")

    Gt, Bt = _PMD._build_bus_shunt_matrices(pm, n, terminals, bus_shunts)

    ungrounded_terminals = [(idx,t) for (idx,t) in enumerate(terminals) if !grounded[idx]]
    
   for (idx, t) in ungrounded_terminals
    # print(sum(crg[g][t] for (g, conns) in bus_gens if t in conns))
    # print([sum(crg[g][t]         for (g, conns) in bus_gens if t in conns)])
    # print([sum(crd[d][t]         for (d, conns) in bus_loads if t in conns)])
    JuMP.@constraint(pm.model,  #, #crs, crsw, crt,
                                      sum(cr[a][t] for (a, conns) in bus_arcs if t in conns)
                                   # + sum(crsw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                   # + sum(crt[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(crg[g][t]         for (g, conns) in bus_gens if t in conns)
                                    #- sum(crs[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(crd[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vr[u] -Bt[idx,jdx]*vi[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
    JuMP.@constraint(pm.model, #[ci, cid, cig,  vi], #cis, cisw, cit,
                                      sum(ci[a][t] for (a, conns) in bus_arcs if t in conns)
                                    # + sum(cisw[a_sw][t] for (a_sw, conns) in bus_arcs_sw if t in conns)
                                    # + sum(cit[a_trans][t] for (a_trans, conns) in bus_arcs_trans if t in conns)
                                    ==
                                      sum(cig[g][t]         for (g, conns) in bus_gens if t in conns)
                                    # - sum(cis[s][t]         for (s, conns) in bus_storage if t in conns)
                                    - sum(cid[d][t]         for (d, conns) in bus_loads if t in conns)
                                    - sum( Gt[idx,jdx]*vi[u] +Bt[idx,jdx]*vr[u] for (jdx,u) in ungrounded_terminals) # shunts
                                    )
    end
end

# galerkin projection constraint
## branch
""
function constraint_mc_gp_branch_series_current_magnitude_squared(pm::AbstractUnbalancedIVRModel, n::Int, i::Int, f_bus::Int, t_bus::Int, f_idx::Tuple{Int64, Int64, Int64}, t_idx::Tuple{Int64, Int64, Int64}, f_connections::Vector{Int}, t_connections::Vector{Int}, T2, T3)
    cmss  = _PMD.var(pm, n, :cmss, i)
    csr = Dict(nw => _PMD.var(pm, nw, :csr, i) for nw in _PMD.nw_ids(pm))
    csi = Dict(nw => _PMD.var(pm, nw, :csi, i) for nw in _PMD.nw_ids(pm))
    for (idx, c) in enumerate(t_connections)
        JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * cmss[c]
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (csr[n1][c] * csr[n2][c] + csi[n1][c] * csi[n2][c]) 
                                    for n1 in _PMD.nw_ids(pm), n2 in _PMD.nw_ids(pm))
                    )
                                    end
end

## generator
""

function constraint_mc_gp_generator_power_wye(pm::IVRUPowerModel, n::Int, id::Int, bus_id::Int, connections::Vector{Int}, pmin::Vector{<:Real}, pmax::Vector{<:Real}, qmin::Vector{<:Real}, qmax::Vector{<:Real}, T2, T3; report::Bool=true, bounded::Bool=true)
    vr = Dict(nw =>_PMD.var(pm, nw, :vr, bus_id) for nw in _PMD.nw_ids(pm))
    vi = Dict(nw =>_PMD.var(pm, nw, :vi, bus_id) for nw in _PMD.nw_ids(pm))
    crg = Dict(nw =>_PMD.var(pm, nw, :crg, id) for nw in _PMD.nw_ids(pm))
    cig = Dict(nw =>_PMD.var(pm, nw, :cig, id) for nw in _PMD.nw_ids(pm))
    # display(_PMD.var(pm, n, :pg, id)) 
    # pg = Vector{JuMP.NonlinearExpression}([])
    # qg = Vector{JuMP.NonlinearExpression}([])
    for (idx, c) in enumerate(connections)
        pg  = _PMD.var(pm, n, :pg, id)[c]
        qg  =  _PMD.var(pm, n, :qg, id)[c]
        JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pg 
                                            == 
                                            sum(T3.get([n1-1,n2-1,n-1]) *
                                                (vr[n1][c]*crg[n2][c]+vi[n1][c]*cig[n2][c])
                                                for n1 in _PMD.nw_ids(pm), n2 in _PMD.nw_ids(pm))
                                            )

        JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qg
                                        == 
                                        sum(T3.get([n1-1,n2-1,n-1]) *
                                        (-vr[n1][c]*cig[n2][c]+vi[n1][c]*crg[n2][c])
                                        for n1 in _PMD.nw_ids(pm), n2 in _PMD.nw_ids(pm))
                                        )
    end

    _PMD.var(pm, n, :crg_bus)[id] = crg[n]
    _PMD.var(pm, n, :cig_bus)[id] = cig[n]
    if report
        sol(pm, n, :gen, id)[:crg_bus] = _PMD.var(pm, n, :crg_bus, id)
        sol(pm, n, :gen, id)[:cig_bus] = _PMD.var(pm, n, :cig_bus, id)

    end
end


## load
""
function constraint_mc_gp_load_power_wye(pm::IVRUPowerModel, n::Int, id::Int, bus_id::Int, connections::Vector{Int}, a::Vector{<:Real}, alpha::Vector{<:Real}, b::Vector{<:Real}, beta::Vector{<:Real}, T2, T3, T4; report::Bool=true)
    vr = Dict(nw => _PMD.var(pm, nw, :vr, bus_id) for nw in _PMD.nw_ids(pm))
    vi = Dict(nw => _PMD.var(pm, nw, :vi, bus_id) for nw in _PMD.nw_ids(pm))
    vms = Dict(nw => _PMD.var(pm, nw, :vms, bus_id) for nw in _PMD.nw_ids(pm))
    
    crd = Dict(nw =>_PMD.var(pm, nw, :crd, id) for nw in _PMD.nw_ids(pm))
    cid = Dict(nw =>_PMD.var(pm, nw, :cid, id) for nw in _PMD.nw_ids(pm))
    if all(alpha.==0) && all(beta.==0)
        pd=a
        qd=b
        display("pd $n bus_id $bus_id: $pd \n")
        display("qd $n bus_id $bus_id: $qd \n")
    end
    for (idx, c) in enumerate(connections)
            print(idx)
            print(c)
            JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pd[idx] 
                                                == 
                                                sum(T3.get([n1-1,n2-1,n-1]) *
                                                    (vr[n1][c]*crd[n2][c]+vi[n1][c]*cid[n2][c])
                                                    for n1 in _PMD.nw_ids(pm), n2 in _PMD.nw_ids(pm))
                                                )


            JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qd[idx]
                                            == 
                                            sum(T3.get([n1-1,n2-1,n-1]) *
                                            (-vr[n1][c]*cid[n2][c]+vi[n1][c]*crd[n2][c])
                                            for n1 in _PMD.nw_ids(pm), n2 in _PMD.nw_ids(pm))
                                            )
        end
       
        _PMD.var(pm, n, :crd_bus)[id] = crd[n]
        _PMD.var(pm, n, :cid_bus)[id] = cid[n]
    if report
        sol(pm, n, :load, id)[:crd_bus] = _PMD.var(pm, n, :crd_bus, id)
        sol(pm, n, :load, id)[:cid_bus] = _PMD.var(pm, n, :cid_bus, id)
    end
end


function constraint_gp_load_power_real(pm::AbstractUnbalancedIVRModel, n::Int, i, l, pd, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))

    crd = Dict(nw => _PM.var(pm, nw, :crd, l) for nw in _PM.nw_ids(pm))
    cid = Dict(nw => _PM.var(pm, nw, :cid, l) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pd
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vr[n1] * crd[n2] + vi[n1] * cid[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end
""
function constraint_gp_load_power_imaginary(pm::AbstractUnbalancedIVRModel, n::Int, i, l, qd, T2, T3)
    vr  = Dict(n => _PM.var(pm, n, :vr, i) for n in _PM.nw_ids(pm))
    vi  = Dict(n => _PM.var(pm, n, :vi, i) for n in _PM.nw_ids(pm))

    crd = Dict(n => _PM.var(pm, n, :crd, l) for n in _PM.nw_ids(pm))
    cid = Dict(n => _PM.var(pm, n, :cid, l) for n in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qd
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crd[n2] - vr[n1] * cid[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

# solution
""
function sol_data_model!(pm::AbstractUnbalancedIVRModel, solution::Dict)
    _PMD.apply_pmd!(_sol_data_model_ivr!, solution)
end

""
function _sol_data_model_ivr!(solution::Dict)
    if haskey(solution, "bus")
        for (i, bus) in solution["bus"]
            if haskey(bus, "vr") && haskey(bus, "vi")
                for c=1:length(bus["vr"])
                    bus["vm"]=fill(1,length(bus["vr"]))
                    # bus["vm"][c] = hypot(bus["vr"][c], bus["vi"][c])
                    # bus["va"][c] = atan(bus["vi"][c], bus["vr"][c])
                end
            end
        end
    end
end