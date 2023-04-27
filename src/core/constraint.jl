################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# general constraints
## bus

""
function constraint_bus_voltage_ref(pm::AbstractACRModel, n::Int, i::Int)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    vn = ifelse(n == 1, 1.0, 0.0)

    JuMP.@constraint(pm.model, vr == vn)
    JuMP.@constraint(pm.model, vi == 0.0)
end

""
function constraint_mc_bus_voltage_ref(pm::AbstractUnbalancedACRModel, n::Int, i::Int, va_ref::Vector{<:Real})
    vr = _PMD.var(pm, n, :vr, i)
    vi = _PMD.var(pm, n, :vi, i)

    for (idx, t) in enumerate(_PMD.ref(pm, n, :bus, i)["terminals"])
        if n==1
            display(va_ref[t])
            display("terminal $t \n")
            # if va_ref[t] == pi/2
            #     JuMP.@constraint(pm.model, vr[t] == 0)
            #     JuMP.@constraint(pm.model, vi[t] >= 0)
            # elseif va_ref[t] == -pi/2
            #     JuMP.@constraint(pm.model, vr[t] == 0)
            #     JuMP.@constraint(pm.model, vi[t] <= 0)
            # elseif va_ref[t] == 0
            #     JuMP.@constraint(pm.model, vr[t] >= 0)
            #     JuMP.@constraint(pm.model, vi[t] == 0)
            # elseif va_ref[t] == pi
            #     JuMP.@constraint(pm.model, vr[t] >= 0)
            #     JuMP.@constraint(pm.model, vi[t] == 0)
            # else
                if t==1
                    JuMP.fix(vr[t], 1; force=true)
                    
                    JuMP.fix(vi[t], 0; force=true)
                    print("\n Printing ref $vi[t]")
                elseif t==2
                    JuMP.fix(vr[t],-0.5; force=true)
                    print("\n Printing ref $vr[t]")
                    JuMP.fix(vi[t], -sqrt(3)/2; force=true)
                else
                    JuMP.fix(vr[t],-0.5; force=true)
                    JuMP.fix(vi[t],sqrt(3)/2; force=true)
                end

                # va_ref also implies a sign for vr, vi
                # if 0<=va_ref[t] && va_ref[t] <= pi
                #     JuMP.@constraint(pm.model, vi[t] >= 0)
                # else
                #     JuMP.@constraint(pm.model, vi[t] <= 0)
                # end
            # end
        else
            JuMP.fix(vr[t], 0; force=true)
            JuMP.fix(vi[t], 0; force=true)
        end
    end
    #     if 0<=va_ref[t] && va_ref[t] <= pi
    #         if t==1
    #             JuMP.@constraint(pm.model, vr[t] == vn)
    #         else
    #             JuMP.@constraint(pm.model, vr[t] == vn*cos(va_ref[t]))
    #         end
    #     else
    #         if t==1
    #             JuMP.@constraint(pm.model, vr[t] == vn)
    #         end
    #             JuMP.@constraint(pm.model, vr[t] == vn*cos(va_ref[t]))
    #     end
    #     JuMP.@constraint(pm.model, vi[t] == vn*sin(va_ref[t]))
        
    # end
end

# galerkin projection constraints
## bus
""
function constraint_gp_bus_voltage_magnitude_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    vms = _PM.var(pm, n, :vms, i)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vms 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * vr[n2] + vi[n1] * vi[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

## bus
""
function constraint_mc_gp_bus_voltage_magnitude_squared(pm::AbstractUnbalancedACRModel, n::Int, i, connections::Vector{Int}, T2, T3)
    vms = _PMD.var(pm, n, :vms, i)
   
    vr  = Dict(nw => _PMD.var(pm, nw, :vr, i) for nw in _PMD.nw_ids(pm))
    vi  = Dict(nw => _PMD.var(pm, nw, :vi, i) for nw in _PMD.nw_ids(pm))
  
    for (idx, c) in enumerate(connections)
        JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * vms[c] 
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1][c] * vr[n2][c] + vi[n1][c] * vi[n2][c]) 
                                    for n1 in _PMD.nw_ids(pm), n2 in _PMD.nw_ids(pm))
                    )
                                    end
end

# chance constraints
## bus
""
function constraint_cc_bus_voltage_magnitude_squared(pm::AbstractACRModel, i, vmin, vmax, λmin, λmax, T2, mop)
    vms  = [_PM.var(pm, n, :vms, i) for n in sorted_nw_ids(pm)]
    
    # bounds on the expectation
    JuMP.@constraint(pm.model, vmin^2 <= _PCE.mean(vms, mop))
    JuMP.@constraint(pm.model, _PCE.mean(vms, mop) <= vmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(vms, T2)
                                <=
                               ((_PCE.mean(vms, mop) - vmin^2) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(vms, T2)
                               <=
                                ((vmax^2 - _PCE.mean(vms, mop)) / λmax)^2
                    )
end

## branch
""
function constraint_cc_branch_series_current_magnitude_squared(pm::AbstractACRModel, b, cmax, λcmax, T2, mop)
    cmss = [_PM.var(pm, nw, :cmss, b) for nw in sorted_nw_ids(pm)]

    # bound on the expectation
    JuMP.@constraint(pm.model,  _PCE.mean(cmss, mop) <= cmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(cmss,T2)
                                <=
                                ((cmax^2 - _PCE.mean(cmss,mop)) / λcmax)^2
                    )
end

## generator
""
function constraint_cc_gen_power_real(pm::AbstractACRModel, g, pmin, pmax, λmin, λmax, T2, mop)
    pg  = [_PM.var(pm, nw, :pg, g) for nw in sorted_nw_ids(pm)]

     # bounds on the expectation 
     JuMP.@constraint(pm.model,  pmin <= _PCE.mean(pg, mop))
     JuMP.@constraint(pm.model,  _PCE.mean(pg, mop) <= pmax)
     # chance constraint bounds
     JuMP.@constraint(pm.model,  _PCE.var(pg, T2)
                                 <=
                                ((_PCE.mean(pg, mop) - pmin) / λmin)^2
                   )
     JuMP.@constraint(pm.model,  _PCE.var(pg, T2)
                                 <=
                                 ((pmax - _PCE.mean(pg, mop)) / λmax)^2
                   )
end
""
function constraint_cc_gen_power_imaginary(pm::AbstractACRModel, g, qmin, qmax, λmin, λmax, T2, mop)
    qg  = [_PM.var(pm, nw, :qg, g) for nw in sorted_nw_ids(pm)]

    # bounds on the expectation 
    JuMP.@constraint(pm.model,  qmin <= _PCE.mean(qg, mop))
    JuMP.@constraint(pm.model,  _PCE.mean(qg, mop) <= qmax)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                                <=
                                ((_PCE.mean(qg,mop) - qmin) / λmin)^2
                    )
    JuMP.@constraint(pm.model,  _PCE.var(qg,T2)
                               <=
                                ((qmax - _PCE.mean(qg,mop)) / λmax)^2
                    )
end

# chance constraints unbalanced
## bus
""
function constraint_mc_cc_bus_voltage_magnitude_squared(pm::AbstractUnbalancedACRModel, i, connections, vmin, vmax, λmin, λmax, T2, mop)
    
    for (idx, c) in enumerate(connections)
        # bounds on the expectation
        vms  = [_PMD.var(pm, n, :vms, i)[c] for n in sorted_nw_ids(pm)]
        # JuMP.@constraint(pm.model, vmin[c]^2 <= _PCE.mean(vms[c], mop))
        JuMP.@constraint(pm.model, _PCE.mean(vms, mop) <= vmax[c]^2)
        # chance constraint bounds
        JuMP.@constraint(pm.model,  _PCE.var(vms, T2)
                                    <=
                                ((_PCE.mean(vms, mop) - vmin[c]^2) / λmin)^2
                        )
        JuMP.@constraint(pm.model,  _PCE.var(vms, T2)
                                <=
                                    ((vmax[c]^2 - _PCE.mean(vms, mop)) / λmax)^2
                        )
    end
end

## branch
""
function constraint_mc_cc_branch_series_current_magnitude_squared(pm::AbstractUnbalancedACRModel, b,connections, cmax, λcmax, T2, mop)
    

    if cmax[1]!=Inf
        for (idx, c) in enumerate(connections)
            cmss = [_PMD.var(pm, nw, :cmss, b)[c] for nw in sorted_nw_ids(pm)]
            display("current $(_PCE.var(cmss,T2))")
        # bound on the expectation
            JuMP.@constraint(pm.model,  _PCE.mean(cmss, mop) <= cmax[c]^2)
            # chance constraint bounds
            JuMP.@constraint(pm.model,  _PCE.var(cmss,T2)
                                        <=
                                        ((cmax[c]^2 - _PCE.mean(cmss,mop)) / λcmax)^2
                            )
        end
    end
end

## generator
""
function constraint_mc_cc_gen_power_real(pm::AbstractUnbalancedACRModel, g, connections, pmin, pmax, λmin, λmax, T2, mop)
    
    display("gen")
    for (idx, c) in enumerate(connections)
        pg  = [_PMD.var(pm, nw, :pg, g)[c] for nw in sorted_nw_ids(pm)]
        # bounds on the expectation 
        display(pg[c])
        JuMP.@constraint(pm.model,  pmin[c] <= _PCE.mean(pg, mop))
        JuMP.@constraint(pm.model,  _PCE.mean(pg, mop) <= pmax[c])
        # chance constraint bounds
        JuMP.@constraint(pm.model,  _PCE.var(pg, T2)
                                    <=
                                    ((_PCE.mean(pg, mop) - pmin[c]) / λmin)^2
                    )
        JuMP.@constraint(pm.model,  _PCE.var(pg[c], T2)
                                    <=
                                    ((pmax[c] - _PCE.mean(pg, mop)) / λmax)^2
                    )
    end
end
""
function constraint_mc_cc_gen_power_imaginary(pm::AbstractUnbalancedACRModel, g, connections, qmin, qmax, λmin, λmax, T2, mop)
    qg  = [_PMD.var(pm, nw, :qg, g) for nw in sorted_nw_ids(pm)]
    for (idx, c) in enumerate(connections)
        # bounds on the expectation 
        JuMP.@constraint(pm.model,  qmin[c] <= _PCE.mean(qg[c], mop))
        JuMP.@constraint(pm.model,  _PCE.mean(qg[c], mop) <= qmax[c])
        # chance constraint bounds
        JuMP.@constraint(pm.model,  _PCE.var(qg[c],T2)
                                    <=
                                    ((_PCE.mean(qg[c],mop) - qmin[c]) / λmin)^2
                        )
        JuMP.@constraint(pm.model,  _PCE.var(qg[c],T2)
                                <=
                                    ((qmax[c] - _PCE.mean(qg[c],mop)) / λmax)^2
                        )

    end
end