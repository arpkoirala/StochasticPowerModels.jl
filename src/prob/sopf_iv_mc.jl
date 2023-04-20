################################################################################
#  2023, Arpan Koirala                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# Multi conductor AC-OPF using IV formulation                                  # 
# NOTE: dc lines are omitted from the current formulation                      #
################################################################################

""
function solve_sopf_iv_mc(data::Dict, model_constructor, optimizer; deg::Int=1, solution_processors=[sol_data_model!], kwargs...)
    @assert _IM.ismultinetwork(data) == false "The data supplied is multinetwork, it should be single-network"
    @assert model_constructor <: _PMD.AbstractUnbalancedIVRModel "This problem type only supports the IVRModel"
    
    sdata = build_stochastic_data_mc(data, deg)
    result = _PMD.solve_mc_model(sdata, model_constructor, optimizer, build_sopf_iv_mc; multinetwork=true, solution_processors=solution_processors, kwargs...)
    result["mop"] = sdata["mop"]
    
    return result
end

""
function solve_sopf_iv_mc(file::String, model_constructor, optimizer; deg::Int=1, solution_processors=[sol_data_model!], kwargs...)
    data = _PM.parse_file(file)
    
    return solve_sopf_iv_mc(data, model_constructor, optimizer; deg=deg, solution_processors=solution_processors, kwargs...)
end

""
function build_sopf_iv_mc(pm::AbstractUnbalancedPowerModel)
    for (n, network) in _PMD.nws(pm) 
        variable_mc_bus_voltage(pm, nw=n)

        variable_mc_branch_current(pm, nw=n)
        # _PMD.variable_mc_switch_current(pm, nw=n)
        # _PMD.variable_mc_transformer_current(pm, nw=n)
        _PMD.variable_mc_generator_power(pm, nw=n, bounded=false)                             # enforcing bounds alters the objective 
        variable_mc_generator_current(pm, nw=n, bounded=false)                           # enforcing bounds makes problem infeasible
        # _PMD.variable_mc_load_power(pm, nw=n, bounded=false)
        # _PMD.variable_mc_storage_power(pm, nw=n, bounded=false)
        variable_mc_load_current(pm, nw=n, bounded=false)
    end

    
    for i in _PMD.ids(pm, :bus, nw=1)
        constraint_mc_cc_bus_voltage_magnitude_squared(pm, i, nw=1)
    end

    for b in _PMD.ids(pm, :branch, nw=1)
        constraint_mc_cc_branch_series_current_magnitude_squared(pm, b, nw=1)
    end

    for g in _PMD.ids(pm, :gen, nw=1)
        constraint_mc_cc_gen_power(pm, g, nw=1)
    end


    for (n, network) in _PMD.nws(pm)
        for i in _PMD.ids(pm, :ref_buses, nw=n)
            constraint_mc_bus_voltage_ref(pm, i, nw=n)
        end

        for g in _PMD.ids(pm, :gen, nw=n)
            constraint_mc_gp_gen_power(pm, g, nw=n; report=true)
        end

        for l in _PMD.ids(pm, :load, nw=n)
            constraint_mc_gp_load_power(pm, l, nw=n; report=true)
        end

        
        for i in _PMD.ids(pm, :bus, nw=n)
            constraint_mc_gp_current_balance(pm, i, nw=n)

            constraint_mc_gp_bus_voltage_magnitude_squared(pm, i, nw=n)
            
        end

        for b in _PMD.ids(pm, :branch, nw=n)
            _PMD.constraint_mc_current_from(pm, b, nw=n)
            _PMD.constraint_mc_current_to(pm, b, nw=n)
            # _PMD.constraint_mc_voltage_angle_difference(pm,b,nw=n)
            _PMD.constraint_mc_bus_voltage_drop(pm, b, nw=n)
             ########
            constraint_mc_gp_branch_series_current_magnitude_squared(pm, b, nw=n)
        end

        

    end

    objective_mc_min_expected_generation_cost(pm)
end