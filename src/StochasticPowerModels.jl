################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels(Distribution).jl for Stochastic (Optimal)#
# Power Flow                                                                   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

module StochasticPowerModels

    # import pkgs
    import CSV
    import DataFrames
    import Distributions
    import Ipopt
    import JuMP
    import PolyChaos
    import PowerModels
    import InfrastructureModels
    import Memento

    # import types
    import PowerModels: AbstractPowerModel, AbstractACRModel, AbstractIVRModel
    import PowerModels: comp_start_value, sol_component_value
    import InfrastructureModels: ids, ref, var, con, sol, nw_ids, nws
    import InfrastructureModels: sol_component_value, sol_component_value_edge

    # pkgs const
    const _DFs = DataFrames
    const _DST = Distributions
    const _PCE = PolyChaos
    const _PMs = PowerModels
    const _IMs = InfrastructureModels
    const _SPM = StochasticPowerModels

    # memento logger
    function __init__()
        global _LOGGER = Memento.getlogger(PowerModels)
    end

    # const 
    const nw_id_default = 1

    # funct
    sorted_nw_ids(pm) = sort(collect(_PMs.nw_ids(pm)))

    # paths
    const BASE_DIR = dirname(@__DIR__)

    # include
    include("core/constraint_template.jl")
    include("core/objective.jl")
    include("core/variable.jl")

    include("form/acr.jl")
    include("form/iv.jl")

    include("prob/sopf_acr.jl")
    include("prob/sopf_iv.jl")
    include("prob/sopf_iv_red.jl")
    include("prob/sopf_iv_itr.jl")
    include("prob/sopf_acr_red.jl")
    include("prob/sopf_acr_itr.jl")

    include("util/util.jl")
    # include("util/plot.jl")
    
    include("server/case_5.jl")
    include("server/case_14.jl")
    include("server/case_30.jl")
    include("server/case_57.jl")
    include("server/case_118.jl")

    # export
    export BASE_DIR

    export run_sopf_iv, run_sopf_acr, run_sopf_iv_reduced
    export run_sopf_iv_itr, run_sopf_acr_itr
end 
