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
    import InfrastructureModels
    import Ipopt
    import JuMP
    import KernelDensity
    import Memento
    import PolyChaos
    import PowerModels
    using CSV
    using DataFrames
    using JSON
    using PowerModelsDistribution
    # import types
    import PowerModels: AbstractPowerModel, AbstractACRModel, AbstractIVRModel

    # pkgs const
    const _IM = InfrastructureModels
    const _KDE = KernelDensity
    const _PCE = PolyChaos
    const _PM = PowerModels
    const _SPM = StochasticPowerModels
    const _PMD = PowerModelsDistribution

    # memento logger
    function __init__()
        global _LOGGER = Memento.getlogger(PowerModels)
    end

    # const 
    const nw_id_default = 1

    # funct
    sorted_nw_ids(pm) = sort(collect(_PMD.nw_ids(pm)))

    # paths
    const BASE_DIR = dirname(@__DIR__)

    # include
    include("core/constraint.jl")
    include("core/constraint_template.jl")
    include("core/objective.jl")
    include("core/variable.jl")

    include("form/acr.jl")
    include("form/iv.jl")
    include("form/iv_mc.jl")

    include("prob/sopf_acr.jl")
    include("prob/sopf_iv.jl")
    include("prob/sopf_iv_mc.jl")

    include("util/data.jl")
    include("util/util.jl")

    # export
    export BASE_DIR

    export solve_sopf_iv, solve_sopf_acr, solve_sopf_iv_mc, solve_sopf_iv_mc_hc

    export build_stochastic_data
    export build_stochastic_data_mc
    export build_stochastic_data_mc_hc
    export extend_matlab_file
    export pce_coeff, sample, density, print_summary
end 
