################################################################################
#  Copyright 2021, Tom Van Acker, Arpan Koirala                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# input data
""
function parse_dst(dst, pa, pb, deg)
    dst == "Beta"    && return _PCE.Beta01OrthoPoly(deg, pa, pb; Nrec=5*deg)
    dst == "Normal"  && return _PCE.GaussOrthoPoly(deg; Nrec=5*deg)
    dst == "Uniform" && return _PCE.Uniform01OrthoPoly(deg; Nrec=5*deg)
end

"""
    StochasticPowerModels.build_stochastic_data(data::Dict{String,Any}, deg::Int)
"""
function build_stochastic_data(data::Dict{String,Any}, deg::Int)
    # add maximum current
    for (nb, branch) in data["branch"]
        f_bus = branch["f_bus"]
        branch["cmax"] = branch["rate_a"] / data["bus"]["$f_bus"]["vmin"]
    end

    # build mop
    opq = [parse_dst(ns[2]["dst"], ns[2]["pa"], ns[2]["pb"], deg) for ns in data["sdata"]]
    mop = _PCE.MultiOrthoPoly(opq, deg)

    # build load matrix
    Nd, Npce = length(data["load"]), mop.dim
    pd, qd = zeros(Nd, Npce), zeros(Nd, Npce)
    for nd in 1:Nd 
        # reactive power
        qd[nd,1] = data["load"]["$nd"]["qd"]
        # active power
        nb = data["load"]["$nd"]["load_bus"]
        ni = data["bus"]["$nb"]["dst_id"]
        if ni == 0
            pd[nd,1] = data["load"]["$nd"]["pd"]
        else
            base = data["baseMVA"]
            μ, σ = data["bus"]["$nb"]["μ"] / base, data["bus"]["$nb"]["σ"] / base
            if mop.uni[ni] isa _PCE.GaussOrthoPoly
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni])
            else
                pd[nd,[1,ni+1]] = _PCE.convert2affinePCE(μ, σ, mop.uni[ni], kind="μσ")
            end
        end
    end

    # replicate the data
    data = _PMs.replicate(data, Npce)

    # add the stochastic data 
    data["T2"] = _PCE.Tensor(2,mop)
    data["T3"] = _PCE.Tensor(3,mop)
    data["T4"] = _PCE.Tensor(4,mop)
    data["mop"] = mop
    for nw in 1:Npce, nd in 1:Nd
        data["nw"]["$nw"]["load"]["$nd"]["pd"] = pd[nd,nw]
        data["nw"]["$nw"]["load"]["$nd"]["qd"] = qd[nd,nw]
    end

    return data
end

# output data
""
function sample(sdata, result, element::String, id::Int, var::String; sample_size::Int=1000)
    coeff = [nw[element]["$id"][var] for nw in values(result["solution"]["nw"])]
    return _PCE.samplePCE(sample_size, coeff, sdata["mop"])
 end
 
 ""
 function density(sdata, result, element::String, id::Int, var::String; sample_size::Int=1000)
    coeff = [nw[element]["$id"][var] for nw in values(result["solution"]["nw"])]
    return _KDE.kde(_PCE.samplePCE(sample_size, coeff, sdata["mop"]))
 end
 