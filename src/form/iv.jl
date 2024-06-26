################################################################################
#  Copyright 2021, Tom Van Acker, Frederik Geth                                #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# variable
## branch 
""
function variable_branch_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_branch_series_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_branch_series_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)

    expression_variable_branch_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    expression_variable_branch_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    
    variable_branch_series_current_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
"variable: `cr[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_variable_branch_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cr = _PM.var(pm, nw)[:cr] = Dict()

    bus = _PM.ref(pm, nw, :bus)
    branch = _PM.ref(pm, nw, :branch)

    for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
        b = branch[l]
        tm = b["tap"]
        tr, ti = _PM.calc_branch_t(b)
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]

        vr_fr = _PM.var(pm, nw, :vr, i)
        vi_fr = _PM.var(pm, nw, :vi, i)
    
        vr_to = _PM.var(pm, nw, :vr, j)
        vi_to = _PM.var(pm, nw, :vi, j)
    
        csr_fr = _PM.var(pm, nw, :csr, l)
        csi_fr = _PM.var(pm, nw, :csi, l)

        cr[(l,i,j)] = (tr * csr_fr - ti * csi_fr + g_sh_fr * vr_fr - b_sh_fr * vi_fr) / tm^2
        cr[(l,j,i)] = -csr_fr + g_sh_to * vr_to - b_sh_to * vi_to

        # ub = Inf
        # if haskey(b, "rate_a")
        #     rate_fr = b["rate_a"]*b["tap"]
        #     rate_to = b["rate_a"]
        #     ub = max(rate_fr/bus[i]["vmin"], rate_to/bus[j]["vmin"])
        # end
        # if haskey(b, "c_rating_a")
        #     ub = b["c_rating_a"]
        # end

        # if !isinf(ub)
        #     JuMP.@constraint(pm.model, cr[(l,i,j)] >= -ub)
        #     JuMP.@constraint(pm.model, cr[(l,i,j)] <= ub)

        #     JuMP.@constraint(pm.model, cr[(l,j,i)] >= -ub)
        #     JuMP.@constraint(pm.model, cr[(l,j,i)] <= ub)
        # end
    end

    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branch, :cr_fr, :cr_to, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), cr)
end
"variable: `ci[l,i,j]` for `(l,i,j)` in `arcs`"
function expression_variable_branch_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    ci = _PM.var(pm, nw)[:ci] = Dict()

    bus = _PM.ref(pm, nw, :bus)
    branch = _PM.ref(pm, nw, :branch)

    for (l,i,j) in _PM.ref(pm, nw, :arcs_from)
        b = branch[l]
        tm = b["tap"]
        tr, ti = _PM.calc_branch_t(b)
        g_sh_fr, b_sh_fr = b["g_fr"], b["b_fr"]
        g_sh_to, b_sh_to = b["g_to"], b["b_to"]

        vr_fr = _PM.var(pm, nw, :vr, i)
        vi_fr = _PM.var(pm, nw, :vi, i)
    
        vr_to = _PM.var(pm, nw, :vr, j)
        vi_to = _PM.var(pm, nw, :vi, j)
    
        csr_fr = _PM.var(pm, nw, :csr, l)
        csi_fr = _PM.var(pm, nw, :csi, l)

        ci[(l,i,j)] = (tr * csi_fr + ti * csr_fr + g_sh_fr * vi_fr + b_sh_fr * vr_fr) / tm^2
        ci[(l,j,i)] = -csi_fr + g_sh_to * vi_to + b_sh_to * vr_to

        # ub = Inf
        # if haskey(b, "rate_a")
        #     rate_fr = b["rate_a"]*b["tap"]
        #     rate_to = b["rate_a"]
        #     ub = max(rate_fr/bus[i]["vmin"], rate_to/bus[j]["vmin"])
        # end
        # if haskey(b, "c_rating_a")
        #     ub = b["c_rating_a"]
        # end

        # if !isinf(ub)
        #     JuMP.@constraint(pm.model, ci[(l,i,j)] >= -ub)
        #     JuMP.@constraint(pm.model, ci[(l,i,j)] <= ub)

        #     JuMP.@constraint(pm.model, ci[(l,j,i)] >= -ub)
        #     JuMP.@constraint(pm.model, ci[(l,j,i)] <= ub)
        # end
    end

    report && _IM.sol_component_value_edge(pm, _PM.pm_it_sym, nw, :branch, :ci_fr, :ci_to, _PM.ref(pm, nw, :arcs_from), _PM.ref(pm, nw, :arcs_to), ci)
end

## load
""
function variable_load_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_load_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_load_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end
"variable: `crd[j]` for `j` in `load`"
function variable_load_current_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    crd = _PM.var(pm, nw)[:crd] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_crd",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "crd_start")
    )

    report && _PM.sol_component_value(pm, nw, :load, :crd, _PM.ids(pm, nw, :load), crd)
end
"variable: `cid[j]` for `j` in `load`"
function variable_load_current_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cid = _PM.var(pm, nw)[:cid] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_cid",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "cid_start")
    )

    report && _PM.sol_component_value(pm, nw, :load, :cid, _PM.ids(pm, nw, :load), cid)
end


## generator
""
function variable_gen_current(pm::AbstractIVRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_gen_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_gen_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

# general constraints
## bus
""
function constraint_current_balance(pm::AbstractIVRModel, n::Int, i, bus_arcs, bus_gens, bus_loads, bus_gs, bus_bs)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    cr = _PM.var(pm, n, :cr)
    ci = _PM.var(pm, n, :ci)

    crd = _PM.var(pm, n, :crd)
    cid = _PM.var(pm, n, :cid)
    crg = _PM.var(pm, n, :crg)
    cig = _PM.var(pm, n, :cig)

    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

# galerkin projection constraint
## branch
""
function constraint_gp_branch_series_current_magnitude_squared(pm::AbstractIVRModel, n::Int, i, T2, T3)
    cmss  = _PM.var(pm, n, :cmss, i)
    csr = Dict(nw => _PM.var(pm, nw, :csr, i) for nw in _PM.nw_ids(pm))
    csi = Dict(nw => _PM.var(pm, nw, :csi, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * cmss
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (csr[n1] * csr[n2] + csi[n1] * csi[n2]) 
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

## generator
""
function constraint_gp_gen_power_real(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))
    
    crg = Dict(nw => _PM.var(pm, nw, :crg, g) for nw in _PM.nw_ids(pm))
    cig = Dict(nw => _PM.var(pm, nw, :cig, g) for nw in _PM.nw_ids(pm))

    pg  = _PM.var(pm, n, :pg, g)
    
    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * pg
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) * 
                                    (vr[n1] * crg[n2] + vi[n1] * cig[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end
""
function constraint_gp_gen_power_imaginary(pm::AbstractIVRModel, n::Int, i, g, T2, T3)
    vr  = Dict(nw => _PM.var(pm, nw, :vr, i) for nw in _PM.nw_ids(pm))
    vi  = Dict(nw => _PM.var(pm, nw, :vi, i) for nw in _PM.nw_ids(pm))
    
    crg = Dict(nw => _PM.var(pm, nw, :crg, g) for nw in _PM.nw_ids(pm))
    cig = Dict(nw => _PM.var(pm, nw, :cig, g) for nw in _PM.nw_ids(pm))

    qg  = _PM.var(pm, n, :qg, g)
    
    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * qg
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (vi[n1] * crg[n2] - vr[n1] * cig[n2])
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_gen_power_n(pm::AbstractPowerModel, n::Int, i, g, T2, T3)

    vms = _PM.var(pm, n, :vms, i)
    
    crg = Dict(nw => _PM.var(pm, nw, :crg, g) for nw in _PM.nw_ids(pm))
    cig = Dict(nw => _PM.var(pm, nw, :cig, g) for nw in _PM.nw_ids(pm))
    cmsg = _PM.var(pm, n, :cmsg, g)

    pg  = _PM.var(pm, n, :pg, g)
    qg  = _PM.var(pm, n, :qg, g)    

    JuMP.@constraint(pm.model, pg^2 + qg^2 >= vms * cmsg)
    # JuMP.@constraint(pm.model, [(vms + cmsg), (vms - cmsg), 2*pg, 2*qg] in JuMP.SecondOrderCone())
    # JuMP.@constraint(pm.model, [vms, cmsg, pg, qg] in JuMP.RotatedSecondOrderCone())
    # JuMP.@constraint(pm.model, [vms*cmsg+1, 2*pg, 2*qg, vms*cmsg-1] in JuMP.SecondOrderCone())

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * cmsg 
                    ==
                        sum(T3.get([n1-1,n2-1,n-1]) * 
                            (crg[n1] * crg[n2] + cig[n1] * cig[n2]) 
                            for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end



## load
""
function constraint_gp_load_power_real(pm::AbstractIVRModel, n::Int, i, l, pd, T2, T3)
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
function constraint_gp_load_power_imaginary(pm::AbstractIVRModel, n::Int, i, l, qd, T2, T3)
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

function constraint_gp_load_power_n(pm::AbstractIVRModel, n::Int, i, l, pd, qd, T2, T3)

    vms = _PM.var(pm, n, :vi, i)

    crd = Dict(nw => _PM.var(pm, nw, :crd, l) for nw in _PM.nw_ids(pm))
    cid = Dict(nw => _PM.var(pm, nw, :cid, l) for nw in _PM.nw_ids(pm))
    cmsd = _PM.var(pm, n, :cmsd, l)

    JuMP.@constraint(pm.model, pd^2 + qd^2 >= vms * cmsd)
    # JuMP.@constraint(pm.model, [(vms + cmsd), (vms - cmsd), 2*pd, 2*qd] in JuMP.SecondOrderCone())

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * cmsd 
    ==
        sum(T3.get([n1-1,n2-1,n-1]) * 
            (crd[n1] * crd[n2] + cid[n1] * cid[n2]) 
            for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
    )

end




# solution
""
function sol_data_model!(pm::AbstractIVRModel, solution::Dict)
    _PM.apply_pm!(_sol_data_model_ivr!, solution)
end

""
function _sol_data_model_ivr!(solution::Dict)
    if haskey(solution, "bus")
        for (i, bus) in solution["bus"]
            if haskey(bus, "vr") && haskey(bus, "vi")
                bus["vm"] = hypot(bus["vr"], bus["vi"])
                bus["va"] = atan(bus["vi"], bus["vr"])
            end
        end
    end
end