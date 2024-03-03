""
function constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    cut = _PM.ref(pm, nw, :branch_flow_cuts, i)
    branch = _PM.ref(pm, nw, :branch, cut.branch_id)

    if haskey(branch, "rate_c")
        constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm, nw, i, cut.bus_injection, branch["rate_c"])
    end
end

""
function constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    cut = _PM.ref(pm, nw, :branch_flow_cuts, i)
    branch = _PM.ref(pm, nw, :branch, cut.branch_id)

    if haskey(branch, "rate_c")
        constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm, nw, i, cut.bus_injection, branch["rate_c"])
    end
end

""
function constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm::_PM.AbstractPowerModel, n::Int, i::Int, cut_map, rate)
    bus_injection = _PM.var(pm, n, :bus_pg)
    bus_withdrawal = _PM.var(pm, n, :bus_wdp)

    JuMP.@constraint(pm.model, sum(weight*(bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut_map) <= rate + _PM.var(pm, n, :branch_cont_flow_vio, i))
end


""
function constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm::_PM.AbstractPowerModel, n::Int, i::Int, cut_map, rate)
    bus_injection = _PM.var(pm, n, :bus_pg)
    bus_withdrawal = _PM.var(pm, n, :bus_wdp)

    JuMP.@constraint(pm.model, -sum(weight*(bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut_map) <= rate + _PM.var(pm, n, :branch_cont_flow_vio, i))
end

function constraint_power_shunt_dispatch(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus_arcs   = _PM.ref(pm, nw, :bus_arcs, i)
    bus_gens   = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads  = _PM.ref(pm, nw, :bus_loads, i)

    bus_shunts_const = _PM.ref(pm, nw, :bus_shunts_const, i)
    bus_shunts_var = _PM.ref(pm, nw, :bus_shunts_var, i)

    bus_pd = Dict(k => _PM.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PM.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs_const = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts_const)
    bus_bs_const = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts_const)

    constraint_power_shunt_dispatch(pm, nw, i, bus_arcs, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
end

function constraint_power_shunt_dispatch(pm::AbstractACRModel, n::Int, i::Int, bus_arcs, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    p    = _PM.get(_PM.var(pm, n),   :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = _PM.get(_PM.var(pm, n),   :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    
    pg   = _PM.get(_PM.var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = _PM.get(_PM.var(pm, n),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")

    bs = _PM.get(_PM.var(pm, n), :bs, Dict()); _PM._check_var_keys(bs, bus_shunts_var, "reactive power", "shunt")

    JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs_const))*(vr^2 + vi^2)
    )
    JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs_const))*(vr^2 + vi^2)
        + sum(bs[s]*(vi^2 + vr^2) for s in bus_shunts_var)
    )
end

function constraint_ne_power_shunt_dispatch(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus_arcs   = _PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_ne = _PM.ref(pm, nw, :ne_bus_arcs, i)
    bus_gens   = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads  = _PM.ref(pm, nw, :bus_loads, i)

    bus_shunts_const = _PM.ref(pm, nw, :bus_shunts_const, i)
    bus_shunts_var = _PM.ref(pm, nw, :bus_shunts_var, i)

    bus_pd = Dict(k => _PM.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PM.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs_const = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts_const)
    bus_bs_const = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts_const)

    constraint_ne_power_shunt_dispatch(pm, nw, i, bus_arcs, bus_arcs_ne, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
end

function constraint_ne_power_shunt_dispatch(pm::AbstractACRModel, n::Int, i::Int, bus_arcs, bus_arcs_ne, bus_gens, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    p    = _PM.get(_PM.var(pm, n),   :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = _PM.get(_PM.var(pm, n),   :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")

    p_ne = _PM.get(_PM.var(pm, n),   :p_ne, Dict()); _PM._check_var_keys(p_ne, bus_arcs_ne, "active power", "ne_branch")
    q_ne = _PM.get(_PM.var(pm, n),   :q_ne, Dict()); _PM._check_var_keys(q_ne, bus_arcs_ne, "reactive power", "ne_branch")
    
    pg   = _PM.get(_PM.var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = _PM.get(_PM.var(pm, n),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")

    bs = _PM.get(_PM.var(pm, n), :bs, Dict()); _PM._check_var_keys(bs, bus_shunts_var, "reactive power", "shunt")

    JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(p_ne[a] for a in bus_arcs_ne)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs_const))*(vr^2 + vi^2)
    )
    JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_ne[a] for a in bus_arcs_ne)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs_const))*(vr^2 + vi^2)
        + sum(bs[s]*(vi^2 + vr^2) for s in bus_shunts_var)
    )
end

function constraint_branch_ne_voltage(pm::AbstractACRModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :ne_branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    
    vbdr_ne  = _PM.var(pm, nw, :vbdr_ne, i)
    vbdi_ne  = _PM.var(pm, nw, :vbdi_ne, i)
    
    vr_fr = _PM.var(pm, nw, :vr, f_bus)
    vr_to = _PM.var(pm, nw, :vr, t_bus)
    vi_fr = _PM.var(pm, nw, :vi, f_bus)
    vi_to = _PM.var(pm, nw, :vi, t_bus)

    JuMP.@constraint(pm.model,  vbdr_ne == (vr_fr-vr_to))
    JuMP.@constraint(pm.model,  vbdi_ne == (vi_fr-vi_to))         
end

function constraint_gp_power_branch_ne_to(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    branch = _PM.ref(pm, nw, :ne_branch, i)
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

    constraint_gp_power_branch_ne_to(pm, nw, i, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, T2, T3)
end

function constraint_gp_power_branch_ne_to(pm::AbstractACRModel, n::Int, i::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm, T2, T3)
    p_to = _PM.var(pm, n, :p_ne, t_idx)
    q_to = _PM.var(pm, n, :q_ne, t_idx)
    z = _PM.var(pm, n, :branch_ne, i)  
    
    vr_fr = Dict(nw => _PM.var(pm, nw, :vr, f_bus) for nw in _PM.nw_ids(pm))
    vr_to = Dict(nw => _PM.var(pm, nw, :vr, t_bus) for nw in _PM.nw_ids(pm))
    vi_fr = Dict(nw => _PM.var(pm, nw, :vi, f_bus) for nw in _PM.nw_ids(pm))
    vi_to = Dict(nw => _PM.var(pm, nw, :vi, t_bus) for nw in _PM.nw_ids(pm))
   
    JuMP.@constraint(pm.model,  p_to * T2.get([n-1,n-1])
                                ==
                                z * sum(T3.get([n1-1,n2-1,n-1]) *
                                    ((g + g_to) * (vr_to[n1] * vr_to[n2] + vi_to[n1] * vi_to[n2]) + 
                                     (-g * tr - b * ti) / tm^2 * (vr_fr[n1] * vr_to[n2] + vi_fr[n1] * vi_to[n2]) + 
                                     (-b * tr + g * ti) / tm^2 * (-(vi_fr[n1] * vr_to[n2] - vr_fr[n1] * vi_to[n2]))
                                    )
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
   
    JuMP.@constraint(pm.model,  q_to * T2.get([n-1,n-1])
                                ==
                                z * sum(T3.get([n1-1,n2-1,n-1]) *
                                    (-(b + b_to) * (vr_to[n1] * vr_to[n2] + vi_to[n1] * vi_to[n2]) - 
                                     (-b * tr + g * ti) / tm^2 * (vr_fr[n1] * vr_to[n2] + vi_fr[n1] * vi_to[n2]) + 
                                     (-g * tr - b * ti) / tm^2 * (-(vi_fr[n1] * vr_to[n2] - vr_fr[n1] * vi_to[n2]))
                                    )
                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_power_branch_ne_from(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
    branch = _PM.ref(pm, nw, :ne_branch, i)
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

    constraint_gp_power_branch_ne_from(pm, nw, i, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
end

function constraint_gp_power_branch_ne_from(pm::AbstractACRModel, n::Int, i::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
    p_fr = _PM.var(pm, n, :p_ne, f_idx)
    q_fr = _PM.var(pm, n, :q_ne, f_idx)
    z = _PM.var(pm, n, :branch_ne, i)

    vr_fr = Dict(nw => _PM.var(pm, nw, :vr, f_bus) for nw in _PM.nw_ids(pm))
    vr_to = Dict(nw => _PM.var(pm, nw, :vr, t_bus) for nw in _PM.nw_ids(pm))
    vi_fr = Dict(nw => _PM.var(pm, nw, :vi, f_bus) for nw in _PM.nw_ids(pm))
    vi_to = Dict(nw => _PM.var(pm, nw, :vi, t_bus) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  p_fr * T2.get([n-1,n-1])
                                ==
                                z * sum(T3.get([n1-1,n2-1,n-1]) *
                                    ((g + g_fr) / tm^2 * (vr_fr[n1] * vr_fr[n2] + vi_fr[n1] * vi_fr[n2]) + 
                                     (-g * tr + b * ti) / tm^2 * (vr_fr[n1] * vr_to[n2] + vi_fr[n1] * vi_to[n2]) + 
                                     (-b * tr - g * ti) / tm^2 * (vi_fr[n1] * vr_to[n2] - vr_fr[n1] * vi_to[n2])
                                    )
                                for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )

    JuMP.@constraint(pm.model,  q_fr * T2.get([n-1,n-1])
                                ==
                                z * sum(T3.get([n1-1,n2-1,n-1]) *
                                    (-(b + b_fr) / tm^2 * (vr_fr[n1] * vr_fr[n2] + vi_fr[n1] * vi_fr[n2]) - 
                                     (-b * tr - g * ti) / tm^2 * (vr_fr[n1] * vr_to[n2] + vi_fr[n1] * vi_to[n2]) + 
                                     (-g * tr + b * ti) / tm^2 * (vi_fr[n1] * vr_to[n2] - vr_fr[n1] * vi_to[n2]) 
                                    )
                                for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_gp_branch_ne_series_current_magnitude_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    T2  = pm.data["T2"]
    T3  = pm.data["T3"]

    constraint_gp_branch_ne_series_current_magnitude_squared(pm, nw, i, T2, T3)
end

function constraint_gp_branch_ne_series_current_magnitude_squared(pm::AbstractACRModel, n::Int, i, T2, T3)
    cmss = _PM.var(pm, n, :cmss_ne, i)
    z = _PM.var(pm, n, :branch_ne, i)
    
    branch = _PM.ref(pm, n, :ne_branch, i)
    g, b  = _PM.calc_branch_y(branch)
    
    vbdr = Dict(nw => _PM.var(pm, nw, :vbdr_ne, i) for nw in _PM.nw_ids(pm))
    vbdi = Dict(nw => _PM.var(pm, nw, :vbdi_ne, i) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  T2.get([n-1,n-1]) * cmss
                                ==
                                z * (g^2 + b^2) * sum(  T3.get([n1-1,n2-1,n-1]) * 
                                                    (vbdr[n1] * vbdr[n2] + vbdi[n1] * vbdi[n2])
                                                    for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

function constraint_cc_branch_ne_series_current_magnitude_squared(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    cmax = _PM.ref(pm, nw, :ne_branch, i, "cmax")
    λmax = _PM.ref(pm, nw, :ne_branch, i, "λcmax")
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    constraint_cc_branch_ne_series_current_magnitude_squared(pm, i, cmax, λmax, T2, mop)
end

function constraint_cc_branch_ne_series_current_magnitude_squared(pm::AbstractACRModel, i, cmax, λcmax, T2, mop)
    cmss = [_PM.var(pm, nw, :cmss_ne, i) for nw in sorted_nw_ids(pm)]
    z = _PM.var(pm, nw=1, :branch_ne, i)

    # bound on the expectation
    JuMP.@constraint(pm.model,  _PCE.mean(cmss, mop) <= z * cmax^2)
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(cmss,T2)
                                <=
                                z * ((cmax^2 - _PCE.mean(cmss,mop)) / λcmax)^2
                    )
end

function constraint_gp_goc_power_branch_from(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    
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

    if branch["transformer"]
        constraint_gp_goc_power_branch_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
    else
        constraint_gp_power_branch_from(pm, nw, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
    end
end


function constraint_gp_goc_power_branch_from(pm::AbstractACRModel, n::Int,f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm, T2, T3)
    p_fr = _PM.var(pm, n, :p, f_idx)
    q_fr = _PM.var(pm, n, :q, f_idx)

    vr_fr = Dict(nw => _PM.var(pm, nw, :vr, f_bus) for nw in _PM.nw_ids(pm))
    vr_to = Dict(nw => _PM.var(pm, nw, :vr, t_bus) for nw in _PM.nw_ids(pm))
    vi_fr = Dict(nw => _PM.var(pm, nw, :vi, f_bus) for nw in _PM.nw_ids(pm))
    vi_to = Dict(nw => _PM.var(pm, nw, :vi, t_bus) for nw in _PM.nw_ids(pm))

    JuMP.@constraint(pm.model,  p_fr * T2.get([n-1,n-1])
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    ((g/tm^2 + g_fr) * (vr_fr[n1] * vr_fr[n2] + vi_fr[n1] * vi_fr[n2]) + 
                                     (-g * tr + b * ti) / tm^2 * (vr_fr[n1] * vr_to[n2] + vi_fr[n1] * vi_to[n2]) + 
                                     (-b * tr - g * ti) / tm^2 * (vi_fr[n1] * vr_to[n2] - vr_fr[n1] * vi_to[n2])
                                    )
                                for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )

    JuMP.@constraint(pm.model,  q_fr * T2.get([n-1,n-1])
                                ==
                                sum(T3.get([n1-1,n2-1,n-1]) *
                                    (-(b/tm^2 + b_fr) * (vr_fr[n1] * vr_fr[n2] + vi_fr[n1] * vi_fr[n2]) - 
                                     (-b * tr - g * ti) / tm^2 * (vr_fr[n1] * vr_to[n2] + vi_fr[n1] * vi_to[n2]) + 
                                     (-g * tr + b * ti) / tm^2 * (vi_fr[n1] * vr_to[n2] - vr_fr[n1] * vi_to[n2]) 
                                    )
                                for n1 in _PM.nw_ids(pm), n2 in _PM.nw_ids(pm))
                    )
end

