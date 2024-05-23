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


function constraint_cc_branch_contingency_ptdf_thermal_limit_from_soft(pm::_PM.AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    cut = _PM.ref(pm, nw, :branch_flow_cuts, i)
    branch = _PM.ref(pm, nw, :branch, cut.branch_id)
    cmax = branch["cmax"]
    λmax = branch["λcmax"]
    
    T2  = pm.data["T2"]
    mop = pm.data["mop"]

    if haskey(branch, "cmax")
        constraint_cc_branch_contingency_ptdf_thermal_limit_from_soft(pm, nw, i, cut.bus_injection, cmax, λmax, T2, mop)
    end
end


function constraint_cc_branch_contingency_ptdf_thermal_limit_from_soft(pm::_PM.AbstractPowerModel, n::Int, i::Int, cut_map, cmax, λcmax, T2, mop)
    bus_injection = _PM.var(pm, n, :bus_pg)
    bus_withdrawal = _PM.var(pm, n, :bus_wdp)
    # cmss = [_PM.var(pm, nw, :cmss, b) for nw in sorted_nw_ids(pm)]
    cms = [sum(weight*(_PM.var(pm, n, :bus_pg)[bus_id] - _PM.var(pm, n, :bus_wdp)[bus_id]) for (bus_id, weight) in cut_map) for n in sorted_nw_ids(pm)]

    # bound on the expectation
    JuMP.@constraint(pm.model,  _PCE.mean(cms, mop) <= cmax + _PM.var(pm, n, :branch_cont_flow_vio, i))
    # chance constraint bounds
    JuMP.@constraint(pm.model,  _PCE.var(cms,T2)
                                    <=
                                    ((cmax - _PCE.mean(cms,mop)) / λcmax)
                    )
    # JuMP.@constraint(pm.model, sum(weight*(bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut_map) <= cmax + _PM.var(pm, n, :branch_cont_flow_vio, i))
end




# sub scopf

function constraint_gen_real_setpoint_link(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    pg_master = _PM.ref(pm, nw, :gen, i)["pg"]

    constraint_gen_real_setpoint_link(pm, nw, i, pg_master)
end
function constraint_gen_real_setpoint_link(pm::_PM.AbstractPowerModel, n::Int, i::Int, pg_master)
    pg = _PM.var(pm, n, :pg, i)

    dual_p = JuMP.@constraint(pm.model,  pg == pg_master)

    _PM.sol(pm, n, :gen, i)[:lm_p] = dual_p
end
function constraint_gen_reactive_setpoint_link(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    qg_master = _PM.ref(pm, nw, :gen, i)["qg"]

    constraint_gen_reactive_setpoint_link(pm, nw, i, qg_master)
end
function constraint_gen_reactive_setpoint_link(pm::_PM.AbstractPowerModel, n::Int, i::Int, qg_master)
    qg = _PM.var(pm, n, :qg, i)

    dual_q = JuMP.@constraint(pm.model,  qg == qg_master)

    _PM.sol(pm, n, :gen, i)[:lm_q] = dual_q 
end

 
function constraint_conv_real_setpoint_link(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    pconv_ac_master = -_PM.ref(pm, nw, :convdc, i)["P_g"]  # P_g

    constraint_conv_real_setpoint_link(pm, nw, i, pconv_ac_master)
end
function constraint_conv_real_setpoint_link(pm::_PM.AbstractPowerModel, n::Int, i::Int, pconv_ac_master)
    pconv_ac = _PM.var(pm, n, :pconv_ac, i)

    dual_ac_p = JuMP.@constraint(pm.model,  pconv_ac == pconv_ac_master)

    _PM.sol(pm, n, :convdc, i)[:lm_ac_p] = dual_ac_p
end
function constraint_conv_reactive_setpoint_link(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    qconv_ac_master = -_PM.ref(pm, nw, :convdc, i)["Q_g"]  # Q_g

    constraint_conv_reactive_setpoint_link(pm, nw, i, qconv_ac_master)
end
function constraint_conv_reactive_setpoint_link(pm::_PM.AbstractPowerModel, n::Int, i::Int, qconv_ac_master)
    qconv_ac = _PM.var(pm, n, :qconv_ac, i)

    dual_ac_q = JuMP.@constraint(pm.model,  qconv_ac == qconv_ac_master)

    _PM.sol(pm, n, :convdc, i)[:lm_ac_q] = dual_ac_q
end
##
function constraint_thermal_limit_from_soft(pm:: _PM.AbstractPowerModel, i::Int; nw::Int= _PM.nw_id_default)
    branch =  _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        constraint_thermal_limit_from_soft(pm, nw, i, f_idx, branch["rate_a"])
    end
end
function constraint_thermal_limit_to_soft(pm:: _PM.AbstractPowerModel, i::Int; nw::Int= _PM.nw_id_default)
    branch =  _PM.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    if haskey(branch, "rate_a")
        constraint_thermal_limit_to_soft(pm, nw, i, t_idx, branch["rate_a"])
    end
end

function constraint_thermal_limit_from_soft(pm::_PM.AbstractPowerModel, n::Int, i::Int, f_idx, rate_a)
    p_fr = _PM.var(pm, n, :p, f_idx)
    q_fr = _PM.var(pm, n, :q, f_idx)
    vio_fr = _PM.var(pm, n, :bf_vio_fr, i)

    JuMP.@constraint(pm.model, p_fr^2 + q_fr^2 <= rate_a^2 + vio_fr)
end

function constraint_thermal_limit_to_soft(pm::_PM.AbstractPowerModel, n::Int, i::Int, t_idx, rate_a)
    p_to = _PM.var(pm, n, :p, t_idx)
    q_to = _PM.var(pm, n, :q, t_idx)
    vio_to = _PM.var(pm, n, :bf_vio_to, i)

    JuMP.@constraint(pm.model, p_to^2 + q_to^2 <= rate_a^2 + vio_to)
end

function constraint_benders_cut(pm:: _PM.AbstractPowerModel, i::Int; nw::Int =_PM.nw_id_default)
    pg_master = Dict(i => _PM.ref(pm, nw, :gen, i)["pg"] for (i,gen) in _PM.ref(pm, nw, :gen))
    qg_master = Dict(i => _PM.ref(pm, nw, :gen, i)["qg"] for (i,gen) in _PM.ref(pm, nw, :gen))
    pconv_ac_master = Dict(i => -_PM.ref(pm, nw, :convdc, i)["P_g"] for (i,convdc) in _PM.ref(pm, nw, :convdc))      
    qconv_ac_master = Dict(i => -_PM.ref(pm, nw, :convdc, i)["Q_g"] for (i,convdc) in _PM.ref(pm, nw, :convdc)) 
    cut = _PM.ref(pm, nw, :cuts, i)
    sub_obj = cut.obj
    gen_lm_p = cut.lm_p
    gen_lm_q = cut.lm_q
    conv_lm_p = cut.lm_ac_p
    conv_lm_q = cut.lm_ac_q

    constraint_benders_cut(pm, nw, i, pg_master, qg_master, pconv_ac_master, qconv_ac_master, sub_obj, gen_lm_p, gen_lm_q, conv_lm_p, conv_lm_q)
end

function constraint_benders_cut(pm:: _PM.AbstractPowerModel, n::Int, i::Int, pg_master, qg_master, pconv_ac_master, qconv_ac_master, sub_obj, gen_lm_p, gen_lm_q, conv_lm_p, conv_lm_q)
    pg = _PM.var(pm, n, :pg)
    qg = _PM.var(pm, n, :qg)
    pconv_ac = _PM.var(pm, n, :pconv_ac)
    qconv_ac = _PM.var(pm, n, :qconv_ac)

    JuMP.@constraint(pm.model, sub_obj + sum(val * (pg[i] - pg_master[i]) for (i,val) in gen_lm_p) + sum(val * (qg[i] - qg_master[i]) for (i,val) in gen_lm_q) + sum(val * (pconv_ac[i] - pconv_ac_master[i]) for (i, val) in conv_lm_p) <= 0)

end

function constraint_benders_cut_soft(pm:: _PM.AbstractPowerModel, i::Int; nw::Int =_PM.nw_id_default)
    pg_master = Dict(i => _PM.ref(pm, nw, :gen, i)["pg"] for (i,gen) in _PM.ref(pm, nw, :gen))
    qg_master = Dict(i => _PM.ref(pm, nw, :gen, i)["qg"] for (i,gen) in _PM.ref(pm, nw, :gen))
    pconv_ac_master = Dict(i => -_PM.ref(pm, nw, :convdc, i)["P_g"] for (i,convdc) in _PM.ref(pm, nw, :convdc))      
    qconv_ac_master = Dict(i => -_PM.ref(pm, nw, :convdc, i)["Q_g"] for (i,convdc) in _PM.ref(pm, nw, :convdc)) 
    cut = _PM.ref(pm, nw, :cuts, i)
    sub_obj = cut.obj
    gen_lm_p = cut.lm_p
    gen_lm_q = cut.lm_q
    conv_lm_p = cut.lm_ac_p
    conv_lm_q = cut.lm_ac_q

    constraint_benders_cut_soft(pm, nw, i, pg_master, qg_master, pconv_ac_master, qconv_ac_master, sub_obj, gen_lm_p, gen_lm_q, conv_lm_p, conv_lm_q)
end

function constraint_benders_cut_soft(pm:: _PM.AbstractPowerModel, n::Int, i::Int, pg_master, qg_master, pconv_ac_master, qconv_ac_master, sub_obj, gen_lm_p, gen_lm_q, conv_lm_p, conv_lm_q)
    pg = _PM.var(pm, n, :pg)
    qg = _PM.var(pm, n, :qg)
    pconv_ac = _PM.var(pm, n, :pconv_ac)
    qconv_ac = _PM.var(pm, n, :qconv_ac)

    pg_cost = _PM.var(pm, :pg_cost)
 
    obj_expr = JuMP.@expression(pm.model, 
        sum( pg_cost[i] for (i, gen) in _PM.ref(pm, :gen) ))

    JuMP.@constraint(pm.model, obj_expr >= sub_obj + sum(val * (pg[i] - pg_master[i]) for (i,val) in gen_lm_p) + sum(val * (qg[i] - qg_master[i]) for (i,val) in gen_lm_q) + sum(val * (pconv_ac[i] - pconv_ac_master[i]) for (i, val) in conv_lm_p))

end




function constraint_power_balance_ac_soft(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = _PM.ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_convs_ac = _PM.ref(pm, nw, :bus_convs_ac, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)

    pd = Dict(k => _PM.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    qd = Dict(k => _PM.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance_ac_soft(pm, nw, i, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_shunts, pd, qd, gs, bs)
end

function constraint_power_balance_ac_soft(pm::_PM.AbstractWModels, n::Int,  i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_shunts, pd, qd, gs, bs)
    w = _PM.var(pm, n, :w, i)
    p = _PM.var(pm, n, :p)
    q = _PM.var(pm, n, :q)
    pg = _PM.var(pm, n, :pg)
    qg = _PM.var(pm, n, :qg)
    pconv_grid_ac = _PM.var(pm, n, :pconv_tf_fr)
    qconv_grid_ac = _PM.var(pm, n, :qconv_tf_fr)

    pb_ac_pos_vio = _PM.var(pm, n, :pb_ac_pos_vio, i)
    qb_ac_pos_vio = _PM.var(pm, n, :qb_ac_pos_vio, i)
    pb_ac_neg_vio = _PM.var(pm, n, :pb_ac_neg_vio, i)
    qb_ac_neg_vio = _PM.var(pm, n, :qb_ac_neg_vio, i)

    cstr_p = JuMP.@constraint(pm.model, pb_ac_pos_vio - pb_ac_neg_vio + sum(p[a] for a in bus_arcs) + sum(pconv_grid_ac[c] for c in bus_convs_ac)  == sum(pg[g] for g in bus_gens)  - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*w)
    cstr_q = JuMP.@constraint(pm.model, qb_ac_pos_vio - qb_ac_neg_vio + sum(q[a] for a in bus_arcs) + sum(qconv_grid_ac[c] for c in bus_convs_ac)  == sum(qg[g] for g in bus_gens)  - sum(qd[d] for d in bus_loads) + sum(bs[s] for s in bus_shunts)*w)

    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

function constraint_power_balance_ac_soft(pm::_PM.AbstractACPModel, n::Int,  i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_shunts, pd, qd, gs, bs)
    vm = _PM.var(pm, n,  :vm, i)
    p = _PM.var(pm, n,  :p)
    q = _PM.var(pm, n,  :q)
    pg = _PM.var(pm, n,  :pg)
    qg = _PM.var(pm, n,  :qg)
    pconv_grid_ac = _PM.var(pm, n,  :pconv_tf_fr)
    qconv_grid_ac = _PM.var(pm, n,  :qconv_tf_fr)

    pb_ac_pos_vio = _PM.var(pm, n, :pb_ac_pos_vio, i)
    qb_ac_pos_vio = _PM.var(pm, n, :qb_ac_pos_vio, i)
    pb_ac_neg_vio = _PM.var(pm, n, :pb_ac_neg_vio, i)
    qb_ac_neg_vio = _PM.var(pm, n, :qb_ac_neg_vio, i)

    cstr_p = JuMP.@NLconstraint(pm.model, pb_ac_pos_vio - pb_ac_neg_vio + sum(p[a] for a in bus_arcs) + sum(pconv_grid_ac[c] for c in bus_convs_ac)  == sum(pg[g] for g in bus_gens)   - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*vm^2)
    cstr_q = JuMP.@NLconstraint(pm.model, qb_ac_pos_vio - qb_ac_neg_vio + sum(q[a] for a in bus_arcs) + sum(qconv_grid_ac[c] for c in bus_convs_ac)  == sum(qg[g] for g in bus_gens)  - sum(qd[d] for d in bus_loads) + sum(bs[s] for s in bus_shunts)*vm^2)

    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end



# deterministic IVRModel

function variable_gen_current_n(pm::_PM.AbstractIVRModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    _PM.variable_gen_current_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_gen_current_imaginary(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    _PM.variable_gen_power(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_bus_voltage_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_gen_current_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_auxillary_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    

    # store active and reactive power expressions for use in objective + post processing

    for (i, bus) in _PM.ref(pm, nw, :bus)
        vr = _PM.var(pm, nw, :vr, i)
        vi = _PM.var(pm, nw, :vi, i)
        vms = _PM.var(pm, nw, :vms, i)
        vmsr = _PM.var(pm, nw, :vmsr, i)
        vmsi = _PM.var(pm, nw, :vmsi, i)

        JuMP.@constraint(pm.model, vmsr >= vr^2)
        JuMP.@constraint(pm.model, vmsi >= vi^2)
        JuMP.@constraint(pm.model, vms == vmsr + vmsi)
    end

    for (i, gen) in _PM.ref(pm, nw, :gen)
        crg = _PM.var(pm, nw, :crg, i)
        cig = _PM.var(pm, nw, :cig, i)
        cmsg = _PM.var(pm, nw, :cmsg, i)
        cmsgr = _PM.var(pm, nw, :cmsgr, i)
        cmsgi = _PM.var(pm, nw, :cmsgi, i)

        JuMP.@constraint(pm.model, cmsgr >= crg^2)
        JuMP.@constraint(pm.model, cmsgi >= cig^2)
        JuMP.@constraint(pm.model, cmsg == cmsgr + cmsgi)
    end

    for (i, gen) in _PM.ref(pm, nw, :gen)
        busid = gen["gen_bus"]
        pg = _PM.var(pm, nw, :pg, i)
        qg = _PM.var(pm, nw, :qg, i)
        vms = _PM.var(pm, nw, :vms, busid)
        cmsg = _PM.var(pm, nw, :cmsg, i)

        JuMP.@constraint(pm.model, pg^2 + qg^2 >= vms * cmsg)
    end

end

function variable_gen_current_magnitude_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    cmsg = _PM.var(pm, nw)[:cmsg] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_cmsg",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "cmsg_start")
    )

    if bounded
        bus = _PM.ref(pm, nw, :bus)
        for (i, g) in _PM.ref(pm, nw, :gen)
            vmin = bus[g["gen_bus"]]["vmin"]
            @assert vmin > 0
            s = sqrt(max(abs(g["pmax"]), abs(g["pmin"]))^2 + max(abs(g["qmax"]), abs(g["qmin"]))^2)
            ub = s/vmin

            JuMP.set_lower_bound(cmsg[i], -ub^2)
            JuMP.set_upper_bound(cmsg[i],  ub^2)
        end
    end

    report && _PM.sol_component_value(pm, nw, :gen, :cmsg, _PM.ids(pm, nw, :gen), cmsg)
end

function variable_load_current_magnitude_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    cmsd = _PM.var(pm, nw)[:cmsd] = JuMP.@variable(pm.model,
        [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_cmsd",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "cmsd_start")
    )

    if bounded
        bus = _PM.ref(pm, nw, :bus)
        for (i, l) in _PM.ref(pm, nw, :load)
            vmin = bus[l["load_bus"]]["vmin"]
            @assert vmin > 0
            s = sqrt(abs(l["pd"])^2 + abs(l["qd"])^2)
            ub = s/vmin

            JuMP.set_lower_bound(cmsd[i], -ub^2)
            JuMP.set_upper_bound(cmsd[i],  ub^2)
        end
    end

    report && _PM.sol_component_value(pm, nw, :load, :cmsd, _PM.ids(pm, nw, :load), cmsd)
end

function variable_auxillary_squared(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    cmsgr = _PM.var(pm, nw)[:cmsgr] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_cmsgr",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "cmsgr_start")
    )

    cmsgi = _PM.var(pm, nw)[:cmsgi] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :gen)], base_name="$(nw)_cmsgi",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :gen, i), "cmsgi_start")
    )

    if bounded
        bus = _PM.ref(pm, nw, :bus)
        for (i, g) in _PM.ref(pm, nw, :gen)
            vmin = bus[g["gen_bus"]]["vmin"]
            @assert vmin > 0
            s = sqrt(max(abs(g["pmax"]), abs(g["pmin"]))^2 + max(abs(g["qmax"]), abs(g["qmin"]))^2)
            ub = s/vmin

            JuMP.@constraint(pm.model, cmsgr[i] + cmsgi[i] <= ub^2)
        end
    end


    vmsr = _PM.var(pm, nw)[:vmsr] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_vmsr",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "vmsr_start")
    )

    vmsi = _PM.var(pm, nw)[:vmsi] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :bus)], base_name="$(nw)_vmsi",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :bus, i), "vmsi_start")
    )

    if bounded
        for (i,b) in _PM.ref(pm, nw, :bus)
            vmin = b["vmin"]
            vmax = b["vmax"]

            JuMP.@constraint(pm.model, vmin^2 <= vmsr[i] + vmsi[i] <= vmax^2)
        end
    end
    cmsdr = _PM.var(pm, nw)[:cmsdr] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_cmsdr",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "cmsdr_start")
    )

    cmsdi = _PM.var(pm, nw)[:cmsdi] = JuMP.@variable(pm.model,
    [i in _PM.ids(pm, nw, :load)], base_name="$(nw)_cmsdi",
    start = _PM.comp_start_value(_PM.ref(pm, nw, :load, i), "cmsdi_start")
    )

    if bounded
        bus = _PM.ref(pm, nw, :bus)
        for (i, l) in _PM.ref(pm, nw, :load)
            vmin = bus[l["load_bus"]]["vmin"]
            @assert vmin > 0
            s = sqrt(abs(l["pd"])^2 + abs(l["qd"])^2)
            ub = s/vmin

            JuMP.@constraint(pm.model, cmsdr[i] + cmsdi[i] <= ub^2)
        end
    end

end

function constraint_current_balance_ac(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
 
     bus = _PM.ref(pm, nw, :bus, i)
     bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
     bus_arcs_dc = _PM.ref(pm, nw, :bus_arcs_dc, i)
     bus_gens = _PM.ref(pm, nw, :bus_gens, i)
     bus_convs_ac = _PM.ref(pm, nw, :bus_convs_ac, i)
     bus_loads = _PM.ref(pm, nw, :bus_loads, i)
     bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)
 
      
 
 
     bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
     bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)
 
     constraint_current_balance_ac(pm, nw, i, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_gs, bus_bs)
 
 end

 function constraint_current_balance_ac(pm::_PM.AbstractIVRModel, n::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_convs_ac, bus_loads, bus_gs, bus_bs)
    vr = _PM.var(pm, n, :vr, i)
    vi = _PM.var(pm, n, :vi, i)

    cr = _PM.var(pm, n, :cr)
    ci = _PM.var(pm, n, :ci)
    cidc = _PM.var(pm, n, :cidc)

    iik_r = _PM.var(pm, n, :iik_r)
    iik_i = _PM.var(pm, n, :iik_i)

    crd = _PM.var(pm, n, :crd)
    cid = _PM.var(pm, n, :cid)
    crg = _PM.var(pm, n, :crg)
    cig = _PM.var(pm, n, :cig)

   


    JuMP.@constraint(pm.model,  sum(cr[a] for a in bus_arcs) + sum(iik_r[c] for c in bus_convs_ac)
                                ==
                                sum(crg[g] for g in bus_gens)
                                - sum(crd[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vr + sum(bs for bs in values(bus_bs))*vi
                                )
    
    JuMP.@constraint(pm.model,  sum(ci[a] for a in bus_arcs) + sum(iik_i[c] for c in bus_convs_ac)
                                + sum(cidc[d] for d in bus_arcs_dc)
                                ==
                                sum(cig[g] for g in bus_gens)
                                - sum(cid[d] for d in bus_loads)
                                - sum(gs for gs in values(bus_gs))*vi - sum(bs for bs in values(bus_bs))*vr
                                )
end

function constraint_load_power(pm::AbstractPowerModel, l::Int; nw::Int=_PM.nw_id_default)
    i   = _PM.ref(pm, nw, :load, l, "load_bus") 

    pd  = _PM.ref(pm, nw, :load, l, "pd")
    qd  = _PM.ref(pm, nw, :load, l, "qd")

    constraint_load_power(pm, nw, i, l, pd, qd)
end

function constraint_load_power(pm::AbstractIVRModel, n::Int, i, l, pd, qd)
    vms  = _PM.var(pm, n, :vms, i) 


    crd = _PM.var(pm, n, :crd, l)
    cid = _PM.var(pm, n, :cid, l) 
    cmsd = _PM.var(pm, n, :cmsd, l) 
    cmsdr = _PM.var(pm, n, :cmsdr, l) 
    cmsdi = _PM.var(pm, n, :cmsdi, l) 


    JuMP.@constraint(pm.model, cmsdr == crd^2)
    JuMP.@constraint(pm.model, cmsdi == cid^2)
    JuMP.@constraint(pm.model, cmsd == cmsdr + cmsdi)
    JuMP.@constraint(pm.model, pd^2 + qd^2 == vms * cmsd)
end
   