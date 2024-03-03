function variable_branch_contigency_power_violation(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    branch_cont_flow_vio = _PM.var(pm, nw)[:branch_cont_flow_vio] = JuMP.@variable(pm.model,
        [i in 1:length(_PM.ref(pm, nw, :branch_flow_cuts))], base_name="$(nw)_branch_cont_flow_vio",
        #start = _PM.comp_start_value(ref(pm, nw, :bus, i), "cont_branch_vio_start")
    )

    if bounded
        for i in 1:length(_PM.ref(pm, nw, :branch_flow_cuts))
            JuMP.set_lower_bound(branch_cont_flow_vio[i], 0.0)
        end
    end

    #report && _PM.sol_component_value(pm, nw, :gen, :pg_delta, ids(pm, nw, :gen), pg_delta)
end


function variable_gen_contigency_power_violation(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    gen_cont_flow_vio = _PM.var(pm, nw)[:gen_cont_flow_vio] = JuMP.@variable(pm.model,
        [i in 1:length(_PM.ref(pm, nw, :gen_flow_cuts))], base_name="$(nw)_gen_cont_flow_vio",
        #start = _PM.comp_start_value(ref(pm, nw, :bus, i), "cont_branch_vio_start")
    )

    if bounded
        for i in 1:length(_PM.ref(pm, nw, :gen_flow_cuts))
            JuMP.set_lower_bound(gen_cont_flow_vio[i], 0.0)
        end
    end

    #report && _PM.sol_component_value(pm, nw, :gen, :pg_delta, ids(pm, nw, :gen), pg_delta)
end


function variable_gen_contigency_capacity_violation(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default, bounded::Bool=true, report::Bool=true)
    gen_cont_cap_vio = _PM.var(pm, nw)[:gen_cont_cap_vio] = JuMP.@variable(pm.model,
        [i in 1:length(_PM.ref(pm, nw, :gen_contingencies))], base_name="$(nw)_gen_cont_cap_vio",
        #start = _PM.comp_start_value(ref(pm, nw, :bus, i), "cont_branch_vio_start")
    )

    if bounded
        for i in 1:length(_PM.ref(pm, nw, :gen_contingencies))
            JuMP.set_lower_bound(gen_cont_cap_vio[i], 0.0)
        end
    end

    #report && _PM.sol_component_value(pm, nw, :gen, :pg_delta, ids(pm, nw, :gen), pg_delta)
end

# tnep
function variable_branch_ne_current(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_branch_ne_series_current_magnitude_squared(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_branch_ne_series_current_magnitude_squared(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    cmss_ne = _PM.var(pm, nw)[:cmss_ne] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :ne_branch)], base_name="$(nw)_cmss_ne",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :ne_branch, l), "cmss_ne_start", 0.0)
    )

    if bounded
        bus = _PM.ref(pm, nw, :bus)
        branch = _PM.ref(pm, nw, :ne_branch)

        for (l,i,j) in _PM.ref(pm, nw, :ne_arcs_from)
            b = branch[l]
            ub = Inf
            if haskey(b, "rate_a")
                rate = b["rate_a"] * b["tap"]
                y_fr = abs(b["g_fr"] + im * b["b_fr"])
                y_to = abs(b["g_to"] + im * b["b_to"])
                shunt_current = max(y_fr * bus[i]["vmax"]^2, y_to * bus[j]["vmax"]^2)
                series_current = max(rate / bus[i]["vmin"], rate / bus[j]["vmin"])
                ub = series_current + shunt_current
            end
            if haskey(b, "c_rating_a")
                total_current = b["c_rating_a"]
                y_fr = abs(b["g_fr"] + im * b["b_fr"])
                y_to = abs(b["g_to"] + im * b["b_to"])
                shunt_current = max(y_fr * bus[i]["vmax"]^2, y_to * bus[j]["vmax"]^2)
                ub = total_current + shunt_current
            end

            if !isinf(ub)
                JuMP.set_lower_bound(cmss_ne[l], -2.0 * ub^2)
                JuMP.set_upper_bound(cmss_ne[l],  2.0 * ub^2)
            end
        end
    end

    report && _PM.sol_component_value(pm, nw, :ne_branch, :cmss_ne, _PM.ids(pm, nw, :ne_branch), cmss_ne)
end

function variable_branch_ne_voltage_drop(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    variable_branch_ne_voltage_drop_real(pm, nw=nw, bounded=bounded, report=report; kwargs...)
    variable_branch_ne_voltage_drop_img(pm, nw=nw, bounded=bounded, report=report; kwargs...)
end

function variable_branch_ne_voltage_drop_real(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    vbdr_ne = _PM.var(pm, nw)[:vbdr_ne] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :ne_branch)], base_name="$(nw)_vbdr_ne",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :ne_branch, l), "vbdr_ne_start", 0.01)
    )
    
    report && _PM.sol_component_value(pm, nw, :ne_branch, :vbdr_ne, _PM.ids(pm, nw, :ne_branch), vbdr_ne)
end

function variable_branch_ne_voltage_drop_img(pm::AbstractACRModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)

    vbdi_ne = _PM.var(pm, nw)[:vbdi_ne] = JuMP.@variable(pm.model,
        [l in _PM.ids(pm, nw, :ne_branch)], base_name="$(nw)_vbdi_ne",
        start = _PM.comp_start_value(_PM.ref(pm, nw, :ne_branch, l), "vbdi_ne_start", 0.01)
    )

    report && _PM.sol_component_value(pm, nw, :ne_branch, :vbdi_ne, _PM.ids(pm, nw, :ne_branch), vbdi_ne)
end

