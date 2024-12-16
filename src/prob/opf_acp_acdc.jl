
function run_sopf_acp_soc_bf(data::Dict{String,Any}, model_type::Type{T}, solver; kwargs...) where T <: _PM.AbstractBFModel
    return _PM.solve_model(data, model_type, solver, build_sopf_acp_soc_bf; ref_extensions = [_PMACDC.add_ref_dcgrid!], kwargs...)
end

function build_sopf_acp_soc_bf(pm::_PM.AbstractPowerModel)

    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_branch_current(pm)

    _PM.constraint_model_current(pm)    # QC

    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in _PM.ids(pm, :bus)
        _PMACDC.constraint_power_balance_ac(pm, i)
    end

    for i in _PM.ids(pm, :branch)
        _PM.constraint_power_losses(pm, i)
        _PM.constraint_voltage_magnitude_difference(pm, i)
        _PM.constraint_voltage_angle_difference(pm, i)
        _PM.constraint_thermal_limit_from(pm, i)    # QC
        _PM.constraint_thermal_limit_to(pm, i)      # QC
    end

    _PM.objective_min_fuel_cost(pm)

end




# for i in _PM.ids(pm, :busdc)
#     _PMACDC.constraint_power_balance_dc(pm, i)
# end

# for i in _PM.ids(pm, :branchdc)
#     _PMACDC.constraint_ohms_dc_branch(pm, i)
#     _PMACDC.constraint_dc_branch_current(pm, i)
# end

# for i in _PM.ids(pm, :convdc)
#     _PMACDC.constraint_converter_losses(pm, i)
#     _PMACDC.constraint_converter_current(pm, i)
#     _PMACDC.constraint_conv_transformer(pm, i)
#     _PMACDC.constraint_conv_reactor(pm, i)
#     _PMACDC.constraint_conv_filter(pm, i)
#     if pm.ref[:it][:pm][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1
#         _PMACDC.constraint_conv_firing_angle(pm, i)
#     end
# end

# _PMACDC.variable_active_dcbranch_flow(pm)
# _PMACDC.variable_dcbranch_current(pm)
# _PMACDC.variable_dc_converter(pm)
# _PMACDC.variable_dcgrid_voltage_magnitude(pm)
# _PMACDC.constraint_voltage_dc(pm)