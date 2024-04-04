# export run_master_scopf_bf


function run_master_scopf_bf(data::Dict{String,Any}, model_type::Type{T}, solver; kwargs...) where T <: _PM.AbstractBFModel
    return _PM.solve_model(data, model_type, solver, build_master_scopf_bf; ref_extensions = [_PMACDC.add_ref_dcgrid!], kwargs...)
end

function build_master_scopf_bf(pm::_PM.AbstractPowerModel)

    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_branch_current(pm)
    _PM.constraint_model_current(pm)

    _PMACDC.variable_active_dcbranch_flow(pm)
    _PMACDC.variable_dcbranch_current(pm)
    _PMACDC.variable_dc_converter(pm)
    _PMACDC.variable_dcgrid_voltage_magnitude(pm)
    _PMACDC.constraint_voltage_dc(pm)


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
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end

    for i in _PM.ids(pm, :busdc)
        _PMACDC.constraint_power_balance_dc(pm, i)
    end

    for i in _PM.ids(pm, :branchdc)
        _PMACDC.constraint_ohms_dc_branch(pm, i)
        _PMACDC.constraint_dc_branch_current(pm, i)
    end

    for i in _PM.ids(pm, :convdc)
        _PMACDC.constraint_converter_losses(pm, i)
        _PMACDC.constraint_converter_current(pm, i)
        _PMACDC.constraint_conv_transformer(pm, i)
        _PMACDC.constraint_conv_reactor(pm, i)
        _PMACDC.constraint_conv_filter(pm, i)
        if pm.ref[:it][:pm][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1
            _PMACDC.constraint_conv_firing_angle(pm, i)
        end
    end

    for (i, cut) in enumerate(_PM.ref(pm, :cuts))
        constraint_benders_cut(pm, i)
    end

    _PM.objective_min_fuel_cost(pm)


end



function run_master_scopf_bf_soft(data::Dict{String,Any}, model_type::Type{T}, solver; kwargs...) where T <: _PM.AbstractBFModel
    return _PM.solve_model(data, model_type, solver, build_master_scopf_bf_soft; ref_extensions = [_PMACDC.add_ref_dcgrid!], kwargs...)
end

function build_master_scopf_bf_soft(pm::_PM.AbstractPowerModel)

    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_branch_current(pm)
    _PM.constraint_model_current(pm)

    _PMACDC.variable_active_dcbranch_flow(pm)
    _PMACDC.variable_dcbranch_current(pm)
    _PMACDC.variable_dc_converter(pm)
    _PMACDC.variable_dcgrid_voltage_magnitude(pm)
    _PMACDC.constraint_voltage_dc(pm)


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
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end

    for i in _PM.ids(pm, :busdc)
        _PMACDC.constraint_power_balance_dc(pm, i)
    end

    for i in _PM.ids(pm, :branchdc)
        _PMACDC.constraint_ohms_dc_branch(pm, i)
        _PMACDC.constraint_dc_branch_current(pm, i)
    end

    for i in _PM.ids(pm, :convdc)
        _PMACDC.constraint_converter_losses(pm, i)
        _PMACDC.constraint_converter_current(pm, i)
        _PMACDC.constraint_conv_transformer(pm, i)
        _PMACDC.constraint_conv_reactor(pm, i)
        _PMACDC.constraint_conv_filter(pm, i)
        if pm.ref[:it][:pm][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1
            _PMACDC.constraint_conv_firing_angle(pm, i)
        end
    end

     # objective
     _PMSC.objective_c1_variable_pg_cost_basecase(pm)
     pg_cost = _PM.var(pm, :pg_cost)
 
     obj_expr = JuMP.@expression(pm.model, 
         sum( pg_cost[i] for (i, gen) in _PM.ref(pm, :gen) ))
     
     JuMP.@objective(pm.model, Min, obj_expr)

    for (i, cut) in enumerate(_PM.ref(pm, :cuts))
        constraint_benders_cut_soft(pm, i)
    end

end




function run_sub_scopf_bf(data::Dict{String,Any}, model_type::Type{T}, solver; kwargs...) where T <: _PM.AbstractBFModel
    return _PM.solve_model(data, model_type, solver, build_sub_scopf_bf; ref_extensions = [_PMACDC.add_ref_dcgrid!], kwargs...)
end

function build_sub_scopf_bf(pm::_PM.AbstractPowerModel)

    _PM.variable_bus_voltage(pm, bounded=false)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_branch_current(pm)
    _PM.constraint_model_current(pm)

    _PMACDC.variable_active_dcbranch_flow(pm)
    _PMACDC.variable_dcbranch_current(pm)
    _PMACDC.variable_dc_converter(pm)
    _PMACDC.variable_dcgrid_voltage_magnitude(pm)
    _PMACDC.constraint_voltage_dc(pm)

    _PMACDCsc.variable_branch_thermal_limit_violation(pm)            

    for i in _PM.ids(pm, :gen)
        constraint_gen_real_setpoint_link(pm,i)
        constraint_gen_reactive_setpoint_link(pm,i)
    end

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
        constraint_thermal_limit_from_soft(pm, i)
        constraint_thermal_limit_to_soft(pm, i)
    end

    for i in _PM.ids(pm, :busdc)
        _PMACDC.constraint_power_balance_dc(pm, i)
    end

    for i in _PM.ids(pm, :branchdc)
        _PMACDC.constraint_ohms_dc_branch(pm, i)
        _PMACDC.constraint_dc_branch_current(pm, i)
    end

    for i in _PM.ids(pm, :convdc)
        _PMACDC.constraint_converter_losses(pm, i)
        _PMACDC.constraint_converter_current(pm, i)
        _PMACDC.constraint_conv_transformer(pm, i)
        _PMACDC.constraint_conv_reactor(pm, i)
        _PMACDC.constraint_conv_filter(pm, i)
        if pm.ref[:it][:pm][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1
            _PMACDC.constraint_conv_firing_angle(pm, i)
        end
        constraint_conv_real_setpoint_link(pm,i)
        constraint_conv_reactive_setpoint_link(pm,i)
    end

    # objective
    _PMSC.objective_c1_variable_pg_cost_basecase(pm)
    pg_cost = _PM.var(pm, :pg_cost)

    obj_expr = JuMP.@expression(pm.model, 
        sum( pg_cost[i] for (i, gen) in _PM.ref(pm, :gen) ) +
        sum(
            sum( 5E5*_PM.var(pm, :bf_vio_fr, i) for i in _PM.ids(pm, :branch) ) +
            sum( 5E5*_PM.var(pm, :bf_vio_to, i) for i in _PM.ids(pm, :branch) )
        ))
    
    JuMP.@objective(pm.model, Min, obj_expr)


end


function run_sub_scopf_bf_soft(data::Dict{String,Any}, model_type::Type{T}, solver; kwargs...) where T <: _PM.AbstractBFModel
    return _PM.solve_model(data, model_type, solver, build_sub_scopf_bf_soft; ref_extensions = [_PMACDC.add_ref_dcgrid!], kwargs...)
end

function build_sub_scopf_bf_soft(pm::_PM.AbstractPowerModel)

    _PM.variable_bus_voltage(pm, bounded=false)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_branch_current(pm)
    _PM.constraint_model_current(pm)

    _PMACDC.variable_active_dcbranch_flow(pm)
    _PMACDC.variable_dcbranch_current(pm)
    _PMACDC.variable_dc_converter(pm)
    _PMACDC.variable_dcgrid_voltage_magnitude(pm)
    _PMACDC.constraint_voltage_dc(pm)

    _PMACDCsc.variable_branch_thermal_limit_violation(pm)
    _PMACDCsc.variable_power_balance_ac_positive_violation(pm)
    _PMACDCsc.variable_power_balance_ac_negative_violation(pm)             

    for i in _PM.ids(pm, :gen)
        constraint_gen_real_setpoint_link(pm,i)
        constraint_gen_reactive_setpoint_link(pm,i)
    end

    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in _PM.ids(pm, :bus)
        constraint_power_balance_ac_soft(pm, i)
    end

    for i in _PM.ids(pm, :branch)
        _PM.constraint_power_losses(pm, i)
        _PM.constraint_voltage_magnitude_difference(pm, i)
        _PM.constraint_voltage_angle_difference(pm, i)
        constraint_thermal_limit_from_soft(pm, i)
        constraint_thermal_limit_to_soft(pm, i)
    end

    for i in _PM.ids(pm, :busdc)
        _PMACDC.constraint_power_balance_dc(pm, i)
    end

    for i in _PM.ids(pm, :branchdc)
        _PMACDC.constraint_ohms_dc_branch(pm, i)
        _PMACDC.constraint_dc_branch_current(pm, i)
    end

    for i in _PM.ids(pm, :convdc)
        _PMACDC.constraint_converter_losses(pm, i)
        _PMACDC.constraint_converter_current(pm, i)
        _PMACDC.constraint_conv_transformer(pm, i)
        _PMACDC.constraint_conv_reactor(pm, i)
        _PMACDC.constraint_conv_filter(pm, i)
        if pm.ref[:it][:pm][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1
            _PMACDC.constraint_conv_firing_angle(pm, i)
        end
        constraint_conv_real_setpoint_link(pm,i)
        constraint_conv_reactive_setpoint_link(pm,i)
    end

    # objective
    _PMSC.objective_c1_variable_pg_cost_basecase(pm)
    pg_cost = _PM.var(pm, :pg_cost)
    JuMP.@objective(pm.model, Min,
    sum( pg_cost[i] for (i, gen) in _PM.ref(pm, :gen) ) +
    sum(
        sum( 5E5*_PM.var(pm, :bf_vio_fr, i) for i in _PM.ids(pm, :branch) ) +
        sum( 5E5*_PM.var(pm, :bf_vio_to, i) for i in _PM.ids(pm, :branch) ) + 
        sum( 5E5*_PM.var(pm, :pb_ac_pos_vio, i) for i in _PM.ids(pm, :bus) ) +
        sum( 5E5*_PM.var(pm, :pb_ac_neg_vio, i) for i in _PM.ids(pm, :bus) ) +
        sum( 5E5*_PM.var(pm, :qb_ac_pos_vio, i) for i in _PM.ids(pm, :bus) ) +
        sum( 5E5*_PM.var(pm, :qb_ac_neg_vio, i) for i in _PM.ids(pm, :bus) )
        )
    )

end




function run_acdcopf_iv(data::Dict{String,Any}, model_type::Type, optimizer; kwargs...)
    return _PM.solve_model(data, model_type, optimizer, build_acdcopf_iv; ref_extensions = [ _PMACDC.add_ref_dcgrid!], kwargs...)
end

function build_acdcopf_iv(pm::_PM.AbstractIVRModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_branch_current(pm)

    variable_gen_current_n(pm)
    variable_load_current(pm, nw = 0)
    variable_load_current_magnitude_squared(pm)

    _PM.variable_dcline_current(pm)

    _PMACDC.variable_active_dcbranch_flow(pm)
    _PMACDC.variable_dcbranch_current(pm)
    _PMACDC.variable_dcgrid_voltage_magnitude(pm)
  
    _PMACDC.variable_dc_converter(pm)

    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in _PM.ids(pm, :bus)
        # _PMACDC.constraint_current_balance_ac(pm, i)
        constraint_current_balance_ac(pm,i)
    end

    for l in _PM.ids(pm, :load)
        constraint_load_power(pm,l)
    end

    for i in _PM.ids(pm, :branch)
        _PM.constraint_current_from(pm, i)
        _PM.constraint_current_to(pm, i)

        _PM.constraint_voltage_drop(pm, i)
        _PM.constraint_voltage_angle_difference(pm, i)

        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end

    for i in _PM.ids(pm, :busdc)
        _PMACDC.constraint_current_balance_dc(pm, i)
    end
    for i in _PM.ids(pm, :branchdc)
        _PMACDC.constraint_ohms_dc_branch(pm, i)
    end
    for i in _PM.ids(pm, :convdc)
        _PMACDC.constraint_converter_limits(pm, i)
        _PMACDC.constraint_converter_losses(pm, i)
        _PMACDC.constraint_converter_current(pm, i)
        _PMACDC.constraint_conv_transformer(pm, i)
        _PMACDC.constraint_conv_reactor(pm, i)
        _PMACDC.constraint_conv_filter(pm, i)
    end

    _PM.objective_min_fuel_and_flow_cost(pm)
end