# export run_master_scopf_bf

# exact master bf soc model

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

# soft master bf soc model

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

    _PMACDCsc.variable_power_balance_ac_positive_violation(pm)
    _PMACDCsc.variable_power_balance_ac_negative_violation(pm)    


    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in _PM.ids(pm, :bus)
        # _PMACDC.constraint_power_balance_ac(pm, i)
        constraint_power_balance_ac_soft(pm, i)
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
         sum( pg_cost[i] for (i, gen) in _PM.ref(pm, :gen) ) +
         sum( 5E5*_PM.var(pm, :pb_ac_pos_vio, i) for i in _PM.ids(pm, :bus) ) +
         sum( 5E5*_PM.var(pm, :pb_ac_neg_vio, i) for i in _PM.ids(pm, :bus) ) +
         sum( 5E5*_PM.var(pm, :qb_ac_pos_vio, i) for i in _PM.ids(pm, :bus) ) +
         sum( 5E5*_PM.var(pm, :qb_ac_neg_vio, i) for i in _PM.ids(pm, :bus) )
         )
     
     JuMP.@objective(pm.model, Min, obj_expr)

    for (i, cut) in enumerate(_PM.ref(pm, :cuts))
        constraint_benders_cut_soft(pm, i)
    end

end


# exact subproblem bf soc model

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

# soft subproblem bf soc model

function run_sub_scopf_bf_soft(data::Dict{String,Any}, model_type::Type{T}, solver; kwargs...) where T <: _PM.AbstractBFModel
    return _PM.solve_model(data, model_type, solver, build_sub_scopf_bf_soft; ref_extensions = [_PMACDC.add_ref_dcgrid!], kwargs...)
end

function build_sub_scopf_bf_soft(pm::_PM.AbstractPowerModel)

    _PM.variable_bus_voltage(pm, bounded=false, report=true)
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


# soft non-convex subproblem bf model

function run_nc_sub_scopf_soft(data::Dict{String,Any}, model_type::Type, solver; kwargs...)
    return _PM.solve_model(data, model_type, solver, build_nc_sub_scopf_soft; ref_extensions = [_PMACDC.add_ref_dcgrid!], kwargs...)
end

""
function build_nc_sub_scopf_soft(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)

    _PMACDC.variable_active_dcbranch_flow(pm)
    _PMACDC.variable_dcbranch_current(pm)
    _PMACDC.variable_dc_converter(pm)
    _PMACDC.variable_dcgrid_voltage_magnitude(pm)

    _PMACDCsc.variable_branch_thermal_limit_violation(pm)
    _PMACDCsc.variable_power_balance_ac_positive_violation(pm)
    _PMACDCsc.variable_power_balance_ac_negative_violation(pm) 

    # _PM.objective_min_fuel_cost(pm)

    _PM.constraint_model_voltage(pm)
    _PMACDC.constraint_voltage_dc(pm)


    for i in _PM.ids(pm, :gen)
        # constraint_gen_real_setpoint_link(pm,i)
        # constraint_gen_reactive_setpoint_link(pm,i)
    end

    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in _PM.ids(pm, :bus)
        # _PMACDC.constraint_power_balance_ac(pm, i)
        constraint_power_balance_ac_soft(pm, i)
    end

    for i in _PM.ids(pm, :branch)
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)
        _PM.constraint_voltage_angle_difference(pm, i) 
        # _PM.constraint_thermal_limit_from(pm, i)
        # _PM.constraint_thermal_limit_to(pm, i)
        constraint_thermal_limit_from_soft(pm, i)
        constraint_thermal_limit_to_soft(pm, i)
    end
    for i in _PM.ids(pm, :busdc)
        _PMACDC.constraint_power_balance_dc(pm, i)
    end
    for i in _PM.ids(pm, :branchdc)
        _PMACDC.constraint_ohms_dc_branch(pm, i)
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
        # constraint_conv_real_setpoint_link(pm,i)
        constraint_conv_reactive_setpoint_link(pm,i)
    end

    # objective
    _PMSC.objective_c1_variable_pg_cost_basecase(pm)
    pg_cost = _PM.var(pm, :pg_cost)
   
    obj_expr = JuMP.@expression(pm.model,
    sum( pg_cost[i] for (i, gen) in _PM.ref(pm, :gen) ) +
    sum(
        sum( 5E5*_PM.var(pm, :bf_vio_fr, i) for i in _PM.ids(pm, :branch) ) +
        sum( 5E5*_PM.var(pm, :bf_vio_to, i) for i in _PM.ids(pm, :branch) ) + 
        sum( 5E5*_PM.var(pm, :pb_ac_pos_vio, i) for i in _PM.ids(pm, :bus) ) +
        sum( 5E5*_PM.var(pm, :pb_ac_neg_vio, i) for i in _PM.ids(pm, :bus) ) +
        sum( 5E5*_PM.var(pm, :qb_ac_pos_vio, i) for i in _PM.ids(pm, :bus) ) +
        sum( 5E5*_PM.var(pm, :qb_ac_neg_vio, i) for i in _PM.ids(pm, :bus) ) #+
        # sum( 5E5*(_PM.var(pm, :pg, i) - _PM.ref(pm, :gen, i, "pg")) for i in _PM.ids(pm, :gen) ) +
        # sum( 5E5*(_PM.var(pm, :qg, i) - _PM.ref(pm, :gen, i, "qg")) for i in _PM.ids(pm, :gen) )
        )
    )
    

    JuMP.@objective(pm.model, Min, obj_expr)

    # Setting a lower bound   
    JuMP.@constraint(pm.model, obj_expr >= _PM.ref(pm, :soc_master_obj))


end