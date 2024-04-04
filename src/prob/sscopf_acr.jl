

""
function solve_sscopf_acr(data::Dict, model_constructor::Type, optimizer; deg::Int=1, solution_processors=[sol_data_model!], kwargs...)
    @assert _IM.ismultinetwork(data) == false "The data supplied is multinetwork, it should be single-network"
    @assert model_constructor <: _PM.AbstractACRModel "This problem type only supports the ACRModel"

    sdata = build_stochastic_data(data, deg)
    result = _PM.solve_model(sdata, model_constructor, optimizer, build_sscopf_acr; multinetwork=true, solution_processors=solution_processors, ref_extensions = [_PMSC.ref_c1!], kwargs...)
    result["mop"] = sdata["mop"]

    return result
end

""
function solve_sscopf_acr(file::String, model_constructor, optimizer; deg::Int=1, solution_processors=[sol_data_model!], kwargs...)
    data = _PM.parse_file(file)

    return solve_sscopf_acr(data, model_constructor, optimizer; deg=deg, solution_processors=solution_processors, kwargs...)
end

""
function build_sscopf_acr(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 
        variable_bus_voltage(pm, nw=n)
        variable_gen_power(pm, nw=n, bounded=false)
        variable_branch_power(pm, nw=n, bounded=false)
        variable_branch_current(pm, nw=n, bounded=false)
        variable_branch_voltage_drop(pm, nw=n) 

        variable_branch_contigency_power_violation(pm, nw=n)
        variable_gen_contigency_power_violation(pm, nw=n)
        variable_gen_contigency_capacity_violation(pm, nw=n)
    end



    for (n, network) in _PM.nws(pm)
        for i in _PM.ids(pm, :bus, nw=n)
            _PMSC.expression_c1_bus_generation(pm, i, nw=n)
            _PMSC.expression_c1_bus_withdrawal(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_power_shunt_dispatch(pm, i, nw=n)
            constraint_gp_bus_voltage_magnitude_squared(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)   
            constraint_branch_voltage(pm, b, nw=n)
            constraint_gp_power_branch_to(pm, b, nw=n)
            constraint_gp_goc_power_branch_from(pm, b, nw=n)
            constraint_gp_branch_series_current_magnitude_squared(pm, b, nw=n)
        end
        
        for (i,cut) in enumerate(_PM.ref(pm, nw=n, :branch_flow_cuts))
            constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm, i, nw=n)
            constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm, i, nw=n)
        end
        
        bus_withdrawal = _PM.var(pm, nw=n, :bus_wdp)
        
        for (i,cut) in enumerate(_PM.ref(pm, nw=n, :gen_flow_cuts))
            branch = _PM.ref(pm, nw=n, :branch, cut.branch_id)
            gen = _PM.ref(pm, nw=n, :gen, cut.gen_id)
            gen_bus = _PM.ref(pm, nw=n, :bus, gen["gen_bus"])
            gen_set = _PM.ref(pm, nw=n, :area_gens)[gen_bus["area"]]
            alpha_total = sum(gen["alpha"] for (i,gen) in _PM.ref(pm, nw=n, :gen) if gen["index"] != cut.gen_id && i in gen_set)
        
            cont_bus_injection = Dict{Int,Any}()
            for (i, bus) in _PM.ref(pm, nw=n, :bus)
                inj = 0.0
                for g in _PM.ref(pm, nw=n, :bus_gens, i)
                    if g != cut.gen_id
                        if g in gen_set
                            inj += _PM.var(pm, nw=n, :pg, g) + gen["alpha"]*_PM.var(pm, nw=n, :pg, cut.gen_id)/alpha_total
                        else
                            inj += _PM.var(pm, nw=n, :pg, g)
                        end
                    end
                end
                cont_bus_injection[i] = inj
            end
        
            # rate = branch["rate_a"]
            rate = branch["rate_c"]
            JuMP.@constraint(pm.model,  sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + _PM.var(pm, nw=n, :gen_cont_flow_vio, i))
            JuMP.@constraint(pm.model, -sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + _PM.var(pm, nw=n, :gen_cont_flow_vio, i))
        end
        
        for (i,gen_cont) in enumerate(_PM.ref(pm, nw=n, :gen_contingencies))
            #println(gen_cont)
            gen = _PM.ref(pm, nw=n, :gen, gen_cont.idx)
            gen_bus = _PM.ref(pm, nw=n, :bus, gen["gen_bus"])
            gen_set = _PM.ref(pm, nw=n, :area_gens)[gen_bus["area"]]
            response_gens = Dict(g => _PM.ref(pm, nw=n, :gen, g) for g in gen_set if g != gen_cont.idx)
        
            # factor of 1.2 accounts for losses in a DC model
            # JuMP.@constraint(pm.model, sum(gen["pmax"] - var(pm, :pg, g) for (g,gen) in response_gens) >= 1.2*var(pm, :pg, gen_cont.idx))
            JuMP.@constraint(pm.model, _PM.var(pm, nw=n, :gen_cont_cap_vio, i) + sum(gen["pmax"] - _PM.var(pm, nw=n, :pg, g) for (g,gen) in response_gens) >= _PM.var(pm, nw=n, :pg, gen_cont.idx))
            # JuMP.@constraint(pm.model, sum(gen["pmin"] - var(pm, :pg, g) for (g,gen) in response_gens) <= var(pm, :pg, gen_cont.idx))
        end
    end

    for i in _PM.ids(pm, :bus,nw=1)
        constraint_cc_bus_voltage_magnitude_squared(pm, i, nw=1)
    end

    for b in _PM.ids(pm, :branch, nw=1)
        constraint_cc_branch_series_current_magnitude_squared(pm, b, nw=1)
    end

    for g in _PM.ids(pm, :gen, nw=1)
        constraint_cc_gen_power(pm, g, nw=1)
    end

    # for (i,cut) in enumerate(_PM.ref(pm, nw=1, :branch_flow_cuts))
    #     constraint_cc_branch_contingency_ptdf_thermal_limit_from_soft(pm, i, nw=1)
    #     # constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm, i, nw=n)
    # end

    objective_min_expected_generation_cost_soft(pm)
end









""
function solve_stnep_acr(data::Dict, model_constructor::Type, optimizer; deg::Int=1, solution_processors=[sol_data_model!], kwargs...)
    @assert _IM.ismultinetwork(data) == false "The data supplied is multinetwork, it should be single-network"
    @assert model_constructor <: _PM.AbstractACRModel "This problem type only supports the ACRModel"

    sdata = build_stochastic_data(data, deg)
    result = _PM.solve_model(sdata, model_constructor, optimizer, build_stnep_acr; multinetwork=true, solution_processors=solution_processors, ref_extensions = [_PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!, _PMSC.ref_c1!], kwargs...)
    result["mop"] = sdata["mop"]

    return result
end

""
function solve_stnep_acr(file::String, model_constructor, optimizer; deg::Int=1, solution_processors=[sol_data_model!], kwargs...)
    data = _PM.parse_file(file)

    return solve_stnep_acr(data, model_constructor, optimizer; deg=deg, solution_processors=solution_processors, kwargs...)
end

""
function build_stnep_acr(pm::AbstractPowerModel)
    for (n, network) in _PM.nws(pm) 
        variable_bus_voltage(pm, nw=n)
        variable_gen_power(pm, nw=n, bounded=false)
        variable_branch_power(pm, nw=n, bounded=false)
        variable_branch_current(pm, nw=n, bounded=false)
        variable_branch_voltage_drop(pm, nw=n) 

        variable_branch_contigency_power_violation(pm, nw=n)
        variable_gen_contigency_power_violation(pm, nw=n)
        variable_gen_contigency_capacity_violation(pm, nw=n)

        _PM.variable_ne_branch_indicator(pm, nw=n)
        _PM.variable_ne_branch_power(pm, nw=n, bounded=false)
        variable_branch_ne_current(pm, nw=n, bounded=false)
        variable_branch_ne_voltage_drop(pm, nw=n) 

        _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=n)
    end

    for (n, network) in _PM.nws(pm)
        for i in _PM.ids(pm, :bus, nw=n)
            _PMSC.expression_c1_bus_generation(pm, i, nw=n)
            _PMSC.expression_c1_bus_withdrawal(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :ref_buses, nw=n)
            constraint_bus_voltage_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_ne_power_shunt_dispatch(pm, i, nw=n) 
            constraint_gp_bus_voltage_magnitude_squared(pm, i, nw=n)
        end

        for b in _PM.ids(pm, :branch, nw=n)   
            constraint_branch_voltage(pm, b, nw=n)
            constraint_gp_power_branch_to(pm, b, nw=n)
            constraint_gp_goc_power_branch_from(pm, b, nw=n)
            constraint_gp_branch_series_current_magnitude_squared(pm, b, nw=n)
        end

        for b_ne in _PM.ids(pm, :ne_branch, nw=n)
            constraint_branch_ne_voltage(pm, b_ne, nw=n)
            constraint_gp_power_branch_ne_to(pm, b_ne, nw=n)
            constraint_gp_power_branch_ne_from(pm, b_ne, nw=n)
            constraint_gp_branch_ne_series_current_magnitude_squared(pm, b_ne, nw=n)
        end
        
        for (i,cut) in enumerate(_PM.ref(pm, nw=n, :branch_flow_cuts))
            constraint_branch_contingency_ptdf_thermal_limit_from_soft(pm, i, nw=n)
            constraint_branch_contingency_ptdf_thermal_limit_to_soft(pm, i, nw=n)
        end
        
        bus_withdrawal = _PM.var(pm, nw=n, :bus_wdp)
        
        for (i,cut) in enumerate(_PM.ref(pm, nw=n, :gen_flow_cuts))
            branch = _PM.ref(pm, nw=n, :branch, cut.branch_id)
            gen = _PM.ref(pm, nw=n, :gen, cut.gen_id)
            gen_bus = _PM.ref(pm, nw=n, :bus, gen["gen_bus"])
            gen_set = _PM.ref(pm, nw=n, :area_gens)[gen_bus["area"]]
            alpha_total = sum(gen["alpha"] for (i,gen) in _PM.ref(pm, nw=n, :gen) if gen["index"] != cut.gen_id && i in gen_set)
        
            cont_bus_injection = Dict{Int,Any}()
            for (i, bus) in _PM.ref(pm, nw=n, :bus)
                inj = 0.0
                for g in _PM.ref(pm, nw=n, :bus_gens, i)
                    if g != cut.gen_id
                        if g in gen_set
                            inj += _PM.var(pm, nw=n, :pg, g) + gen["alpha"]*_PM.var(pm, nw=n, :pg, cut.gen_id)/alpha_total
                        else
                            inj += _PM.var(pm, nw=n, :pg, g)
                        end
                    end
                end
                cont_bus_injection[i] = inj
            end
        
            # rate = branch["rate_a"]
            rate = branch["rate_c"]
            JuMP.@constraint(pm.model,  sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + _PM.var(pm, nw=n, :gen_cont_flow_vio, i))
            JuMP.@constraint(pm.model, -sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + _PM.var(pm, nw=n, :gen_cont_flow_vio, i))
        end
        
        for (i,gen_cont) in enumerate(_PM.ref(pm, nw=n, :gen_contingencies))
            #println(gen_cont)
            gen = _PM.ref(pm, nw=n, :gen, gen_cont.idx)
            gen_bus = _PM.ref(pm, nw=n, :bus, gen["gen_bus"])
            gen_set = _PM.ref(pm, nw=n, :area_gens)[gen_bus["area"]]
            response_gens = Dict(g => _PM.ref(pm, nw=n, :gen, g) for g in gen_set if g != gen_cont.idx)
        
            # factor of 1.2 accounts for losses in a DC model
            # JuMP.@constraint(pm.model, sum(gen["pmax"] - var(pm, :pg, g) for (g,gen) in response_gens) >= 1.2*var(pm, :pg, gen_cont.idx))
            JuMP.@constraint(pm.model, _PM.var(pm, nw=n, :gen_cont_cap_vio, i) + sum(gen["pmax"] - _PM.var(pm, nw=n, :pg, g) for (g,gen) in response_gens) >= _PM.var(pm, nw=n, :pg, gen_cont.idx))
            # JuMP.@constraint(pm.model, sum(gen["pmin"] - var(pm, :pg, g) for (g,gen) in response_gens) <= var(pm, :pg, gen_cont.idx))
        end
    end

    for i in _PM.ids(pm, :bus,nw=1)
        constraint_cc_bus_voltage_magnitude_squared(pm, i, nw=1)
    end

    for b in _PM.ids(pm, :branch, nw=1)
        constraint_cc_branch_series_current_magnitude_squared(pm, b, nw=1)
    end

    for b_ne in _PM.ids(pm, :ne_branch, nw=1)
        constraint_cc_branch_ne_series_current_magnitude_squared(pm, b_ne, nw=1)
    end

    for g in _PM.ids(pm, :gen, nw=1)
        constraint_cc_gen_power(pm, g, nw=1)
    end

    objective_min_expected_tnep_generation_cost_soft(pm)
end

