

# function run_scopf_cuts(ini_file::String, model_type::Type, optimizer; scenario_id::String="", kwargs...)
#     goc_data = parse_c1_case(ini_file, scenario_id=scenario_id)
#     network = build_c1_pm_model(goc_data)
#     return run_scopf_cuts!(network, model_type, optimizer; kwargs...)
# end

"""
Solves a SCOPF problem by iteratively checking for violated branch flow
constraints in contingencies and resolving until a fixed-point is reached.

The base-case model is formulation agnostic.  The flow cuts are based on PTDF
and utilize the DC Power Flow assumption.
"""
function run_sscopf_cuts(network::Dict{String,<:Any}, model_type::Type, optimizer; deg::Int=1, max_iter::Int=100, time_limit::Float64=Inf)
    if _IM.ismultinetwork(network)
        Memento.error(_LOGGER, "run_c1_scopf_ptdf_cuts can only be used on single networks")
    end

    time_start = time()

    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = _PM.solve_opf(network, model_type, optimizer)
    if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
        Memento.error(_LOGGER, "base-case OPF solve failed in run_c1_scopf_ptdf_cuts, status $(result["termination_status"])")
    end
    Memento.info(_LOGGER, "objective: $(result["objective"])")
    _PM.update_data!(network, result["solution"])

    result["iterations"] = 0

    iteration = 1
    cuts_found = 1
    while cuts_found > 0
        time_start_iteration = time()

        cuts = _PMSC.check_c1_contingencies_branch_power(network, total_cut_limit=iteration, gen_flow_cuts=[], branch_flow_cuts=[])

        cuts_found = 0

        for cut in cuts.gen_cuts
            if !(cut in network["gen_flow_cuts"])
                append!(network["gen_flow_cuts"], cuts.gen_cuts)
                cuts_found +=1
            else
                Memento.warn(_LOGGER, "flow cut on cont $(cut.cont_label) branch $(cut.branch_id) is active but not secure")
            end
        end
        for cut in cuts.branch_cuts
            if !(cut in network["branch_flow_cuts"])
                append!(network["branch_flow_cuts"], cuts.branch_cuts)
                cuts_found +=1
            else
                Memento.warn(_LOGGER, "flow cut on cont $(cut.cont_label) branch $(cut.branch_id) is active but not secure")
            end
        end

        # cuts_found = length(cuts.gen_cuts)+length(cuts.branch_cuts)
        if cuts_found <= 0
            Memento.info(_LOGGER, "no violated cuts found scopf fixed-point reached")
            break
        else
            Memento.info(_LOGGER, "found $(cuts_found) branch flow violations")
        end

        # append!(network["gen_flow_cuts"], cuts.gen_cuts)
        # append!(network["branch_flow_cuts"], cuts.branch_cuts)

        Memento.info(_LOGGER, "active cuts: gen $(length(network["gen_flow_cuts"])), branch $(length(network["branch_flow_cuts"]))")

        time_solve_start = time()
        result = solve_sscopf_acr(network, model_type, optimizer, deg=deg)
        if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
            Memento.warn(_LOGGER, "scopf solve failed with status $(result["termination_status"]), terminating fixed-point early")
            break
        end
        Memento.info(_LOGGER, "objective: $(result["objective"])")
        solution = result["solution"]["nw"]["1"]
        solution["per_unit"] = result["solution"]["per_unit"]
        _PM.update_data!(network, solution)

        time_iteration = time() - time_start_iteration
        time_remaining = time_limit - (time() - time_start)
        if time_remaining < time_iteration
            Memento.warn(_LOGGER, "insufficent time for next iteration, time remaining $(time_remaining), estimated iteration time $(time_iteration)")
            break
        end
        iteration += 1
    end

    result["iterations"] = iteration
    return result
end


function run_ssctnep_cuts(network::Dict{String,<:Any}, model_type::Type, optimizer; deg::Int=1, max_iter::Int=100, time_limit::Float64=Inf)
    if _IM.ismultinetwork(network)
        Memento.error(_LOGGER, "run_c1_scopf_ptdf_cuts can only be used on single networks")
    end

    time_start = time()

    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = _PM.solve_opf(network, model_type, optimizer)
    if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
        Memento.error(_LOGGER, "base-case OPF solve failed in run_c1_scopf_ptdf_cuts, status $(result["termination_status"])")
    end
    Memento.info(_LOGGER, "objective: $(result["objective"])")
    _PM.update_data!(network, result["solution"])

    result["iterations"] = 0

    iteration = 1
    cuts_found = 1
    while cuts_found > 0
        time_start_iteration = time()

        cuts = _PMSC.check_c1_contingencies_branch_power(network, total_cut_limit=iteration, gen_flow_cuts=[], branch_flow_cuts=[])

        cuts_found = 0

        for cut in cuts.gen_cuts
            if !(cut in network["gen_flow_cuts"])
                append!(network["gen_flow_cuts"], cuts.gen_cuts)
                cuts_found +=1
            else
                Memento.warn(_LOGGER, "flow cut on cont $(cut.cont_label) branch $(cut.branch_id) is active but not secure")
            end
        end
        for cut in cuts.branch_cuts
            if !(cut in network["branch_flow_cuts"])
                append!(network["branch_flow_cuts"], cuts.branch_cuts)
                cuts_found +=1
            else
                Memento.warn(_LOGGER, "flow cut on cont $(cut.cont_label) branch $(cut.branch_id) is active but not secure")
            end
        end

        # cuts_found = length(cuts.gen_cuts)+length(cuts.branch_cuts)
        if cuts_found <= 0
            Memento.info(_LOGGER, "no violated cuts found scopf fixed-point reached")
            break
        else
            Memento.info(_LOGGER, "found $(cuts_found) branch flow violations")
        end

        # append!(network["gen_flow_cuts"], cuts.gen_cuts)
        # append!(network["branch_flow_cuts"], cuts.branch_cuts)

        Memento.info(_LOGGER, "active cuts: gen $(length(network["gen_flow_cuts"])), branch $(length(network["branch_flow_cuts"]))")

        for (i,branch) in network["branch"]
            if haskey(branch, "merge")
                delete!(network["branch"], i)
                Memento.info(_LOGGER, "branch $i unmerged")
            end
        end

        time_solve_start = time()
        result = solve_stnep_acr(network, model_type, optimizer, deg=deg)
        if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
            Memento.warn(_LOGGER, "scopf solve failed with status $(result["termination_status"]), terminating fixed-point early")
            break
        end
        Memento.info(_LOGGER, "objective: $(result["objective"])")
        solution = result["solution"]["nw"]["1"]
        solution["per_unit"] = result["solution"]["per_unit"]
        _PM.update_data!(network, solution)

        for (i, branch) in result["solution"]["nw"]["1"]["ne_branch"]
            b = length(network["branch"])
            if isapprox(branch["built"], 1.0; atol = 0.1)
                network["branch"]["$(b+1)"] = network["ne_branch"][i]
                network["branch"]["$(b+1)"]["pf"]  = branch["pf"]
                network["branch"]["$(b+1)"]["pt"]  = branch["pt"]
                network["branch"]["$(b+1)"]["qf"]  = branch["qf"]
                network["branch"]["$(b+1)"]["qt"]  = branch["qt"]
                network["branch"]["$(b+1)"]["cmss"]  = branch["cmss_ne"]
                network["branch"]["$(b+1)"]["vbdi"]  = branch["vbdi_ne"]
                network["branch"]["$(b+1)"]["vbdr"]  = branch["vbdr_ne"]
                network["branch"]["$(b+1)"]["merge"]  = 1
                Memento.info(_LOGGER, "ne_branch $i is merged and updated as branch $(b+1) for contingency checks")
            end
        end

        time_iteration = time() - time_start_iteration
        time_remaining = time_limit - (time() - time_start)
        if time_remaining < time_iteration
            Memento.warn(_LOGGER, "insufficent time for next iteration, time remaining $(time_remaining), estimated iteration time $(time_iteration)")
            break
        end
        iteration += 1
    end

    result["iterations"] = iteration
    return result
end


function run_scopf_cuts(network::Dict{String,<:Any}, model_type::Type, optimizer; max_iter::Int=100, time_limit::Float64=Inf)
    if _IM.ismultinetwork(network)
        Memento.error(_LOGGER, "run_c1_scopf_ptdf_cuts can only be used on single networks")
    end

    time_start = time()

    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = _PM.solve_opf(network, model_type, optimizer)
    if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
        Memento.error(_LOGGER, "base-case OPF solve failed in run_c1_scopf_ptdf_cuts, status $(result["termination_status"])")
    end
    Memento.info(_LOGGER, "objective: $(result["objective"])")
    _PM.update_data!(network, result["solution"])

    result["iterations"] = 0

    iteration = 1
    cuts_found = 1
    while cuts_found > 0
        time_start_iteration = time()

        cuts = _PMSC.check_c1_contingencies_branch_power(network, total_cut_limit=iteration, gen_flow_cuts=[], branch_flow_cuts=[])

        cuts_found = 0

        for cut in cuts.gen_cuts
            if !(cut in network["gen_flow_cuts"])
                append!(network["gen_flow_cuts"], cuts.gen_cuts)
                cuts_found +=1
            else
                Memento.warn(_LOGGER, "flow cut on cont $(cut.cont_label) branch $(cut.branch_id) is active but not secure")
            end
        end
        for cut in cuts.branch_cuts
            if !(cut in network["branch_flow_cuts"])
                append!(network["branch_flow_cuts"], cuts.branch_cuts)
                cuts_found +=1
            else
                Memento.warn(_LOGGER, "flow cut on cont $(cut.cont_label) branch $(cut.branch_id) is active but not secure")
            end
        end

        # cuts_found = length(cuts.gen_cuts)+length(cuts.branch_cuts)
        if cuts_found <= 0
            Memento.info(_LOGGER, "no violated cuts found scopf fixed-point reached")
            break
        else
            Memento.info(_LOGGER, "found $(cuts_found) branch flow violations")
        end

        Memento.info(_LOGGER, "active cuts: gen $(length(network["gen_flow_cuts"])), branch $(length(network["branch_flow_cuts"]))")

        time_solve_start = time()
        result = run_scopf_cuts_soft(network, model_type, optimizer)
        if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
            Memento.warn(_LOGGER, "scopf solve failed with status $(result["termination_status"]), terminating fixed-point early")
            break
        end
        Memento.info(_LOGGER, "objective: $(result["objective"])")
        _PM.update_data!(network, result["solution"])

        time_iteration = time() - time_start_iteration
        time_remaining = time_limit - (time() - time_start)
        if time_remaining < time_iteration
            Memento.warn(_LOGGER, "insufficent time for next iteration, time remaining $(time_remaining), estimated iteration time $(time_iteration)")
            break
        end
        iteration += 1
    end

    result["iterations"] = iteration
    return result
end

function run_sctnep_cuts(network::Dict{String,<:Any}, model_type::Type, optimizer; max_iter::Int=100, time_limit::Float64=Inf)
    if _IM.ismultinetwork(network)
        Memento.error(_LOGGER, "run_c1_scopf_ptdf_cuts can only be used on single networks")
    end

    time_start = time()

    network["gen_flow_cuts"] = []
    network["branch_flow_cuts"] = []

    result = _PM.solve_opf(network, model_type, optimizer)
    if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
        Memento.error(_LOGGER, "base-case OPF solve failed in run_c1_scopf_ptdf_cuts, status $(result["termination_status"])")
    end
    Memento.info(_LOGGER, "objective: $(result["objective"])")
    _PM.update_data!(network, result["solution"])

    result["iterations"] = 0

    iteration = 1
    cuts_found = 1
    while cuts_found > 0
        time_start_iteration = time()

        cuts = _PMSC.check_c1_contingencies_branch_power(network, total_cut_limit=iteration, gen_flow_cuts=[], branch_flow_cuts=[])

        cuts_found = 0

        for cut in cuts.gen_cuts
            if !(cut in network["gen_flow_cuts"])
                append!(network["gen_flow_cuts"], cuts.gen_cuts)
                cuts_found +=1
            else
                Memento.warn(_LOGGER, "flow cut on cont $(cut.cont_label) branch $(cut.branch_id) is active but not secure")
            end
        end
        for cut in cuts.branch_cuts
            if !(cut in network["branch_flow_cuts"])
                append!(network["branch_flow_cuts"], cuts.branch_cuts)
                cuts_found +=1
            else
                Memento.warn(_LOGGER, "flow cut on cont $(cut.cont_label) branch $(cut.branch_id) is active but not secure")
            end
        end

        # cuts_found = length(cuts.gen_cuts)+length(cuts.branch_cuts)
        if cuts_found <= 0
            Memento.info(_LOGGER, "no violated cuts found scopf fixed-point reached")
            break
        else
            Memento.info(_LOGGER, "found $(cuts_found) branch flow violations")
        end

        Memento.info(_LOGGER, "active cuts: gen $(length(network["gen_flow_cuts"])), branch $(length(network["branch_flow_cuts"]))")

        time_solve_start = time()
        result = run_tnep_scopf_cuts_soft(network, model_type, optimizer)
        if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
            Memento.warn(_LOGGER, "scopf solve failed with status $(result["termination_status"]), terminating fixed-point early")
            break
        end
        Memento.info(_LOGGER, "objective: $(result["objective"])")
        _PM.update_data!(network, result["solution"])

        time_iteration = time() - time_start_iteration
        time_remaining = time_limit - (time() - time_start)
        if time_remaining < time_iteration
            Memento.warn(_LOGGER, "insufficent time for next iteration, time remaining $(time_remaining), estimated iteration time $(time_iteration)")
            break
        end
        iteration += 1
    end

    result["iterations"] = iteration
    return result
end


function run_scopf_cuts_soft(file, model_constructor, solver; kwargs...)
    return _PM.solve_model(file, model_constructor, solver, build_scopf_cuts_soft; ref_extensions=[_PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!, _PMSC.ref_c1!], kwargs...)
end

""
function build_scopf_cuts_soft(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)

    _PMSC.variable_c1_shunt_admittance_imaginary(pm)

    _PMSC.variable_c1_branch_contigency_power_violation(pm)
    _PMSC.variable_c1_gen_contigency_power_violation(pm)
    _PMSC.variable_c1_gen_contigency_capacity_violation(pm)

    for i in _PM.ids(pm, :bus)
        _PMSC.expression_c1_bus_generation(pm, i)
        _PMSC.expression_c1_bus_withdrawal(pm, i)
    end


    _PM.constraint_model_voltage(pm)

    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in _PM.ids(pm, :bus)
        _PMSC.constraint_c1_power_balance_shunt_dispatch(pm, i)
    end

    for i in _PM.ids(pm, :branch)
        _PMSC.constraint_goc_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)

        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end

    for (i,cut) in enumerate(_PM.ref(pm, :branch_flow_cuts))
        _PMSC.constraint_c1_branch_contingency_ptdf_thermal_limit_from_soft(pm, i)
        _PMSC.constraint_c1_branch_contingency_ptdf_thermal_limit_to_soft(pm, i)
    end

    bus_withdrawal = _PM.var(pm, :bus_wdp)

    for (i,cut) in enumerate(_PM.ref(pm, :gen_flow_cuts))
        branch = _PM.ref(pm, :branch, cut.branch_id)
        gen = _PM.ref(pm, :gen, cut.gen_id)
        gen_bus = _PM.ref(pm, :bus, gen["gen_bus"])
        gen_set = _PM.ref(pm, :area_gens)[gen_bus["area"]]
        alpha_total = sum(gen["alpha"] for (i,gen) in _PM.ref(pm, :gen) if gen["index"] != cut.gen_id && i in gen_set)

        cont_bus_injection = Dict{Int,Any}()
        for (i, bus) in _PM.ref(pm, :bus)
            inj = 0.0
            for g in _PM.ref(pm, :bus_gens, i)
                if g != cut.gen_id
                    if g in gen_set
                        inj += _PM.var(pm, :pg, g) + gen["alpha"]*_PM.var(pm, :pg, cut.gen_id)/alpha_total
                    else
                        inj += _PM.var(pm, :pg, g)
                    end
                end
            end
            cont_bus_injection[i] = inj
        end

        #rate = branch["rate_a"]
        rate = branch["rate_c"]
        JuMP.@constraint(pm.model,  sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + _PM.var(pm, :gen_cont_flow_vio, i))
        JuMP.@constraint(pm.model, -sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + _PM.var(pm, :gen_cont_flow_vio, i))
    end

    for (i,gen_cont) in enumerate(_PM.ref(pm, :gen_contingencies))
        #println(gen_cont)
        gen = _PM.ref(pm, :gen, gen_cont.idx)
        gen_bus = _PM.ref(pm, :bus, gen["gen_bus"])
        gen_set = _PM.ref(pm, :area_gens)[gen_bus["area"]]
        response_gens = Dict(g => _PM.ref(pm, :gen, g) for g in gen_set if g != gen_cont.idx)

        # factor of 1.2 accounts for losses in a DC model
        #@constraint(pm.model, sum(gen["pmax"] - var(pm, :pg, g) for (g,gen) in response_gens) >= 1.2*var(pm, :pg, gen_cont.idx))
        JuMP.@constraint(pm.model, _PM.var(pm, :gen_cont_cap_vio, i) + sum(gen["pmax"] - _PM.var(pm, :pg, g) for (g,gen) in response_gens) >= _PM.var(pm, :pg, gen_cont.idx))
        #@constraint(pm.model, sum(gen["pmin"] - var(pm, :pg, g) for (g,gen) in response_gens) <= var(pm, :pg, gen_cont.idx))
    end

    ##### Setup Objective #####
    _PMSC.objective_c1_variable_pg_cost(pm)
    # explicit network id needed because of conductor-less
    pg_cost = _PM.var(pm, :pg_cost)
    branch_cont_flow_vio = _PM.var(pm, :branch_cont_flow_vio)
    gen_cont_flow_vio = _PM.var(pm, :gen_cont_flow_vio)
    gen_cont_cap_vio = _PM.var(pm, :gen_cont_cap_vio)

    JuMP.@objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in _PM.ref(pm, :gen) ) +
        sum( 5e5*branch_cont_flow_vio[i] for i in 1:length(_PM.ref(pm, :branch_flow_cuts)) ) +
        sum( 5e5*gen_cont_flow_vio[i] for i in 1:length(_PM.ref(pm, :gen_flow_cuts)) ) + 
        sum( 5e5*gen_cont_cap_vio[i] for i in 1:length(_PM.ref(pm, :gen_contingencies)) )
    )
end


function run_tnep_scopf_cuts_soft(file, model_constructor, solver; kwargs...)
    return _PM.solve_model(file, model_constructor, solver, build_tnep_scopf_cuts_soft; ref_extensions=[_PM.ref_add_on_off_va_bounds!, _PM.ref_add_ne_branch!, _PMSC.ref_c1!], kwargs...)
end

""
function build_tnep_scopf_cuts_soft(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)

    _PM.variable_ne_branch_indicator(pm)
    _PM.variable_ne_branch_power(pm)
    _PM.variable_ne_branch_voltage(pm)

    _PMSC.variable_c1_shunt_admittance_imaginary(pm)

    _PMSC.variable_c1_branch_contigency_power_violation(pm)
    _PMSC.variable_c1_gen_contigency_power_violation(pm)
    _PMSC.variable_c1_gen_contigency_capacity_violation(pm)

    for i in _PM.ids(pm, :bus)
        _PMSC.expression_c1_bus_generation(pm, i)
        _PMSC.expression_c1_bus_withdrawal(pm, i)
    end


    _PM.constraint_model_voltage(pm)
    _PM.constraint_ne_model_voltage(pm)

    for i in _PM.ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in _PM.ids(pm, :bus)
        constraint_ne_power_balance_shunt_dispatch(pm, i)
    end

    for i in _PM.ids(pm, :branch)
        _PMSC.constraint_goc_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)

        _PM.constraint_voltage_angle_difference(pm, i)

        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
    end

    for i in _PM.ids(pm, :ne_branch)
        _PM.constraint_ne_ohms_yt_from(pm, i)
        _PM.constraint_ne_ohms_yt_to(pm, i)

        _PM.constraint_ne_voltage_angle_difference(pm, i)

        _PM.constraint_ne_thermal_limit_from(pm, i)
        _PM.constraint_ne_thermal_limit_to(pm, i)
    end


    for (i,cut) in enumerate(_PM.ref(pm, :branch_flow_cuts))
        _PMSC.constraint_c1_branch_contingency_ptdf_thermal_limit_from_soft(pm, i)
        _PMSC.constraint_c1_branch_contingency_ptdf_thermal_limit_to_soft(pm, i)
    end

    bus_withdrawal = _PM.var(pm, :bus_wdp)

    for (i,cut) in enumerate(_PM.ref(pm, :gen_flow_cuts))
        branch = _PM.ref(pm, :branch, cut.branch_id)
        gen = _PM.ref(pm, :gen, cut.gen_id)
        gen_bus = _PM.ref(pm, :bus, gen["gen_bus"])
        gen_set = _PM.ref(pm, :area_gens)[gen_bus["area"]]
        alpha_total = sum(gen["alpha"] for (i,gen) in _PM.ref(pm, :gen) if gen["index"] != cut.gen_id && i in gen_set)

        cont_bus_injection = Dict{Int,Any}()
        for (i, bus) in _PM.ref(pm, :bus)
            inj = 0.0
            for g in _PM.ref(pm, :bus_gens, i)
                if g != cut.gen_id
                    if g in gen_set
                        inj += _PM.var(pm, :pg, g) + gen["alpha"]*_PM.var(pm, :pg, cut.gen_id)/alpha_total
                    else
                        inj += _PM.var(pm, :pg, g)
                    end
                end
            end
            cont_bus_injection[i] = inj
        end

        #rate = branch["rate_a"]
        rate = branch["rate_c"]
        JuMP.@constraint(pm.model,  sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + _PM.var(pm, :gen_cont_flow_vio, i))
        JuMP.@constraint(pm.model, -sum( weight*(cont_bus_injection[bus_id] - bus_withdrawal[bus_id]) for (bus_id, weight) in cut.bus_injection) <= rate + _PM.var(pm, :gen_cont_flow_vio, i))
    end

    for (i,gen_cont) in enumerate(_PM.ref(pm, :gen_contingencies))
        #println(gen_cont)
        gen = _PM.ref(pm, :gen, gen_cont.idx)
        gen_bus = _PM.ref(pm, :bus, gen["gen_bus"])
        gen_set = _PM.ref(pm, :area_gens)[gen_bus["area"]]
        response_gens = Dict(g => _PM.ref(pm, :gen, g) for g in gen_set if g != gen_cont.idx)

        # factor of 1.2 accounts for losses in a DC model
        #@constraint(pm.model, sum(gen["pmax"] - var(pm, :pg, g) for (g,gen) in response_gens) >= 1.2*var(pm, :pg, gen_cont.idx))
        JuMP.@constraint(pm.model, _PM.var(pm, :gen_cont_cap_vio, i) + sum(gen["pmax"] - _PM.var(pm, :pg, g) for (g,gen) in response_gens) >= _PM.var(pm, :pg, gen_cont.idx))
        #@constraint(pm.model, sum(gen["pmin"] - var(pm, :pg, g) for (g,gen) in response_gens) <= var(pm, :pg, gen_cont.idx))
    end

    ##### Setup Objective #####
    _PMSC.objective_c1_variable_pg_cost(pm)
    # explicit network id needed because of conductor-less
    pg_cost = _PM.var(pm, :pg_cost)
    branch_cont_flow_vio = _PM.var(pm, :branch_cont_flow_vio)
    gen_cont_flow_vio = _PM.var(pm, :gen_cont_flow_vio)
    gen_cont_cap_vio = _PM.var(pm, :gen_cont_cap_vio)

    JuMP.@objective(pm.model, Min,
        sum( pg_cost[i] for (i,gen) in _PM.ref(pm, :gen) ) +
        sum( 5e5*branch_cont_flow_vio[i] for i in 1:length(_PM.ref(pm, :branch_flow_cuts)) ) +
        sum( 5e5*gen_cont_flow_vio[i] for i in 1:length(_PM.ref(pm, :gen_flow_cuts)) ) + 
        sum( 5e5*gen_cont_cap_vio[i] for i in 1:length(_PM.ref(pm, :gen_contingencies)) ) + 
        sum( branch["construction_cost"]*_PM.var(pm, :branch_ne, i) for (i,branch) in _PM.ref(pm, :ne_branch) )
    )
end

function constraint_ne_power_balance_shunt_dispatch(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_ne = _PM.ref(pm, nw, :ne_bus_arcs, i)
    bus_arcs_dc = _PM.ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = _PM.ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_storage = _PM.ref(pm, nw, :bus_storage, i)

    bus_shunts_const = _PM.ref(pm, :bus_shunts_const, i)
    bus_shunts_var = _PM.ref(pm, :bus_shunts_var, i)

    bus_pd = Dict(k => _PM.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PM.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs_const = Dict(k => _PM.ref(pm, :shunt, k, "gs") for k in bus_shunts_const)
    bus_bs_const = Dict(k => _PM.ref(pm, :shunt, k, "bs") for k in bus_shunts_const)

    constraint_ne_power_balance_shunt_dispatch(pm, nw, i, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
end

function constraint_ne_power_balance_shunt_dispatch(pm::_PM.AbstractACRModel, n::Int, i::Int, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vi = _PM.var(pm, n, :vi, i)
    vr = _PM.var(pm, n, :vr, i)
    p    = get(_PM.var(pm, n),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PM.var(pm, n),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PM.var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PM.var(pm, n),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PM.var(pm, n),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PM.var(pm, n),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PM.var(pm, n),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PM.var(pm, n),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(_PM.var(pm, n), :p_dc, Dict()); _PM._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(_PM.var(pm, n), :q_dc, Dict()); _PM._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    bs = get(_PM.var(pm, n), :bs, Dict()); _PM._check_var_keys(bs, bus_shunts_var, "reactive power", "shunt")

    p_ne = get(_PM.var(pm, n), :p_ne, Dict()); _PM._check_var_keys(p_ne, bus_arcs_ne, "active power", "ne_branch")
    q_ne = get(_PM.var(pm, n), :q_ne, Dict()); _PM._check_var_keys(q_ne, bus_arcs_ne, "reactive power", "ne_branch")


    # possibly can save 2x in function eval, but no the dominant runtime at this moment
    #vm_sqr = @variable(pm.model, start=1.0, base_name="$(0)_vm_sqr_$(i)")

    #JuMP.@constraint(pm.model, vm_sqr == vi^2 + vr^2)
    #cstr_p = JuMP.@constraint(pm.model, 0 == - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vm_sqr)
    #cstr_q = JuMP.@constraint(pm.model, 0 == - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vm_sqr + sum(bs[s]*vm_sqr for s in bus_shunts_var))

    cstr_p = JuMP.@constraint(pm.model, 0 == - sum(p_ne[a] for a in bus_arcs_ne) - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*(vi^2 + vr^2))
    cstr_q = JuMP.@constraint(pm.model, 0 == - sum(q_ne[a] for a in bus_arcs_ne) - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*(vi^2 + vr^2) + sum(bs[s]*(vi^2 + vr^2) for s in bus_shunts_var))

    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

function constraint_ne_power_balance_shunt_dispatch(pm::_PM.AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_ne, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_shunts_var, bus_pd, bus_qd, bus_gs_const, bus_bs_const)
    vm   = _PM.var(pm, n, :vm, i)
    p    = get(_PM.var(pm, n),    :p, Dict()); _PM._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(_PM.var(pm, n),    :q, Dict()); _PM._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(_PM.var(pm, n),   :pg, Dict()); _PM._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(_PM.var(pm, n),   :qg, Dict()); _PM._check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(_PM.var(pm, n),   :ps, Dict()); _PM._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(_PM.var(pm, n),   :qs, Dict()); _PM._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(_PM.var(pm, n),  :psw, Dict()); _PM._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(_PM.var(pm, n),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    qsw  = get(_PM.var(pm, n),  :qsw, Dict()); _PM._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(_PM.var(pm, n), :p_dc, Dict()); _PM._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(_PM.var(pm, n), :q_dc, Dict()); _PM._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    bs = get(_PM.var(pm, n), :bs, Dict()); _PM._check_var_keys(bs, bus_shunts_var, "reactive power", "shunt")

    p_ne = get(_PM.var(pm, n), :p_ne, Dict()); _PM._check_var_keys(p_ne, bus_arcs_ne, "active power", "ne_branch")
    q_ne = get(_PM.var(pm, n), :q_ne, Dict()); _PM._check_var_keys(q_ne, bus_arcs_ne, "reactive power", "ne_branch")

    cstr_p = JuMP.@constraint(pm.model, 0 == - sum(p_ne[a] for a in bus_arcs_ne) - sum(p[a] for a in bus_arcs) + sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs_const))*vm^2)
    cstr_q = JuMP.@constraint(pm.model, 0 == - sum(q_ne[a] for a in bus_arcs_ne) - sum(q[a] for a in bus_arcs) + sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs_const))*vm^2 + sum(bs[s]*vm^2 for s in bus_shunts_var))

    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end




function run_scopf_acdc_benders_cuts(network::Dict{String,<:Any}, model_type::Type, optimizer, setting; max_iter::Int=1000, time_limit::Float64=Inf)
    results = Dict{String, Any}()
    time_start = time()
 
    # initiallize solve
    result = _SPM.run_master_scopf_bf(network, model_type, optimizer, setting=setting);
    if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
        Memento.warn(_LOGGER, "initial master problem solve failed in run_scopf_acdc_benders_cuts with status $(result["termination_status"])")
    end
    Memento.info(_LOGGER, "objective: $(result["objective"])")
    _PM.update_data!(network, result["solution"]);
    _PMACDCsc.update_data_converter_setpoints!(network, result["solution"]);

    # iterations
    iteration = 1
    cuts_found = 1
    while cuts_found > 0
        time_start_iteration = time()

        cuts = _SPM.check_acdc_subproblems(network, model_type, optimizer, setting, cuts_limit=iteration)
        network["secured_contingencies"] = cuts.secured_contingencies
        cuts_found = length(cuts.benders_cuts)

        if cuts_found <= 0
            Memento.info(_LOGGER, "no benders cuts found acdc scopf fixed-point reached")
            break
        else
            Memento.info(_LOGGER, "found $(cuts_found) benders cuts")
        end

        append!(network["cuts"], cuts.benders_cuts)
        Memento.info(_LOGGER, "active benders cuts: $(length(network["cuts"]))")

        time_solve_start = time()
    
        # master solve
        result = _SPM.run_master_scopf_bf(network, model_type, optimizer, setting=setting);
        if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
            Memento.warn(_LOGGER, "master problem solve failed in run_scopf_acdc_benders_cuts with status $(result["termination_status"]), terminating fixed-point early")
            break
        end
        Memento.info(_LOGGER, "objective: $(result["objective"])")
        _PM.update_data!(network, result["solution"]);
        _PMACDCsc.update_data_converter_setpoints!(network, result["solution"]);

        time_iteration = time() - time_start_iteration
        time_remaining = time_limit - (time() - time_start)
        if time_remaining < time_iteration
            Memento.warn(_LOGGER, "insufficent time for next iteration, time remaining $(time_remaining), estimated iteration time $(time_iteration)")
            break
        end
        iteration += 1
    end
    results["iterations"] = iteration
    results["final"] = result 

    return results
end

function run_scopf_acdc_benders_cuts_soft(network::Dict{String,<:Any}, model_type::Type, optimizer, setting; max_iter::Int=1000, time_limit::Float64=Inf)
    results = Dict{String, Any}()
    time_start = time()
 
    mobj = 0.0
    # initiallize solve
    result = _SPM.run_master_scopf_bf(network, model_type, optimizer, setting=setting);
    if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
        Memento.warn(_LOGGER, "initial master problem solve failed in run_scopf_acdc_benders_cuts with status $(result["termination_status"])")
    end
    Memento.info(_LOGGER, "initial master problem objective: $(result["objective"])")
    _PM.update_data!(network, result["solution"]);
    _PMACDCsc.update_data_converter_setpoints!(network, result["solution"]);
    mobj = result["objective"]

    # iterations
    iteration = 1
    cuts_found = 1
    while cuts_found > 0
        time_start_iteration = time()

        cuts = _SPM.check_acdc_subproblems_soft(network, model_type, optimizer, setting, cuts_limit=iteration)

        # convergence check
        if sum( ((cuts.sub_obj[val] - mobj)/cuts.sub_obj[val]) for val in eachindex(cuts.sub_obj)) < 0.001
            Memento.info(_LOGGER, "modified benders decoposition converged, fixed-point reached")
            break
        end
        time_iteration = time() - time_start_iteration
        time_remaining = time_limit - (time() - time_start)
        if time_remaining < time_iteration
            Memento.warn(_LOGGER, "insufficent time for next iteration, time remaining $(time_remaining), estimated iteration time $(time_iteration)")
            break
        end
        if iteration >= max_iter
            Memento.warn(_LOGGER, "maximum iterations reached $(max_iter), terminating fixed-point early")
            break
        end

        cuts_found = length(cuts.benders_cuts)
        if cuts_found <= 0
            Memento.info(_LOGGER, "no benders cuts found acdc scopf fixed-point reached")
            break
        else
            Memento.info(_LOGGER, "found $(cuts_found) benders cuts")
        end

        append!(network["cuts"], cuts.benders_cuts)
        Memento.info(_LOGGER, "active benders cuts: $(length(network["cuts"]))")

            
        # master solve
        result = _SPM.run_master_scopf_bf_soft(network, model_type, optimizer, setting=setting);
        if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
            Memento.warn(_LOGGER, "master problem solve failed in run_scopf_acdc_benders_cuts with status $(result["termination_status"]), terminating fixed-point early")
            break
        end
        Memento.info(_LOGGER, "master problem objective: $(result["objective"])")
        _PM.update_data!(network, result["solution"]);
        _PMACDCsc.update_data_converter_setpoints!(network, result["solution"]);
        mobj = result["objective"]


        iteration += 1
    end
    results["iterations"] = iteration
    results["final"] = result 

    return results
end

function run_scopf_acdc_nc_benders_cuts_soft(network::Dict{String,<:Any}, model_type::Type, optimizer, setting; max_iter::Int=1000, time_limit::Float64=Inf)
    results = Dict{String, Any}()
    time_start = time()
 
    mobj = 0.0
    # initiallize solve
    result = _SPM.run_master_scopf_bf(network, model_type, optimizer, setting=setting);
    if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
        Memento.warn(_LOGGER, "initial master problem solve failed in run_scopf_acdc_benders_cuts with status $(result["termination_status"])")
    end
    Memento.info(_LOGGER, "initial master problem objective: $(result["objective"])")
    _PM.update_data!(network, result["solution"]);
    _PMACDCsc.update_data_converter_setpoints!(network, result["solution"]);
    mobj = result["objective"]

    outer_iteration = 1
    while outer_iteration > 0

        # inner iterations
        inner_iteration = 1
        cuts_found = 1
        while cuts_found > 0
            time_start_iteration = time()

            cuts = _SPM.check_acdc_subproblems_soft(network, model_type, optimizer, setting, cuts_limit=inner_iteration)

            # convergence check
            if sum( ((cuts.sub_obj[val] - mobj)/cuts.sub_obj[val]) for val in eachindex(cuts.sub_obj)) < 0.001
                Memento.info(_LOGGER, "modified benders decoposition converged, fixed-point reached")
                break
            end
            time_iteration = time() - time_start_iteration
            time_remaining = time_limit - (time() - time_start)
            if time_remaining < time_iteration
                Memento.warn(_LOGGER, "insufficent time for next iteration, time remaining $(time_remaining), estimated iteration time $(time_iteration)")
                break
            end
            if inner_iteration >= max_iter
                Memento.warn(_LOGGER, "maximum iterations reached $(max_iter), terminating fixed-point early")
                break
            end

            cuts_found = length(cuts.benders_cuts)
            if cuts_found <= 0
                Memento.info(_LOGGER, "no benders cuts found acdc scopf fixed-point reached")
                break
            else
                Memento.info(_LOGGER, "found $(cuts_found) benders cuts")
            end

            append!(network["cuts"], cuts.benders_cuts)
            Memento.info(_LOGGER, "active benders cuts: $(length(network["cuts"]))")

                
            # master solve
            result = _SPM.run_master_scopf_bf_soft(network, model_type, optimizer, setting=setting);
            if !(result["termination_status"] == _PM.OPTIMAL || result["termination_status"] == _PM.LOCALLY_SOLVED || result["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
                Memento.warn(_LOGGER, "master problem solve failed in run_scopf_acdc_benders_cuts with status $(result["termination_status"]), terminating fixed-point early")
                break
            end
            Memento.info(_LOGGER, "master problem objective: $(result["objective"])")
            _PM.update_data!(network, result["solution"]);
            _PMACDCsc.update_data_converter_setpoints!(network, result["solution"]);
            mobj = result["objective"]


            inner_iteration += 1
        end
        results["inner_iterations"] = inner_iteration
        results["final"] = result 

        # check nonlinear subproblems
        _PM.update_data!(network, results["final"]["solution"]);
        _PMACDCsc.update_data_converter_setpoints!(network, results["final"]["solution"]);
        network["soc_master_obj"] = results["final"]["objective"]

        result_nc_sub = _SPM.check_nc_acdc_subproblems_soft(network, _PM.ACPPowerModel, optimizer, setting)
        argmin_result = select_argmin_solution(result_nc_sub);
        _PM.update_data!(network, argmin_result["solution"]);
        _PMACDCsc.update_data_converter_setpoints!(network, argmin_result["solution"]);
        mobj = argmin_result["objective"]

        # convergence check
        if sum( ((r["objective"] - results["final"]["objective"])/r["objective"]) for r in result_nc_sub) < 0.1
            Memento.info(_LOGGER, "nonconvex benders decoposition converged, fixed-point reached")
            break
        end

        outer_iteration += 1

    end
    results["outer_iterations"] = outer_iteration
    return results
end






function check_acdc_subproblems(network, model_type, optimizer, setting; cuts_limit=typemax(Int64), gen_eval_limit=typemax(Int64),
    branch_eval_limit=typemax(Int64), branchdc_eval_limit=typemax(Int64), convdc_eval_limit=typemax(Int64))    

    time_start_subproblems = time()
        
    gen_contingencies = _PMSC.calc_c1_gen_contingency_subset(network, gen_eval_limit=gen_eval_limit)
    branch_contingencies = _PMSC.calc_c1_branch_contingency_subset(network, branch_eval_limit=branch_eval_limit)
    branchdc_contingencies = _PMACDCsc.calc_branchdc_contingency_subset(network, branchdc_eval_limit=branchdc_eval_limit)   
    convdc_contingencies = _PMACDCsc.calc_convdc_contingency_subset(network, convdc_eval_limit=convdc_eval_limit)

    benders_cuts = [];
    for (i,cont) in enumerate(gen_contingencies)
        if !(cont in network["secured_contingencies"])
            # if length(benders_cuts) >= cuts_limit
            #     Memento.info(_LOGGER, "hit total cut limit $(cuts_limit)")           
            #     break
            # end

            cont_gen = network["gen"]["$(cont.idx)"]
            pg_lost = cont_gen["pg"]
            cont_gen["gen_status"] = 0
            cont_gen["pg"] = 0.0

            result_sub = _SPM.run_sub_scopf_bf(network, model_type, optimizer, setting=setting);
            if (result_sub["termination_status"] == _PM.OPTIMAL || result_sub["termination_status"] == _PM.LOCALLY_SOLVED || result_sub["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
                Memento.info(_LOGGER, "subproblem solve is feasible in check_acdc_subproblems, status $(result_sub["termination_status"]): contingency $(cont.label) is secured")
                push!(network["secured_contingencies"], cont)
            end
            Memento.info(_LOGGER, "objective: $(result_sub["objective"])")
            
            # cuts
            if !(cont in network["secured_contingencies"])
                lm_p = Dict(parse(Int, i) => gen["lm_p"] for (i, gen) in result_sub["solution"]["gen"]);
                lm_q = Dict(parse(Int, i) => gen["lm_q"] for (i, gen) in result_sub["solution"]["gen"]);
                lm_ac_p = Dict(parse(Int, i) => conv["lm_ac_p"] for (i, conv) in result_sub["solution"]["convdc"]);
                lm_ac_q = Dict(parse(Int, i) => conv["lm_ac_q"] for (i, conv) in result_sub["solution"]["convdc"]);
                obj = result_sub["objective"];

                cut = (obj=obj, lm_p=lm_p, lm_q=lm_q, lm_ac_p=lm_ac_p, lm_ac_q=lm_ac_q);
                
                push!(benders_cuts, cut);
            end

            cont_gen["gen_status"] = 1
            cont_gen["pg"] = pg_lost
        end
    end

    for (i,cont) in enumerate(branch_contingencies)
        if !(cont in network["secured_contingencies"])
            # if length(benders_cuts) >= cuts_limit
            #     Memento.info(_LOGGER, "hit total cut limit $(cuts_limit)")                                    
            #     break
            # end

            cont_branch = network["branch"]["$(cont.idx)"]
            cont_branch["br_status"] = 0
            _PMACDC.fix_data!(network)

            result_sub = _SPM.run_sub_scopf_bf(network, model_type, optimizer, setting=setting);
            if (result_sub["termination_status"] == _PM.OPTIMAL || result_sub["termination_status"] == _PM.LOCALLY_SOLVED || result_sub["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
                Memento.info(_LOGGER, "subproblem solve is feasible in check_acdc_subproblems, status $(result_sub["termination_status"]): contingency $(cont.label) is secured")
                push!(network["secured_contingencies"], cont)
            end
            Memento.info(_LOGGER, "objective: $(result_sub["objective"])")

            # cuts
            if !(cont in network["secured_contingencies"])
                lm_p = Dict(parse(Int, i) => gen["lm_p"] for (i, gen) in result_sub["solution"]["gen"]);
                lm_q = Dict(parse(Int, i) => gen["lm_q"] for (i, gen) in result_sub["solution"]["gen"]);
                lm_ac_p = Dict(parse(Int, i) => conv["lm_ac_p"] for (i, conv) in result_sub["solution"]["convdc"]);
                lm_ac_q = Dict(parse(Int, i) => conv["lm_ac_q"] for (i, conv) in result_sub["solution"]["convdc"]);
                obj = result_sub["objective"];

                cut = (obj=obj, lm_p=lm_p, lm_q=lm_q, lm_ac_p=lm_ac_p, lm_ac_q=lm_ac_q);
                
                push!(benders_cuts, cut);
            end

            cont_branch["br_status"] = 1
        end
    end

    time_contingencies = time() - time_start_subproblems
    Memento.info(_LOGGER, "subproblems eval time: $(time_contingencies)")            

    return (benders_cuts=benders_cuts, secured_contingencies=network["secured_contingencies"])
end



function check_acdc_subproblems_soft(network, model_type, optimizer, setting; cuts_limit=typemax(Int64), gen_eval_limit=typemax(Int64),
    branch_eval_limit=typemax(Int64), branchdc_eval_limit=typemax(Int64), convdc_eval_limit=typemax(Int64))    

    time_start_subproblems = time()
        
    gen_contingencies = _PMSC.calc_c1_gen_contingency_subset(network, gen_eval_limit=gen_eval_limit)
    branch_contingencies = _PMSC.calc_c1_branch_contingency_subset(network, branch_eval_limit=branch_eval_limit)
    branchdc_contingencies = _PMACDCsc.calc_branchdc_contingency_subset(network, branchdc_eval_limit=branchdc_eval_limit)   
    convdc_contingencies = _PMACDCsc.calc_convdc_contingency_subset(network, convdc_eval_limit=convdc_eval_limit)

    benders_cuts = [];
    sub_obj = [];
    for (i,cont) in enumerate(gen_contingencies)
     
        cont_gen = network["gen"]["$(cont.idx)"]
        pg_lost = cont_gen["pg"]
        cont_gen["gen_status"] = 0
        cont_gen["pg"] = 0.0

        result_sub = _SPM.run_sub_scopf_bf_soft(network, model_type, optimizer, setting=setting);
        if !(result_sub["termination_status"] == _PM.OPTIMAL || result_sub["termination_status"] == _PM.LOCALLY_SOLVED || result_sub["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
            Memento.info(_LOGGER, "subproblem solve is infeasible in check_acdc_subproblems, status $(result_sub["termination_status"])")
        end
        Memento.info(_LOGGER, "subproblem objective: $(result_sub["objective"])")
        push!(sub_obj, result_sub["objective"])
            
        # cuts
        lm_p = Dict(parse(Int, i) => gen["lm_p"] for (i, gen) in result_sub["solution"]["gen"]);
        lm_q = Dict(parse(Int, i) => gen["lm_q"] for (i, gen) in result_sub["solution"]["gen"]);
        lm_ac_p = Dict(parse(Int, i) => conv["lm_ac_p"] for (i, conv) in result_sub["solution"]["convdc"]);
        lm_ac_q = Dict(parse(Int, i) => conv["lm_ac_q"] for (i, conv) in result_sub["solution"]["convdc"]);
        obj = result_sub["objective"];
                 
        cut = (cont_id = cont.idx, obj=obj, lm_p=lm_p, lm_q=lm_q, lm_ac_p=lm_ac_p, lm_ac_q=lm_ac_q);
                
        push!(benders_cuts, cut);

        cont_gen["gen_status"] = 1
        cont_gen["pg"] = pg_lost
    end

    for (i,cont) in enumerate(branch_contingencies)
        if !(cont in network["secured_contingencies"])
        
            cont_branch = network["branch"]["$(cont.idx)"]
            cont_branch["br_status"] = 0
            _PMACDC.fix_data!(network)

            result_sub = _SPM.run_sub_scopf_bf_soft(network, model_type, optimizer, setting=setting);
            if !(result_sub["termination_status"] == _PM.OPTIMAL || result_sub["termination_status"] == _PM.LOCALLY_SOLVED || result_sub["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
                Memento.info(_LOGGER, "subproblem solve is infeasible in check_acdc_subproblems, status $(result_sub["termination_status"])")
            end
            Memento.info(_LOGGER, "subproblem objective: $(result_sub["objective"])")
            push!(sub_obj, result_sub["objective"])

            # cuts
            if !(cont in network["secured_contingencies"])
                lm_p = Dict(parse(Int, i) => gen["lm_p"] for (i, gen) in result_sub["solution"]["gen"]);
                lm_q = Dict(parse(Int, i) => gen["lm_q"] for (i, gen) in result_sub["solution"]["gen"]);
                lm_ac_p = Dict(parse(Int, i) => conv["lm_ac_p"] for (i, conv) in result_sub["solution"]["convdc"]);
                lm_ac_q = Dict(parse(Int, i) => conv["lm_ac_q"] for (i, conv) in result_sub["solution"]["convdc"]);
                obj = result_sub["objective"];

                cut = (cont_id = cont.idx, obj=obj, lm_p=lm_p, lm_q=lm_q, lm_ac_p=lm_ac_p, lm_ac_q=lm_ac_q);
                        
                push!(benders_cuts, cut);
            end
            
            cont_branch["br_status"] = 1
        end
    end

    time_contingencies = time() - time_start_subproblems
    Memento.info(_LOGGER, "subproblems eval time: $(time_contingencies)")            

    return (benders_cuts=benders_cuts, sub_obj=sub_obj)  
end


function check_nc_acdc_subproblems_soft(network, model_type, optimizer, setting; cuts_limit=typemax(Int64), gen_eval_limit=typemax(Int64),
    branch_eval_limit=typemax(Int64), branchdc_eval_limit=typemax(Int64), convdc_eval_limit=typemax(Int64))    

    time_start_subproblems = time()
        
    gen_contingencies = _PMSC.calc_c1_gen_contingency_subset(network, gen_eval_limit=gen_eval_limit)
    branch_contingencies = _PMSC.calc_c1_branch_contingency_subset(network, branch_eval_limit=branch_eval_limit)
    branchdc_contingencies = _PMACDCsc.calc_branchdc_contingency_subset(network, branchdc_eval_limit=branchdc_eval_limit)   
    convdc_contingencies = _PMACDCsc.calc_convdc_contingency_subset(network, convdc_eval_limit=convdc_eval_limit)

    
    result = [];
    for (i,cont) in enumerate(gen_contingencies)
     
        cont_gen = network["gen"]["$(cont.idx)"]
        pg_lost = cont_gen["pg"]
        cont_gen["gen_status"] = 0
        cont_gen["pg"] = 0.0

        result_sub = _SPM.run_nc_sub_scopf_soft(network, model_type, optimizer, setting=setting);
        if !(result_sub["termination_status"] == _PM.OPTIMAL || result_sub["termination_status"] == _PM.LOCALLY_SOLVED || result_sub["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
            Memento.info(_LOGGER, "subproblem solve is infeasible in check_acdc_subproblems, status $(result_sub["termination_status"])")
        end
        Memento.info(_LOGGER, "nonconvex subproblem objective: $(result_sub["objective"])")
        push!(result, result_sub)


        cont_gen["gen_status"] = 1
        cont_gen["pg"] = pg_lost
    end

    for (i,cont) in enumerate(branch_contingencies)
        if !(cont in network["secured_contingencies"])
        
            cont_branch = network["branch"]["$(cont.idx)"]
            cont_branch["br_status"] = 0
            _PMACDC.fix_data!(network)

            result_sub = _SPM.run_nc_sub_scopf_soft(network, model_type, optimizer, setting=setting);
            if !(result_sub["termination_status"] == _PM.OPTIMAL || result_sub["termination_status"] == _PM.LOCALLY_SOLVED || result_sub["termination_status"] == _PM.ALMOST_LOCALLY_SOLVED)
                Memento.info(_LOGGER, "subproblem solve is infeasible in check_acdc_subproblems, status $(result_sub["termination_status"])")
            end
            Memento.info(_LOGGER, "nonconvex subproblem objective: $(result_sub["objective"])")
            push!(result, result_sub)

            
            cont_branch["br_status"] = 1
        end
    end

    time_contingencies = time() - time_start_subproblems
    Memento.info(_LOGGER, "subproblems eval time: $(time_contingencies)")            

    return (result=result)  
end



function select_argmin_solution(result)
    X = []
    Y = []
    for r in result
        push!(X, r["objective"])
        push!(Y, r)
    end
    X = sort(X)
    perm = sortperm(X)
    Y[perm]
    return Y[1]
end