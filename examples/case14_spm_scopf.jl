
using JuMP
using Ipopt
using Gurobi
using Juniper
using PowerModels
using PowerModelsACDC
using StochasticPowerModels
using PowerModelsSecurityConstrained

# constants 
const _PM = PowerModels
const _PMACDC = PowerModelsACDC
const _SPM = StochasticPowerModels
const _PMSC = PowerModelsSecurityConstrained

ENV["GUROBI_HOME"] = "C:\\gurobi952\\win64"

ENV["GRB_LICENSE_FILE"] = "C:\\gurobi952lic\\gurobi.lic"

# solvers
nlp_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-4, "print_level"=>0), "log_levels"=>[])
gurobi_solver = optimizer_with_attributes(Gurobi.Optimizer)
# input
deg  = 2
# case = "case14_ieee_spm.m"
# case = "case5_spm.m"
case = "case5_spm_tnep.m"

# data
file  = joinpath(BASE_DIR, "test/data/matpower", case)
data  = _PM.parse_file(file)

for (b,branch) in data["branch"]
    branch["rate_c"] = branch["rate_c"]*0.6
end

data["gen_contingencies"] = Vector{Any}(undef, 4)

data["gen_contingencies"][1] = (idx = 1, label = "GEN-1-1", type = "gen")
data["gen_contingencies"][2] = (idx = 2, label = "GEN-2-1", type = "gen")
data["gen_contingencies"][3] = (idx = 1, label = "GEN-3-1", type = "gen")
data["gen_contingencies"][4] = (idx = 2, label = "GEN-4-1", type = "gen")

# data["gen_contingencies"][1] = (idx = 3, label = "GEN-3-1", type = "gen")
# data["gen_contingencies"][2] = (idx = 4, label = "GEN-4-1", type = "gen")
# data["gen_contingencies"][3] = (idx = 5, label = "GEN-5-1", type = "gen")
# data["gen_contingencies"][4] = (idx = 2, label = "GEN-2-1", type = "gen")

data["branch_contingencies"] = Vector{Any}(undef, 6)

data["branch_contingencies"][1] = (idx = 1, label = "LINE-1-2-BL", type = "branch")
data["branch_contingencies"][2] = (idx = 2, label = "LINE-1-4-BL", type = "branch")
data["branch_contingencies"][3] = (idx = 3, label = "LINE-1-5-BL", type = "branch")
data["branch_contingencies"][4] = (idx = 4, label = "LINE-2-3-BL", type = "branch")
data["branch_contingencies"][5] = (idx = 5, label = "LINE-3-4-BL", type = "branch")
data["branch_contingencies"][6] = (idx = 6, label = "LINE-4-5-BL", type = "branch")

# data["branch_contingencies"][1] = (idx = 12, label = "LINE-6-12-BL", type = "branch")
# data["branch_contingencies"][2] = (idx = 2, label = "LINE-1-5-BL", type = "branch")
# data["branch_contingencies"][3] = (idx = 3, label = "LINE-2-3-BL", type = "branch")
# data["branch_contingencies"][4] = (idx = 4, label = "LINE-2-4-BL", type = "branch")
# data["branch_contingencies"][2] = (idx = 1, label = "LINE-1-2-BL", type = "branch")

data["gen_contingencies_active"] = []
data["branch_contingencies_active"] = []

data["area_gens"] = Dict{Int64, Set{Int64}}()

data["area_gens"][1] = Set{Int64}([4, 2, 3, 1])

# data["area_gens"][1] = Set{Int64}([5, 4, 2, 3, 1])

data["gen"]["1"]["alpha"] = 5.0
data["gen"]["2"]["alpha"] = 19.0
data["gen"]["3"]["alpha"] = 49.25
data["gen"]["4"]["alpha"] = 38.75

data["ne_branch"]["1"]["cmax"] = 4.44444
data["ne_branch"]["2"]["cmax"] = 4.73333
data["ne_branch"]["3"]["cmax"] = 0.01111

data["ne_branch"]["1"]["λcmax"] = 1.03643
data["ne_branch"]["2"]["λcmax"] = 1.03643
data["ne_branch"]["3"]["λcmax"] = 1.03643

# data["gen"]["1"]["alpha"] = 5.0
# data["gen"]["2"]["alpha"] = 19.0
# data["gen"]["3"]["alpha"] = 49.25
# data["gen"]["4"]["alpha"] = 38.75
# data["gen"]["5"]["alpha"] = 3.0

#-----------------------------------
# run the convenience functions for stochastic SCOPF for ACR

# result_sscopf = _SPM.run_sscopf_cuts(data, _PM.ACRPowerModel, nlp_solver, deg=deg)

# result_ssctnep = _SPM.run_ssctnep_cuts(data, _PM.ACRPowerModel, minlp_solver, deg=deg)

result_scopf = _SPM.run_scopf_cuts(data, _PM.ACPPowerModel, nlp_solver)

# result_sctnep = _SPM.run_sctnep_cuts(data, _PM.ACPPowerModel, minlp_solver)


## soc
case = "case5_ACDC.m"
file  = joinpath(BASE_DIR, "test/data/matpower", case)
data  = _PM.parse_file(file)
_PMACDC.process_additional_data!(data)
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true) 

result = _PMACDC.run_acdcopf(data, _PM.SOCWRPowerModel, nlp_solver, setting=s)
result = _PMACDC.run_acdcopf_bf(data, _PM.SOCBFPowerModel, nlp_solver, setting=s)
result = _SPM.run_master_bf(data, _PM.SOCBFPowerModel, nlp_solver, setting=s)



##


solution = result["solution"]["nw"]["1"]
solution["per_unit"] = result["solution"]["per_unit"]
_PM.update_data!(network, solution)

network["gen"]["2"]["gen_status"] = 0

result_pf = _PM.solve_pf(network, _PM.ACRPowerModel, ipopt_solver)

for (b,bus) in result_pf["solution"]["bus"]
    result_pf["solution"]["bus"][b]["vm"] = sqrt((bus["vr"])^2 + (bus["vi"])^2)
    result_pf["solution"]["bus"][b]["va"] = atan((bus["vi"])/(bus["vr"]))
end

vm = Dict(parse(Int, i) => bus["vm"] for (i,bus) in result_pf["solution"]["bus"])
va = Dict(parse(Int, i) => bus["va"] for (i,bus) in result_pf["solution"]["bus"])

flows = Dict{String,Any}()
for (i,branch) in data["branch"]
    if branch["br_status"] != 0
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]

        g, b = calc_branch_y(branch)
        tr, ti = calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]

        tm = branch["tap"]

        vm_fr = vm[f_bus]
        vm_to = vm[t_bus]
        va_fr = va[f_bus]
        va_to = va[t_bus]

        p_fr =  (g+g_fr)/tm^2*vm_fr^2 + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))
        q_fr = -(b+b_fr)/tm^2*vm_fr^2 - (-b*tr-g*ti)/tm^2*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm^2*(vm_fr*vm_to*sin(va_fr-va_to))

        p_to =  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))
        q_to = -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm^2*(vm_to*vm_fr*cos(va_to-va_fr)) + (-g*tr-b*ti)/tm^2*(vm_to*vm_fr*sin(va_to-va_fr))
    else
        p_fr = NaN
        q_fr = NaN

        p_to = NaN
        q_to = NaN
    end
    flows[i] = Dict(
        "pf" => p_fr,
        "qf" => q_fr,
        "pt" => p_to,
        "qt" => q_to
    )
end
flows = Dict{String,Any}("branch" => flows)

using Plots
pf_scopf = [branch["pf"] for (b, branch) in result["solution"]["nw"]["1"]["branch"]]
pf_pf = [branch["pf"] for (b, branch) in flows["branch"]]
rate_c = [branch["rate_c"] for (b, branch) in data["branch"]]

scatter(pf_scopf)
scatter!(pf_pf)
plot!(rate_c)

result_sopf = _SPM.solve_sopf_acr(data, _PM.ACRPowerModel, ipopt_solver, deg=deg)







pg_opf = [gen["pg"] for (g, gen) in result_sopf["solution"]["nw"]["1"]["gen"]]
pg_scopf = [gen["pg"] for (g, gen) in result["solution"]["nw"]["1"]["gen"]]

pf_opf = [branch["pf"] for (b, branch) in result_sopf["solution"]["nw"]["1"]["branch"]]
pf_scopf = [branch["pf"] for (b, branch) in result["solution"]["nw"]["1"]["branch"]]
rate_c = [branch["rate_c"] for (b, branch) in data["branch"]]

scatter(pg_opf)
scatter!(pg_scopf)

scatter(pf_opf)
scatter!(pf_scopf)
plot!(rate_c)