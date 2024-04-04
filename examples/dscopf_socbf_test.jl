
using JuMP
using Ipopt
using PowerModels
using PowerModelsACDC
using StochasticPowerModels
using PowerModelsSecurityConstrained
using PowerModelsACDCsecurityconstrained
using SCS
using Gurobi


# constants 
const _PM = PowerModels;
const _PMACDC = PowerModelsACDC;
const _SPM = StochasticPowerModels;
const _PMSC = PowerModelsSecurityConstrained;
const _PMACDCsc = PowerModelsACDCsecurityconstrained;

ENV["GUROBI_HOME"] = "C:\\gurobi952\\win64"
ENV["GRB_LICENSE_FILE"] = "C:\\gurobi952lic\\gurobi.lic"

# solvers
nlp_solver = optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=> 100000, "print_level"=>0);
sdp_solver = JuMP.optimizer_with_attributes(SCS.Optimizer, "verbose"=>false);
gurobi_solver = optimizer_with_attributes(Gurobi.Optimizer);

_PMACDCsc.silence()

# data
case = "case5_ACDC.m";
file  = joinpath(BASE_DIR, "test/data/matpower", case);
data  = _PM.parse_file(file);
_PMACDC.process_additional_data!(data);
data["gen_contingencies"] = Vector{Any}(undef, 2)
data["gen_contingencies"][1] = (idx = 2, label = "gen2_bus2", type = "gen")
data["gen_contingencies"][2] = (idx = 1, label = "gen1_bus1", type = "gen")

data["branchdc_contingencies"] = []
data["convdc_contingencies"] = []
data["branch_contingencies"] = Vector{Any}(undef, 5)
data["branch_contingencies"][1] = (idx = 2, label = "branch_bus1_bus3", type = "branch")
data["branch_contingencies"][2] = (idx = 4, label = "branch_bus2_bus4", type = "branch")
data["branch_contingencies"][3] = (idx = 6, label = "branch_bus3_bus4", type = "branch")
data["branch_contingencies"][4] = (idx = 7, label = "branch_bus4_bus5", type = "branch")
data["branch_contingencies"][5] = (idx = 5, label = "branch_bus2_bus5", type = "branch")


s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true);
data["cuts"] = [];
data["secured_contingencies"] = []; # not for soft

data["RES"] = []
pen_level = 0
p_size = 0

result = _SPM.run_acdcopf_iv(data, _PM.IVRPowerModel, nlp_solver, setting = s)
result = _PMACDC.run_acdcopf_iv(data, _PM.IVRPowerModel, nlp_solver, setting = s)


deg  = 2
case = "case5_spm.m"

# data
file  = joinpath(BASE_DIR, "test/data/matpower", case)
result = _SPM.solve_sopf_iv_n(file, _PM.IVRPowerModel, gurobi_solver, deg=deg)
result = _SPM.solve_sopf_iv(file, _PM.IVRPowerModel, nlp_solver, deg=deg)


result = _SPM.run_scopf_acdc_benders_cuts(data, _PM.SOCBFPowerModel, nlp_solver, s)

result = _SPM.run_scopf_acdc_benders_cuts_soft(data, _PM.SOCBFPowerModel, nlp_solver, s)


data["contingencies"] = []
_PMACDCsc.fix_scopf_data_case5_acdc!(data)
resultsc = _PMACDCsc.run_scopf_acdc_cuts(data, _PM.ACPPowerModel, _PM.ACPPowerModel, _PMACDCsc.run_acdc_scopf_cuts_soft, nlp_solver, nlp_solver, s)