################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# StochasticPowerModels.jl                                                     #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################

# using pkgs
using JuMP
using Ipopt
using PowerModels
using StochasticPowerModels

# constants 
const _PMs = PowerModels
const _SPM = StochasticPowerModels

# data
deg  = 1
aux  = false
path = joinpath(_SPM.BASE_DIR,"test/data/matpower/case30_spm_muhlpfordt.m")
data = _PMs.parse_file(path)

# initialize solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer)

# solve problem
result_stc = run_sopf_iv(data, _PMs.IVRPowerModel, solver, aux=aux, deg=deg)

# solve reduced problem
result_stc = run_sopf_iv(data, _PMs.IVRPowerModel, solver, aux=aux, deg=deg)

# solve problem iteratively
(result_dtr, result_itr) = run_sopf_iv_itr(data, _PMs.IVRPowerModel, solver, aux=aux, deg=deg);

# assert
@assert isapprox(result_stc["objective"], result_itr["objective"], rtol=1e-5)