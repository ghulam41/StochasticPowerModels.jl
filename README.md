# StochasticPowerModels

<a href="https://github.com/timmyfaraday/StochasticPowerModels.jl/actions?query=workflow%3ACI"><img src="https://github.com/timmyfaraday/StochasticPowerModels.jl/workflows/CI/badge.svg"></img></a>
<a href="https://codecov.io/gh/timmyfaraday/StochasticPowerModels.jl"><img src="https://img.shields.io/codecov/c/github/timmyfaraday/StochasticPowerModels.jl?logo=Codecov"></img></a>

StochasticPowerModels.jl is an extension package of PowerModels.jl for 
Stochastic (Optimal) Power Flow. It is designed to enable inclusion of 
uncertainty in Steady-State Power Network Optimization. 

## Core Problem Specification

- Stochastic Optimal Power Flow (sOPF)

## Core Network Formulation

- Exact
    - ACR
    - IVR (reduced)

## Core Stochastic Specification

- Polynomial Chaos Expansion
    - with/without auxiliary variables/constraints

## Network Data Formats

- Matpower ".m" files

## Stochastic Data Extension

- Matpower ".m" files, extended to include:
    - stochastic germ: `mpc.sdata`,
    - stochastic bus data: `mpc.bus_sdata`, including: `dst_id`, `μ`, `σ`, `λvmin` and `λvmax`,
    - stochastic gen data: `mpc.gen_sdata`, including: `λpmin`, `λpmax`, `λqmin` and `λqmax`, and
    - stochastic branch data: `mpc.branch_sdata`, including: `λcmax`.

## Installation

The latest stable release of StochasticPowerModels can be installed using the 
Julia package manager:

```
] add https://github.com/timmyfaraday/StochasticPowerModels.jl.git
```

In order to test whether the package works, run:

```
] test StochasticPowerModels
```

## Acknowledgements

The primary developer is Tom Van Acker ([@timmyfaraday](https://github.com/timmyfaraday)), 
with support from the following contributors:
- Arpan Koirala ([@arpkoirala](https://github.com/arpkoirala)), KU Leuven, ACR formulation
- Frederik Geth ([@frederikgeth](https://github.com/frederikgeth)), CSIRO, reduced IVR formulation

## License

This code is provided under a BSD license.
