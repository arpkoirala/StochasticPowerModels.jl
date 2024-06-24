# Dynamic Operating Envelop using Polynomial chaos

StochasticPowerModels.jl is an extension package of PowerModels.jl for 
Stochastic (Optimal) Power Flow. It is designed to enable inclusion of 
uncertainty in Steady-State Power Network Optimization. 

In this project file, the base code of StochasticPowerModels.jl is utilized to obtain Stochastic PV hosting capacity of LVDS. Activating the environment should be enough using:
 ```
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
```
## Core Problem Specification

- Stochastic multi-conductor Optimal Power Flow (sOPF-MC) HC
- Hosting capacity for multi-conductor (HC-mc)
- Dynamic operating Envelop (DoE)


## Varients of DoE

- Efficiency measures
- Equality measures
- Alpha fairmness

## Core Stochastic Specification
For now, we only support Polynomial Chaos Expansion based stochastic HC and deterministic HC in IV formulation.

- Polynomial Chaos Expansion
    - with auxiliary variables/constraints

## Example in DoE calculation
- For an example in different type of HC calculation:
    - `examples/case_test_doe.jl` is pointed out where the codes for DoE paper is present:
    
    Koirala, Arpan; Geth, Frederik, Van Aacker, Tom (2023): Determining Dynamic Operating Envelopes Using Stochastic Unbalanced Optimal Power Flow 

The Folder /text/data/Spanish/Pola has test network on folder in JSON format for both balanced and unbalanced assumption.
- Normal_pm_2021.csv has the gaussinan PV forecast used in the numerical illustration
- beta_lm_2016_8_6_60min.csv has the Beta distribution of diferent load type used in the numerical illustration.
- the parsers are available at \src\util section

The folder also contains load distribution functions, irradiance forecast PDF and a parser (`build_stochastic_data_mc(data,deg)`).
    
## Network Data with Stochastic Data Extension
The Folder /text/data/Spanish/All_feeder has Spanish network on folder in JSON format. 

The file `test/data/Spanish/CreatePMDDictionary.jl` is the parser file to convert the JSON file into `PowerModels.jl` or `PowerModelsDistribution.jl` format.

The original dataset consists of a full Low voltage network of sub-urban region with 30 transformers, 160 feeders, 10290 nodes and 8087 consumers, with load profiles of 20 days from actual smart-meter.
The paper is available in https://www.sciencedirect.com/science/article/pii/S0142061519318836
and the link for full data-set: https://data.mendeley.com/datasets/685vgp64sm/1

For the original networks, the line impedance is specified 4x4 matrice without mutual impedance and the load from smart meter data for 20 days.

However, in this part the load and irradiance are defined as Beta distribution in folders `beta_lm_2016_8_6.csv` and `beta_pm_2016_8_6.csv` respectively for a high irradiance day in spring. 

For each feeder there are 4 JSON file describing the feeder topology:
	
- *_configuration.json, 
- *_branches.json, 
- *_buses.json, and 
- *_devices.json 
and 1 csv file:
- *.csv- which is the linking file between devices and the load uncertainty. 

    
## Installation
For this application the installation of StochasticPowerModels.jl is not suggested as this application might not be updated along with the changes in the main package and will remain as an standalone application.

## Acknowledgements

The primary developer of StochasticPowerModels.jl is Tom Van Acker ([@timmyfaraday](https://github.com/timmyfaraday)), 
with support from the following contributors:
- Arpan Koirala ([@arpkoirala](https://github.com/arpkoirala)), KU Leuven, ACR formulation
- Frederik Geth ([@frederikgeth](https://github.com/frederikgeth)), CSIRO, reduced IVR formulation

The threephase hosting capacity extension and dynamic operating envelope was developed by Arpan Koirala.

The latest stable release of StochasticPowerModels can be obtained at:

```
https://github.com/Electa-Git/StochasticPowerModels.jl
```
## License

This code is provided under a BSD license.

```

## Acknowledgements

The primary developer is Tom Van Acker ([@timmyfaraday](https://github.com/timmyfaraday)), 
with support from the following contributors:
- Arpan Koirala ([@arpkoirala](https://github.com/arpkoirala)), KU Leuven, ACR formulation
- Frederik Geth ([@frederikgeth](https://github.com/frederikgeth)), CSIRO, reduced IVR formulation

## License

This code is provided under a BSD license.
