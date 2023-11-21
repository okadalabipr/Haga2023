# Mathematical modeling of skin aging

We used [`biomass==0.5.2`](https://github.com/biomass-dev/biomass) for mathematical modeling.
The package can be installed via pip:

```
$ pip install biomass==0.5.2
```

This requires Python 3.7 or later.

## Description

We generated two mathematical models, [`TGFb_model/`](./TGFb_model/) and [`TGFb_VEGF_model/`](./TGFb_VEGF_model/), in this study. Each related files can be found under each model directory.

| Name                            | Content                                                                                                                     |
| --------------------------------| ----------------------------------------------------------------------------------------------------------------------------|
| [`MethodS.md`]                  | Differential equation, parameters and initial condition                                                                     |
| [`name2idx/`]                   | Names of model parameters and species                                                                                       |
| [`out/`]                        | Parameter values that are estimated from experimental data                                                                  |
| [`reaction_network.py`]         | Reaction indices grouped according to biological processes                                                                  |
| [`set_model.py`]                | Differential equation, parameters and initial condition to run biomass                                                      |
| [`observalbe.py`]               | Observables, simulations and experimental data                                                                              |
| [`viz.py`]                      | Plotting parameters for customizing figure properties                                                                       |
| [`set_search_param.py`]         | Lower and upper bounds of model parameters to be estimated                                                                  |
| [`fitness.py`]                  | An objective function to be minimized, i.e., the distance between model simulation and experimental data                    |
| [`optimize.py`]                 | Run to estimate parameter values from experimental data                                                                     |
| [`simulate.py`]                 | Run for simulation, visualization of estimated parameter sets and objective function trajectories, and sensitivity analysis |