from biomass import Model
from biomass.models import TGFb_model

#Define model
model = Model(TGFb_model.__package__).create(show_info=True)

#Run simulation
from biomass import run_simulation

run_simulation(model, viz_type="average", stdev=True)

#Estimated parameters
from biomass.result import OptimizationResults

res = OptimizationResults(model)
# Export estimated parameters in CSV format
res.to_csv()

# Visualize estimated parameter sets
res.savefig(figsize=(32,10), boxplot_kws={"orient": "v"})

# Visualize objective function traces for different optimization runs
res.trace_obj()

#Run sensitivity analysis
from biomass import run_analysis

run_analysis(model, target='reaction', metric='integral', style='heatmap')


