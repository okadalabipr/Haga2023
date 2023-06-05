#Biomass Version
import biomass
print('biomass version:', biomass.__version__)

from biomass import Model
from biomass.models import TGFb_VEGF_model

#Define model
model = Model(TGFb_VEGF_model.__package__).create(show_info=True)

#Run optimization
from biomass import optimize

optimize(
    model, x_id=range(1, 31), options={
        "popsize": 5,
        "max_generation": 100,
        "allowable_error": 0.5,
        "local_search_method": "DE",
        "maxiter": 50,
        "workers": 1,
        
    }
)