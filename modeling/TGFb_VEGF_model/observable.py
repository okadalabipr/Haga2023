import numpy as np

from biomass.dynamics.solver import get_steady_state, solve_ode

from .name2idx import C, V
from .set_model import DifferentialEquation


class Observable(DifferentialEquation):
    """
    Correlating model simulations and experimental measurements.

    Attributes
    ----------
    obs_names : list of strings
        Names of self.obs_names.

    t : range
        Simulation time span.

    conditions : list of strings
        Expetimental conditions.

    simulations : numpy.ndarray
        The numpy array to store simulation results.

    normalization : nested dict
        * 'timepoint' : Optional[int]
            The time point at which simulated values are normalized.
            If None, the maximum value will be used for normalization.

        * 'condition' : list of strings
            The experimental conditions to use for normalization.
            If empty, all conditions defined in sim.conditions will be used.

    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """
    def __init__(self):
        super(Observable, self).__init__(perturbation={})
        self.obs_names = [
            'THBS1',
            'Phosphorylated SMAD2',
            'Phosphorylated SMAD3',
            'cFOS',
            'FMOD',
            'Phosphorylated Akt',
        ]
        self.t: range = range(3000) #time span
        self.conditions: list = ["Control", "TGFβ1"]
        self.simulations: np.ndarray = np.empty(
            (len(self.obs_names), len(self.t), len(self.conditions))
        )
        self.normalization: dict = {}
        for observable in self.obs_names:
            self.normalization[observable] = {"timepoint": None, "condition": []}
        self.experiments: list = [None] * len(self.obs_names)
        self.error_bars: list = [None] * len(self.obs_names)

    def simulate(self, x, y0, _perturbation=None):
        if _perturbation is not None:
            self.perturbation = _perturbation
        # add ligand
        for i, condition in enumerate(self.conditions):
            if condition == 'Control':
                y0[V.TGFb_inact] = 0.013
                y0[V.TGFb_act] =0
            elif condition == 'TGFβ1':
                y0[V.TGFb_inact] = 0.013
                y0[V.TGFb_act] =0.0902 #(nM)
            #Define initial states
            y0[V.S2] = x[C.S2tot]
            y0[V.S3] = x[C.S3tot]
            y0[V.S4] = x[C.S4tot]

            sol = solve_ode(self.diffeq, y0, self.t, tuple(x))
            
            if sol is None:
                return False
            else:
                self.simulations[self.obs_names.index("THBS1"), :, i] = sol.y[
                    V.THBS1, :
                ]
                self.simulations[self.obs_names.index("Phosphorylated SMAD2"), :, i] = sol.y[
                    V.ppS2, :
                ]
                self.simulations[self.obs_names.index("Phosphorylated SMAD3"), :, i] = sol.y[
                    V.ppS3, :
                ]
                self.simulations[self.obs_names.index("cFOS"), :, i] = sol.y[
                    V.cFOS, :
                ]
                self.simulations[self.obs_names.index("FMOD"), :, i] = sol.y[
                    V.FMOD, :
                ]
                self.simulations[self.obs_names.index("Phosphorylated Akt"), :, i] = sol.y[
                    V.Aktstar, :
                ]

    def set_data(self):
        self.experiments[self.obs_names.index("THBS1")] = {
            "Control": [0, 0.045, 0.022, 0.059, 0.088, 0.527, 0.257, 0.179],
            "TGFβ1": [0.0239, 0.0243, 0.0045, 0.0627, 0.0911, 0.5217, 0.9951, 1.0000],
        }
        self.error_bars[self.obs_names.index("THBS1")] = {
            "Control": [
                sd / np.sqrt(3) for sd in [0.036, 0.024, 0.031, 0.021, 0.026, 0.197, 0.111, 0.079]
            ],
            "TGFβ1": [
                sd / np.sqrt(3) for sd in [0.0149, 0.0207, 0.0162, 0.0727, 0.0611, 0.1757,0.4279, 0.3416]
            ],
        }
        
        self.experiments[self.obs_names.index("Phosphorylated SMAD2")] = {
            "Control": [0, 0.0517, 0.0837, 0.1092, 0.0961, 0.1125, 0.0902, 0.1092],
            "TGFβ1": [0.0956, 0.6675, 0.9530, 1.0000, 0.5097, 0.2398, 0.1842, 0.2203],
        }
        self.error_bars[self.obs_names.index("Phosphorylated SMAD2")] = {
            "Control": [
                sd / np.sqrt(3) for sd in [0.0224, 0.0444, 0.0625, 0.0619, 0.0763, 0.0888, 0.0543, 0.0833]
            ],
            "TGFβ1": [
                sd / np.sqrt(3) for sd in [0.0498, 0.0483, 0.1089, 0.1552, 0.0895, 0.1364, 0.1064, 0.1538]
            ],
        }
        
        self.experiments[self.obs_names.index("Phosphorylated SMAD3")] = {
            "Control": [0, 0.1296, 0.1353, 0.1528, 0.1115, 0.0722, 0.0229, 0],
            "TGFβ1": [0, 0.9295, 0.9221, 1.0000, 0.7695, 0.4764, 0.2184, 0.1707],
        }
        self.error_bars[self.obs_names.index("Phosphorylated SMAD3")] = {
            "Control": [
                sd / np.sqrt(3) for sd in [0, 0.0200, 0.0254, 0.0226, 0.0071, 0.0044, 0.0397, 0]
            ],
            "TGFβ1": [
                sd / np.sqrt(3) for sd in [0, 0.0607, 0.0601, 0.0417, 0.1941, 0.1036, 0.0373, 0.0452]
            ],
        }
        self.experiments[self.obs_names.index("cFOS")] = {
            "Control": [0, 0.0018, 0.0980, 0.5733, 0.3194, 0.0088, 0, 0],
            "TGFβ1": [0, 0, 0.1229, 1.000, 0.5687, 0.0315, 0, 0],
        }
        self.error_bars[self.obs_names.index("cFOS")] = {
            "Control": [
                sd / np.sqrt(3) for sd in [0, 0.0031, 0.0141, 0.1146, 0.1596, 0.0056, 0, 0]
            ],
            "TGFβ1": [
                sd / np.sqrt(3) for sd in [0 , 0, 0.0258, 0.1256, 0.0757, 0.0122, 0, 0]
            ],
        }
        self.experiments[self.obs_names.index("FMOD")] = {
            "Control": [0, 0, 0, 0, 0.0098, 0.1214, 0.3547, 1.0000],
            "TGFβ1": [0, 0, 0, 0.0021, 0.0186, 0.1381, 0.2935, 0.5292],
        }
        self.error_bars[self.obs_names.index("FMOD")] = {
            "Control": [
                sd / np.sqrt(3) for sd in [0, 0, 0, 0, 0.0169, 0.0249, 0.0580, 0.0429]
            ],
            "TGFβ1": [
                sd / np.sqrt(3) for sd in [0, 0, 0, 0.0037, 0.0169, 0.0118, 0.0181,0.06378]
            ],
        }
        self.experiments[self.obs_names.index("Phosphorylated Akt")] = {
            "Control": [0, 1.0000, 0.6968, 0.8812, 0.7410, 0.5109, 0.2646, 0.0775],
            "TGFβ1": [0.038, 0.8809, 0.7806, 0.9843, 0.9740, 0.8676, 0.8998, 0.4728],
        }
        self.error_bars[self.obs_names.index("Phosphorylated Akt")] = {
            "Control": [
                sd / np.sqrt(3) for sd in [0.0256, 0.3775, 0.4171, 0.3822, 0.1888, 0.1886, 0.1603, 0.1214]
            ],
            "TGFβ1": [
                sd / np.sqrt(3) for sd in [0.0185, 0.1250, 0.0727, 0.1070, 0.1541, 0.0948, 0.1107, 0.1071]
            ],
        }

    @staticmethod
    def get_timepoint(obs_name):
        # 0, 15min, 30 min, 60 min, 120 min,8 h, 24 h, 48 h(Unit: min)
        return [0, 15, 30, 60, 120, 480, 1440, 2880]  
