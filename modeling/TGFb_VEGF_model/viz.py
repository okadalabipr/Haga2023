from matplotlib import pyplot as plt

from biomass.plotting import *

from .observable import Observable


class Visualization(Observable):
    """
    Plotting parameters for customizing figure properties.

    Attributes
    ----------
    cm : matplotlib.colors.ListedColormap (default: `plt.cm.get_cmap('tab10')`)
        Choosing colormaps for `cmap`.
    single_observable_options : list of SingleObservable
        Visualization options for time-course simulation (single-observable).
    multiple_observables_options : MultipleObservables
        Visualization options for time-course simulation (multi-observables).
    sensitivity_options : SensitivityOptions
        Visualization options for sensitivity analysis results.
    """

    def __init__(self):
        super().__init__()

        self.cm = plt.cm.get_cmap("tab10")
        self.single_observable_options = [
            SingleObservable(self.cm, obs_name) for obs_name in self.obs_names
        ]
        self.multiple_observables_options = MultipleObservables(self.cm)
        self.sensitivity_options = SensitivityOptions(self.cm)

    def get_single_observable_options(self):

        for i, _ in enumerate(self.obs_names):
            self.single_observable_options[i].divided_by = 60  # min. -> h.
            self.single_observable_options[i].xlim = (-5, 55)
            self.single_observable_options[i].xticks = [10 * i for i in range(6)]
            self.single_observable_options[i].xlabel = 'Time (h)'
            self.single_observable_options[i].ylim = (-0.1, 1.3)
            self.single_observable_options[i].yticks = [0.0, 0.3, 0.6, 0.9, 1.2]
            self.single_observable_options[i].cmap = ['#4F81BD', '#C00000']
            self.single_observable_options[i].shape = ['D', 's']
            self.single_observable_options[i].dont_show = []
            self.single_observable_options[i].legend_kws = dict(loc="best", frameon=False)
        return self.single_observable_options

    def get_multiple_observables_options(self):

        return self.multiple_observables_options
    
    def get_sensitivity_options(self):
        for i, _ in enumerate(self.obs_names):
            self.single_observable_options[i].cmap = ['#4F81BD', '#C00000']
        return self.sensitivity_options

    @staticmethod
    def set_timecourse_rcParams():
        """figure/simulation"""
        plt.rcParams["figure.figsize"] = (8, 6)
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 1.5
        plt.rcParams["xtick.major.width"] = 1.5
        plt.rcParams["ytick.major.width"] = 1.5
        plt.rcParams["lines.linewidth"] = 1.8
        plt.rcParams["lines.markersize"] = 12
        plt.rcParams["savefig.bbox"] = "tight"
        # plt.rcParams["savefig.format"] = "pdf"
        # plt.rcParams['font.family'] = 'Arial'
        # plt.rcParams['mathtext.fontset'] = 'custom'
        # plt.rcParams['mathtext.it'] = 'Arial:italic'

    @staticmethod
    def set_sensitivity_rcParams():
        """figure/sensitivity"""
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 1.2
        plt.rcParams["xtick.major.width"] = 1.2
        plt.rcParams["ytick.major.width"] = 1.2
        plt.rcParams["savefig.bbox"] = "tight"
        plt.rcParams["savefig.format"] = "png"
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def convert_species_name(name):
        """figure/sensitivity/initial_condition
        - Sensitivity for species with nonzero initial conditions
        """
        return name