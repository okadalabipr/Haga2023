import numpy as np

from biomass.estimation import convert_scale, initialize_search_param

from .name2idx import C, V
from .set_model import initial_values, param_values


class SearchParam(object):
    """ Specify model parameters and/or initial values to optimize"""

    def __init__(self):
        # parameters
        self.idx_params = [
        # ----------
        C.kf_1_TGFbeta,
        C.Kmf_1_TGFbeta,
        C.Rec_act,
        C.pRec_debind,
        # ----------TGFbeta pathway
        C.kf_2_TGFbeta,
        C.Kmf_2_TGFbeta,
        C.k_inhibit_TGF,
        C.S_dephosphos,
        C.S_dephos,
        C.kf_3_TGFbeta,
        C.Kmf_3_TGFbeta,
        C.k_on_ppS2_ppS3,
        C.k_on_ppS2_ppS3_S4,
        C.k_off_ppS2_ppS3_S4,
        # ----------SMAD7
        C.prod_mS7,
        C.Km_1_TGF,
        C.mS7_turn,
        C.prod_S7,
        # ----------cFOS
        C.prod_mcFOS,
        C.Km_2_TGF,
        C.mcFOS_turn,
        C.prod_cFOS,
        C.k_on_ppS2_ppS3_S4_cFOS,
        C.k_off_ppS2_ppS3_S4_cFOS,
        # ----------THBS1
        C.prod_mTHBS1,
        C.Km_4_TGF,
        C.THBS1_turn,
        C.prod_THBS1,
        # weight for initial value
        #C.w_TGF pathway
        C.w_TGFBR1,
        C.w_TGFBR2,
        C.w_THBS1,
        C.w_SMAD7,
        C.w_cFOS,
        C.degrad_cFOS,
    ]

        # initial values
        self.idx_initials = [
        ]

    def get_region(self):
        x = param_values()
        y0 = initial_values()

        search_param = initialize_search_param(
            parameters=C.NAMES,
            species=V.NAMES,
            param_values=x,
            initial_values=y0,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        search_rgn = np.zeros((2, len(x) + len(y0)))
        # Default: 0.1 ~ 10
        for i, j in enumerate(self.idx_params):
            search_rgn[0, j] = search_param[i] * 0.1  # lower bound
            search_rgn[1, j] = search_param[i] * 10  # upper bound
        # Default: 0.5 ~ 2
        for i, j in enumerate(self.idx_initials):
            search_rgn[0, j + len(x)] = search_param[i + len(self.idx_params)] * 0.01  # lower bound
            search_rgn[1, j + len(x)] = search_param[i + len(self.idx_params)] * 100  # upper bound

        # search_rgn[:,C.parameter] = [lower_bound,upper_bound]
        # search_rgn[:,V.specie+len(x)] = [lower_bound,upper_bound]
        search_rgn[:, C.prod_mcFOS] = [1, 1000.0]
        search_rgn[:, C.Km_2_TGF] = [0.0001, 10.0]
        search_rgn[:, C.pRec_debind] = [0.0001, 10.0]
        search_rgn[:, C.prod_mTHBS1] = [1, 1000.0]
        # ----------Turn
        search_rgn[:, C.THBS1_turn] = [0.0001, 10.0]
        search_rgn[:, C.mS7_turn] = [0.0001, 10.0]
        search_rgn[:, C.mcFOS_turn] = [0.0001, 10.0]        
        search_rgn = convert_scale(
            region=search_rgn,
            parameters=C.NAMES,
            species=V.NAMES,
            estimated_params=self.idx_params,
            estimated_initials=self.idx_initials,
        )

        return search_rgn

    def update(self, indiv):
        x = param_values()
        y0 = initial_values()

        for i, j in enumerate(self.idx_params):
            x[j] = indiv[i]
        for i, j in enumerate(self.idx_initials):
            y0[j] = indiv[i + len(self.idx_params)]
        
        # --------------------------------------------------------------------------
        # TPM value used are caculated using below
        # y0 = x[C.w_gene]*TPM value(in vitro(PDL24), 0 min)
        # --------------------------------------------------------------------------
        y0[V.TGFBR_inact] = x[C.w_TGFBR1]*64.76 + x[C.w_TGFBR2]*95.81
        y0[V.mS7] = x[C.w_SMAD7]*20.62
        y0[V.mcFOS] = x[C.w_cFOS]*1.17
        y0[V.mTHBS1] = x[C.w_THBS1]*4522.39
        return x, y0

    def gene2val(self, indiv_gene):
        search_rgn = self.get_region()
        indiv = 10 ** (indiv_gene * (search_rgn[1, :] - search_rgn[0, :]) + search_rgn[0, :])

        return indiv

    def val2gene(self, indiv):
        search_rgn = self.get_region()
        indiv_gene = (np.log10(indiv) - search_rgn[0, :]) / (search_rgn[1, :] - search_rgn[0, :])

        return indiv_gene