import numpy as np

from biomass.estimation import convert_scale, initialize_search_param

from .name2idx import C, V
from .set_model import initial_values, param_values


class SearchParam(object):
    """ Specify model parameters and/or initial values to optimize"""

    def __init__(self):
        # parameters
        self.idx_params = [
        # ----------FMOD complex formation
        C.k_on_FMOD,
        C.k_off_FMOD,
        # ----------Imoto et al
        C.VmaxPY,
        C.KmPY,
        C.kdeg,
        C.kf47,
        C.Vmaxr47,
        C.Kmf47,
        C.Kmr47,
        C.kf48,
        C.Kmf48,
        C.Kmr48,
        C.kf49,
        C.kr49,
        C.Kmf49,
        C.Kmr49,
        C.Kmr49b,
        C.kr49b,
        C.kf51,
        C.Vmaxr51,
        C.Kmf51,
        C.Kmrb51,
        C.kf52,
        C.Vmaxr52,
        C.Kmf52,
        C.Kmr52,
        C.kf54,
        C.Vmaxr54,
        C.Kmf54,
        C.Kmr54,
        C.kf55,
        C.Vmaxr55,
        C.Kmf55,
        C.Kmr55,
        C.kf38,
        C.kf39,
        C.kf50,
        C.a98,
        C.b98,
        C.koff46,
        C.EGF_off,
        C.koff4,
        C.koff16,
        C.koff17,
        C.koff18,
        C.koff40,
        C.koff41,
        C.koff42,
        C.koff43,
        C.koff44,
        C.koff45,
        C.koff57,
        C.koff58,
        C.koff59,
        C.koff60,
        C.kPTP10,
        C.koff73,
        C.kPTP38,
        C.kPTP39,
        C.koff88,
        C.kPTP50,
        C.kf81,
        C.Vmaxr81,
        C.Kmf81,
        C.Kmr81,
        C.kf84,
        C.Vmaxr84,
        C.Kmf84,
        C.Kmr84,
        C.kcon49,
        C.kon1,
        C.kon86,
        C.kon4,
        C.kf10,
        C.kon16,
        C.kon17,
        C.kon18,
        C.kon73,
        C.kon40,
        C.kon41,
        C.kon42,
        C.kon43,
        C.kon44,
        C.kon45,
        C.kon88,
        C.kon46,
        C.kon57,
        C.kon58,
        C.kon59,
        C.kon60,
        # Nakakuki et al., Cell (2010)
        C.V1,
        C.Km1,
        C.V5,
        C.Km5,
        # ----------PI3K-Akt pathway
        C.act_PI3K,
        C.inact_PI3K,
        # ----------VEGF pathway
        C.kf_16_vegf,
        C.Kmf_16_vegf,
        C.kr_17_vegf,
        C.Kmr_17_vegf,
        C.prod_mFMOD,
        C.Km_18_vegf,
        C.degrad_mFMOD,
        C.prod_FMOD,
        # weight for initial value
        #w_E1
        C.w_VEGFR1,
        C.w_VEGFR2,
        #GRB2
        C.w_G,
        #C.w_S
        C.w_SHC1,
        C.w_SHC2,
        C.w_SHC3,
        C.w_SHC4,
        #C.w_I,PI3K
        C.w_PIK3CA,
        C.w_PIK3CB,
        C.w_PIK3CD,
        C.w_PIK3CG,
        C.w_PTEN,
        #C.w_R
        C.w_RASA1,
        C.w_RASA2,
        C.w_RASA3,
        #C.w_O
        C.w_SOS1,
        C.w_SOS2,
        #GAB1
        C.w_A,
        #C.w_Akt
        C.w_AKT1,
        C.w_AKT2,
        #C.w_RsD
        C.w_HRAS,
        C.w_KRAS,
        C.w_NRAS,
        #C.w_Raf
        C.w_ARAF,
        C.w_BRAF,
        C.w_RAF1,
        #C.w_MEK
        C.w_MAP2K1,
        C.w_MAP2K2,
        C.w_T, #PTP-1B(Negative feedback)
        #C.w_ERK
        C.w_MAPK1,
        C.w_MAPK3,
        C.w_FMOD,
    ]

        # initial values
        self.idx_initials = [
        V.P2,
        V.TF_inact,
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
            search_rgn[0, j + len(x)] = search_param[i + len(self.idx_params)] * 1  # lower bound
            search_rgn[1, j + len(x)] = search_param[i + len(self.idx_params)] * 1000  # upper bound

        # search_rgn[:,C.parameter] = [lower_bound,upper_bound]
        # search_rgn[:,V.specie+len(x)] = [lower_bound,upper_bound]
        #n; Hill equations for activate gene
        search_rgn[:, C.k_on_FMOD] = [0.001, 1000]
        search_rgn[:, C.k_off_FMOD] = [0.001, 1000]
        search_rgn[:, C.act_PI3K] = [0.001, 1000]
        search_rgn[:, C.inact_PI3K] = [0.001, 1000]
        search_rgn[:, C.kf_16_vegf] = [0.001, 1000]
        search_rgn[:, C.Kmf_16_vegf] = [0.001, 1000]
        search_rgn[:, C.kr_17_vegf] = [0.001, 1000]
        search_rgn[:, C.Kmr_17_vegf] = [0.001, 1000]
        search_rgn[:, C.prod_mFMOD] = [0.001, 1000]
        search_rgn[:, C.Km_18_vegf] = [0.001, 1000]
        search_rgn[:, C.degrad_mFMOD] = [0.001, 1000]
        search_rgn[:, C.prod_FMOD] = [0.001, 1000]
        #Gene weight
        search_rgn[:, C.w_FMOD] = [0.1, 100]
        search_rgn[:, C.w_VEGFR1] = [0.1, 100.0]
        search_rgn[:, C.w_VEGFR2] = [0.1, 100.0]
        search_rgn[:, C.w_G] = [0.1, 100.0]
        #search_rgn[:, C.w_S] = [0.1, 100.0]
        search_rgn[:, C.w_SHC1] = [0.1, 100.0]
        search_rgn[:, C.w_SHC2] = [0.1, 100.0]
        search_rgn[:, C.w_SHC3] = [0.1, 100.0]
        search_rgn[:, C.w_SHC4] = [0.1, 100.0]
        #search_rgn[:, C.w_I] = [0.1, 100.0]
        search_rgn[:, C.w_PIK3CA] = [0.1, 100.0]
        search_rgn[:, C.w_PIK3CB] = [0.1, 100.0]
        search_rgn[:, C.w_PIK3CD] = [0.1, 100.0]
        search_rgn[:, C.w_PIK3CG] = [0.1, 100.0]
        # ver.5
        search_rgn[:, C.w_PTEN] = [0.1, 100.0]
        #search_rgn[:, C.w_R] = [0.1, 100.0]
        search_rgn[:, C.w_RASA1] = [0.1, 100.0]
        search_rgn[:, C.w_RASA2] = [0.1, 100.0]
        search_rgn[:, C.w_RASA3] = [0.1, 100.0]
        #search_rgn[:, C.w_O] = [0.1, 100.0]
        search_rgn[:, C.w_SOS1] = [0.1, 100.0]
        search_rgn[:, C.w_SOS2] = [0.1, 100.0]
        search_rgn[:, C.w_A] = [0.1, 100.0]
        #search_rgn[:, C.w_Akt] = [0.1, 100.0]
        search_rgn[:, C.w_AKT1] = [0.1, 100.0]
        search_rgn[:, C.w_AKT2] = [0.1, 100.0]
        #search_rgn[:, C.w_RsD] = [0.1, 100.0]
        search_rgn[:, C.w_HRAS] = [0.1, 100.0]
        search_rgn[:, C.w_KRAS] = [0.1, 100.0]
        search_rgn[:, C.w_NRAS] = [0.1, 100.0]
        #search_rgn[:, C.w_Raf] = [0.1, 100.0]
        search_rgn[:, C.w_ARAF] = [0.1, 100.0]
        search_rgn[:, C.w_BRAF] = [0.1, 100.0]
        search_rgn[:, C.w_RAF1] = [0.1, 100.0]
        #search_rgn[:, C.w_MEK] = [0.1, 100.0]
        search_rgn[:, C.w_MAP2K1] = [0.1, 100.0]
        search_rgn[:, C.w_MAP2K2] = [0.1, 100.0]
        search_rgn[:, C.w_T] = [0.1, 100.0]
        #search_rgn[:, C.w_ERKc] = [0.1, 100.0]
        search_rgn[:, C.w_MAPK1] = [0.1, 100.0]
        search_rgn[:, C.w_MAPK3] = [0.1, 100.0]

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
        
        # ----------TGFbeta1 activation
        # TPM value used are caculated using below
        # y0 = x[C.w_gene]*TPM value(in vitro(PDL24), 0 min)
        y0[V.TGFBR_inact] = x[C.w_TGFBR1]*64.76 + x[C.w_TGFBR2]*95.81
        y0[V.mS7] = x[C.w_SMAD7]*20.62
        y0[V.mcFOS] = x[C.w_cFOS]*1.17
        y0[V.mTHBS1] = x[C.w_THBS1]*4522.39
        # --------------------------------------------------------------------------
        y0[V.mFMOD] = x[C.w_FMOD]*20.88
        y0[V.E1] = x[C.w_VEGFR1]*101.53 + x[C.w_VEGFR2]*0.02
        y0[V.G] = x[C.w_G]*128.45 #GRB2
        y0[V.S] = x[C.w_SHC1]*401.13 + x[C.w_SHC2]*1.16 + x[C.w_SHC3]*4.28 + x[C.w_SHC4]*1.26  
        y0[V.I] = x[C.w_PIK3CA]*16.36 + x[C.w_PIK3CB]*15.10 + x[C.w_PIK3CD]*32.27 + x[C.w_PIK3CG]*0
        x[C.PTEN] = x[C.w_PTEN]*51.20
        y0[V.R] = x[C.w_RASA1]*85.47 + x[C.w_RASA2]*26.93 + x[C.w_RASA3]*92.99
        y0[V.A] = x[C.w_A]*7.09 #GAB1
        y0[V.O] = x[C.w_SOS1]*25.34 + x[C.w_SOS2]*20.01
        y0[V.Akt] = x[C.w_AKT1]*84.62 + x[C.w_AKT2]*41.27
        y0[V.RsD] = x[C.w_HRAS]*101.55 + x[C.w_KRAS]*42.43 + x[C.w_NRAS]*113.71
        y0[V.Raf] = x[C.w_ARAF]*85.60 + x[C.w_BRAF]*4.73 + x[C.w_RAF1]*79.97
        y0[V.MEK] = x[C.w_MAP2K1]*73.15 + x[C.w_MAP2K2]*114.17 #MEK1/2
        y0[V.T] = x[C.w_T]*84.13 #PTPN1
        y0[V.ERKc] = x[C.w_MAPK1]*78.20 + x[C.w_MAPK3]*75.96 #ERK1/2
        # constraints --------------------------------------------------------------
        x[C.V6] = x[C.V5]
        x[C.Km6] = x[C.Km5]
        # --------------------------------------------------------------------------


        return x, y0

    def gene2val(self, indiv_gene):
        search_rgn = self.get_region()
        indiv = 10 ** (indiv_gene * (search_rgn[1, :] - search_rgn[0, :]) + search_rgn[0, :])

        return indiv

    def val2gene(self, indiv):
        search_rgn = self.get_region()
        indiv_gene = (np.log10(indiv) - search_rgn[0, :]) / (search_rgn[1, :] - search_rgn[0, :])

        return indiv_gene