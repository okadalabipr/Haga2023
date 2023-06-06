from .name2idx import C, V


class DifferentialEquation(object):
    def __init__(self, perturbation):
        super(DifferentialEquation, self).__init__()
        self.perturbation = perturbation

    # Refined Model
    def diffeq(self, t, y, *x):
        v = {}
    # TGFbeta1 activation(Original reaction)
        v[1] = x[C.kf_1_TGFbeta]*y[V.THBS1]*y[V.TGFb_inact]/(x[C.Kmf_1_TGFbeta] + y[V.TGFb_inact])
        v[2] = x[C.k_on_FMOD]*y[V.FMOD]*y[V.TGFb_act]
        v[3] = x[C.Rec_act]*y[V.TGFBR_inact]*y[V.TGFb_act] #TGFR activation
        v[4] = x[C.pRec_debind]*y[V.TGFBR_act]  #TGFR deactivation
    # THBS1 regulation
        v[5] = x[C.kf_2_TGFbeta]*y[V.S2]*y[V.TGFBR_act]/(x[C.Kmf_2_TGFbeta]*(1+y[V.S7]/x[C.k_inhibit_TGF])+y[V.S2]) #SMAD2 activation by TGFR
        v[6] = x[C.S_dephosphos]*y[V.ppS2]
        v[7] = x[C.S_dephos]*y[V.pS2]
        v[8] = x[C.kf_3_TGFbeta]*y[V.S3]*y[V.TGFBR_act]/(x[C.Kmf_3_TGFbeta]*(1+y[V.S7]/x[C.k_inhibit_TGF])+y[V.S3]) #SMAD3 activation by TGFR
        v[9] = x[C.S_dephosphos]*y[V.ppS3]
        v[10] = x[C.S_dephos]*y[V.pS3]
        v[11] = x[C.k_on_ppS2_ppS3]*y[V.ppS3]*y[V.ppS2] #ppS2-ppS3 formation
        v[12] = x[C.S_dephosphos]*y[V.ppS2_ppS3]
        v[13] = x[C.k_on_ppS2_ppS3_S4]*y[V.ppS2_ppS3]*y[V.S4]
        v[14] = x[C.prod_mS7]*y[V.ppS2_ppS3_S4]**x[C.n1_TGF]/(x[C.Km_1_TGF]**x[C.n1_TGF]+y[V.ppS2_ppS3_S4]**x[C.n1_TGF]) # SMAD7
        v[15] = x[C.mS7_turn]*y[V.mS7]
        v[16] = x[C.prod_S7]*y[V.mS7] #Translation of SMAD7
        v[17] = x[C.prod_mcFOS]*y[V.TGFBR_act]**x[C.n2_TGF]/(x[C.Km_2_TGF]**x[C.n2_TGF]+y[V.TGFBR_act]**x[C.n2_TGF]) # cFOS
        v[18] = x[C.mcFOS_turn]*y[V.mcFOS]
        v[19] = x[C.prod_cFOS]*y[V.mcFOS] #Translation of cFOS
        v[20] = x[C.k_off_FMOD]*y[V.FMOD_complex]
        v[21] = x[C.k_off_ppS2_ppS3_S4]*y[V.ppS2_ppS3_S4]
        v[22] = x[C.k_off_ppS2_ppS3_S4_cFOS]*y[V.ppS2_ppS3_S4_cFOS]
        v[23] = x[C.k_on_ppS2_ppS3_S4_cFOS]*y[V.ppS2_ppS3_S4]*y[V.cFOS] #Epigenetic complex formation
        v[24] = x[C.prod_mTHBS1]*y[V.ppS2_ppS3_S4_cFOS]**x[C.n4_TGF]/(x[C.Km_4_TGF]**x[C.n4_TGF]+y[V.ppS2_ppS3_S4_cFOS]**x[C.n4_TGF]) # THBS1
        v[25] = x[C.THBS1_turn]*y[V.mTHBS1]
        v[26] = x[C.prod_THBS1]*y[V.mTHBS1] #Translation of THBS1
        v[27] = x[C.degrad_cFOS]*y[V.cFOS] #Degradation of cFos proteins

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        dydt = [0] * V.NUM
        # ----------TGFbeta activation (Binding stage)
        dydt[V.TGFb_inact] = - v[1]
        dydt[V.TGFb_act] = + v[1] - v[2] - v[3] + v[20]
        dydt[V.FMOD] = - v[2] + v[20]
        dydt[V.FMOD_complex] = + v[2] - v[20]
        dydt[V.TGFBR_inact] = - v[3] + v[4]
        dydt[V.TGFBR_act] = + v[3] - v[4]
        # ----------TGFbeta pathway
        dydt[V.S2] = - v[5] + v[7]
        dydt[V.S3] = - v[8] + v[10]
        dydt[V.S4] = - v[13] + v[21]
        dydt[V.pS2] = + v[6] - v[7] + v[12]
        dydt[V.pS3] = + v[9] - v[10] + v[12]
        dydt[V.ppS2] = + v[5] - v[6] - v[11]
        dydt[V.ppS3] = + v[8] - v[9] - v[11]
        dydt[V.ppS2_ppS3] = + v[11] - v[12] - v[13] + v[21]
        dydt[V.ppS2_ppS3_S4] = + v[13] - v[23] - v[21] + v[22]
        dydt[V.mS7] = + v[14] - v[15]
        dydt[V.S7] = + v[16] 
        dydt[V.mcFOS] = + v[17] - v[18]
        dydt[V.cFOS] = + v[19] - v[27] - v[23] + v[22] 
        dydt[V.ppS2_ppS3_S4_cFOS] = + v[23] - v[22]
        dydt[V.mTHBS1] = + v[24] - v[25]
        dydt[V.THBS1] = + v[26]
        
        return dydt

def param_values():

    x = [0] * C.NUM
	
	#If the value is 0, it will show error
    # ----------TGFbeta1 activation
    x[C.kf_1_TGFbeta] = 1
    x[C.Kmf_1_TGFbeta] = 1
    x[C.k_on_FMOD] =1
    x[C.k_off_FMOD] = 1
    x[C.Rec_act] = 1
    x[C.pRec_debind] = 1
    # ----------TGFbeta pathway
    x[C.S2tot] = 0.060
    x[C.S3tot] = 0.38
    x[C.S4tot] = 0.0044
    x[C.kf_2_TGFbeta] = 1
    x[C.Kmf_2_TGFbeta] = 1
    x[C.k_inhibit_TGF] = 1
    x[C.S_dephosphos] = 1
    x[C.S_dephos] = 1
    x[C.kf_3_TGFbeta] = 1
    x[C.Kmf_3_TGFbeta] = 1
    x[C.k_on_ppS2_ppS3] = 1
    x[C.k_on_ppS2_ppS3_S4] = 1
    x[C.k_off_ppS2_ppS3_S4] = 1
    # ----------SMAD7
    x[C.prod_mS7] = 1
    x[C.n1_TGF] = 1
    x[C.Km_1_TGF] = 1
    x[C.mS7_turn] = 1
    x[C.prod_S7] = 1
    # ----------cFOS
    x[C.prod_mcFOS] = 1
    x[C.n2_TGF] = 1
    x[C.Km_2_TGF] = 1
    x[C.mcFOS_turn] = 1
    x[C.prod_cFOS] = 1
    x[C.k_on_ppS2_ppS3_S4_cFOS] = 1
    x[C.k_off_ppS2_ppS3_S4_cFOS] = 1
    # ----------THBS1
    x[C.prod_mTHBS1] = 1
    x[C.n4_TGF] = 1
    x[C.Km_4_TGF] = 1 
    x[C.THBS1_turn] = 1
    x[C.prod_THBS1] = 1
	# ----------Weight genes
    x[C.w_THBS1] = 1
    #TGFbeta pathway
    x[C.w_TGFBR1] = 1
    x[C.w_TGFBR2] = 1
    x[C.w_SMAD7] = 1
    x[C.w_cFOS] = 1
	# ----------degrad
    x[C.degrad_cFOS] = 1
    return x


def initial_values():

    y0 = [0] * V.NUM
    # ----------TGFbeta1 activation
    y0[V.TGFb_inact] = 1
    y0[V.TGFb_act] =0 #Set by condition
    y0[V.FMOD] = 0
    y0[V.FMOD_complex] = 0
    y0[V.TGFBR_inact] = 1 
    y0[V.TGFBR_act] = 0
    # ----------TGFbeta pathway
    y0[V.S2] = 0.060
    y0[V.S3] = 0.38
    y0[V.S4] = 0.0044
    y0[V.pS2] = 0
    y0[V.pS3] = 0
    y0[V.ppS2] = 0
    y0[V.ppS3] = 0
    y0[V.ppS2_ppS3] = 0
    y0[V.ppS2_ppS3_S4] = 0
    y0[V.mS7] = 1 #estimate by weight*TPM
    y0[V.S7] = 0
    y0[V.mcFOS] = 1 #estimate by weight*TPM
    y0[V.cFOS] = 0
    y0[V.ppS2_ppS3_S4_cFOS] = 0
    y0[V.mTHBS1] = 1 
    y0[V.THBS1] = 0

    return y0