from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    # ----------TGFbeta1 activation
    'kf_1_TGFbeta',
    'Kmf_1_TGFbeta',
    'k_on_FMOD',
    'k_off_FMOD',
    'Rec_act',	
    'pRec_debind',
    # ----------TGFbeta pathway
	'S2tot',
	'S3tot',
	'S4tot',
	'kf_2_TGFbeta',
	'Kmf_2_TGFbeta',
	'k_inhibit_TGF',
	'S_dephosphos',
	'S_dephos',
	'kf_3_TGFbeta',
	'Kmf_3_TGFbeta',
	'k_on_ppS2_ppS3',
	'k_on_ppS2_ppS3_S4',
	'k_off_ppS2_ppS3_S4',
	# ----------SMAD7
	'prod_mS7',
	'n1_TGF',
	'Km_1_TGF',
	'mS7_turn',
	'prod_S7',
	# ----------cFOS
	'prod_mcFOS',
	'n2_TGF',
	'Km_2_TGF',
	'mcFOS_turn',
	'prod_cFOS',
	'k_on_ppS2_ppS3_S4_cFOS',
	'k_off_ppS2_ppS3_S4_cFOS',
	# ----------THBS1
	'prod_mTHBS1',
	'n4_TGF',
	'Km_4_TGF',
	'THBS1_turn',
	'prod_THBS1',
	# ----------Weight genes
    'w_THBS1',
    'w_TGFBR1',
    'w_TGFBR2',
    'w_SMAD7',
    'w_cFOS',
	# ----------degrad
    'degrad_cFOS',     
]

NUM: int = len(NAMES)

Parameters = make_dataclass(
    cls_name="Parameters",
    fields=[(name, int) for name in NAMES],
    namespace={"NAMES": NAMES, "NUM": NUM},
    frozen=True,
)

name2idx: Dict[str, int] = {k: v for v, k in enumerate(NAMES)}

C = Parameters(**name2idx)

del name2idx
