from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    # ----------TGFbeta1 activation
    'TGFb_inact',
    'TGFb_act',
    'FMOD',
    'FMOD_complex',    
    'TGFBR_inact',
    'TGFBR_act',
    # ----------TGFbeta pathway
    #S; SMAD
    'S2',    
    'S3',
    'S4',
    'pS2',
    'pS3',
    'ppS2',
    'ppS3',
    'ppS2_ppS3',
    'ppS2_ppS3_S4',
    'mS7',
    'S7',
    'mcFOS',
    'cFOS',
    'ppS2_ppS3_S4_cFOS',
    'mTHBS1',
    'THBS1',
]

NUM: int = len(NAMES)

Species = make_dataclass(
    cls_name="Species",
    fields=[(name, int) for name in NAMES],
    namespace={"NAMES": NAMES, "NUM": NUM},
    frozen=True,
)

name2idx: Dict[str, int] = {k: v for v, k in enumerate(NAMES)}

V = Species(**name2idx)

del name2idx