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
    # Birtwistle et al., MSB (2007)
    'E',#VEGF Ligand
    'E1',#VEGFR free receptor
    'E_E1', #VEGF_VEGFR
    'E11',#VEGF_VEGFR_dimer
    'E11P',#VEGF_pVEGFR_dimer
    'G',#Free Grb2
    'S',#Free Shc
    'I',#Free PI3K
    'R',#Free RasGAP
    'T',#Free PTP-1B
    'O',#Free SOS
    'A',#Free GAB1
    'E11G',#E11 dimer bound to Grb2
    'E11S',#E11 dimer bound to Shc
    'E11R',#E11 dimer bound to RasGAP
    'sigmaG',#Membrane localized Grb2
    'sigmaS',#Membrane localized Shc
    'sigmaI',#Membrane localized PI3K
    'sigmaR',#Membrane localized RasGAP
    'sigmaA',#Membrane localized Gab1
    'sigmaSP',#Phosphorylated membrane localized Shc
    'sigmaAP',#Phosphorylated membrane localized Gab1
    'sigmaG_O',#SOS bound to membrane bound Grb2
    'sigmaG_A',#Gab1 bound to membrane bound Grb2
    'sigmaSP_G',#Grb2 bound to membrane localized, phosphorylated Shc
    'sigmaAP_S',#Shc bound to membrane localized, phosphorylated Gab1
    'sigmaAP_I',#PI-3K bound to membrane localized, phosphorylated Gab1
    'sigmaAP_R',#RasGAP bound to membrane localized, phosphorylated Gab1
    'P3_A',#Gab1 bound to PIP3
    'P2',#Free PIP2
    'P3',#Free PIP3
    'Akt',#Free, inactive Akt
    'Aktstar',#Free, active Akt
    'RsD',#Ras-GDP
    'RsT',#Ras-GTP
    'sigmaRP',#Membrane localized, phosphorylated RasGAP
    'Raf',#Free, inactive Raf
    'Rafstar',#Free, active Raf
    'MEK',#Free, inactive MEK
    'ppMEKc',#Free, active MEK
    'OP',#Inactive, phosphorylated SOS
    'AP',#Inactive, phosphorylated Gab1
    'A_sigmaG_O',#Membrane bound Grb2 bound to Gab1 and SOS
    'sigmaA_G',#Gab1 bound to PIP3 and Grb2
    'sigmaA_G_O',#Gab1 bound to PIP3 and Grb2-SOS
    'sigmaO',#Membrane localized SOS
    'E11T',#PTP-1B bound to VEGFR dimers
    'sigmaT',#Membrane localized PTP-1B
    'sigmaAP_T',#PTP-1B bound to membrane-localized Gab1
    'E1_PT',#Threonine Phosphorylated VEGFR
    'E_E1_PT',#Threonine Phosphorylated bond to VEGFR
    'fint',
    # Nakakuki et al., Cell (2010)
    'ERKc',
    'ERKn',
    'pERKc',
    'pERKn',
    'ppERKc',
    'ppERKn',
    # ----------VEGF pathway
    'TF_inact',
    'TF_act',
    'mFMOD',
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