Maximize
Z: R("R_BiomassEcoli")

Subject To
M10fthfc: -R("R_AICART") + R("R_AICART_Rev") - R("R_FTHFD") -
 R("R_GARFT") + R("R_GARFT_Rev") + R("R_MTHFC") - R("R_MTHFC_Rev")  = 0.

M12dgrECc: -0.02 R("R_DAGK_EC") + 0.02 R("R_PAPA_EC") + 0.02
 R("R_PEPT_EC") - 0.02 R("R_PEPT_EC_Rev")  = 0.

M12ppdSb: R("M12ppdSb_Transp") + R("R_EX_12ppd_S_e_") -
 R("R_EX_12ppd_S_e__Rev")  = 0.

M12ppdSc: R("R_12PPDt") - R("R_12PPDt_Rev") + R("R_LCARS") -
 R("R_LCARS_Rev")  = 0.

M12ppdSe: -R("R_12PPDt") + R("R_12PPDt_Rev") - R("R_EX_12ppd_S_e_")
 + R("R_EX_12ppd_S_e__Rev")  = 0.

M13dpgc: R("R_GAPD") - R("R_GAPD_Rev") + R("R_PGK") - R("R_PGK_Rev")
  = 0.

M15dapb: R("M15dapb_Transp") + R("R_EX_15dap_e_") -
 R("R_EX_15dap_e__Rev")  = 0.

M15dapc: -R("R_CADVt") + R("R_LYSDC")  = 0.

M15dape: R("R_CADVt") - R("R_EX_15dap_e_") + R("R_EX_15dap_e__Rev") 
 = 0.

M1pyr5cc: R("R_G5SADs") - R("R_P5CD") - R("R_P5CR") + R("R_PROD2") 
 = 0.

M23ddhbc: -R("R_DHBD") + R("R_DHBD_Rev") + R("R_ICHORT")  = 0.

M23dhbac: R("R_DHBSr") - R("R_DHBSr_Rev") - 3 R("R_ENTCS")  = 0.

M23dhbc: R("R_DHBD") - R("R_DHBD_Rev") - R("R_DHBSr") +
 R("R_DHBSr_Rev")  = 0.

M23dhdpc: -R("R_DHDPRy") + R("R_DHDPS")  = 0.

M23dhmbc: -R("R_DHAD1") + R("R_KARA1i")  = 0.

M23dhmpc: -R("R_DHAD2") + R("R_KARA2i")  = 0.

M23dogulnc: -R("R_DOGULNR")  = 0.

M25aicsc: -R("R_ADSL2r") + R("R_ADSL2r_Rev") + R("R_PRASCS") -
 R("R_PRASCS_Rev")  = 0.

M25dkglcnc: -R("R_DKGLCNR1") - R("R_DKGLCNR2x") - R("R_DKGLCNR2y") 
 = 0.

M25drappc: -R("R_DHPPDA2") + R("R_GTPCII2")  = 0.

M26dapLLc: -R("R_DAPE") + R("R_DAPE_Rev") + R("R_SDPDS")  = 0.

M26dapMb: R("M26dapMb_Transp") + R("R_EX_26dap_M_e_") -
 R("R_EX_26dap_M_e__Rev")  = 0.

M26dapMc: R("R_DAPabc") - R("R_DAPDC") + R("R_DAPE") -
 R("R_DAPE_Rev") - R("R_UAAGDS")  = 0.

M26dapMe: -R("R_DAPabc") - R("R_EX_26dap_M_e_") +
 R("R_EX_26dap_M_e__Rev")  = 0.

M2ahbutc: R("R_ACHBS") - R("R_KARA2i")  = 0.

M2aobutc: -R("R_GLYATi") + R("R_THRD")  = 0.

M2cpr5pc: -R("R_IGPS") + R("R_PRAIi")  = 0.

M2dda7pc: R("R_DDPA") - R("R_DHQS")  = 0.

M2ddg6pc: R("R_DDGLK") - R("R_EDA") + R("R_EDD")  = 0.

M2ddglcnb: R("M2ddglcnb_Transp") + R("R_EX_2ddglcn_e_") -
 R("R_EX_2ddglcn_e__Rev")  = 0.

M2ddglcnc: R("R_ALTRH") + R("R_DDGLCNt2r") - R("R_DDGLCNt2r_Rev") -
 R("R_DDGLK") + R("R_MNNH")  = 0.

M2ddglcne: -R("R_DDGLCNt2r") + R("R_DDGLCNt2r_Rev") -
 R("R_EX_2ddglcn_e_") + R("R_EX_2ddglcn_e__Rev")  = 0.

M2dh3dgal6pc: R("R_DDGALK") - R("R_DDPGALA") + R("R_DDPGALA_Rev")  =
 0.

M2dh3dgalc: -R("R_DDGALK") + R("R_GALCTND")  = 0.

M2dhglcnc: -R("R_2DGLCNRx") - R("R_2DGLCNRy")  = 0.

M2dhgulnc: -R("R_2DGULRx") - R("R_2DGULRy") + R("R_DKGLCNR1")  = 0.

M2dhpc: -R("R_DPR") + R("R_MOHMT")  = 0.

M2dmmq8c: -R("R_AMMQT8_2") + R("R_DHNAOT") + R("R_DMSOR2") +
 R("R_DMSOR2e") + R("R_FRD3") - R("R_G3PD7") - R("R_GLYCTO4") - R("R_HYD3") -
 R("R_NADH8") - R("R_NADH9") + R("R_TMAOR2") + R("R_TMAOR2e")  = 0.

M2dmmql8c: -R("R_DMSOR2") - R("R_DMSOR2e") - R("R_FRD3") +
 R("R_G3PD7") + R("R_GLYCTO4") + R("R_HYD3") + R("R_NADH8") + R("R_NADH9") -
 R("R_TMAOR2") - R("R_TMAOR2e")  = 0.

M2dr1pc: R("R_DURIPP") - R("R_DURIPP_Rev") - R("R_PPM2") +
 R("R_PPM2_Rev") + R("R_PUNP2") - R("R_PUNP2_Rev") + R("R_PUNP4") -
 R("R_PUNP4_Rev") + R("R_PUNP6") - R("R_PUNP6_Rev") + R("R_TMDPP") -
 R("R_TMDPP_Rev")  = 0.

M2dr5pc: -R("R_DRPA") + R("R_PPM2") - R("R_PPM2_Rev")  = 0.

M2h3oppanc: R("R_GLCRAL") + R("R_GLXCL") + R("R_HPYRI") -
 R("R_HPYRI_Rev") - R("R_TRSAR")  = 0.

M2ippmc: R("R_IPPMIa") - R("R_IPPMIa_Rev") - R("R_IPPMIb") +
 R("R_IPPMIb_Rev")  = 0.

M2kmbc: R("R_DKMPPD") + R("R_DKMPPD2") - R("R_UNK3")  = 0.

M2mahmpc: R("R_PMPK") - R("R_TMPPP")  = 0.

M2mcacnc: R("R_MCITD") - R("R_MICITD")  = 0.

M2mcitc: -R("R_MCITD") + R("R_MCITS")  = 0.

M2me4pc: R("R_DXPRIi") - R("R_MEPCT")  = 0.

M2mecdpc: -R("R_MECDPDH") + R("R_MECDPS")  = 0.

M2obutc: -R("R_ACHBS") - R("R_OBTFL") + R("R_THRD_L")  = 0.

M2ohphc: -R("R_OHPHM") + R("R_OPHHX")  = 0.

M2ombzlc: -R("R_OMBZLM") + R("R_OMPHHX")  = 0.

M2omhmblc: -R("R_DMQMT") + R("R_OMMBLHX")  = 0.

M2ommblc: R("R_OMBZLM") - R("R_OMMBLHX")  = 0.

M2omphc: R("R_OHPHM") - R("R_OMPHHX")  = 0.

M2ophc: R("R_OPHBDC") - R("R_OPHHX")  = 0.

M2p4c2mec: R("R_CDPMEK") - R("R_MECDPS")  = 0.

M2pgc: -R("R_ENO") + R("R_ENO_Rev") - R("R_PGM") + R("R_PGM_Rev")  =
 0.

M2pglycc: -R("R_PGLYCP")  = 0.

M2shchcc: R("R_SHCHCS2") - R("R_SUCBZS")  = 0.

M34hppc: R("R_PPND") + R("R_TYRTA") - R("R_TYRTA_Rev")  = 0.

M3c2hmpc: -R("R_IPMD") - R("R_IPPMIa") + R("R_IPPMIa_Rev")  = 0.

M3c3hmpc: R("R_IPPMIb") - R("R_IPPMIb_Rev") + R("R_IPPS")  = 0.

M3c4mopc: R("R_IPMD") - R("R_OMCDC")  = 0.

M3dhgulnc: R("R_DOGULNR")  = 0.

M3dhgulnpc: -R("R_KG6PDC")  = 0.

M3dhqc: -R("R_DHQD") + R("R_DHQD_Rev") + R("R_DHQS")  = 0.

M3dhskc: R("R_DHQD") - R("R_DHQD_Rev") - R("R_SHK3Dr") +
 R("R_SHK3Dr_Rev")  = 0.

M3hcinnmb: R("M3hcinnmb_Transp") + R("R_EX_3hcinnm_e_") -
 R("R_EX_3hcinnm_e__Rev")  = 0.

M3hcinnmc: -R("R_3HCINNMH") + R("R_HCINNMt2r") -
 R("R_HCINNMt2r_Rev")  = 0.

M3hcinnme: -R("R_EX_3hcinnm_e_") + R("R_EX_3hcinnm_e__Rev") -
 R("R_HCINNMt2r") + R("R_HCINNMt2r_Rev")  = 0.

M3hmrsACPc: R("R_KAS16") - R("R_U23GAAT") - R("R_UAGAAT") +
 R("R_UAGAAT_Rev")  = 0.

M3hpppnb: R("M3hpppnb_Transp") + R("R_EX_3hpppn_e_") -
 R("R_EX_3hpppn_e__Rev")  = 0.

M3hpppnc: -R("R_3HPPPNH") + R("R_HPPPNt2r") - R("R_HPPPNt2r_Rev")  =
 0.

M3hpppne: -R("R_EX_3hpppn_e_") + R("R_EX_3hpppn_e__Rev") -
 R("R_HPPPNt2r") + R("R_HPPPNt2r_Rev")  = 0.

M3ig3pc: R("R_IGPS") - R("R_TRPS1") - R("R_TRPS3")  = 0.

M3mobc: R("R_DHAD1") - R("R_IPPS") - R("R_MOHMT") + R("R_VALTA") -
 R("R_VALTA_Rev") - R("R_VPAMT")  = 0.

M3mopc: R("R_DHAD2") + R("R_ILETA") - R("R_ILETA_Rev")  = 0.

M3ophbc: R("R_HBZOPT") - R("R_OPHBDC")  = 0.

M3pgc: R("R_GLYCK") - R("R_PGCD") - R("R_PGK") + R("R_PGK_Rev") +
 R("R_PGM") - R("R_PGM_Rev")  = 0.

M3phpc: R("R_PGCD") - R("R_PSERT")  = 0.

M3psmec: -R("R_CHORS") + R("R_PSCVT") - R("R_PSCVT_Rev")  = 0.

M4abutb: R("M4abutb_Transp") + R("R_EX_4abut_e_") -
 R("R_EX_4abut_e__Rev")  = 0.

M4abutc: -R("R_ABTA") + R("R_ABUTD") + R("R_ABUTt2") -
 R("R_GLUABUTt7") + R("R_GLUABUTt7_Rev") + R("R_GLUDC")  = 0.

M4abute: -R("R_ABUTt2") - R("R_EX_4abut_e_") +
 R("R_EX_4abut_e__Rev") + R("R_GLUABUTt7") - R("R_GLUABUTt7_Rev")  = 0.

M4abutnc: -R("R_ABUTD") + R("R_PTRCTA")  = 0.

M4abzc: R("R_ADCL") - R("R_DHPS2")  = 0.

M4adchoc: -R("R_ADCL") + R("R_ADCS")  = 0.

M4ahmmpc: -R("R_HMPK1")  = 0.

M4ampmc: R("R_AMPMS") + R("R_HMPK1") - R("R_PMPK")  = 0.

M4c2mec: -R("R_CDPMEK") + R("R_MEPCT")  = 0.

M4h2opntnc: -R("R_HOPNTAL") + R("R_OP4ENH")  = 0.

M4hbac: R("R_THZPSN")  = 0.

M4hbzc: R("R_CHRPL") - R("R_HBZOPT")  = 0.

M4hthrc: R("R_4HTHRS")  = 0.

M4mhetzc: -R("R_HETZK")  = 0.

M4mopc: -R("R_LEUTAi") + R("R_OMCDC")  = 0.

M4mpetzc: R("R_HETZK") + R("R_THZPSN") - R("R_TMPPP")  = 0.

M4paspc: R("R_ASAD") - R("R_ASAD_Rev") + R("R_ASPK") -
 R("R_ASPK_Rev")  = 0.

M4perc: R("R_E4PD") - R("R_E4PD_Rev") - R("R_PERD") +
 R("R_PERD_Rev")  = 0.

M4ppanc: R("R_PNTK") - R("R_PPNCL2")  = 0.

M4ppcysc: -R("R_PPCDC") + R("R_PPNCL2")  = 0.

M4r5auc: R("R_PMDPHT") - R("R_RBFSa") + R("R_RBFSb")  = 0.

M5aizcc: -R("R_AIRC3") + R("R_AIRC3_Rev") - R("R_PRASCS") +
 R("R_PRASCS_Rev")  = 0.

M5aopc: R("R_G1SATi") - 2 R("R_PPBNGS")  = 0.

M5aprbuc: R("R_APRAUR") - R("R_PMDPHT")  = 0.

M5apruc: -R("R_APRAUR") + R("R_DHPPDA2")  = 0.

M5caizc: R("R_AIRC2") + R("R_AIRC3") - R("R_AIRC3_Rev")  = 0.

M5dglcnc: -R("R_5DGLCNR") + R("R_5DGLCNR_Rev") + R("R_DKGLCNR2x") +
 R("R_DKGLCNR2y") - R("R_IDOND") - R("R_IDOND2") + R("R_IDOND_Rev")  =
 0.

M5dh4dglcc: R("R_GALCTD") - R("R_GLCRAL") + R("R_GLCRD")  = 0.

M5mdr1pc: -R("R_MTRI") + R("R_MTRI_Rev") + R("R_MTRK")  = 0.

M5mdru1pc: -R("R_MDRPD") + R("R_MTRI") - R("R_MTRI_Rev")  = 0.

M5mtac: -R("R_MTAN") + R("R_SPMS")  = 0.

M5mthfc: -0.05 R("R_BiomassEcoli") - R("R_METS") + R("R_MTHFR2")  =
 0.

M5mtrc: R("R_MTAN") - R("R_MTRK")  = 0.

M5prdmbzc: R("R_NNDMBRT") - R("R_RZ5PP")  = 0.

M6hmhptc: R("R_DHNPA2") - R("R_HPPK2")  = 0.

M6hmhptppc: -R("R_DHPS2") + R("R_HPPK2")  = 0.

M6pgcc: -R("R_EDD") - R("R_GND") + R("R_GNK") + R("R_PGL")  = 0.

M6pglc: R("R_G6PDH2r") - R("R_G6PDH2r_Rev") - R("R_PGL")  = 0.

M8aonnc: -R("R_AMAOTr") + R("R_AMAOTr_Rev") + R("R_AOXSr") -
 R("R_AOXSr_Rev")  = 0.

Maacaldc: R("R_TAUDO")  = 0.

Maacoac: R("R_ACACCT") + R("R_ACACT1r") - R("R_ACACT1r_Rev") +
 R("R_FAO4")  = 0.

Macacb: R("Macacb_Transp") + R("R_EX_acac_e_") -
 R("R_EX_acac_e__Rev")  = 0.

Macacc: -R("R_ACACCT") + R("R_ACACt2") - R("R_ACACt2_Rev")  = 0.

Macace: -R("R_ACACt2") + R("R_ACACt2_Rev") - R("R_EX_acac_e_") +
 R("R_EX_acac_e__Rev")  = 0.

MacACPc: R("R_ACOATA") - R("R_ACOATA_Rev") - R("R_KAS14") +
 R("R_MACPD")  = 0.

Macaldb: R("Macaldb_Transp") + R("R_EX_acald_e_") -
 R("R_EX_acald_e__Rev")  = 0.

Macaldc: -R("R_ACALDi") + R("R_ACALDt") - R("R_ACALDt_Rev") -
 R("R_ALDD2x") + R("R_DRPA") + R("R_ETHAAL") + R("R_HOPNTAL") + R("R_THRAr")
 - R("R_THRAr_Rev")  = 0.

Macalde: -R("R_ACALDt") + R("R_ACALDt_Rev") - R("R_EX_acald_e_") +
 R("R_EX_acald_e__Rev")  = 0.

Macb: R("Macb_Transp") + R("R_EX_ac_e_") - R("R_EX_ac_e__Rev")  = 0.

Macc: R("R_ACACCT") - R("R_ACKr") + R("R_ACKr_Rev") + R("R_ACODA") -
 R("R_ACS") + R("R_ACt2r") - R("R_ACt2r_Rev") + R("R_AGDC") +
 R("R_ALDD2x") + R("R_BUTCT") + R("R_CITL") + R("R_CYSS") + R("R_NACODA") +
 R("R_POX") + R("R_UHGADA")  = 0.

Maccoac: -R("R_ACACCT") - 2 R("R_ACACT1r") + 2 R("R_ACACT1r_Rev") +
 R("R_ACALDi") - R("R_ACCOACr") + R("R_ACCOACr_Rev") - R("R_ACGS") -
 R("R_ACOATA") + R("R_ACOATA_Rev") + R("R_ACS") - R("R_ADHEr") +
 R("R_ADHEr_Rev") - 0.00005 R("R_BiomassEcoli") - R("R_BUTCT") - R("R_CS") + 7
 R("R_FAO1") + 8 R("R_FAO2") + 9 R("R_FAO3") - R("R_G1PACT") +
 R("R_GLYATi") - R("R_IPPS") - R("R_KAS15") - R("R_MALS") + R("R_PDH") +
 R("R_PFL") - R("R_PTAr") + R("R_PTAr_Rev") - R("R_SERAT") + R("R_SERAT_Rev")
 - R("R_SPMDAT1") - R("R_SPMDAT2") - R("R_TDPADGAT")  = 0.

Mace: -R("R_ACt2r") + R("R_ACt2r_Rev") - R("R_EX_ac_e_") +
 R("R_EX_ac_e__Rev")  = 0.

Macg5pc: R("R_ACGK") + R("R_AGPR") - R("R_AGPR_Rev")  = 0.

Macg5sac: R("R_ACOTA") - R("R_ACOTA_Rev") - R("R_AGPR") +
 R("R_AGPR_Rev") - R("R_NACODA")  = 0.

Macgam1pc: R("R_G1PACT") - R("R_UAGDP")  = 0.

Macgam6pc: R("R_ACGApts") - R("R_AGDC") + R("R_AMANAPE")  = 0.

Macgamb: R("Macgamb_Transp") + R("R_EX_acgam_e_") -
 R("R_EX_acgam_e__Rev")  = 0.

Macgame: -R("R_ACGApts") - R("R_EX_acgam_e_") +
 R("R_EX_acgam_e__Rev")  = 0.

Macgluc: -R("R_ACGK") + R("R_ACGS")  = 0.

Macmanab: R("Macmanab_Transp") + R("R_EX_acmana_e_") -
 R("R_EX_acmana_e__Rev")  = 0.

Macmanac: R("R_ACNML") - R("R_AMANK")  = 0.

Macmanae: -R("R_ACMANApts") - R("R_EX_acmana_e_") +
 R("R_EX_acmana_e__Rev")  = 0.

Macmanapc: R("R_ACMANApts") - R("R_AMANAPE") + R("R_AMANK")  = 0.

Macnamb: R("Macnamb_Transp") + R("R_EX_acnam_e_") -
 R("R_EX_acnam_e__Rev")  = 0.

Macnamc: R("R_ACNAMt2") - R("R_ACNML")  = 0.

Macname: -R("R_ACNAMt2") - R("R_EX_acnam_e_") +
 R("R_EX_acnam_e__Rev")  = 0.

Maconmc: R("R_ACONMT")  = 0.

MaconTc: -R("R_ACONMT")  = 0.

Macornc: -R("R_ACODA") - R("R_ACOTA") + R("R_ACOTA_Rev")  = 0.

MACPc: -R("R_AACPS1") - R("R_AACPS2") - R("R_AACPS3") -
 R("R_AACPS4") - R("R_AACPS5") - R("R_ACOATA") + R("R_ACOATA_Rev") + R("R_ACPS1") +
 4 R("R_C120SN") + 5 R("R_C140SN") + 5 R("R_C141SN") + 6 R("R_C160SN")
 + 6 R("R_C161SN") + 7 R("R_C181SN") + R("R_EDTXS1") + R("R_EDTXS2") +
 R("R_EDTXS3") + R("R_EDTXS4") + R("R_KAS14") + R("R_KAS16") -
 R("R_MCOATA") + R("R_MCOATA_Rev") + 2 R("R_PASYN_EC") + R("R_U23GAAT") +
 R("R_UAGAAT") - R("R_UAGAAT_Rev")  = 0.

Macserc: -R("R_CYSS") + R("R_SERAT") - R("R_SERAT_Rev")  = 0.

MactACPc: -R("R_C120SN") - R("R_C140SN") - R("R_C141SN") -
 R("R_C160SN") - R("R_C161SN") - R("R_C181SN") + R("R_KAS14") + R("R_KAS15")  =
 0.

Mactpc: R("R_ACKr") - R("R_ACKr_Rev") + R("R_PTAr") -
 R("R_PTAr_Rev")  = 0.

Madeb: R("Madeb_Transp") + R("R_EX_ade_e_") - R("R_EX_ade_e__Rev") 
 = 0.

Madec: -R("R_ADD") + R("R_ADEt2r") - R("R_ADEt2r_Rev") - R("R_ADPT")
 + R("R_AHCYSNS") + R("R_AMPN") + R("R_MTAN") + R("R_PUNP1") -
 R("R_PUNP1_Rev") + R("R_PUNP2") - R("R_PUNP2_Rev")  = 0.

Madee: -R("R_ADEt2r") + R("R_ADEt2r_Rev") - R("R_EX_ade_e_") +
 R("R_EX_ade_e__Rev")  = 0.

Madnb: R("Madnb_Transp") + R("R_EX_adn_e_") - R("R_EX_adn_e__Rev") 
 = 0.

Madnc: -R("R_ADA") - R("R_ADNK1") + R("R_ADNt2") + R("R_ADNt2r") -
 R("R_ADNt2r_Rev") + R("R_AHC") - R("R_AHC_Rev") + R("R_NTD7") -
 R("R_PUNP1") + R("R_PUNP1_Rev")  = 0.

Madne: -R("R_ADNt2") - R("R_ADNt2r") + R("R_ADNt2r_Rev") -
 R("R_EX_adn_e_") + R("R_EX_adn_e__Rev")  = 0.

Madocbic: -R("R_ADOCBIK") + R("R_CBIAT_DELETE") -
 R("R_CBIAT_DELETE_Rev")  = 0.

Madocbipc: -R("R_ACBIPGT") + R("R_ADOCBIK")  = 0.

Madocblc: R("R_ADOCBLS") + R("R_CBLAT_DELETE") -
 R("R_CBLAT_DELETE_Rev")  = 0.

Madpc: R("R_ACCOACr") - R("R_ACCOACr_Rev") + R("R_ACCOAL") +
 R("R_ACGK") + R("R_ACKr") - R("R_ACKr_Rev") + 2 R("R_ADK1") - 2
 R("R_ADK1_Rev") + R("R_ADK3") - R("R_ADK3_Rev") + R("R_ADK4") - R("R_ADK4_Rev") +
 R("R_ADNK1") + R("R_ADOCBIK") + R("R_ADSK") + R("R_AIRC2") +
 R("R_ALAabc") + R("R_ALAALAr") - R("R_ALAALAr_Rev") + R("R_AMANK") + 2
 R("R_AP4AH") + R("R_AP5AH") + R("R_ARBabc") + R("R_ARGabc") + R("R_ASNabc") +
 R("R_ASPabc") + R("R_ASPK") - R("R_ASPK_Rev") + R("R_ATPM") -
 R("R_ATPS4r") + R("R_ATPS4r_Rev") + 45.5608 R("R_BiomassEcoli") +
 R("R_CBL1abc") + R("R_CBMK") + 2 R("R_CBPS") + R("R_CDPMEK") + R("R_CHLabc") +
 R("R_CTPS2") + R("R_CYSabc") + R("R_CYTK1") - R("R_CYTK1_Rev") +
 R("R_CYTK2") - R("R_CYTK2_Rev") + R("R_DADK") - R("R_DADK_Rev") +
 R("R_DAGK_EC") + R("R_DAPabc") + R("R_DBTSr") - R("R_DBTSr_Rev") + R("R_DDGALK") +
 R("R_DDGLK") + R("R_DGK1") - R("R_DGK1_Rev") + R("R_DHFS") +
 R("R_DPCOAK") + R("R_DTMPK") - R("R_DTMPK_Rev") + R("R_DURIK1") + R("R_DXYLK")
 + R("R_FCLK") + R("R_FE2abc") + R("R_FRUK") + R("R_GALabc") +
 R("R_GALKr") - R("R_GALKr_Rev") + R("R_GART") + R("R_GK1") - R("R_GK1_Rev") +
 R("R_GLCS1") + R("R_GLNabc") + R("R_GLNS") + R("R_GLU5K") +
 R("R_GLUabc") + R("R_GLUCYS") + R("R_GLYBabc") + R("R_GLYC3Pabc") + R("R_GLYCK")
 + R("R_GLYK") + R("R_GMHEPK") + R("R_GNK") + R("R_GSNK") +
 R("R_GSPMDS") + R("R_GTHS") + R("R_HETZK") + R("R_HEX1") + R("R_HISabc") +
 R("R_HMPK1") + R("R_HSK") + R("R_ILEabc") + R("R_INSK") + R("R_Kabc") +
 R("R_LEUabc") + 3 R("R_LPSSYN_EC") + R("R_LYSabc") + R("R_MALTabc") +
 R("R_MALTHXabc") + R("R_MALTPTabc") + R("R_MALTTRabc") +
 R("R_MALTTTRabc") + R("R_METabc") + R("R_METDabc") + R("R_MTRK") + R("R_NADK") +
 R("R_NAMNPP") + R("R_NDPK1") - R("R_NDPK1_Rev") + R("R_NDPK2") -
 R("R_NDPK2_Rev") + R("R_NDPK3") - R("R_NDPK3_Rev") + R("R_NDPK4") -
 R("R_NDPK4_Rev") + R("R_NDPK5") - R("R_NDPK5_Rev") + R("R_NDPK6") -
 R("R_NDPK6_Rev") + R("R_NDPK7") - R("R_NDPK7_Rev") + R("R_NDPK8") - R("R_NDPK8_Rev")
 + R("R_ORNabc") + R("R_PFK") + R("R_PFK_2") + R("R_PGK") -
 R("R_PGK_Rev") + R("R_PIabc") + R("R_PMPK") + R("R_PNTK") - R("R_PPAKr") +
 R("R_PPAKr_Rev") + R("R_PPCK") + R("R_PRAGSr") - R("R_PRAGSr_Rev") +
 R("R_PRAIS") + R("R_PRASCS") - R("R_PRASCS_Rev") + R("R_PRFGS") +
 R("R_PROabc") + R("R_PTRCabc") + R("R_PYDAMK") + R("R_PYDXK") + R("R_PYDXNK") -
 R("R_PYK") + R("R_RBFK") + R("R_RBK") + R("R_RBK_L1") + R("R_RIBabc")
 + R("R_RMK") - R("R_RNDR1") + R("R_SHKK") + R("R_SPMDabc") +
 R("R_SUCCabc") + R("R_SUCOAS") - R("R_SUCOAS_Rev") + R("R_SULabc") +
 R("R_TAURabc") + R("R_TDSK") + R("R_THMabc") + R("R_THRabc") + R("R_TMDK1") +
 R("R_TMKr") - R("R_TMKr_Rev") + R("R_TMPKr") - R("R_TMPKr_Rev") +
 R("R_TSULabc") + R("R_UAAGDS") + R("R_UAMAGS") + R("R_UAMAS") +
 R("R_UGMDDS") + R("R_UMPK") - R("R_UMPK_Rev") + R("R_URIDK2r") -
 R("R_URIDK2r_Rev") + R("R_VALabc") + R("R_XYLabc") + R("R_XYLK")  = 0.

Madpglcc: -R("R_GLCS1") + R("R_GLGC")  = 0.

MadphepDDc: -R("R_AGMHE") + R("R_GMHEPAT")  = 0.

MadphepLDc: R("R_AGMHE") - 3 R("R_LPSSYN_EC")  = 0.

Magdpcbic: R("R_ACBIPGT") - R("R_ADOCBLS")  = 0.

Magmc: -R("R_AGMT") + R("R_ARGDC")  = 0.

MagpcECc: -0.02 R("R_LPLIPA3") - 0.02 R("R_LPLIPA5") + 0.02
 R("R_PLIPA3")  = 0.

MagpeECc: -0.02 R("R_LPLIPA2") - 0.02 R("R_LPLIPA4") + 0.02
 R("R_PLIPA2")  = 0.

MagpgECc: -0.02 R("R_LPLIPA1") - 0.02 R("R_LPLIPA6") + 0.02
 R("R_PLIPA1")  = 0.

Mahcysc: R("R_ACONMT") - R("R_AHC") + R("R_AHC_Rev") -
 R("R_AHCYSNS") + R("R_AMMQT8_2") + R("R_DMQMT") + R("R_OHPHM") + R("R_OMBZLM") + 2
 R("R_UPP3MT")  = 0.

Mahdtc: -R("R_DNTPPA") + R("R_GTPCI")  = 0.

Maicarc: R("R_ADSL2r") - R("R_ADSL2r_Rev") - R("R_AICART") +
 R("R_AICART_Rev") + R("R_IG3PS")  = 0.

Mairc: -R("R_AIRC2") - R("R_AMPMS") + R("R_PRAIS")  = 0.

Makgb: R("Makgb_Transp") + R("R_EX_akg_e_") - R("R_EX_akg_e__Rev") 
 = 0.

Makgc: -R("R_ABTA") - R("R_ACOTA") + R("R_ACOTA_Rev") - R("R_AKGDH")
 + R("R_AKGt2r") - R("R_AKGt2r_Rev") - R("R_ALATA_L") +
 R("R_ALATA_L_Rev") - R("R_ASPTA") + R("R_ASPTA_Rev") + R("R_GLUDy") -
 R("R_GLUDy_Rev") - R("R_GLUSy") + R("R_HSTPT") + R("R_ICDHyr") - R("R_ICDHyr_Rev") -
 R("R_ILETA") + R("R_ILETA_Rev") + R("R_LEUTAi") + R("R_OHPBAT") -
 R("R_OHPBAT_Rev") - R("R_ORNTA") - R("R_OXGDC2") - R("R_PHETA1") +
 R("R_PHETA1_Rev") + R("R_PSERT") - R("R_PTRCTA") - R("R_SDPTA") +
 R("R_SDPTA_Rev") - R("R_SOTA") - R("R_TAUDO") + R("R_TDPAGTA") - R("R_TYRTA") +
 R("R_TYRTA_Rev") + R("R_UNK3") - R("R_VALTA") + R("R_VALTA_Rev")  =
 0.

Makge: -R("R_AKGt2r") + R("R_AKGt2r_Rev") - R("R_EX_akg_e_") +
 R("R_EX_akg_e__Rev")  = 0.

Malaalac: R("R_ALAALAr") - R("R_ALAALAr_Rev") - R("R_UGMDDS")  = 0.

MalaBc: R("R_ASP1DC") - R("R_PANTS")  = 0.

MalacSc: R("R_ACLS") - R("R_KARA1i")  = 0.

MalaDb: R("MalaDb_Transp") + R("R_EX_ala_D_e_") -
 R("R_EX_ala_D_e__Rev")  = 0.

MalaDc: -2 R("R_ALAALAr") + 2 R("R_ALAALAr_Rev") + R("R_ALAR") +
 R("R_ALARi") - R("R_ALAR_Rev") - R("R_ALATA_D2") - R("R_DAAD") +
 R("R_DALAt2r") - R("R_DALAt2r_Rev")  = 0.

MalaDe: -R("R_DALAt2r") + R("R_DALAt2r_Rev") - R("R_EX_ala_D_e_") +
 R("R_EX_ala_D_e__Rev")  = 0.

MalaLb: R("MalaLb_Transp") + R("R_EX_ala_L_e_") -
 R("R_EX_ala_L_e__Rev")  = 0.

MalaLc: R("R_ALAabc") - R("R_ALAR") - R("R_ALARi") + R("R_ALAR_Rev")
 + R("R_ALAt2r") - R("R_ALAt2r_Rev") - R("R_ALATA_L") -
 R("R_ALATA_L2") + R("R_ALATA_L_Rev") - R("R_AOXSr") + R("R_AOXSr_Rev") - 0.488
 R("R_BiomassEcoli") + R("R_BTS2") - R("R_BTS2_Rev") + R("R_THZPSN") -
 R("R_UAMAS") - R("R_VPAMT")  = 0.

MalaLe: -R("R_ALAabc") - R("R_ALAt2r") + R("R_ALAt2r_Rev") -
 R("R_EX_ala_L_e_") + R("R_EX_ala_L_e__Rev")  = 0.

Malltnb: R("Malltnb_Transp") + R("R_EX_alltn_e_") -
 R("R_EX_alltn_e__Rev")  = 0.

Malltnc: -R("R_ALLTN") + R("R_ALLTNt2r") - R("R_ALLTNt2r_Rev")  = 0.

Malltne: -R("R_ALLTNt2r") + R("R_ALLTNt2r_Rev") - R("R_EX_alltn_e_")
 + R("R_EX_alltn_e__Rev")  = 0.

Mallttc: -R("R_ALLTAH") + R("R_ALLTN")  = 0.

Maltrnc: -R("R_ALTRH") - R("R_TAGURr") + R("R_TAGURr_Rev")  = 0.

Mametamc: R("R_ADMDCr") - R("R_ADMDCr_Rev") - R("R_SPMS")  = 0.

Mametc: -R("R_ACONMT") - R("R_ADMDCr") + R("R_ADMDCr_Rev") -
 R("R_AMAOTr") + R("R_AMAOTr_Rev") - R("R_AMMQT8_2") - R("R_DMQMT") +
 R("R_METAT") - R("R_OHPHM") - R("R_OMBZLM") - 2 R("R_UPP3MT")  = 0.

Mamobc: R("R_AMAOTr") - R("R_AMAOTr_Rev")  = 0.

Mampb: R("Mampb_Transp") + R("R_EX_amp_e_") - R("R_EX_amp_e__Rev") 
 = 0.

Mampc: R("R_AACPS1") + R("R_AACPS2") + R("R_AACPS3") + R("R_AACPS4")
 + R("R_AACPS5") + R("R_ACS") - R("R_ADK1") + R("R_ADK1_Rev") -
 R("R_ADK3") + R("R_ADK3_Rev") - R("R_ADK4") + R("R_ADK4_Rev") + R("R_ADNK1")
 + R("R_ADPT") + R("R_ADSL1r") - R("R_ADSL1r_Rev") - R("R_AMPN") +
 R("R_ARGSS") + R("R_ASNS1") + R("R_ASNS2") - 0.001 R("R_BiomassEcoli") +
 R("R_BPNT") + 6 R("R_ENTCS") + R("R_FAO1") + R("R_FAO2") + R("R_FAO3")
 + R("R_GLUTRS") + R("R_GMPS2") + R("R_HPPK2") + R("R_NADDP") +
 R("R_NADS1") - R("R_NTD7") + R("R_NTPP6") + R("R_PACCOAL") + R("R_PANTS") +
 R("R_PPS") + R("R_PRPPS") - R("R_PRPPS_Rev") + R("R_SELNPS") +
 R("R_SUCBZL") + R("R_THZPSN")  = 0.

Mampe: -R("R_EX_amp_e_") + R("R_EX_amp_e__Rev") + R("R_NADDPe")  =
 0.

Manthc: -R("R_ANPRT") + R("R_ANS")  = 0.

Map4ac: -R("R_AP4AH")  = 0.

Map5ac: -R("R_AP5AH")  = 0.

MapgECc: 0.02 R("R_LPLIPA4") + 0.02 R("R_LPLIPA5") + 0.02
 R("R_LPLIPA6")  = 0.

MapoACPc: -R("R_ACPS1")  = 0.

Mapsc: -R("R_ADSK") + R("R_SADT2")  = 0.

Mara5pc: R("R_A5PISO") - R("R_A5PISO_Rev") - R("R_KDOPS")  = 0.

MarabLb: R("MarabLb_Transp") + R("R_EX_arab_L_e_") -
 R("R_EX_arab_L_e__Rev")  = 0.

MarabLc: -R("R_ARAI") + R("R_ARAI_Rev") + R("R_ARBabc") +
 R("R_ARBt2r") - R("R_ARBt2r_Rev")  = 0.

MarabLe: -R("R_ARBabc") - R("R_ARBt2r") + R("R_ARBt2r_Rev") -
 R("R_EX_arab_L_e_") + R("R_EX_arab_L_e__Rev")  = 0.

Marbt6pc: -R("R_AB6PGH")  = 0.

MargLb: R("MargLb_Transp") + R("R_EX_arg_L_e_") -
 R("R_EX_arg_L_e__Rev")  = 0.

MargLc: R("R_ARGabc") - R("R_ARGDC") + R("R_ARGORNt7") -
 R("R_ARGORNt7_Rev") + R("R_ARGSL") - R("R_ARGSL_Rev") - R("R_AST") - 0.281
 R("R_BiomassEcoli")  = 0.

MargLe: -R("R_ARGabc") - R("R_ARGORNt7") + R("R_ARGORNt7_Rev") -
 R("R_EX_arg_L_e_") + R("R_EX_arg_L_e__Rev")  = 0.

Margsucc: -R("R_ARGSL") + R("R_ARGSL_Rev") + R("R_ARGSS")  = 0.

MasnLb: R("MasnLb_Transp") + R("R_EX_asn_L_e_") -
 R("R_EX_asn_L_e__Rev")  = 0.

MasnLc: R("R_ASNabc") - R("R_ASNN") + R("R_ASNS1") + R("R_ASNS2") +
 R("R_ASNt2r") - R("R_ASNt2r_Rev") - 0.229 R("R_BiomassEcoli")  = 0.

MasnLe: -R("R_ASNabc") - R("R_ASNt2r") + R("R_ASNt2r_Rev") -
 R("R_EX_asn_L_e_") + R("R_EX_asn_L_e__Rev")  = 0.

MaspLb: R("MaspLb_Transp") + R("R_EX_asp_L_e_") -
 R("R_EX_asp_L_e__Rev")  = 0.

MaspLc: -R("R_ADSS") - R("R_ARGSS") + R("R_ASNN") - R("R_ASNS1") -
 R("R_ASNS2") - R("R_ASP1DC") + R("R_ASPabc") - R("R_ASPCT") -
 R("R_ASPK") + R("R_ASPK_Rev") - R("R_ASPO3") - R("R_ASPO4") - R("R_ASPO5") -
 R("R_ASPO6") - R("R_ASPT") + R("R_ASPt2") + R("R_ASPt2_2") +
 R("R_ASPt2_3") - R("R_ASPTA") + R("R_ASPTA_Rev") - 0.229 R("R_BiomassEcoli") -
 R("R_PRASCS") + R("R_PRASCS_Rev")  = 0.

MaspLe: -R("R_ASPabc") - R("R_ASPt2") - R("R_ASPt2_2") -
 R("R_ASPt2_3") - R("R_EX_asp_L_e_") + R("R_EX_asp_L_e__Rev")  = 0.

Maspsac: -R("R_ASAD") + R("R_ASAD_Rev") - R("R_DHDPS") + R("R_HSDy")
 - R("R_HSDy_Rev")  = 0.

Matpc: -R("R_AACPS1") - R("R_AACPS2") - R("R_AACPS3") -
 R("R_AACPS4") - R("R_AACPS5") - R("R_ACCOACr") + R("R_ACCOACr_Rev") -
 R("R_ACCOAL") - R("R_ACGK") - R("R_ACKr") + R("R_ACKr_Rev") - R("R_ACS") -
 R("R_ADK1") + R("R_ADK1_Rev") - R("R_ADNCYC") - R("R_ADNK1") - R("R_ADOCBIK")
 - R("R_ADSK") - R("R_AIRC2") - R("R_ALAabc") - R("R_ALAALAr") +
 R("R_ALAALAr_Rev") - R("R_AMANK") + R("R_AP5AH") - R("R_ARBabc") -
 R("R_ARGabc") - R("R_ARGSS") - R("R_ASNabc") - R("R_ASNS1") - R("R_ASNS2") -
 R("R_ASPabc") - R("R_ASPK") + R("R_ASPK_Rev") - R("R_ATPM") -
 R("R_ATPPRT") + R("R_ATPS4r") - R("R_ATPS4r_Rev") - 45.7318 R("R_BiomassEcoli")
 - R("R_CBIAT_DELETE") + R("R_CBIAT_DELETE_Rev") - R("R_CBL1abc") -
 R("R_CBLAT_DELETE") + R("R_CBLAT_DELETE_Rev") - R("R_CBMK") - 2
 R("R_CBPS") - R("R_CDPMEK") - R("R_CHLabc") - R("R_CTPS2") - R("R_CYSabc") -
 R("R_CYTK1") + R("R_CYTK1_Rev") - R("R_CYTK2") + R("R_CYTK2_Rev") -
 R("R_DADK") + R("R_DADK_Rev") - R("R_DAGK_EC") - R("R_DAPabc") -
 R("R_DBTSr") + R("R_DBTSr_Rev") - R("R_DDGALK") - R("R_DDGLK") - R("R_DGK1") +
 R("R_DGK1_Rev") - R("R_DHBSr") + R("R_DHBSr_Rev") - R("R_DHFS") -
 R("R_DPCOAK") - R("R_DTMPK") + R("R_DTMPK_Rev") - R("R_DURIK1") -
 R("R_DXYLK") - R("R_FAO1") - R("R_FAO2") - R("R_FAO3") - R("R_FCLK") -
 R("R_FE2abc") - R("R_FMNAT") - R("R_FRUK") - R("R_GALabc") - R("R_GALKr") +
 R("R_GALKr_Rev") - R("R_GART") - R("R_GK1") + R("R_GK1_Rev") -
 R("R_GLGC") - R("R_GLNabc") - R("R_GLNS") - R("R_GLU5K") - R("R_GLUabc") -
 R("R_GLUCYS") - R("R_GLUTRS") - R("R_GLYBabc") - R("R_GLYC3Pabc") -
 R("R_GLYCK") - R("R_GLYK") - R("R_GMHEPAT") - R("R_GMHEPK") - R("R_GMPS2")
 - R("R_GNK") - R("R_GSNK") - R("R_GSPMDS") - R("R_GTHS") -
 R("R_HETZK") - R("R_HEX1") - R("R_HISabc") - R("R_HMPK1") - R("R_HPPK2") -
 R("R_HSK") - R("R_ILEabc") - R("R_INSK") - R("R_Kabc") - R("R_LEUabc") -
 R("R_LYSabc") - R("R_MALTabc") - R("R_MALTHXabc") - R("R_MALTPTabc") -
 R("R_MALTTRabc") - R("R_MALTTTRabc") - R("R_METabc") - R("R_METAT") -
 R("R_METDabc") - R("R_MTRK") - R("R_NADK") - R("R_NADS1") -
 R("R_NAMNPP") - R("R_NDPK1") + R("R_NDPK1_Rev") - R("R_NDPK2") + R("R_NDPK2_Rev")
 - R("R_NDPK3") + R("R_NDPK3_Rev") - R("R_NDPK4") + R("R_NDPK4_Rev") -
 R("R_NDPK5") + R("R_NDPK5_Rev") - R("R_NDPK6") + R("R_NDPK6_Rev") -
 R("R_NDPK7") + R("R_NDPK7_Rev") - R("R_NDPK8") + R("R_NDPK8_Rev") -
 R("R_NMNAT") - R("R_NNAT") - R("R_NTPP6") - R("R_ORNabc") -
 R("R_PACCOAL") - R("R_PANTS") - R("R_PFK") - R("R_PFK_2") - R("R_PGK") +
 R("R_PGK_Rev") - R("R_PIabc") - R("R_PMPK") - R("R_PNTK") + R("R_PPAKr") -
 R("R_PPAKr_Rev") - R("R_PPCK") - R("R_PPS") - R("R_PRAGSr") +
 R("R_PRAGSr_Rev") - R("R_PRAIS") - R("R_PRASCS") + R("R_PRASCS_Rev") - R("R_PRFGS")
 - R("R_PROabc") - R("R_PRPPS") + R("R_PRPPS_Rev") - R("R_PTPATi") -
 R("R_PTRCabc") - R("R_PYDAMK") - R("R_PYDXK") - R("R_PYDXNK") +
 R("R_PYK") - R("R_RBFK") - R("R_RBK") - R("R_RBK_L1") - R("R_RIBabc") -
 R("R_RMK") - R("R_RNTR1") - R("R_SADT2") - R("R_SELNPS") - R("R_SERASr") +
 R("R_SERASr_Rev") - R("R_SHKK") - R("R_SPMDabc") - R("R_SUCBZL") -
 R("R_SUCCabc") - R("R_SUCOAS") + R("R_SUCOAS_Rev") - R("R_SULabc") -
 R("R_TAURabc") - R("R_TDSK") - R("R_THMabc") - R("R_THRabc") -
 R("R_THZPSN") - R("R_TMDK1") - R("R_TMKr") + R("R_TMKr_Rev") - R("R_TMPKr") +
 R("R_TMPKr_Rev") - R("R_TSULabc") - R("R_UAAGDS") - R("R_UAMAGS") -
 R("R_UAMAS") - R("R_UGMDDS") - R("R_UMPK") + R("R_UMPK_Rev") -
 R("R_URIDK2r") + R("R_URIDK2r_Rev") - R("R_VALabc") - R("R_XYLabc") - R("R_XYLK") 
 = 0.

Mbbtcoac: -R("R_CRNBTCT") + R("R_CRNBTCT_Rev")  = 0.

Mbetaldc: -R("R_BETALDHx") - R("R_BETALDHy")  = 0.

Mbtcoac: R("R_BUTCT") - R("R_FAO4")  = 0.

Mbtnc: R("R_BSORx") + R("R_BSORy") + R("R_BTS2") - R("R_BTS2_Rev") 
 = 0.

Mbtnsoc: -R("R_BSORx") - R("R_BSORy")  = 0.

Mbutb: R("Mbutb_Transp") + R("R_EX_but_e_") - R("R_EX_but_e__Rev") 
 = 0.

Mbutc: -R("R_BUTCT") + R("R_BUTt2r") - R("R_BUTt2r_Rev")  = 0.

Mbute: -R("R_BUTt2r") + R("R_BUTt2r_Rev") - R("R_EX_but_e_") +
 R("R_EX_but_e__Rev")  = 0.

Mcampc: R("R_ADNCYC")  = 0.

Mcbaspc: R("R_ASPCT") + R("R_DHORTS") - R("R_DHORTS_Rev")  = 0.

Mcbic: -R("R_CBIAT_DELETE") + R("R_CBIAT_DELETE_Rev")  = 0.

Mcbl1b: R("Mcbl1b_Transp") + R("R_EX_cbl1_e_") -
 R("R_EX_cbl1_e__Rev")  = 0.

Mcbl1c: R("R_CBL1abc") - R("R_CBLAT_DELETE") +
 R("R_CBLAT_DELETE_Rev")  = 0.

Mcbl1e: -R("R_CBL1abc") - R("R_EX_cbl1_e_") + R("R_EX_cbl1_e__Rev") 
 = 0.

Mcbpc: -R("R_ASPCT") + R("R_CBMK") + R("R_CBPS") - R("R_OCBT") +
 R("R_OCBT_Rev")  = 0.

Mcdpc: R("R_CYTK1") - R("R_CYTK1_Rev") + 2 R("R_LPSSYN_EC") -
 R("R_NDPK3") + R("R_NDPK3_Rev") - R("R_RNDR3")  = 0.

McdpdagECc: -0.02 R("R_CDAPPA_EC") + 0.02 R("R_DASYN_EC") - 0.02
 R("R_DASYN_EC_Rev") - 0.02 R("R_PGSA_EC") + 0.02 R("R_PGSA_EC_Rev") -
 0.02 R("R_PSSA_EC") + 0.02 R("R_PSSA_EC_Rev")  = 0.

Mcdpeac: -2 R("R_LPSSYN_EC") + R("R_PEPT_EC") - R("R_PEPT_EC_Rev") 
 = 0.

Mcechdddc: -R("R_DHPPD") + R("R_PPPNDO")  = 0.

Mcenchdddc: R("R_CINNDO") - R("R_DHCIND")  = 0.

Mcholb: R("Mcholb_Transp") + R("R_EX_chol_e_") -
 R("R_EX_chol_e__Rev")  = 0.

Mcholc: R("R_CHLabc") + R("R_CHLt2r") - R("R_CHLt2r_Rev") +
 R("R_GPDDA1")  = 0.

Mchole: -R("R_CHLabc") - R("R_CHLt2r") + R("R_CHLt2r_Rev") -
 R("R_EX_chol_e_") + R("R_EX_chol_e__Rev")  = 0.

Mchorc: -R("R_ADCS") - R("R_ANS") - R("R_CHORM") + R("R_CHORS") -
 R("R_CHRPL") - R("R_ICHORSi")  = 0.

Mcinnmc: -R("R_CINNDO")  = 0.

Mcitb: R("Mcitb_Transp") + R("R_EX_cit_e_") - R("R_EX_cit_e__Rev") 
 = 0.

Mcitc: -R("R_ACONT") + R("R_ACONT_Rev") - R("R_CITL") + R("R_CITt7")
 + R("R_CS")  = 0.

Mcite: -R("R_CITt7") - R("R_EX_cit_e_") + R("R_EX_cit_e__Rev")  = 0.

McitrLc: -R("R_ARGSS") + R("R_OCBT") - R("R_OCBT_Rev")  = 0.

Mckdoc: R("R_KDOCT2") - 3 R("R_LPSSYN_EC") - R("R_MOAT") -
 R("R_MOAT2")  = 0.

MclpnECc: -0.000129 R("R_BiomassEcoli") + 0.02 R("R_CLPNS_EC") -
 0.02 R("R_CLPNS_EC_Rev")  = 0.

Mcmpc: R("R_CDAPPA_EC") - R("R_CMPN") + R("R_CYTDK2") - R("R_CYTK1")
 + R("R_CYTK1_Rev") + 3 R("R_LPSSYN_EC") + R("R_MECDPS") + R("R_MOAT")
 + R("R_MOAT2") - R("R_NTD4") + R("R_NTPP4") - R("R_PEPT_EC") +
 R("R_PEPT_EC_Rev") + R("R_PGSA_EC") - R("R_PGSA_EC_Rev") + R("R_PPNCL2") +
 R("R_PSSA_EC") - R("R_PSSA_EC_Rev")  = 0.

Mco2b: R("Mco2b_Transp") + R("R_EX_co2_e_") - R("R_EX_co2_e__Rev") 
 = 0.

Mco2c: R("R_ACHBS") + R("R_ACLS") + R("R_ADMDCr") -
 R("R_ADMDCr_Rev") + R("R_AKGDH") + R("R_AOXSr") - R("R_AOXSr_Rev") + R("R_ARGDC") +
 R("R_ASP1DC") + 4 R("R_C120SN") + 5 R("R_C140SN") + 5 R("R_C141SN") + 6
 R("R_C160SN") + 6 R("R_C161SN") + 7 R("R_C181SN") - R("R_CBMK") +
 R("R_CO2t") - R("R_CO2t_Rev") + 2 R("R_CPPPGO") + 2 R("R_CYNTAH") +
 R("R_DAPDC") - R("R_DBTSr") + R("R_DBTSr_Rev") + R("R_DHNAOT") + R("R_DXPS")
 + R("R_FDH2") + R("R_FDH3") + R("R_FHL") + R("R_GLUDC") +
 R("R_GLXCL") + R("R_GLYCL") + R("R_GND") - R("R_HCO3E") + R("R_HCO3E_Rev") +
 R("R_ICDHyr") - R("R_ICDHyr_Rev") + R("R_IGPS") + R("R_KAS14") +
 R("R_KAS15") + R("R_KAS16") + R("R_KG6PDC") + R("R_LYSDC") + R("R_MACPD") +
 R("R_ME1") + R("R_ME2") + R("R_MMCD") + R("R_NNDPR") + R("R_OMCDC") +
 R("R_OMPDC") + R("R_OPHBDC") + R("R_ORNDC") + R("R_OXGDC2") + R("R_PDH")
 + R("R_PDX5PS") + R("R_POX") - R("R_PPC") + R("R_PPCDC") + R("R_PPCK")
 + R("R_PPND") + R("R_PPNDH") + R("R_PSD_EC") + R("R_SADH") +
 R("R_TAUDO") + R("R_THZPSN") + R("R_UGLYCH") + 4 R("R_UPPDC1")  = 0.

Mco2e: -R("R_CO2t") + R("R_CO2t_Rev") - R("R_EX_co2_e_") +
 R("R_EX_co2_e__Rev")  = 0.

Mcoac: R("R_ACACT1r") - R("R_ACACT1r_Rev") - R("R_ACALDi") -
 R("R_ACCOAL") + R("R_ACGS") + R("R_ACOATA") - R("R_ACOATA_Rev") -
 R("R_ACPS1") - R("R_ACS") + R("R_ADHEr") - R("R_ADHEr_Rev") - R("R_AKGDH") +
 R("R_AOXSr") - R("R_AOXSr_Rev") + R("R_AST") - 6.e-6 R("R_BiomassEcoli") +
 R("R_CS") + R("R_DPCOAK") - 7 R("R_FAO1") - 8 R("R_FAO2") - 9
 R("R_FAO3") + R("R_G1PACT") - R("R_GLYATi") + R("R_HSST") + R("R_IPPS") +
 R("R_KAS15") + R("R_MALS") + R("R_MCITS") + R("R_MCOATA") -
 R("R_MCOATA_Rev") + R("R_NPHS") - R("R_OBTFL") - R("R_PACCOAL") - R("R_PDH") -
 R("R_PFL") + R("R_PTA2") + R("R_PTAr") - R("R_PTAr_Rev") + R("R_SERAT") -
 R("R_SERAT_Rev") + R("R_SPMDAT1") + R("R_SPMDAT2") - R("R_SUCBZL") -
 R("R_SUCOAS") + R("R_SUCOAS_Rev") + R("R_TDPADGAT") + R("R_THDPS")  =
 0.

Mcpppg3c: -R("R_CPPPGO") + R("R_UPPDC1")  = 0.

Mcrnb: R("Mcrnb_Transp") + R("R_EX_crn_e_") - R("R_EX_crn_e__Rev") 
 = 0.

Mcrnc: -R("R_CRNBTCT") + R("R_CRNBTCT_Rev") - R("R_CRNCBCT") +
 R("R_CRNCBCT_Rev") + R("R_CRNt7")  = 0.

Mcrncoac: R("R_CRNBTCT") - R("R_CRNBTCT_Rev") + R("R_CRNCBCT") -
 R("R_CRNCBCT_Rev") - R("R_CRNCDH") + R("R_CRNCDH_Rev")  = 0.

Mcrne: -R("R_CRNt7") - R("R_EX_crn_e_") + R("R_EX_crn_e__Rev")  = 0.

Mcsnb: R("Mcsnb_Transp") + R("R_EX_csn_e_") - R("R_EX_csn_e__Rev") 
 = 0.

Mcsnc: R("R_CMPN") - R("R_CSND") + R("R_CSNt2")  = 0.

Mcsne: -R("R_CSNt2") - R("R_EX_csn_e_") + R("R_EX_csn_e__Rev")  = 0.

Mctbtc: R("R_CRNCBCT") - R("R_CRNCBCT_Rev")  = 0.

Mctbtcoac: -R("R_CRNCBCT") + R("R_CRNCBCT_Rev") + R("R_CRNCDH") -
 R("R_CRNCDH_Rev")  = 0.

Mctpc: -0.126 R("R_BiomassEcoli") + R("R_CTPS2") - R("R_DASYN_EC") +
 R("R_DASYN_EC_Rev") - R("R_KDOCT2") - R("R_MEPCT") + R("R_NDPK3") -
 R("R_NDPK3_Rev") - R("R_NTPP4") - R("R_PPNCL2") - R("R_RNTR3")  = 0.

Mcyanc: -R("R_CYANST")  = 0.

Mcyntb: R("Mcyntb_Transp") + R("R_EX_cynt_e_") -
 R("R_EX_cynt_e__Rev")  = 0.

Mcyntc: -R("R_CYNTAH") + R("R_CYNTt2")  = 0.

Mcynte: -R("R_CYNTt2") - R("R_EX_cynt_e_") + R("R_EX_cynt_e__Rev") 
 = 0.

McysLb: R("McysLb_Transp") + R("R_EX_cys_L_e_") -
 R("R_EX_cys_L_e__Rev")  = 0.

McysLc: -0.087 R("R_BiomassEcoli") - R("R_BTS2") + R("R_BTS2_Rev") +
 R("R_CYSabc") - R("R_CYSDS") + R("R_CYSS") - R("R_GLUCYS") -
 R("R_PPNCL2") - R("R_SHSL1") - R("R_THZPSN")  = 0.

McysLe: -R("R_CYSabc") - R("R_EX_cys_L_e_") + R("R_EX_cys_L_e__Rev")
  = 0.

McystLc: -R("R_CYSTL") + R("R_SHSL1")  = 0.

Mcytdb: R("Mcytdb_Transp") + R("R_EX_cytd_e_") -
 R("R_EX_cytd_e__Rev")  = 0.

Mcytdc: -R("R_CYTD") - R("R_CYTDK2") + R("R_CYTDt2") +
 R("R_CYTDt2r") - R("R_CYTDt2r_Rev") + R("R_NTD4")  = 0.

Mcytde: -R("R_CYTDt2") - R("R_CYTDt2r") + R("R_CYTDt2r_Rev") -
 R("R_EX_cytd_e_") + R("R_EX_cytd_e__Rev")  = 0.

Mdad2b: R("Mdad2b_Transp") + R("R_EX_dad_2_e_") -
 R("R_EX_dad_2_e__Rev")  = 0.

Mdad2c: -R("R_DADA") + R("R_DADNt2") + R("R_NTD6") - R("R_PUNP2") +
 R("R_PUNP2_Rev")  = 0.

Mdad2e: -R("R_DADNt2") - R("R_EX_dad_2_e_") + R("R_EX_dad_2_e__Rev")
  = 0.

Mdadpc: R("R_DADK") - R("R_DADK_Rev") - R("R_NDPK8") +
 R("R_NDPK8_Rev") + R("R_RNDR1")  = 0.

Mdampc: -R("R_DADK") + R("R_DADK_Rev") - R("R_NTD6") + R("R_NTPP5") 
 = 0.

Mdannc: R("R_AMAOTr") - R("R_AMAOTr_Rev") - R("R_DBTSr") +
 R("R_DBTSr_Rev")  = 0.

Mdatpc: -0.0247 R("R_BiomassEcoli") + R("R_NDPK8") -
 R("R_NDPK8_Rev") - R("R_NTPP5") + R("R_RNTR1")  = 0.

Mdb4pc: R("R_DB4PS") - R("R_RBFSa")  = 0.

Mdcampc: -R("R_ADSL1r") + R("R_ADSL1r_Rev") + R("R_ADSS")  = 0.

Mdcdpc: R("R_CYTK2") - R("R_CYTK2_Rev") - R("R_NDPK7") +
 R("R_NDPK7_Rev") + R("R_RNDR3")  = 0.

Mdcmpc: -R("R_CYTK2") + R("R_CYTK2_Rev") - R("R_NTD3") +
 R("R_NTPP3")  = 0.

Mdctpc: -0.0254 R("R_BiomassEcoli") - R("R_DCTPD") + R("R_NDPK7") -
 R("R_NDPK7_Rev") - R("R_NTPP3") + R("R_RNTR3")  = 0.

Mdcytb: R("Mdcytb_Transp") + R("R_EX_dcyt_e_") -
 R("R_EX_dcyt_e__Rev")  = 0.

Mdcytc: -R("R_DCYTD") + R("R_DCYTt2") + R("R_NTD3")  = 0.

Mdcyte: -R("R_DCYTt2") - R("R_EX_dcyt_e_") + R("R_EX_dcyt_e__Rev") 
 = 0.

MddcaACPc: R("R_C120SN") - R("R_EDTXS1") - R("R_KAS16")  = 0.

Mdgdpc: R("R_DGK1") - R("R_DGK1_Rev") - R("R_NDPK5") +
 R("R_NDPK5_Rev") + R("R_RNDR2")  = 0.

Mdgmpc: -R("R_DGK1") + R("R_DGK1_Rev") - R("R_NTD8") + R("R_NTPP1") 
 = 0.

Mdgsnb: R("Mdgsnb_Transp") + R("R_EX_dgsn_e_") -
 R("R_EX_dgsn_e__Rev")  = 0.

Mdgsnc: R("R_DGSNt2") + R("R_NTD8") + R("R_NTPTP1") - R("R_PUNP4") +
 R("R_PUNP4_Rev")  = 0.

Mdgsne: -R("R_DGSNt2") - R("R_EX_dgsn_e_") + R("R_EX_dgsn_e__Rev") 
 = 0.

Mdgtpc: -0.0254 R("R_BiomassEcoli") + R("R_NDPK5") -
 R("R_NDPK5_Rev") - R("R_NTPP1") - R("R_NTPTP1") + R("R_RNTR2")  = 0.

Mdhab: R("Mdhab_Transp") + R("R_EX_dha_e_") - R("R_EX_dha_e__Rev") 
 = 0.

Mdhac: -R("R_DHAPT") + R("R_DHAt") - R("R_DHAt_Rev") + R("R_F6PA") -
 R("R_F6PA_Rev") + R("R_GLYCDx")  = 0.

Mdhae: -R("R_DHAt") + R("R_DHAt_Rev") - R("R_EX_dha_e_") +
 R("R_EX_dha_e__Rev")  = 0.

Mdhapc: R("R_DHAPT") + R("R_FBA") - R("R_FBA_Rev") + R("R_FCLPA") -
 R("R_FCLPA_Rev") + R("R_G3PD2") - R("R_G3PD2_Rev") + R("R_G3PD5") +
 R("R_G3PD6") + R("R_G3PD7") - R("R_MGSA") - R("R_QULNS") + R("R_RMPA") -
 R("R_RMPA_Rev") + R("R_TGBPA") - R("R_TGBPA_Rev") - R("R_TPI") +
 R("R_TPI_Rev")  = 0.

Mdhcinnmc: R("R_3HCINNMH") + R("R_DHCIND") - R("R_DHCINDO")  = 0.

Mdhfc: -R("R_DHFR") + R("R_DHFR_Rev") + R("R_DHFS") + R("R_TMDS")  =
 0.

Mdhnac: -R("R_DHNAOT") + R("R_NPHS")  = 0.

Mdhnptc: -R("R_DHNPA2") + R("R_DNMPPA")  = 0.

MdhorSc: -R("R_DHORD2") - R("R_DHORD5") - R("R_DHORTS") +
 R("R_DHORTS_Rev")  = 0.

Mdhpmpc: -R("R_DNMPPA") + R("R_DNTPPA")  = 0.

Mdhpppnc: R("R_3HPPPNH") + R("R_DHPPD") - R("R_HPPPNDO")  = 0.

Mdhptc: -R("R_DHFS") + R("R_DHPS2")  = 0.

Mdhptdc: -R("R_DHPTDCs") + R("R_RHCCE")  = 0.

Mdinb: R("Mdinb_Transp") + R("R_EX_din_e_") - R("R_EX_din_e__Rev") 
 = 0.

Mdinc: R("R_DADA") + R("R_DINSt2") - R("R_PUNP6") + R("R_PUNP6_Rev")
  = 0.

Mdine: -R("R_DINSt2") - R("R_EX_din_e_") + R("R_EX_din_e__Rev")  =
 0.

Mdkmppc: -R("R_DKMPPD") - R("R_DKMPPD2") + R("R_MDRPD")  = 0.

Mdmbzidc: -R("R_NNDMBRT")  = 0.

Mdmlzc: R("R_RBFSa") - 2 R("R_RBFSb")  = 0.

Mdmppc: -R("R_DMATT") + R("R_DMPPS") + R("R_IPDDIi")  = 0.

Mdmsb: R("Mdmsb_Transp") + R("R_EX_dms_e_") - R("R_EX_dms_e__Rev") 
 = 0.

Mdmsc: R("R_DMSOR1") + R("R_DMSOR2")  = 0.

Mdmse: R("R_DMSOR1e") + R("R_DMSOR2e") - R("R_EX_dms_e_") +
 R("R_EX_dms_e__Rev")  = 0.

Mdmsob: R("Mdmsob_Transp") + R("R_EX_dmso_e_") -
 R("R_EX_dmso_e__Rev")  = 0.

Mdmsoc: -R("R_DMSOR1") - R("R_DMSOR2")  = 0.

Mdmsoe: -R("R_DMSOR1e") - R("R_DMSOR2e") - R("R_EX_dmso_e_") +
 R("R_EX_dmso_e__Rev")  = 0.

Mdnadc: -R("R_NADS1") + R("R_NNAT")  = 0.

Mdpcoac: -R("R_DPCOAK") + R("R_PTPATi")  = 0.

Mdsclc: -R("R_SHCHD2") + R("R_UPP3MT")  = 0.

Mdtbtc: -R("R_BTS2") + R("R_BTS2_Rev") + R("R_DBTSr") -
 R("R_DBTSr_Rev")  = 0.

Mdtdp4aaddgc: -R("R_AADDGT") + R("R_TDPADGAT")  = 0.

Mdtdp4addgc: -R("R_TDPADGAT") + R("R_TDPAGTA")  = 0.

Mdtdp4d6dgc: -R("R_TDPAGTA") - R("R_TDPDRE") + R("R_TDPGDH")  = 0.

Mdtdp4d6dmc: R("R_TDPDRE") - R("R_TDPDRR")  = 0.

Mdtdpc: R("R_AADDGT") + R("R_DTMPK") - R("R_DTMPK_Rev") -
 R("R_NDPK4") + R("R_NDPK4_Rev")  = 0.

Mdtdpgluc: R("R_G1PTT") - R("R_TDPGDH")  = 0.

Mdtdprmnc: R("R_TDPDRR")  = 0.

Mdtmpc: -R("R_DTMPK") + R("R_DTMPK_Rev") - R("R_NTD5") +
 R("R_NTPP7") + R("R_TMDK1") + R("R_TMDS")  = 0.

Mdttpc: -0.0247 R("R_BiomassEcoli") - R("R_G1PTT") + R("R_NDPK4") -
 R("R_NDPK4_Rev") - R("R_NTPP7")  = 0.

Mdudpc: -R("R_NDPK6") + R("R_NDPK6_Rev") + R("R_RNDR4") +
 R("R_URIDK2r") - R("R_URIDK2r_Rev")  = 0.

Mdumpc: R("R_DURIK1") + R("R_DUTPDP") - R("R_NTD1") - R("R_TMDS") -
 R("R_URIDK2r") + R("R_URIDK2r_Rev")  = 0.

Mdurib: R("Mdurib_Transp") + R("R_EX_duri_e_") -
 R("R_EX_duri_e__Rev")  = 0.

Mduric: R("R_DCYTD") - R("R_DURIK1") - R("R_DURIPP") +
 R("R_DURIPP_Rev") + R("R_DURIt2") + R("R_NTD1")  = 0.

Mdurie: -R("R_DURIt2") - R("R_EX_duri_e_") + R("R_EX_duri_e__Rev") 
 = 0.

Mdutpc: R("R_DCTPD") - R("R_DUTPDP") + R("R_NDPK6") -
 R("R_NDPK6_Rev") + R("R_RNTR4")  = 0.

Mdxyl5pc: -R("R_DXPRIi") + R("R_DXPS") + R("R_DXYLK") -
 R("R_PDX5PS") - R("R_THZPSN")  = 0.

Mdxylc: -R("R_DXYLK")  = 0.

Me4pc: -R("R_DDPA") - R("R_E4PD") + R("R_E4PD_Rev") + R("R_TALA") -
 R("R_TALA_Rev") - R("R_TKT2") + R("R_TKT2_Rev")  = 0.

MecaECc: R("R_ECAP_EC")  = 0.

Meig3pc: R("R_IG3PS") - R("R_IGPDH")  = 0.

Menterc: R("R_ENTCS")  = 0.

Methac: -R("R_ETHAAL") + R("R_GPDDA2")  = 0.

Metohb: R("Metohb_Transp") + R("R_EX_etoh_e_") -
 R("R_EX_etoh_e__Rev")  = 0.

Metohc: R("R_ADHEr") - R("R_ADHEr_Rev") + R("R_ETOHt2r") -
 R("R_ETOHt2r_Rev")  = 0.

Metohe: -R("R_ETOHt2r") + R("R_ETOHt2r_Rev") - R("R_EX_etoh_e_") +
 R("R_EX_etoh_e__Rev")  = 0.

Mf1pc: -R("R_FRUK") + R("R_FRUpts")  = 0.

Mf6pc: -R("R_F6PA") + R("R_F6PA_Rev") + R("R_FBP") + R("R_FRUpts2")
 + R("R_G6PDA") - R("R_GF6PTA") + R("R_M1PD") - R("R_M1PD_Rev") +
 R("R_MAN6PI") - R("R_MAN6PI_Rev") - R("R_PFK") + R("R_PGI") -
 R("R_PGI_Rev") + R("R_SBTPD") - R("R_SBTPD_Rev") + R("R_TALA") - R("R_TALA_Rev") +
 R("R_TKT2") - R("R_TKT2_Rev")  = 0.

Mfadc: -0.00001 R("R_BiomassEcoli") - R("R_DAAD") - 6 R("R_FAO1") -
 7 R("R_FAO2") - 8 R("R_FAO3") - R("R_FAO4") + R("R_FMNAT") -
 R("R_PROD2") - R("R_SUCD1i") + R("R_SUCD4") - R("R_SUCD4_Rev")  = 0.

Mfadh2c: R("R_DAAD") + 6 R("R_FAO1") + 7 R("R_FAO2") + 8 R("R_FAO3")
 + R("R_FAO4") + R("R_PROD2") + R("R_SUCD1i") - R("R_SUCD4") +
 R("R_SUCD4_Rev")  = 0.

Mfc1pc: R("R_FCLK") - R("R_FCLPA") + R("R_FCLPA_Rev")  = 0.

MfclLc: R("R_FCI") - R("R_FCI_Rev") - R("R_FCLK")  = 0.

Mfdpc: -R("R_FBA") + R("R_FBA_Rev") - R("R_FBP") + R("R_FRUK") +
 R("R_PFK")  = 0.

Mfe2b: R("Mfe2b_Transp") + R("R_EX_fe2_e_") - R("R_EX_fe2_e__Rev") 
 = 0.

Mfe2c: -R("R_FCLT") + R("R_FE2abc") - R("R_SHCHF")  = 0.

Mfe2e: -R("R_EX_fe2_e_") + R("R_EX_fe2_e__Rev") - R("R_FE2abc")  =
 0.

Mfgamc: R("R_GARFT") - R("R_GARFT_Rev") + R("R_GART") - R("R_PRFGS")
  = 0.

Mfmnc: -R("R_FMNAT") + R("R_RBFK")  = 0.

Mforb: R("Mforb_Transp") + R("R_EX_for_e_") - R("R_EX_for_e__Rev") 
 = 0.

Mforc: 2 R("R_AMPMS") + R("R_DB4PS") + R("R_DKMPPD") +
 R("R_DKMPPD2") - R("R_FDH2") - R("R_FDH3") - R("R_FHL") + R("R_FORt") -
 R("R_FORt_Rev") + R("R_FTHFD") - R("R_GART") + R("R_GTPCI") + R("R_GTPCII2") +
 R("R_OBTFL") + R("R_PFL")  = 0.

Mfore: -R("R_EX_for_e_") + R("R_EX_for_e__Rev") - R("R_FORt") +
 R("R_FORt_Rev")  = 0.

Mfpramc: -R("R_PRAIS") + R("R_PRFGS")  = 0.

Mfpricac: R("R_AICART") - R("R_AICART_Rev") + R("R_IMPC") -
 R("R_IMPC_Rev")  = 0.

Mfrdpc: R("R_GRTT") - R("R_HEMEOS") - R("R_OCTDPS") - R("R_UDCPDPS")
  = 0.

Mfrub: R("Mfrub_Transp") + R("R_EX_fru_e_") - R("R_EX_fru_e__Rev") 
 = 0.

Mfruc: R("R_FFSD") - R("R_XYLI2i")  = 0.

Mfrue: -R("R_EX_fru_e_") + R("R_EX_fru_e__Rev") - R("R_FRUpts") -
 R("R_FRUpts2")  = 0.

Mfruurc: R("R_GUI1") - R("R_GUI1_Rev") + R("R_MANAO") -
 R("R_MANAO_Rev")  = 0.

Mfuc1pLb: R("Mfuc1pLb_Transp") + R("R_EX_fuc1p_L_e_") -
 R("R_EX_fuc1p_L_e__Rev")  = 0.

Mfuc1pLc: R("R_FUCPt6_2")  = 0.

Mfuc1pLe: -R("R_EX_fuc1p_L_e_") + R("R_EX_fuc1p_L_e__Rev") -
 R("R_FUCPt6_2")  = 0.

MfucLb: R("MfucLb_Transp") + R("R_EX_fuc_L_e_") -
 R("R_EX_fuc_L_e__Rev")  = 0.

MfucLc: -R("R_FCI") + R("R_FCI_Rev") + R("R_FUCt") - R("R_FUCt_Rev")
  = 0.

MfucLe: -R("R_EX_fuc_L_e_") + R("R_EX_fuc_L_e__Rev") - R("R_FUCt") +
 R("R_FUCt_Rev")  = 0.

Mfumb: R("Mfumb_Transp") + R("R_EX_fum_e_") - R("R_EX_fum_e__Rev") 
 = 0.

Mfumc: R("R_ADSL1r") - R("R_ADSL1r_Rev") + R("R_ADSL2r") -
 R("R_ADSL2r_Rev") + R("R_ARGSL") - R("R_ARGSL_Rev") - R("R_ASPO5") +
 R("R_ASPT") - R("R_FRD2") - R("R_FRD3") - R("R_FUM") + R("R_FUM_Rev") +
 R("R_FUMt2_2") + R("R_FUMt2_3") + R("R_HKNTDH") + R("R_SUCD1i") +
 R("R_SUCFUMt") - R("R_SUCFUMt_Rev")  = 0.

Mfume: -R("R_EX_fum_e_") + R("R_EX_fum_e__Rev") - R("R_FUMt2_2") -
 R("R_FUMt2_3") - R("R_SUCFUMt") + R("R_SUCFUMt_Rev")  = 0.

Mg1pc: -R("R_G1PP") - R("R_G1PTT") - R("R_GALU") - R("R_GALUi") +
 R("R_GALU_Rev") + R("R_GLCP") - R("R_GLGC") + R("R_MLTP1") -
 R("R_MLTP1_Rev") + R("R_MLTP2") - R("R_MLTP2_Rev") + R("R_MLTP3") -
 R("R_MLTP3_Rev") - R("R_PGMT") + R("R_PGMT_Rev") + R("R_UGLT") - R("R_UGLT_Rev")  =
 0.

Mg3pc: R("R_DDPGALA") - R("R_DDPGALA_Rev") + R("R_DRPA") -
 R("R_DXPS") + R("R_EDA") + R("R_F6PA") - R("R_F6PA_Rev") + R("R_FBA") -
 R("R_FBA_Rev") - R("R_GAPD") + R("R_GAPD_Rev") - R("R_TALA") +
 R("R_TALA_Rev") + R("R_TGBPA") - R("R_TGBPA_Rev") + R("R_TKT1") - R("R_TKT1_Rev") +
 R("R_TKT2") - R("R_TKT2_Rev") + R("R_TPI") - R("R_TPI_Rev") +
 R("R_TRPS1") + R("R_TRPS3")  = 0.

Mg3pcc: -R("R_GPDDA1") + R("R_LPLIPA3") + R("R_LPLIPA5")  = 0.

Mg3pec: -R("R_GPDDA2") + R("R_LPLIPA2") + R("R_LPLIPA4")  = 0.

Mg3pgc: -R("R_GPDDA4") + R("R_LPLIPA1") + R("R_LPLIPA6")  = 0.

Mg3pic: -R("R_GPDDA5")  = 0.

Mg3psc: -R("R_GPDDA3")  = 0.

Mg6pb: R("Mg6pb_Transp") + R("R_EX_g6p_e_") - R("R_EX_g6p_e__Rev") 
 = 0.

Mg6pc: R("R_AB6PGH") + R("R_FFSD") - R("R_G6PDH2r") +
 R("R_G6PDH2r_Rev") + R("R_G6Pt6_2") + R("R_GLCpts") + R("R_HEX1") - R("R_PGI") +
 R("R_PGI_Rev") + R("R_PGMT") - R("R_PGMT_Rev") + R("R_TRE6PH") -
 R("R_TRE6PS")  = 0.

Mg6pe: -R("R_EX_g6p_e_") + R("R_EX_g6p_e__Rev") - R("R_G6Pt6_2")  =
 0.

Mgal1pc: R("R_GALKr") - R("R_GALKr_Rev") - R("R_UGLT") +
 R("R_UGLT_Rev")  = 0.

Mgalb: R("Mgalb_Transp") + R("R_EX_gal_e_") - R("R_EX_gal_e__Rev") 
 = 0.

Mgalc: R("R_GALabc") - R("R_GALKr") + R("R_GALKr_Rev") +
 R("R_GALS3") + R("R_GALt2") + R("R_LACZ")  = 0.

MgalctDb: R("MgalctDb_Transp") + R("R_EX_galct_D_e_") -
 R("R_EX_galct_D_e__Rev")  = 0.

MgalctDc: -R("R_GALCTD") + R("R_GALCTt2r") - R("R_GALCTt2r_Rev")  =
 0.

MgalctDe: -R("R_EX_galct_D_e_") + R("R_EX_galct_D_e__Rev") -
 R("R_GALCTt2r") + R("R_GALCTt2r_Rev")  = 0.

MgalctnDb: R("MgalctnDb_Transp") + R("R_EX_galctn_D_e_") -
 R("R_EX_galctn_D_e__Rev")  = 0.

MgalctnDc: -R("R_GALCTND") + R("R_GALCTNt2r") - R("R_GALCTNt2r_Rev")
  = 0.

MgalctnDe: -R("R_EX_galctn_D_e_") + R("R_EX_galctn_D_e__Rev") -
 R("R_GALCTNt2r") + R("R_GALCTNt2r_Rev")  = 0.

Mgale: -R("R_EX_gal_e_") + R("R_EX_gal_e__Rev") - R("R_GALabc") -
 R("R_GALt2")  = 0.

Mgalt1pc: R("R_GALTpts") - R("R_GLTPD") + R("R_GLTPD_Rev")  = 0.

Mgaltb: R("Mgaltb_Transp") + R("R_EX_galt_e_") -
 R("R_EX_galt_e__Rev")  = 0.

Mgalte: -R("R_EX_galt_e_") + R("R_EX_galt_e__Rev") - R("R_GALTpts") 
 = 0.

Mgalurb: R("Mgalurb_Transp") + R("R_EX_galur_e_") -
 R("R_EX_galur_e__Rev")  = 0.

Mgalurc: R("R_GALURt2r") - R("R_GALURt2r_Rev") - R("R_GUI2") +
 R("R_GUI2_Rev")  = 0.

Mgalure: -R("R_EX_galur_e_") + R("R_EX_galur_e__Rev") -
 R("R_GALURt2r") + R("R_GALURt2r_Rev")  = 0.

Mgam1pc: -R("R_G1PACT") - R("R_PGAMT") + R("R_PGAMT_Rev")  = 0.

Mgam6pc: R("R_AGDC") - R("R_G6PDA") + R("R_GAMpts") + R("R_GF6PTA")
 + R("R_PGAMT") - R("R_PGAMT_Rev")  = 0.

Mgamb: R("Mgamb_Transp") + R("R_EX_gam_e_") - R("R_EX_gam_e__Rev") 
 = 0.

Mgame: -R("R_EX_gam_e_") + R("R_EX_gam_e__Rev") - R("R_GAMpts")  =
 0.

Mgarc: -R("R_GARFT") + R("R_GARFT_Rev") - R("R_GART") +
 R("R_PRAGSr") - R("R_PRAGSr_Rev")  = 0.

Mgbbtnb: R("Mgbbtnb_Transp") + R("R_EX_gbbtn_e_") -
 R("R_EX_gbbtn_e__Rev")  = 0.

Mgbbtnc: R("R_CRNBTCT") - R("R_CRNBTCT_Rev") - R("R_CRNt7")  = 0.

Mgbbtne: R("R_CRNt7") - R("R_EX_gbbtn_e_") + R("R_EX_gbbtn_e__Rev") 
 = 0.

Mgcaldc: R("R_DHNPA2") - R("R_GCALDD")  = 0.

Mgdpc: R("R_ADK3") - R("R_ADK3_Rev") + R("R_ADSS") + R("R_CYTDK2") +
 R("R_GK1") - R("R_GK1_Rev") + 2 R("R_GP4GH") - R("R_MAN1PT2") -
 R("R_NDPK1") + R("R_NDPK1_Rev") - R("R_RNDR2") + R("R_SADT2") +
 R("R_URIK2")  = 0.

Mgdpddmanc: -R("R_GDMANE") + R("R_GMAND")  = 0.

Mgdpfucc: R("R_GOFUCR")  = 0.

Mgdpmannc: -R("R_GMAND") + R("R_MAN1PT2")  = 0.

Mgdpofucc: R("R_GDMANE") - R("R_GOFUCR")  = 0.

MglcDb: R("MglcDb_Transp") + R("R_EX_glc_e_") - R("R_EX_glc_e__Rev")
  = 0.

MglcDc: R("R_AMALT1") + R("R_AMALT2") + R("R_AMALT3") +
 R("R_AMALT4") + R("R_G1PP") + R("R_GALS3") + R("R_GLCt2") - R("R_HEX1") +
 R("R_LACZ") + R("R_MLTG1") + R("R_MLTG2") + R("R_MLTG3") + R("R_MLTG4") +
 R("R_MLTG5") + R("R_TRE6PH") + 2 R("R_TREH") + R("R_XYLI2i")  = 0.

MglcDe: -R("R_EX_glc_e_") + R("R_EX_glc_e__Rev") - R("R_GLCDe") -
 R("R_GLCpts") - R("R_GLCt2") + 2 R("R_TREHe")  = 0.

Mglcnb: R("Mglcnb_Transp") + R("R_EX_glcn_e_") -
 R("R_EX_glcn_e__Rev")  = 0.

Mglcnc: R("R_2DGLCNRx") + R("R_2DGLCNRy") + R("R_5DGLCNR") -
 R("R_5DGLCNR_Rev") + R("R_GLCNt2r") - R("R_GLCNt2r_Rev") - R("R_GNK")  = 0.

Mglcne: -R("R_EX_glcn_e_") + R("R_EX_glcn_e__Rev") + R("R_GLCDe") -
 R("R_GLCNt2r") + R("R_GLCNt2r_Rev")  = 0.

Mglcrb: R("Mglcrb_Transp") + R("R_EX_glcr_e_") -
 R("R_EX_glcr_e__Rev")  = 0.

Mglcrc: -R("R_GLCRD") + R("R_GLCRt2r") - R("R_GLCRt2r_Rev")  = 0.

Mglcre: -R("R_EX_glcr_e_") + R("R_EX_glcr_e__Rev") - R("R_GLCRt2r")
 + R("R_GLCRt2r_Rev")  = 0.

Mglcurb: R("Mglcurb_Transp") + R("R_EX_glcur_e_") -
 R("R_EX_glcur_e__Rev")  = 0.

Mglcurc: R("R_GLCURt2r") - R("R_GLCURt2r_Rev") - R("R_GUI1") +
 R("R_GUI1_Rev")  = 0.

Mglcure: -R("R_EX_glcur_e_") + R("R_EX_glcur_e__Rev") -
 R("R_GLCURt2r") + R("R_GLCURt2r_Rev")  = 0.

MglnLb: R("MglnLb_Transp") + R("R_EX_gln_L_e_") -
 R("R_EX_gln_L_e__Rev")  = 0.

MglnLc: -R("R_ADCS") - R("R_ANS") - R("R_ASNS1") - 0.25
 R("R_BiomassEcoli") - R("R_CBPS") - R("R_CTPS2") - R("R_GF6PTA") + R("R_GLNabc") +
 R("R_GLNS") - R("R_GLUN") - R("R_GLUPRT") - R("R_GLUSy") -
 R("R_GMPS2") - R("R_IG3PS") - R("R_PRFGS")  = 0.

MglnLe: -R("R_EX_gln_L_e_") + R("R_EX_gln_L_e__Rev") - R("R_GLNabc")
  = 0.

Mglu1sac: -R("R_G1SATi") + R("R_GLUTRR")  = 0.

Mglu5pc: -R("R_G5SD") + R("R_GLU5K")  = 0.

Mglu5sac: -R("R_G5SADs") + R("R_G5SD") + R("R_NACODA") +
 R("R_ORNTA")  = 0.

Mglucysc: R("R_GLUCYS") - R("R_GTHS")  = 0.

MgluDc: -R("R_GLUR") + R("R_GLUR_Rev") - R("R_UAMAGS")  = 0.

MgluLb: R("MgluLb_Transp") + R("R_EX_glu_L_e_") -
 R("R_EX_glu_L_e__Rev")  = 0.

MgluLc: R("R_ABTA") - R("R_ACGS") + R("R_ACOTA") - R("R_ACOTA_Rev")
 + R("R_ADCS") + R("R_ALATA_L") - R("R_ALATA_L_Rev") + R("R_ANS") +
 R("R_ASNS1") + R("R_ASPTA") - R("R_ASPTA_Rev") - 0.25 R("R_BiomassEcoli")
 + R("R_CBPS") + R("R_CTPS2") - R("R_DHFS") + R("R_GF6PTA") -
 R("R_GLNS") - R("R_GLU5K") + R("R_GLUabc") + R("R_GLUABUTt7") -
 R("R_GLUABUTt7_Rev") - R("R_GLUCYS") - R("R_GLUDC") - R("R_GLUDy") +
 R("R_GLUDy_Rev") + R("R_GLUN") + R("R_GLUPRT") + R("R_GLUR") - R("R_GLUR_Rev") + 2
 R("R_GLUSy") + R("R_GLUt2r") - R("R_GLUt2r_Rev") + R("R_GLUt4") -
 R("R_GLUTRS") + R("R_GMPS2") - R("R_HSTPT") + R("R_IG3PS") + R("R_ILETA") -
 R("R_ILETA_Rev") - R("R_LEUTAi") - R("R_OHPBAT") + R("R_OHPBAT_Rev") +
 R("R_ORNTA") + R("R_P5CD") + R("R_PHETA1") - R("R_PHETA1_Rev") +
 R("R_PRFGS") - R("R_PSERT") + R("R_PTRCTA") + R("R_SDPTA") -
 R("R_SDPTA_Rev") + R("R_SGDS") + R("R_SOTA") - R("R_TDPAGTA") + R("R_TYRTA") -
 R("R_TYRTA_Rev") - R("R_UNK3") + R("R_VALTA") - R("R_VALTA_Rev")  = 0.

MgluLe: -R("R_EX_glu_L_e_") + R("R_EX_glu_L_e__Rev") - R("R_GLUabc")
 - R("R_GLUABUTt7") + R("R_GLUABUTt7_Rev") - R("R_GLUt2r") +
 R("R_GLUt2r_Rev") - R("R_GLUt4")  = 0.

Mglutrnac: -R("R_GLUTRR") + R("R_GLUTRS")  = 0.

Mglxc: -2 R("R_GLXCL") - R("R_GLYCLTDx") - R("R_GLYCLTDy") +
 R("R_GLYCTO2") + R("R_GLYCTO3") + R("R_GLYCTO4") + R("R_ICL") - R("R_MALS") +
 R("R_UGLYCH")  = 0.

Mglyaldb: R("Mglyaldb_Transp") + R("R_EX_glyald_e_") -
 R("R_EX_glyald_e__Rev")  = 0.

Mglyaldc: -R("R_ALCD19") + R("R_ALCD19_Rev") + R("R_GLYALDt") -
 R("R_GLYALDt_Rev")  = 0.

Mglyalde: -R("R_EX_glyald_e_") + R("R_EX_glyald_e__Rev") -
 R("R_GLYALDt") + R("R_GLYALDt_Rev")  = 0.

Mglyb: R("Mglyb_Transp") + R("R_EX_gly_e_") - R("R_EX_gly_e__Rev") 
 = 0.

Mglybb: R("Mglybb_Transp") + R("R_EX_glyb_e_") -
 R("R_EX_glyb_e__Rev")  = 0.

Mglybc: R("R_BETALDHx") + R("R_BETALDHy") + R("R_GLYBabc") +
 R("R_GLYBt2r") - R("R_GLYBt2r_Rev")  = 0.

Mglybe: -R("R_EX_glyb_e_") + R("R_EX_glyb_e__Rev") - R("R_GLYBabc")
 - R("R_GLYBt2r") + R("R_GLYBt2r_Rev")  = 0.

Mglyc: -0.582 R("R_BiomassEcoli") + R("R_GHMT2") + R("R_GLYATi") -
 R("R_GLYCL") + R("R_GLYt2r") - R("R_GLYt2r_Rev") - R("R_GTHS") -
 R("R_PRAGSr") + R("R_PRAGSr_Rev") + R("R_THRAr") - R("R_THRAr_Rev")  = 0.

Mglyc3pb: R("Mglyc3pb_Transp") + R("R_EX_glyc3p_e_") -
 R("R_EX_glyc3p_e__Rev")  = 0.

Mglyc3pc: -R("R_G3PD2") + R("R_G3PD2_Rev") - R("R_G3PD5") -
 R("R_G3PD6") - R("R_G3PD7") + R("R_GLYC3Pabc") + R("R_GLYC3Pt6") + R("R_GLYK")
 + R("R_GPDDA1") + R("R_GPDDA2") + R("R_GPDDA3") + R("R_GPDDA4") +
 R("R_GPDDA5") - R("R_PASYN_EC") - R("R_PGSA_EC") + R("R_PGSA_EC_Rev")  =
 0.

Mglyc3pe: -R("R_EX_glyc3p_e_") + R("R_EX_glyc3p_e__Rev") -
 R("R_GLYC3Pabc") - R("R_GLYC3Pt6")  = 0.

Mglycb: R("Mglycb_Transp") + R("R_EX_glyc_e_") -
 R("R_EX_glyc_e__Rev")  = 0.

Mglycc: R("R_ALCD19") - R("R_ALCD19_Rev") + R("R_CLPNS_EC") -
 R("R_CLPNS_EC_Rev") - R("R_GLYCDx") - R("R_GLYCt") + R("R_GLYCt_Rev") -
 R("R_GLYK") + R("R_GPDDA4")  = 0.

Mglyce: -R("R_EX_glyc_e_") + R("R_EX_glyc_e__Rev") + R("R_GLYCt") -
 R("R_GLYCt_Rev")  = 0.

Mglycltb: R("Mglycltb_Transp") + R("R_EX_glyclt_e_") -
 R("R_EX_glyclt_e__Rev")  = 0.

Mglycltc: R("R_GCALDD") + R("R_GLYCLTDx") + R("R_GLYCLTDy") +
 R("R_GLYCLTt2r") - R("R_GLYCLTt2r_Rev") - R("R_GLYCTO2") - R("R_GLYCTO3") -
 R("R_GLYCTO4") + R("R_PGLYCP")  = 0.

Mglyclte: -R("R_EX_glyclt_e_") + R("R_EX_glyclt_e__Rev") -
 R("R_GLYCLTt2r") + R("R_GLYCLTt2r_Rev")  = 0.

Mglycogenc: -0.154 R("R_BiomassEcoli") - R("R_GLCP") + R("R_GLCS1") 
 = 0.

MglycRc: -R("R_GLYCK") + R("R_HPYRRx") + R("R_HPYRRy") +
 R("R_TRSAR")  = 0.

Mglye: -R("R_EX_gly_e_") + R("R_EX_gly_e__Rev") - R("R_GLYt2r") +
 R("R_GLYt2r_Rev")  = 0.

Mgmhep17bpc: R("R_GMHEPK") - R("R_GMHEPPA")  = 0.

Mgmhep1pc: -R("R_GMHEPAT") + R("R_GMHEPPA")  = 0.

Mgmhep7pc: -R("R_GMHEPK") + R("R_S7PI")  = 0.

Mgmpc: R("R_ADOCBLS") - R("R_GK1") + R("R_GK1_Rev") - R("R_GMPR") +
 R("R_GMPS2") + R("R_GSNK") + R("R_GUAPRT") - R("R_NTD9") +
 R("R_NTPP2")  = 0.

Mgp4gc: -R("R_GP4GH")  = 0.

Mgrdpc: R("R_DMATT") - R("R_GRTT")  = 0.

Mgsnb: R("Mgsnb_Transp") + R("R_EX_gsn_e_") - R("R_EX_gsn_e__Rev") 
 = 0.

Mgsnc: -R("R_GSNK") + R("R_GSNt2") + R("R_NTD9") + R("R_NTPTP2") -
 R("R_PUNP3") + R("R_PUNP3_Rev")  = 0.

Mgsne: -R("R_EX_gsn_e_") + R("R_EX_gsn_e__Rev") - R("R_GSNt2")  = 0.

Mgthoxc: -R("R_GTHOr") + R("R_GTHOr_Rev")  = 0.

Mgthrdc: R("R_GLYOX") + R("R_GSPMDA") - R("R_GSPMDS") + 2
 R("R_GTHOr") - 2 R("R_GTHOr_Rev") + R("R_GTHS") - R("R_LGTHL")  = 0.

Mgtpc: -R("R_ACBIPGT") - R("R_ADK3") + R("R_ADK3_Rev") - R("R_ADSS")
 - 0.203 R("R_BiomassEcoli") - R("R_CYTDK2") - R("R_GTPCI") -
 R("R_GTPCII2") + R("R_NDPK1") - R("R_NDPK1_Rev") - R("R_NTPP2") -
 R("R_NTPTP2") - R("R_RNTR2") - R("R_SADT2") - R("R_URIK2")  = 0.

Mgtspmdc: -R("R_GSPMDA") + R("R_GSPMDS")  = 0.

Mguab: R("Mguab_Transp") + R("R_EX_gua_e_") - R("R_EX_gua_e__Rev") 
 = 0.

Mguac: -R("R_GUAD") - R("R_GUAPRT") + R("R_GUAt") + R("R_GUAt2") -
 R("R_GUAt_Rev") + R("R_PUNP3") - R("R_PUNP3_Rev") + R("R_PUNP4") -
 R("R_PUNP4_Rev")  = 0.

Mguae: -R("R_EX_gua_e_") + R("R_EX_gua_e__Rev") - R("R_GUAt") -
 R("R_GUAt2") + R("R_GUAt_Rev")  = 0.

Mh2c: R("R_FHL") - R("R_HYD1") - R("R_HYD2") - R("R_HYD3")  = 0.

Mh2mb4pc: -R("R_DMPPS") - R("R_IPDPS") + R("R_MECDPDH")  = 0.

Mh2o2c: R("R_ASPO6") - 2 R("R_CAT") + R("R_PDX5PO") -
 R("R_PDX5PO_Rev") + R("R_PEAMNO") + R("R_PYAM5PO") + R("R_SPODM")  = 0.

Mh2ob: R("Mh2ob_Transp") + R("R_EX_h2o_e_") - R("R_EX_h2o_e__Rev") 
 = 0.

Mh2oc: R("R_3HCINNMH") + R("R_3HPPPNH") - R("R_4HTHRS") -
 R("R_AB6PGH") - R("R_ABUTD") - R("R_ACODA") - R("R_ADA") - R("R_ADD") -
 R("R_AGDC") - R("R_AGMT") - R("R_AHC") + R("R_AHC_Rev") - R("R_AHCYSNS") -
 R("R_ALAabc") - R("R_ALDD19x") - R("R_ALDD2x") - R("R_ALLTAH") -
 R("R_ALLTN") + R("R_ALTRH") - R("R_AMPMS") - R("R_AMPN") - R("R_AP4AH") -
 R("R_AP5AH") - R("R_ARBabc") - R("R_ARGabc") - R("R_ASNabc") - R("R_ASNN")
 - R("R_ASNS1") - R("R_ASPabc") - R("R_ATPM") + R("R_ATPS4r") -
 R("R_ATPS4r_Rev") - R("R_BETALDHx") - R("R_BETALDHy") - 45.5608
 R("R_BiomassEcoli") - R("R_BPNT") + R("R_BSORx") + R("R_BSORy") + 5 R("R_C120SN")
 + 6 R("R_C140SN") + 6 R("R_C141SN") + 7 R("R_C160SN") + 7
 R("R_C161SN") + 8 R("R_C181SN") + 2 R("R_CAT") - R("R_CBIAT_DELETE") +
 R("R_CBIAT_DELETE_Rev") - R("R_CBL1abc") - R("R_CBLAT_DELETE") +
 R("R_CBLAT_DELETE_Rev") - R("R_CBPS") - R("R_CDAPPA_EC") - R("R_CHLabc") - R("R_CMPN")
 + 2 R("R_CPPPGO") + R("R_CRNCDH") - R("R_CRNCDH_Rev") - R("R_CS") -
 R("R_CSND") - R("R_CTPS2") - R("R_CYSabc") - R("R_CYSDS") -
 R("R_CYSTL") + R("R_CYTBD") + R("R_CYTBO3") - R("R_CYTD") - R("R_DAAD") -
 R("R_DADA") - R("R_DAPabc") - R("R_DCTPD") - R("R_DCYTD") - R("R_DDPA") +
 R("R_DHAD1") + R("R_DHAD2") + 2 R("R_DHDPS") - R("R_DHORTS") +
 R("R_DHORTS_Rev") - R("R_DHPPDA2") + R("R_DHPTDCs") + R("R_DHQD") -
 R("R_DHQD_Rev") - R("R_DKMPPD") - 3 R("R_DKMPPD2") + R("R_DMPPS") + R("R_DMSOR1")
 + R("R_DMSOR2") - R("R_DNMPPA") - R("R_DNTPPA") - R("R_DUTPDP") -
 R("R_E4PD") + R("R_E4PD_Rev") + R("R_EDD") + R("R_ENO") - R("R_ENO_Rev") -
 6 R("R_FAO1") - 7 R("R_FAO2") - 8 R("R_FAO3") - R("R_FAO4") -
 R("R_FBP") - R("R_FE2abc") - R("R_FFSD") - R("R_FTHFD") - R("R_FUM") +
 R("R_FUM_Rev") - R("R_G1PP") + R("R_G5SADs") - R("R_G6PDA") - R("R_GALabc")
 + R("R_GALCTD") + R("R_GALCTND") - R("R_GALS3") - R("R_GCALDD") +
 R("R_GHMT2") + R("R_GLCRD") - R("R_GLNabc") - R("R_GLUabc") - R("R_GLUDy")
 + R("R_GLUDy_Rev") - R("R_GLUN") - R("R_GLUPRT") - R("R_GLYBabc") -
 R("R_GLYC3Pabc") - R("R_GLYOX") + R("R_GMAND") - R("R_GMHEPPA") -
 R("R_GMPS2") - R("R_GP4GH") - R("R_GPDDA1") - R("R_GPDDA2") - R("R_GPDDA3")
 - R("R_GPDDA4") - R("R_GPDDA5") - R("R_GSPMDA") - R("R_GTPCI") - 3
 R("R_GTPCII2") - R("R_GUAD") + R("R_H2Ot") - R("R_H2Ot_Rev") -
 R("R_HCO3E") + R("R_HCO3E_Rev") - R("R_HEMEOS") - R("R_HISabc") - R("R_HISTD")
 - R("R_HISTP") - R("R_HKNDDH") - R("R_HKNTDH") - R("R_HMBS") -
 R("R_HYPOE") - R("R_ICHORT") + R("R_IGPDH") + R("R_IGPS") - R("R_ILEabc") -
 R("R_IMPC") + R("R_IMPC_Rev") - R("R_IMPD") + R("R_IPDPS") +
 R("R_IPPMIa") - R("R_IPPMIa_Rev") - R("R_IPPMIb") + R("R_IPPMIb_Rev") -
 R("R_IPPS") - R("R_Kabc") - R("R_KDOPP") - R("R_KDOPS") - R("R_LACZ") -
 R("R_LCAD") - R("R_LCADi") + R("R_LCAD_Rev") - R("R_LEUabc") -
 R("R_LPLIPA1") - R("R_LPLIPA2") - R("R_LPLIPA3") - R("R_LYSabc") - R("R_MALS") -
 R("R_MALTabc") - R("R_MALTHXabc") - R("R_MALTPTabc") - R("R_MALTTRabc")
 - R("R_MALTTTRabc") + R("R_MCITD") - R("R_MCITS") + R("R_MDRPD") +
 R("R_MECDPDH") - R("R_METabc") - R("R_METAT") - R("R_METDabc") -
 R("R_MI1PP") - R("R_MICITD") - R("R_MLTG1") - R("R_MLTG2") - R("R_MLTG3") -
 R("R_MLTG4") - R("R_MLTG5") + R("R_MNNH") - R("R_MOHMT") - R("R_MTAN") -
 R("R_MTHFC") + R("R_MTHFC_Rev") - R("R_NACODA") - R("R_NADDP") -
 R("R_NADPPPS") - R("R_NAMNPP") - R("R_NMNDA") - R("R_NMNN") - R("R_NMNt7")
 - R("R_NNAM") + R("R_NO3R1") + R("R_NO3R2") - R("R_NTD1") -
 R("R_NTD10") - R("R_NTD11") - R("R_NTD2") - R("R_NTD3") - R("R_NTD4") -
 R("R_NTD5") - R("R_NTD6") - R("R_NTD7") - R("R_NTD8") - R("R_NTD9") -
 R("R_NTPP1") - R("R_NTPP2") - R("R_NTPP3") - R("R_NTPP4") - R("R_NTPP5") -
 R("R_NTPP6") - R("R_NTPP7") - R("R_NTPP8") - R("R_NTPTP1") -
 R("R_NTPTP2") + 2 R("R_NTRIR2x") - R("R_OP4ENH") - R("R_ORNabc") - 2 R("R_P5CD")
 - R("R_PAPA_EC") + 2 R("R_PDX5PS") - R("R_PDXPP") - R("R_PEAMNO") -
 R("R_PGL") - R("R_PGLYCP") - R("R_PGPP_EC") - R("R_PIabc") -
 R("R_PLIPA1") - R("R_PLIPA2") - R("R_PLIPA3") - R("R_PMDPHT") - R("R_POX") -
 R("R_PPA") + 2 R("R_PPBNGS") - R("R_PPC") + R("R_PPNDH") + 3 R("R_PPPGO")
 - R("R_PPS") - R("R_PRAMPC") - R("R_PRATPP") - R("R_PRFGS") -
 R("R_PROabc") - R("R_PSP_L") - R("R_PTRCabc") - R("R_PYAM5PO") - R("R_PYDXPP")
 + 2 R("R_QULNS") + 2 R("R_RBFSa") - R("R_RIBabc") + R("R_RNDR1") +
 R("R_RNDR2") + R("R_RNDR3") + R("R_RNDR4") + R("R_RNTR1") + R("R_RNTR2")
 + R("R_RNTR3") + R("R_RNTR4") - R("R_RZ5PP") - 2 R("R_SADH") -
 R("R_SADT2") - R("R_SDPDS") - R("R_SELNPS") - R("R_SGDS") - R("R_SGSAD") -
 R("R_SPMDabc") - R("R_SSALx") - R("R_SSALy") + R("R_SUCBZS") -
 R("R_SUCCabc") - R("R_SULabc") - 3 R("R_SULR") + 3 R("R_SULR_Rev") +
 R("R_TARTD") - R("R_TAURabc") + R("R_TDPGDH") - R("R_THDPS") - R("R_THMabc") -
 R("R_THRabc") - R("R_THRS") + R("R_THZPSN") + R("R_TMAOR1") +
 R("R_TMAOR2") - R("R_TRE6PH") - R("R_TRE6PP") - R("R_TREH") - R("R_TRPAS2") +
 R("R_TRPAS2_Rev") + R("R_TRPS1") + R("R_TRPS2") - R("R_TSULabc") -
 R("R_UACMAMO") - R("R_UDCPDP") - R("R_UDPGD") - R("R_UGLYCH") -
 R("R_UHGADA") + R("R_UPP3S") - R("R_USHD") - R("R_VALabc") - R("R_XYLabc")  =
 0.

Mh2oe: R("R_DMSOR1e") + R("R_DMSOR2e") - R("R_EX_h2o_e_") +
 R("R_EX_h2o_e__Rev") - R("R_GLCDe") - R("R_H2Ot") + R("R_H2Ot_Rev") -
 R("R_NADDPe") + R("R_TMAOR1e") + R("R_TMAOR2e") - R("R_TREHe")  = 0.

Mh2sc: R("R_CYSDS") - R("R_CYSS") - R("R_SULR") + R("R_SULR_Rev")  =
 0.

Mhb: R("Mhb_Transp") + R("R_EX_h_e_") - R("R_EX_h_e__Rev")  = 0.

Mhc: -R("R_2DGLCNRx") - R("R_2DGLCNRy") - R("R_2DGULRx") -
 R("R_2DGULRy") - R("R_3HCINNMH") - R("R_3HPPPNH") - R("R_5DGLCNR") +
 R("R_5DGLCNR_Rev") + R("R_AADDGT") + 2 R("R_ABUTD") + R("R_ABUTt2") +
 R("R_ACACt2") - R("R_ACACt2_Rev") + R("R_ACALDi") - R("R_ACBIPGT") +
 R("R_ACCOACr") - R("R_ACCOACr_Rev") + R("R_ACGS") - R("R_ACHBS") - R("R_ACLS") +
 R("R_ACMAMUT") + R("R_ACNAMt2") + R("R_ACPS1") + R("R_ACt2r") -
 R("R_ACt2r_Rev") - R("R_ADA") + R("R_ADCL") - R("R_ADD") + R("R_ADEt2r") -
 R("R_ADEt2r_Rev") - 2 R("R_ADHEr") + 2 R("R_ADHEr_Rev") - R("R_ADMDCr")
 + R("R_ADMDCr_Rev") + R("R_ADNK1") + R("R_ADNt2") + R("R_ADNt2r") -
 R("R_ADNt2r_Rev") + R("R_ADOCBIK") + R("R_ADOCBLS") + R("R_ADSK") + 2
 R("R_ADSS") + R("R_AGPR") - R("R_AGPR_Rev") + R("R_AIRC2") +
 R("R_AKGt2r") - R("R_AKGt2r_Rev") + R("R_ALAabc") + R("R_ALAALAr") -
 R("R_ALAALAr_Rev") + R("R_ALAt2r") - R("R_ALAt2r_Rev") - R("R_ALCD19") +
 R("R_ALCD19_Rev") + 2 R("R_ALDD19x") + 2 R("R_ALDD2x") + R("R_ALLTN") +
 R("R_ALLTNt2r") - R("R_ALLTNt2r_Rev") + R("R_AMANK") + R("R_AMMQT8_2") + 4
 R("R_AMPMS") + R("R_ANS") - R("R_AOXSr") + R("R_AOXSr_Rev") + 2
 R("R_AP4AH") + 2 R("R_AP5AH") - R("R_APRAUR") + R("R_ARBabc") + R("R_ARBt2r")
 - R("R_ARBt2r_Rev") + R("R_ARGabc") - R("R_ARGDC") + R("R_ARGSS") +
 R("R_ASAD") - R("R_ASAD_Rev") + R("R_ASNabc") + R("R_ASNS1") +
 R("R_ASNS2") + R("R_ASNt2r") - R("R_ASNt2r_Rev") - R("R_ASP1DC") +
 R("R_ASPabc") + R("R_ASPCT") + R("R_ASPt2") + 2 R("R_ASPt2_2") + 3 R("R_ASPt2_3")
 + R("R_AST") + R("R_ATPM") + 3 R("R_ATPS4r") - 3 R("R_ATPS4r_Rev") +
 2 R("R_BETALDHx") + 2 R("R_BETALDHy") + 45.56035 R("R_BiomassEcoli") -
 R("R_BSORx") - R("R_BSORy") + 2 R("R_BTS2") - 2 R("R_BTS2_Rev") +
 R("R_BUTt2r") - R("R_BUTt2r_Rev") - 14 R("R_C120SN") - 17 R("R_C140SN") -
 16 R("R_C141SN") - 20 R("R_C160SN") - 19 R("R_C161SN") - 22
 R("R_C181SN") + R("R_CADVt") + R("R_CBL1abc") + 2 R("R_CBMK") + 2 R("R_CBPS") +
 2 R("R_CDAPPA_EC") + R("R_CDPMEK") + R("R_CHLabc") + R("R_CHLt2r") -
 R("R_CHLt2r_Rev") - R("R_CINNDO") - 2 R("R_CPPPGO") + R("R_CS") -
 R("R_CSND") + R("R_CSNt2") + 2 R("R_CTPS2") + R("R_CYANST") - 3
 R("R_CYNTAH") + R("R_CYNTt2") + R("R_CYSabc") + R("R_CYSS") - 2 R("R_CYTBD") -
 2.5 R("R_CYTBO3") - R("R_CYTD") + R("R_CYTDK2") + R("R_CYTDt2") +
 R("R_CYTDt2r") - R("R_CYTDt2r_Rev") - R("R_DADA") + R("R_DADNt2") +
 R("R_DAGK_EC") + R("R_DALAt2r") - R("R_DALAt2r_Rev") + R("R_DAPabc") -
 R("R_DAPDC") - R("R_DASYN_EC") + R("R_DASYN_EC_Rev") + R("R_DB4PS") + 3
 R("R_DBTSr") - 3 R("R_DBTSr_Rev") - R("R_DCTPD") - R("R_DCYTD") +
 R("R_DCYTt2") + R("R_DDGALK") + R("R_DDGLCNt2r") - R("R_DDGLCNt2r_Rev") +
 R("R_DDGLK") + R("R_DGSNt2") + R("R_DHBD") - R("R_DHBD_Rev") +
 R("R_DHCIND") - R("R_DHDPRy") + R("R_DHDPS") - R("R_DHFR") + R("R_DHFR_Rev") +
 R("R_DHNAOT") + R("R_DHORTS") - R("R_DHORTS_Rev") + R("R_DHPPD") -
 R("R_DHPPDA2") + R("R_DHPS2") + R("R_DINSt2") - R("R_DKGLCNR1") -
 R("R_DKGLCNR2x") - R("R_DKGLCNR2y") + 2 R("R_DKMPPD") + 6 R("R_DKMPPD2") +
 R("R_D_LACt2") - R("R_D_LACt2_Rev") - R("R_DMPPS") + R("R_DMQMT") +
 R("R_DNTPPA") - R("R_DOGULNR") + R("R_DPCOAK") - R("R_DPR") + R("R_DSERt2r")
 - R("R_DSERt2r_Rev") + R("R_DURIK1") + R("R_DURIt2") + R("R_DUTPDP")
 - R("R_DXPRIi") - R("R_DXPS") + R("R_DXYLK") + 2 R("R_E4PD") - 2
 R("R_E4PD_Rev") + R("R_ECAP_EC") + 6 R("R_ENTCS") + R("R_ETOHt2r") -
 R("R_ETOHt2r_Rev") + 6 R("R_FAO1") + 7 R("R_FAO2") + 8 R("R_FAO3") +
 R("R_FAO4") + R("R_FCLK") + 2 R("R_FCLT") - 3 R("R_FDH2") - 3 R("R_FDH3") +
 R("R_FE2abc") - R("R_FHL") - R("R_FMNAT") + R("R_FRUK") + R("R_FTHFD")
 + R("R_FUCt") - R("R_FUCt_Rev") + 2 R("R_FUMt2_2") + 3 R("R_FUMt2_3")
 + R("R_G1PACT") - R("R_G1PTT") + R("R_G3PD2") - R("R_G3PD2_Rev") +
 R("R_G5SADs") - R("R_G5SD") + R("R_G6PDH2r") - R("R_G6PDH2r_Rev") +
 R("R_GALabc") + R("R_GALCTNt2r") - R("R_GALCTNt2r_Rev") + R("R_GALCTt2r") -
 R("R_GALCTt2r_Rev") + R("R_GALKr") - R("R_GALKr_Rev") + R("R_GALt2")
 - R("R_GALU") - R("R_GALUi") + R("R_GALU_Rev") + R("R_GALURt2r") -
 R("R_GALURt2r_Rev") + R("R_GAPD") - R("R_GAPD_Rev") + R("R_GARFT") -
 R("R_GARFT_Rev") + R("R_GART") + 2 R("R_GCALDD") + R("R_GLCNt2r") -
 R("R_GLCNt2r_Rev") + R("R_GLCRt2r") - R("R_GLCRt2r_Rev") + R("R_GLCS1") +
 R("R_GLCt2") + R("R_GLCURt2r") - R("R_GLCURt2r_Rev") - R("R_GLGC") +
 R("R_GLNabc") + R("R_GLNS") + R("R_GLTPD") - R("R_GLTPD_Rev") +
 R("R_GLUabc") + R("R_GLUCYS") - R("R_GLUDC") + R("R_GLUDy") - R("R_GLUDy_Rev")
 - R("R_GLUSy") + R("R_GLUt2r") - R("R_GLUt2r_Rev") - R("R_GLUTRR") -
 R("R_GLXCL") + R("R_GLYBabc") + R("R_GLYBt2r") - R("R_GLYBt2r_Rev") +
 R("R_GLYC3Pabc") + R("R_GLYCDx") + R("R_GLYCK") - R("R_GLYCLTDx") -
 R("R_GLYCLTDy") + R("R_GLYCLTt2r") - R("R_GLYCLTt2r_Rev") + R("R_GLYK") +
 R("R_GLYOX") + R("R_GLYt2r") - R("R_GLYt2r_Rev") - R("R_GMHEPAT") +
 R("R_GMHEPK") - 2 R("R_GMPR") + 2 R("R_GMPS2") + R("R_GNK") -
 R("R_GOFUCR") + 2 R("R_GP4GH") + R("R_GPDDA1") + R("R_GPDDA2") + R("R_GPDDA3")
 + R("R_GPDDA4") + R("R_GPDDA5") + R("R_GSNK") + R("R_GSNt2") +
 R("R_GSPMDS") - R("R_GTHOr") + R("R_GTHOr_Rev") + R("R_GTHS") + 2
 R("R_GTPCII2") - R("R_GUAD") + R("R_GUAt2") + R("R_HCINNMt2r") -
 R("R_HCINNMt2r_Rev") + R("R_HCO3E") - R("R_HCO3E_Rev") + R("R_HDCAt2") + R("R_HETZK")
 + R("R_HEX1") + R("R_HISabc") + R("R_HISt2r") - R("R_HISt2r_Rev") + 3
 R("R_HISTD") + 2 R("R_HKNDDH") + 2 R("R_HKNTDH") + R("R_HMPK1") +
 R("R_HPPK2") + R("R_HPPPNt2r") - R("R_HPPPNt2r_Rev") - R("R_HPYRRx") -
 R("R_HPYRRy") + R("R_HSDy") - R("R_HSDy_Rev") + R("R_HSK") - 2
 R("R_HYD1") - 2 R("R_HYD2") - 2 R("R_HYD3") - R("R_IDOND") - R("R_IDOND2") +
 R("R_IDOND_Rev") + R("R_IDONt2r") - R("R_IDONt2r_Rev") + R("R_IG3PS") -
 R("R_IGPS") + R("R_ILEabc") + R("R_ILEt2r") - R("R_ILEt2r_Rev") +
 R("R_IMPD") + R("R_INDOLEt2r") - R("R_INDOLEt2r_Rev") + R("R_INSK") +
 R("R_INSt2") + R("R_INSt2r") - R("R_INSt2r_Rev") - R("R_IPDPS") +
 R("R_IPMD") + R("R_IPPS") + R("R_Kabc") - R("R_KARA1i") - R("R_KARA2i") -
 R("R_KAS14") - R("R_KAS15") - 2 R("R_KAS16") - R("R_KG6PDC") + R("R_Kt2r")
 - R("R_Kt2r_Rev") + 2 R("R_LCAD") + 2 R("R_LCADi") - 2 R("R_LCAD_Rev")
 - R("R_LCARS") + R("R_LCARS_Rev") + R("R_LCTSt") - R("R_LCTSt_Rev") +
 R("R_LDH_D") - R("R_LDH_D_Rev") + R("R_LEUabc") + R("R_LEUt2r") -
 R("R_LEUt2r_Rev") + R("R_L_LACt2r") - R("R_L_LACt2r_Rev") + R("R_LPADSS")
 + R("R_LPLIPA1") + R("R_LPLIPA2") + R("R_LPLIPA3") + 10
 R("R_LPSSYN_EC") + R("R_LYSabc") - R("R_LYSDC") + R("R_LYSt2r") - R("R_LYSt2r_Rev")
 + R("R_M1PD") - R("R_M1PD_Rev") - R("R_MACPD") + R("R_MALS") + 2
 R("R_MALt2_2") + 3 R("R_MALt2_3") + R("R_MALTabc") + R("R_MALTHXabc") +
 R("R_MALTPTabc") + R("R_MALTTRabc") + R("R_MALTTTRabc") - R("R_MAN1PT2")
 + R("R_MANAO") - R("R_MANAO_Rev") + R("R_MCITS") + R("R_MDH") -
 R("R_MDH_Rev") - R("R_MECDPDH") + R("R_MELIBt2") - R("R_MEPCT") +
 R("R_METabc") + R("R_METDabc") - R("R_MMCD") + R("R_MOAT") + R("R_MOAT2") +
 R("R_MTHFD") - R("R_MTHFD_Rev") - R("R_MTHFR2") + R("R_MTRK") + 2
 R("R_NADDP") - R("R_NADH10") - R("R_NADH5") - 4.5 R("R_NADH6") - 3
 R("R_NADH7") - 3.8 R("R_NADH8") - R("R_NADH9") + R("R_NADK") + R("R_NADS1") +
 R("R_NAt3_1") + 3 R("R_NAt3_15") - R("R_NAt3_1_Rev") + 2 R("R_NAt3_2") -
 R("R_NMNAT") + R("R_NMNN") + R("R_NMNt7") - R("R_NNAT") +
 R("R_NNDMBRT") - 2 R("R_NNDPR") + R("R_NO2t2r") - R("R_NO2t2r_Rev") - 2
 R("R_NO3R1") - 2 R("R_NO3R2") + R("R_NTPP1") + R("R_NTPP2") + R("R_NTPP3") +
 R("R_NTPP4") + R("R_NTPP5") + R("R_NTPP6") + R("R_NTPP7") + R("R_NTPP8")
 - 5 R("R_NTRIR2x") + R("R_OCBT") - R("R_OCBT_Rev") + R("R_OCDCAt2") +
 R("R_OHPHM") + R("R_OMBZLM") - R("R_OMCDC") - R("R_OMPDC") -
 R("R_OPHBDC") + R("R_ORNabc") - R("R_ORNDC") - R("R_OXGDC2") + R("R_P5CD") - 2
 R("R_P5CR") + R("R_PANTS") + 2 R("R_PAPSR") + R("R_PDX5PS") -
 R("R_PEPT_EC") + R("R_PEPT_EC_Rev") + R("R_PERD") - R("R_PERD_Rev") +
 R("R_PFK") + R("R_PFK_2") + R("R_PGCD") + R("R_PGL") + R("R_PGSA_EC") -
 R("R_PGSA_EC_Rev") + R("R_PHEt2r") - R("R_PHEt2r_Rev") + R("R_PIabc") +
 R("R_PIt2r") - R("R_PIt2r_Rev") + R("R_PLIPA1") + R("R_PLIPA2") +
 R("R_PLIPA3") + R("R_PNTK") + R("R_PPA") + R("R_PPBNGS") + R("R_PPC") -
 R("R_PPCDC") + R("R_PPNCL2") - R("R_PPNDH") - R("R_PPPNDO") +
 R("R_PPPNt2r") - R("R_PPPNt2r_Rev") + 2 R("R_PPS") + R("R_PPTGS") + R("R_PRAGSr") -
 R("R_PRAGSr_Rev") + 2 R("R_PRAIS") + R("R_PRASCS") -
 R("R_PRASCS_Rev") + R("R_PRATPP") + R("R_PRFGS") + R("R_PROabc") + R("R_PROD2") +
 R("R_PROt2r") - R("R_PROt2r_Rev") + R("R_PRPPS") - R("R_PRPPS_Rev") -
 R("R_PSD_EC") + R("R_PSSA_EC") - R("R_PSSA_EC_Rev") - R("R_PTPATi") +
 R("R_PTRCabc") + R("R_PTRCt2r") - R("R_PTRCt2r_Rev") + R("R_PYDAMK") +
 R("R_PYDXK") + R("R_PYDXNK") - R("R_PYK") + R("R_PYRt2r") -
 R("R_PYRt2r_Rev") + R("R_QULNS") + R("R_RBFK") + R("R_RBK") + R("R_RBK_L1") +
 R("R_RIBabc") + R("R_RMK") + R("R_RMNt") - 2 R("R_SADH") + R("R_SBTPD") -
 R("R_SBTPD_Rev") - R("R_SERASr") + R("R_SERASr_Rev") + R("R_SERt2r") -
 R("R_SERt2r_Rev") + 2 R("R_SGSAD") + R("R_SHCHD2") + 3 R("R_SHCHF") -
 R("R_SHK3Dr") + R("R_SHK3Dr_Rev") + R("R_SHKK") + R("R_SHSL1") +
 R("R_SPMDabc") + R("R_SPMDAT1") + R("R_SPMDAT2") + R("R_SPMS") - 2
 R("R_SPODM") + 2 R("R_SSALx") + 2 R("R_SSALy") + R("R_SUCCabc") + 2
 R("R_SUCCt2_2") + 3 R("R_SUCCt2_3") - R("R_SUCCt2b") + R("R_SULabc") + 5
 R("R_SULR") - 5 R("R_SULR_Rev") + R("R_TAGURr") - R("R_TAGURr_Rev") +
 R("R_TAUDO") + R("R_TAURabc") + R("R_TDPADGAT") - R("R_TDPDRR") + R("R_TDSK")
 + 2 R("R_THD2") + R("R_THMabc") + R("R_THMDt2") + R("R_THMDt2r") -
 R("R_THMDt2r_Rev") + R("R_THRabc") + R("R_THRD") + R("R_THRt2r") -
 R("R_THRt2r_Rev") + R("R_THZPSN") - R("R_TMAOR1") - R("R_TMAOR2") +
 R("R_TMDK1") + R("R_TMKr") - R("R_TMKr_Rev") - R("R_TMPPP") - R("R_TRDR") +
 R("R_TRE6PS") + R("R_TRPt2r") - R("R_TRPt2r_Rev") - R("R_TRSAR") +
 R("R_TSULabc") + R("R_TTDCAt2") + R("R_TYRt2r") - R("R_TYRt2r_Rev") +
 R("R_U23GAAT") + R("R_UAAGDS") + 3 R("R_UACMAMO") - R("R_UAGDP") +
 R("R_UAGPT3") + R("R_UAMAGS") + R("R_UAMAS") - R("R_UAPGR") + R("R_UDCPDP") +
 3 R("R_UDPGD") - 2 R("R_UGLYCH") + R("R_UGMDDS") + R("R_UPP3MT") - 4
 R("R_UPPDC1") + R("R_URAt2") + R("R_URAt2r") - R("R_URAt2r_Rev") +
 R("R_URIK2") + R("R_URIt2") + R("R_URIt2r") - R("R_URIt2r_Rev") + 2
 R("R_USHD") + R("R_VALabc") + R("R_VALt2r") - R("R_VALt2r_Rev") +
 R("R_XANt2") + R("R_XTSNt2r") - R("R_XTSNt2r_Rev") + R("R_XYLabc") +
 R("R_XYLK") + R("R_XYLt2")  = 0.

Mhco3c: -R("R_ACCOACr") + R("R_ACCOACr_Rev") - R("R_AIRC2") -
 R("R_CBPS") - R("R_CYNTAH") + R("R_HCO3E") - R("R_HCO3E_Rev")  = 0.

MhcysLc: R("R_AHC") - R("R_AHC_Rev") + R("R_CYSTL") - R("R_METS") +
 R("R_RHCCE")  = 0.

Mhdcab: R("Mhdcab_Transp") + R("R_EX_hdca_e_") -
 R("R_EX_hdca_e__Rev")  = 0.

Mhdcac: -R("R_AACPS3") - R("R_FAO2") + R("R_HDCAt2") + 0.36
 R("R_LPLIPA1") + 0.36 R("R_LPLIPA2") + 0.36 R("R_LPLIPA3") + 0.36
 R("R_PLIPA1") + 0.36 R("R_PLIPA2") + 0.36 R("R_PLIPA3")  = 0.

Mhdcae: -R("R_EX_hdca_e_") + R("R_EX_hdca_e__Rev") - R("R_HDCAt2") 
 = 0.

Mhdceac: -R("R_AACPS4") + 0.07 R("R_LPLIPA1") + 0.07 R("R_LPLIPA2")
 + 0.07 R("R_LPLIPA3") + 0.07 R("R_PLIPA1") + 0.07 R("R_PLIPA2") + 0.07
 R("R_PLIPA3")  = 0.

MhdeACPc: R("R_AACPS4") + R("R_C161SN") - R("R_EDTXS3") - 0.14
 R("R_PASYN_EC")  = 0.

Mhe: -R("R_ABUTt2") - R("R_ACACt2") + R("R_ACACt2_Rev") -
 R("R_ACNAMt2") - R("R_ACt2r") + R("R_ACt2r_Rev") - R("R_ADEt2r") +
 R("R_ADEt2r_Rev") - R("R_ADNt2") - R("R_ADNt2r") + R("R_ADNt2r_Rev") -
 R("R_AKGt2r") + R("R_AKGt2r_Rev") - R("R_ALAt2r") + R("R_ALAt2r_Rev") -
 R("R_ALLTNt2r") + R("R_ALLTNt2r_Rev") - R("R_ARBt2r") + R("R_ARBt2r_Rev") -
 R("R_ASNt2r") + R("R_ASNt2r_Rev") - R("R_ASPt2") - 2 R("R_ASPt2_2") - 3
 R("R_ASPt2_3") - 4 R("R_ATPS4r") + 4 R("R_ATPS4r_Rev") - R("R_BUTt2r") +
 R("R_BUTt2r_Rev") - R("R_CADVt") - R("R_CHLt2r") + R("R_CHLt2r_Rev")
 - R("R_CSNt2") - R("R_CYNTt2") + 2 R("R_CYTBD") + 2.5 R("R_CYTBO3") -
 R("R_CYTDt2") - R("R_CYTDt2r") + R("R_CYTDt2r_Rev") - R("R_DADNt2") -
 R("R_DALAt2r") + R("R_DALAt2r_Rev") - R("R_DCYTt2") - R("R_DDGLCNt2r")
 + R("R_DDGLCNt2r_Rev") - R("R_DGSNt2") - R("R_DINSt2") -
 R("R_D_LACt2") + R("R_D_LACt2_Rev") - R("R_DSERt2r") + R("R_DSERt2r_Rev") -
 R("R_DURIt2") - R("R_ETOHt2r") + R("R_ETOHt2r_Rev") - R("R_EX_h_e_") +
 R("R_EX_h_e__Rev") + 2 R("R_FDH2") + 2 R("R_FDH3") - R("R_FUCt") +
 R("R_FUCt_Rev") - 2 R("R_FUMt2_2") - 3 R("R_FUMt2_3") - R("R_GALCTNt2r") +
 R("R_GALCTNt2r_Rev") - R("R_GALCTt2r") + R("R_GALCTt2r_Rev") -
 R("R_GALt2") - R("R_GALURt2r") + R("R_GALURt2r_Rev") + R("R_GLCDe") -
 R("R_GLCNt2r") + R("R_GLCNt2r_Rev") - R("R_GLCRt2r") + R("R_GLCRt2r_Rev") -
 R("R_GLCt2") - R("R_GLCURt2r") + R("R_GLCURt2r_Rev") - R("R_GLUt2r") +
 R("R_GLUt2r_Rev") - R("R_GLYBt2r") + R("R_GLYBt2r_Rev") -
 R("R_GLYCLTt2r") + R("R_GLYCLTt2r_Rev") - R("R_GLYt2r") + R("R_GLYt2r_Rev") -
 R("R_GSNt2") - R("R_GUAt2") - R("R_HCINNMt2r") + R("R_HCINNMt2r_Rev") -
 R("R_HDCAt2") - R("R_HISt2r") + R("R_HISt2r_Rev") - R("R_HPPPNt2r") +
 R("R_HPPPNt2r_Rev") + 2 R("R_HYD1") + 2 R("R_HYD2") + 2 R("R_HYD3") -
 R("R_IDONt2r") + R("R_IDONt2r_Rev") - R("R_ILEt2r") + R("R_ILEt2r_Rev") -
 R("R_INDOLEt2r") + R("R_INDOLEt2r_Rev") - R("R_INSt2") - R("R_INSt2r")
 + R("R_INSt2r_Rev") - R("R_Kt2r") + R("R_Kt2r_Rev") - R("R_LCTSt") +
 R("R_LCTSt_Rev") - R("R_LEUt2r") + R("R_LEUt2r_Rev") - R("R_L_LACt2r")
 + R("R_L_LACt2r_Rev") - R("R_LYSt2r") + R("R_LYSt2r_Rev") - 2
 R("R_MALt2_2") - 3 R("R_MALt2_3") - R("R_MELIBt2") + 2 R("R_NADDPe") + 3.5
 R("R_NADH6") + 2 R("R_NADH7") + 2.8 R("R_NADH8") - R("R_NAt3_1") - 3
 R("R_NAt3_15") + R("R_NAt3_1_Rev") - 2 R("R_NAt3_2") - R("R_NO2t2r") +
 R("R_NO2t2r_Rev") + 2 R("R_NO3R1") + 2 R("R_NO3R2") - R("R_OCDCAt2") -
 R("R_PHEt2r") + R("R_PHEt2r_Rev") - R("R_PIt2r") + R("R_PIt2r_Rev") -
 R("R_PPPNt2r") + R("R_PPPNt2r_Rev") - R("R_PROt2r") + R("R_PROt2r_Rev")
 - R("R_PTRCt2r") + R("R_PTRCt2r_Rev") - R("R_PYRt2r") +
 R("R_PYRt2r_Rev") - R("R_RMNt") - R("R_SERt2r") + R("R_SERt2r_Rev") - 2
 R("R_SUCCt2_2") - 3 R("R_SUCCt2_3") + R("R_SUCCt2b") - 2 R("R_THD2") -
 R("R_THMDt2") - R("R_THMDt2r") + R("R_THMDt2r_Rev") - R("R_THRt2r") +
 R("R_THRt2r_Rev") - R("R_TMAOR1e") - R("R_TMAOR2e") - R("R_TRPt2r") +
 R("R_TRPt2r_Rev") - R("R_TTDCAt2") - R("R_TYRt2r") + R("R_TYRt2r_Rev") -
 R("R_URAt2") - R("R_URAt2r") + R("R_URAt2r_Rev") - R("R_URIt2") -
 R("R_URIt2r") + R("R_URIt2r_Rev") - R("R_VALt2r") + R("R_VALt2r_Rev") -
 R("R_XANt2") - R("R_XTSNt2r") + R("R_XTSNt2r_Rev") - R("R_XYLt2")  = 0.

MhemeOc: R("R_HEMEOS")  = 0.

MhisLb: R("MhisLb_Transp") + R("R_EX_his_L_e_") -
 R("R_EX_his_L_e__Rev")  = 0.

MhisLc: -0.09 R("R_BiomassEcoli") + R("R_HISabc") + R("R_HISt2r") -
 R("R_HISt2r_Rev") + R("R_HISTD")  = 0.

MhisLe: -R("R_EX_his_L_e_") + R("R_EX_his_L_e__Rev") - R("R_HISabc")
 - R("R_HISt2r") + R("R_HISt2r_Rev")  = 0.

Mhispc: -R("R_HISTP") + R("R_HSTPT")  = 0.

Mhistdc: -R("R_HISTD") + R("R_HISTP")  = 0.

Mhknddc: -R("R_HKNDDH") + R("R_HPPPNDO")  = 0.

Mhkntdc: R("R_DHCINDO") - R("R_HKNTDH")  = 0.

Mhmbilc: R("R_HMBS") - R("R_UPP3S")  = 0.

Mhmfurnc: R("R_DHPTDCs")  = 0.

MhomLc: -R("R_HSDy") + R("R_HSDy_Rev") - R("R_HSK") - R("R_HSST")  =
 0.

Mhpyrc: -R("R_HPYRI") + R("R_HPYRI_Rev") - R("R_HPYRRx") -
 R("R_HPYRRy")  = 0.

Mhqnc: R("R_AB6PGH")  = 0.

Mhxanb: R("Mhxanb_Transp") + R("R_EX_hxan_e_") -
 R("R_EX_hxan_e__Rev")  = 0.

Mhxanc: R("R_ADD") - R("R_HXPRT") + R("R_HYXNt") - R("R_HYXNt_Rev")
 + R("R_PUNP5") - R("R_PUNP5_Rev") + R("R_PUNP6") - R("R_PUNP6_Rev")  =
 0.

Mhxane: -R("R_EX_hxan_e_") + R("R_EX_hxan_e__Rev") - R("R_HYXNt") +
 R("R_HYXNt_Rev")  = 0.

Miaspc: R("R_ASPO3") + R("R_ASPO4") + R("R_ASPO5") + R("R_ASPO6") -
 R("R_QULNS")  = 0.

Michorc: R("R_ICHORSi") - R("R_ICHORT") - R("R_SHCHCS2")  = 0.

Micitc: R("R_ACONT") - R("R_ACONT_Rev") - R("R_ICDHyr") +
 R("R_ICDHyr_Rev") - R("R_ICL")  = 0.

MidonLb: R("MidonLb_Transp") + R("R_EX_idon_L_e_") -
 R("R_EX_idon_L_e__Rev")  = 0.

MidonLc: R("R_2DGULRx") + R("R_2DGULRy") + R("R_IDOND") +
 R("R_IDOND2") - R("R_IDOND_Rev") + R("R_IDONt2r") - R("R_IDONt2r_Rev")  = 0.

MidonLe: -R("R_EX_idon_L_e_") + R("R_EX_idon_L_e__Rev") -
 R("R_IDONt2r") + R("R_IDONt2r_Rev")  = 0.

Midpc: R("R_ADK4") - R("R_ADK4_Rev")  = 0.

MileLb: R("MileLb_Transp") + R("R_EX_ile_L_e_") -
 R("R_EX_ile_L_e__Rev")  = 0.

MileLc: -0.276 R("R_BiomassEcoli") + R("R_ILEabc") + R("R_ILEt2r") -
 R("R_ILEt2r_Rev") - R("R_ILETA") + R("R_ILETA_Rev")  = 0.

MileLe: -R("R_EX_ile_L_e_") + R("R_EX_ile_L_e__Rev") - R("R_ILEabc")
 - R("R_ILEt2r") + R("R_ILEt2r_Rev")  = 0.

Mimacpc: -R("R_HSTPT") + R("R_IGPDH")  = 0.

Mimpc: -R("R_ADSS") + R("R_GMPR") + R("R_HXPRT") - R("R_IMPC") +
 R("R_IMPC_Rev") - R("R_IMPD") + R("R_INSK") - R("R_NTD11")  = 0.

Mindoleb: R("Mindoleb_Transp") + R("R_EX_indole_e_") -
 R("R_EX_indole_e__Rev")  = 0.

Mindolec: R("R_INDOLEt2r") - R("R_INDOLEt2r_Rev") + R("R_TRPAS2") -
 R("R_TRPAS2_Rev") - R("R_TRPS2") + R("R_TRPS3")  = 0.

Mindolee: -R("R_EX_indole_e_") + R("R_EX_indole_e__Rev") -
 R("R_INDOLEt2r") + R("R_INDOLEt2r_Rev")  = 0.

Minostc: R("R_GPDDA5") + R("R_MI1PP")  = 0.

Minsb: R("Minsb_Transp") + R("R_EX_ins_e_") - R("R_EX_ins_e__Rev") 
 = 0.

Minsc: R("R_ADA") - R("R_INSK") + R("R_INSt2") + R("R_INSt2r") -
 R("R_INSt2r_Rev") + R("R_NTD11") - R("R_PUNP5") + R("R_PUNP5_Rev")  = 0.

Minse: -R("R_EX_ins_e_") + R("R_EX_ins_e__Rev") - R("R_INSt2") -
 R("R_INSt2r") + R("R_INSt2r_Rev")  = 0.

Mipdpc: -R("R_DMATT") - R("R_GRTT") - R("R_IPDDIi") + R("R_IPDPS") -
 5 R("R_OCTDPS") - 8 R("R_UDCPDPS")  = 0.

Mitpc: -R("R_ADK4") + R("R_ADK4_Rev")  = 0.

Mkb: R("Mkb_Transp") + R("R_EX_k_e_") - R("R_EX_k_e__Rev")  = 0.

Mkc: R("R_Kabc") + R("R_Kt2r") - R("R_Kt2r_Rev")  = 0.

Mkdo2lipid4c: -R("R_EDTXS1") - R("R_EDTXS3") + R("R_MOAT2")  = 0.

Mkdo2lipid4Lc: R("R_EDTXS1") - R("R_EDTXS2")  = 0.

Mkdo2lipid4pc: R("R_EDTXS3") - R("R_EDTXS4")  = 0.

Mkdo8pc: -R("R_KDOPP") + R("R_KDOPS")  = 0.

Mkdoc: -R("R_KDOCT2") + R("R_KDOPP")  = 0.

Mkdolipid4c: R("R_MOAT") - R("R_MOAT2")  = 0.

Mke: -R("R_EX_k_e_") + R("R_EX_k_e__Rev") - R("R_Kabc") -
 R("R_Kt2r") + R("R_Kt2r_Rev")  = 0.

MlacDb: R("MlacDb_Transp") + R("R_EX_lac_D_e_") -
 R("R_EX_lac_D_e__Rev")  = 0.

MlacDc: R("R_D_LACt2") - R("R_D_LACt2_Rev") + R("R_GLYOX") -
 R("R_LDH_D") - R("R_LDH_D2") + R("R_LDH_D_Rev")  = 0.

MlacDe: -R("R_D_LACt2") + R("R_D_LACt2_Rev") - R("R_EX_lac_D_e_") +
 R("R_EX_lac_D_e__Rev")  = 0.

MlacLb: R("MlacLb_Transp") + R("R_EX_lac_L_e_") -
 R("R_EX_lac_L_e__Rev")  = 0.

MlacLc: R("R_LCAD") + R("R_LCADi") - R("R_LCAD_Rev") -
 R("R_L_LACD2") - R("R_L_LACD3") + R("R_L_LACt2r") - R("R_L_LACt2r_Rev")  = 0.

MlacLe: -R("R_EX_lac_L_e_") + R("R_EX_lac_L_e__Rev") -
 R("R_L_LACt2r") + R("R_L_LACt2r_Rev")  = 0.

MlaldLc: R("R_FCLPA") - R("R_FCLPA_Rev") - R("R_LCAD") -
 R("R_LCADi") + R("R_LCAD_Rev") - R("R_LCARS") + R("R_LCARS_Rev") + R("R_RMPA") -
 R("R_RMPA_Rev")  = 0.

Mlctsb: R("Mlctsb_Transp") + R("R_EX_lcts_e_") -
 R("R_EX_lcts_e__Rev")  = 0.

Mlctsc: -R("R_LACZ") + R("R_LCTSt") - R("R_LCTSt_Rev")  = 0.

Mlctse: -R("R_EX_lcts_e_") + R("R_EX_lcts_e__Rev") - R("R_LCTSt") +
 R("R_LCTSt_Rev")  = 0.

MleuLb: R("MleuLb_Transp") + R("R_EX_leu_L_e_") -
 R("R_EX_leu_L_e__Rev")  = 0.

MleuLc: -0.428 R("R_BiomassEcoli") + R("R_LEUabc") + R("R_LEUt2r") -
 R("R_LEUt2r_Rev") + R("R_LEUTAi")  = 0.

MleuLe: -R("R_EX_leu_L_e_") + R("R_EX_leu_L_e__Rev") - R("R_LEUabc")
 - R("R_LEUt2r") + R("R_LEUt2r_Rev")  = 0.

MlgtSc: -R("R_GLYOX") + R("R_LGTHL")  = 0.

Mlipac: R("R_EDTXS2") - R("R_LPSSYN_EC")  = 0.

Mlipacoldc: R("R_EDTXS4")  = 0.

MlipidAc: -R("R_MOAT") + R("R_TDSK")  = 0.

MlipidAdsc: R("R_LPADSS") - R("R_TDSK")  = 0.

MlipidXc: -R("R_LPADSS") + R("R_USHD")  = 0.

MlpsECc: -0.0084 R("R_BiomassEcoli") + R("R_LPSSYN_EC")  = 0.

MlysLb: R("MlysLb_Transp") + R("R_EX_lys_L_e_") -
 R("R_EX_lys_L_e__Rev")  = 0.

MlysLc: -0.326 R("R_BiomassEcoli") + R("R_CADVt") + R("R_DAPDC") +
 R("R_LYSabc") - R("R_LYSDC") + R("R_LYSt2r") - R("R_LYSt2r_Rev")  = 0.

MlysLe: -R("R_CADVt") - R("R_EX_lys_L_e_") + R("R_EX_lys_L_e__Rev")
 - R("R_LYSabc") - R("R_LYSt2r") + R("R_LYSt2r_Rev")  = 0.

MmalACPc: -4 R("R_C120SN") - 5 R("R_C140SN") - 5 R("R_C141SN") - 6
 R("R_C160SN") - 6 R("R_C161SN") - 7 R("R_C181SN") - R("R_KAS14") -
 R("R_KAS15") - R("R_KAS16") - R("R_MACPD") + R("R_MCOATA") -
 R("R_MCOATA_Rev")  = 0.

Mmalcoac: R("R_ACCOACr") - R("R_ACCOACr_Rev") - R("R_MCOATA") +
 R("R_MCOATA_Rev")  = 0.

MmalLb: R("MmalLb_Transp") + R("R_EX_mal_L_e_") -
 R("R_EX_mal_L_e__Rev")  = 0.

MmalLc: R("R_FUM") - R("R_FUM_Rev") + R("R_MALS") + R("R_MALt2_2") +
 R("R_MALt2_3") - R("R_MDH") - R("R_MDH2") - R("R_MDH3") +
 R("R_MDH_Rev") - R("R_ME1") - R("R_ME2")  = 0.

MmalLe: -R("R_EX_mal_L_e_") + R("R_EX_mal_L_e__Rev") -
 R("R_MALt2_2") - R("R_MALt2_3")  = 0.

Mmalt6pc: R("R_MALTpts")  = 0.

Mmaltb: R("Mmaltb_Transp") + R("R_EX_malt_e_") -
 R("R_EX_malt_e__Rev")  = 0.

Mmaltc: -R("R_AMALT1") - R("R_AMALT2") - R("R_AMALT3") -
 R("R_AMALT4") + R("R_MALTabc") + R("R_MLTG1")  = 0.

Mmalte: -R("R_EX_malt_e_") + R("R_EX_malt_e__Rev") - R("R_MALTabc")
 - R("R_MALTpts")  = 0.

Mmalthpc: R("R_AMALT4") - R("R_MLTG5") - R("R_MLTP3") +
 R("R_MLTP3_Rev")  = 0.

Mmalthxb: R("Mmalthxb_Transp") + R("R_EX_malthx_e_") -
 R("R_EX_malthx_e__Rev")  = 0.

Mmalthxc: R("R_AMALT3") - R("R_AMALT4") + R("R_MALTHXabc") -
 R("R_MLTG4") + R("R_MLTG5") - R("R_MLTP2") + R("R_MLTP2_Rev") + R("R_MLTP3")
 - R("R_MLTP3_Rev")  = 0.

Mmalthxe: -R("R_EX_malthx_e_") + R("R_EX_malthx_e__Rev") -
 R("R_MALTHXabc")  = 0.

Mmaltptb: R("Mmaltptb_Transp") + R("R_EX_maltpt_e_") -
 R("R_EX_maltpt_e__Rev")  = 0.

Mmaltptc: R("R_AMALT2") - R("R_AMALT3") + R("R_MALTPTabc") -
 R("R_MLTG3") + R("R_MLTG4") - R("R_MLTP1") + R("R_MLTP1_Rev") + R("R_MLTP2")
 - R("R_MLTP2_Rev")  = 0.

Mmaltpte: -R("R_EX_maltpt_e_") + R("R_EX_maltpt_e__Rev") -
 R("R_MALTPTabc")  = 0.

Mmalttrb: R("Mmalttrb_Transp") + R("R_EX_malttr_e_") -
 R("R_EX_malttr_e__Rev")  = 0.

Mmalttrc: -R("R_AMALT1") + R("R_MALTTRabc") - R("R_MLTG1") +
 R("R_MLTG2")  = 0.

Mmalttre: -R("R_EX_malttr_e_") + R("R_EX_malttr_e__Rev") -
 R("R_MALTTRabc")  = 0.

Mmaltttrb: R("Mmaltttrb_Transp") + R("R_EX_maltttr_e_") -
 R("R_EX_maltttr_e__Rev")  = 0.

Mmaltttrc: R("R_AMALT1") - R("R_AMALT2") + R("R_MALTTTRabc") -
 R("R_MLTG2") + R("R_MLTG3") + R("R_MLTP1") - R("R_MLTP1_Rev")  = 0.

Mmaltttre: -R("R_EX_maltttr_e_") + R("R_EX_maltttr_e__Rev") -
 R("R_MALTTTRabc")  = 0.

Mman1pc: -R("R_MAN1PT2") - R("R_PMANM") + R("R_PMANM_Rev")  = 0.

Mman6pb: R("Mman6pb_Transp") + R("R_EX_man6p_e_") -
 R("R_EX_man6p_e__Rev")  = 0.

Mman6pc: -R("R_MAN6PI") + R("R_MAN6PI_Rev") + R("R_MAN6Pt6_2") +
 R("R_MANpts") + R("R_PMANM") - R("R_PMANM_Rev")  = 0.

Mman6pe: -R("R_EX_man6p_e_") + R("R_EX_man6p_e__Rev") -
 R("R_MAN6Pt6_2")  = 0.

Mmanac: -R("R_MANAO") + R("R_MANAO_Rev") - R("R_MNNH")  = 0.

Mmanb: R("Mmanb_Transp") + R("R_EX_man_e_") - R("R_EX_man_e__Rev") 
 = 0.

Mmane: -R("R_EX_man_e_") + R("R_EX_man_e__Rev") - R("R_MANpts")  =
 0.

Mmelibb: R("Mmelibb_Transp") + R("R_EX_melib_e_") -
 R("R_EX_melib_e__Rev")  = 0.

Mmelibc: -R("R_GALS3") + R("R_MELIBt2")  = 0.

Mmelibe: -R("R_EX_melib_e_") + R("R_EX_melib_e__Rev") -
 R("R_MELIBt2")  = 0.

MmetDb: R("MmetDb_Transp") + R("R_EX_met_D_e_") -
 R("R_EX_met_D_e__Rev")  = 0.

MmetDc: R("R_METDabc")  = 0.

MmetDe: -R("R_EX_met_D_e_") + R("R_EX_met_D_e__Rev") -
 R("R_METDabc")  = 0.

Mmethfc: -R("R_MTHFC") + R("R_MTHFC_Rev") + R("R_MTHFD") -
 R("R_MTHFD_Rev")  = 0.

MmetLb: R("MmetLb_Transp") + R("R_EX_met_L_e_") -
 R("R_EX_met_L_e__Rev")  = 0.

MmetLc: -0.146 R("R_BiomassEcoli") + R("R_METabc") - R("R_METAT") +
 R("R_METS") + R("R_UNK3")  = 0.

MmetLe: -R("R_EX_met_L_e_") + R("R_EX_met_L_e__Rev") - R("R_METabc")
  = 0.

Mmi1pDc: -R("R_MI1PP")  = 0.

Mmicitc: -R("R_MCITL2") + R("R_MCITL2_Rev") + R("R_MICITD")  = 0.

Mmlthfc: R("R_GHMT2") + R("R_GLYCL") - R("R_MOHMT") - R("R_MTHFD") +
 R("R_MTHFD_Rev") - R("R_MTHFR2") - R("R_TMDS")  = 0.

MmmcoaRc: -R("R_MME") + R("R_MME_Rev") + R("R_MMM2")  = 0.

MmmcoaSc: -R("R_MMCD") + R("R_MME") - R("R_MME_Rev")  = 0.

Mmnl1pc: -R("R_M1PD") + R("R_M1PD_Rev") + R("R_MNLpts")  = 0.

Mmnlb: R("Mmnlb_Transp") + R("R_EX_mnl_e_") - R("R_EX_mnl_e__Rev") 
 = 0.

Mmnle: -R("R_EX_mnl_e_") + R("R_EX_mnl_e__Rev") - R("R_MNLpts")  =
 0.

Mmql8c: R("R_ASPO4") + R("R_DHORD5") - R("R_DMSOR1") -
 R("R_DMSOR1e") + R("R_FDH3") - R("R_FRD2") + R("R_G3PD6") + R("R_GLYCTO3") +
 R("R_HYD2") + R("R_L_LACD3") + R("R_MDH3") + R("R_NADH10") + R("R_NADH7") -
 R("R_NO3R2") - R("R_TMAOR1") - R("R_TMAOR1e")  = 0.

Mmqn8c: R("R_AMMQT8_2") - R("R_ASPO4") - R("R_DHORD5") +
 R("R_DMSOR1") + R("R_DMSOR1e") - R("R_FDH3") + R("R_FRD2") - R("R_G3PD6") -
 R("R_GLYCTO3") - R("R_HYD2") - R("R_L_LACD3") - R("R_MDH3") - R("R_NADH10")
 - R("R_NADH7") + R("R_NO3R2") + R("R_TMAOR1") + R("R_TMAOR1e")  = 0.

Mmthgxlc: -R("R_LGTHL") + R("R_MGSA")  = 0.

MmyrsACPc: R("R_AACPS1") + R("R_C140SN") - R("R_EDTXS2") -
 R("R_EDTXS4") - 0.04 R("R_PASYN_EC")  = 0.

MN1aspmdc: R("R_SPMDAT1")  = 0.

Mn8aspmdc: R("R_SPMDAT2")  = 0.

Mna1b: R("Mna1b_Transp") + R("R_EX_na1_e_") - R("R_EX_na1_e__Rev") 
 = 0.

Mna1c: R("R_GLUt4") - R("R_NAt3_1") - 2 R("R_NAt3_15") +
 R("R_NAt3_1_Rev") - R("R_NAt3_2") + R("R_PNTOt4") + R("R_PROt4") + R("R_SERt4") +
 R("R_THRt4")  = 0.

Mna1e: -R("R_EX_na1_e_") + R("R_EX_na1_e__Rev") - R("R_GLUt4") +
 R("R_NAt3_1") + 2 R("R_NAt3_15") - R("R_NAt3_1_Rev") + R("R_NAt3_2") -
 R("R_PNTOt4") - R("R_PROt4") - R("R_SERt4") - R("R_THRt4")  = 0.

Mnacb: R("Mnacb_Transp") + R("R_EX_nac_e_") - R("R_EX_nac_e__Rev") 
 = 0.

Mnacc: R("R_NACUP") - R("R_NAMNPP") + R("R_NNAM") + R("R_NNDMBRT") 
 = 0.

Mnace: -R("R_EX_nac_e_") + R("R_EX_nac_e__Rev") - R("R_NACUP")  = 0.

Mnadb: R("Mnadb_Transp") + R("R_EX_nad_e_") - R("R_EX_nad_e__Rev") 
 = 0.

Mnadc: R("R_2DGLCNRx") + R("R_2DGULRx") + R("R_3HCINNMH") +
 R("R_3HPPPNH") - R("R_ABUTD") - R("R_ACALDi") + 2 R("R_ADHEr") - 2
 R("R_ADHEr_Rev") - R("R_AKGDH") + R("R_ALCD19") - R("R_ALCD19_Rev") -
 R("R_ALDD19x") - R("R_ALDD2x") - R("R_BETALDHx") - 0.00215 R("R_BiomassEcoli") +
 R("R_BSORx") + R("R_CINNDO") - R("R_DHBD") + R("R_DHBD_Rev") -
 R("R_DHCIND") - R("R_DHPPD") + R("R_DKGLCNR2x") + R("R_DMPPS") +
 R("R_DOGULNR") - R("R_E4PD") + R("R_E4PD_Rev") - 6 R("R_FAO1") - 7 R("R_FAO2") - 8
 R("R_FAO3") - R("R_FAO4") - R("R_GAPD") + R("R_GAPD_Rev") -
 R("R_GCALDD") - R("R_GLTPD") + R("R_GLTPD_Rev") - R("R_GLYCDx") - R("R_GLYCL")
 + R("R_GLYCLTDx") - 2 R("R_HISTD") + R("R_HPYRRx") + R("R_IDOND") -
 R("R_IDOND_Rev") - R("R_IMPD") + R("R_IPDPS") - R("R_IPMD") -
 R("R_LCAD") - R("R_LCADi") + R("R_LCAD_Rev") + R("R_LCARS") - R("R_LCARS_Rev") -
 R("R_LDH_D") + R("R_LDH_D_Rev") - R("R_M1PD") + R("R_M1PD_Rev") -
 R("R_MANAO") + R("R_MANAO_Rev") - R("R_MDH") + R("R_MDH_Rev") -
 R("R_ME1") + R("R_MTHFR2") - R("R_NADDP") + R("R_NADH10") + R("R_NADH5") +
 R("R_NADH6") + R("R_NADH7") + R("R_NADH8") + R("R_NADH9") - R("R_NADK") +
 R("R_NADPPPS") + R("R_NADS1") - R("R_NADTRHD") + R("R_NMNAT") + 3
 R("R_NTRIR2x") - R("R_P5CD") - R("R_PDH") - R("R_PDX5PS") - R("R_PERD") +
 R("R_PERD_Rev") - R("R_PGCD") - R("R_PPND") + R("R_PPPNDO") -
 R("R_SBTPD") + R("R_SBTPD_Rev") - R("R_SGSAD") - R("R_SHCHD2") - R("R_SSALx")
 - R("R_TAGURr") + R("R_TAGURr_Rev") + R("R_THD2") - R("R_THRD") +
 R("R_TRSAR") - 2 R("R_UACMAMO") - 2 R("R_UDPGD")  = 0.

Mnade: -R("R_EX_nad_e_") + R("R_EX_nad_e__Rev") - R("R_NADDPe")  =
 0.

Mnadhc: -R("R_2DGLCNRx") - R("R_2DGULRx") - R("R_3HCINNMH") -
 R("R_3HPPPNH") + R("R_ABUTD") + R("R_ACALDi") - 2 R("R_ADHEr") + 2
 R("R_ADHEr_Rev") + R("R_AKGDH") - R("R_ALCD19") + R("R_ALCD19_Rev") +
 R("R_ALDD19x") + R("R_ALDD2x") + R("R_BETALDHx") - 0.00005 R("R_BiomassEcoli")
 - R("R_BSORx") - R("R_CINNDO") + R("R_DHBD") - R("R_DHBD_Rev") +
 R("R_DHCIND") + R("R_DHPPD") - R("R_DKGLCNR2x") - R("R_DMPPS") -
 R("R_DOGULNR") + R("R_E4PD") - R("R_E4PD_Rev") + 6 R("R_FAO1") + 7 R("R_FAO2") +
 8 R("R_FAO3") + R("R_FAO4") + R("R_GAPD") - R("R_GAPD_Rev") +
 R("R_GCALDD") + R("R_GLTPD") - R("R_GLTPD_Rev") + R("R_GLYCDx") +
 R("R_GLYCL") - R("R_GLYCLTDx") + 2 R("R_HISTD") - R("R_HPYRRx") - R("R_IDOND") +
 R("R_IDOND_Rev") + R("R_IMPD") - R("R_IPDPS") + R("R_IPMD") +
 R("R_LCAD") + R("R_LCADi") - R("R_LCAD_Rev") - R("R_LCARS") + R("R_LCARS_Rev")
 + R("R_LDH_D") - R("R_LDH_D_Rev") + R("R_M1PD") - R("R_M1PD_Rev") +
 R("R_MANAO") - R("R_MANAO_Rev") + R("R_MDH") - R("R_MDH_Rev") +
 R("R_ME1") - R("R_MTHFR2") - R("R_NADH10") - R("R_NADH5") - R("R_NADH6") -
 R("R_NADH7") - R("R_NADH8") - R("R_NADH9") + R("R_NADTRHD") - 3
 R("R_NTRIR2x") + R("R_P5CD") + R("R_PDH") + R("R_PDX5PS") + R("R_PERD") -
 R("R_PERD_Rev") + R("R_PGCD") + R("R_PPND") - R("R_PPPNDO") + R("R_SBTPD")
 - R("R_SBTPD_Rev") + R("R_SGSAD") + R("R_SHCHD2") + R("R_SSALx") +
 R("R_TAGURr") - R("R_TAGURr_Rev") - R("R_THD2") + R("R_THRD") -
 R("R_TRSAR") + 2 R("R_UACMAMO") + 2 R("R_UDPGD")  = 0.

Mnadpc: R("R_2DGLCNRy") + R("R_2DGULRy") + R("R_5DGLCNR") -
 R("R_5DGLCNR_Rev") - R("R_AGPR") + R("R_AGPR_Rev") + R("R_APRAUR") -
 R("R_ASAD") + R("R_ASAD_Rev") - R("R_BETALDHy") - 0.00013 R("R_BiomassEcoli") +
 R("R_BSORy") + 10 R("R_C120SN") + 12 R("R_C140SN") + 11 R("R_C141SN")
 + 14 R("R_C160SN") + 13 R("R_C161SN") + 15 R("R_C181SN") +
 R("R_DHDPRy") + R("R_DHFR") - R("R_DHFR_Rev") + R("R_DKGLCNR1") +
 R("R_DKGLCNR2y") + R("R_DPR") + R("R_DXPRIi") - R("R_G3PD2") + R("R_G3PD2_Rev") +
 R("R_G5SD") - R("R_G6PDH2r") + R("R_G6PDH2r_Rev") - R("R_GLUDy") +
 R("R_GLUDy_Rev") + R("R_GLUSy") + R("R_GLUTRR") + R("R_GLYCLTDy") +
 R("R_GMPR") - R("R_GND") + R("R_GOFUCR") + R("R_GTHOr") - R("R_GTHOr_Rev") +
 R("R_HPYRRy") - R("R_HSDy") + R("R_HSDy_Rev") - R("R_ICDHyr") +
 R("R_ICDHyr_Rev") + R("R_IDOND2") + R("R_KARA1i") + R("R_KARA2i") +
 R("R_KAS16") - R("R_ME2") - R("R_MTHFD") + R("R_MTHFD_Rev") + R("R_NADK") -
 R("R_NADPPPS") + R("R_NADTRHD") + R("R_P5CR") + R("R_SHK3Dr") -
 R("R_SHK3Dr_Rev") - R("R_SSALy") - 3 R("R_SULR") + 3 R("R_SULR_Rev") +
 R("R_TDPDRR") - R("R_THD2") + R("R_TRDR") + R("R_UAPGR")  = 0.

Mnadphc: -R("R_2DGLCNRy") - R("R_2DGULRy") - R("R_5DGLCNR") +
 R("R_5DGLCNR_Rev") + R("R_AGPR") - R("R_AGPR_Rev") - R("R_APRAUR") +
 R("R_ASAD") - R("R_ASAD_Rev") + R("R_BETALDHy") - 0.0004 R("R_BiomassEcoli")
 - R("R_BSORy") - 10 R("R_C120SN") - 12 R("R_C140SN") - 11
 R("R_C141SN") - 14 R("R_C160SN") - 13 R("R_C161SN") - 15 R("R_C181SN") -
 R("R_DHDPRy") - R("R_DHFR") + R("R_DHFR_Rev") - R("R_DKGLCNR1") -
 R("R_DKGLCNR2y") - R("R_DPR") - R("R_DXPRIi") + R("R_G3PD2") - R("R_G3PD2_Rev") -
 R("R_G5SD") + R("R_G6PDH2r") - R("R_G6PDH2r_Rev") + R("R_GLUDy") -
 R("R_GLUDy_Rev") - R("R_GLUSy") - R("R_GLUTRR") - R("R_GLYCLTDy") -
 R("R_GMPR") + R("R_GND") - R("R_GOFUCR") - R("R_GTHOr") + R("R_GTHOr_Rev") -
 R("R_HPYRRy") + R("R_HSDy") - R("R_HSDy_Rev") + R("R_ICDHyr") -
 R("R_ICDHyr_Rev") - R("R_IDOND2") - R("R_KARA1i") - R("R_KARA2i") -
 R("R_KAS16") + R("R_ME2") + R("R_MTHFD") - R("R_MTHFD_Rev") - R("R_NADTRHD")
 - R("R_P5CR") - R("R_SHK3Dr") + R("R_SHK3Dr_Rev") + R("R_SSALy") + 3
 R("R_SULR") - 3 R("R_SULR_Rev") - R("R_TDPDRR") + R("R_THD2") -
 R("R_TRDR") - R("R_UAPGR")  = 0.

Mncamc: R("R_NMNN") + R("R_NMNt7") - R("R_NNAM")  = 0.

Mnh4b: R("Mnh4b_Transp") + R("R_EX_nh4_e_") - R("R_EX_nh4_e__Rev") 
 = 0.

Mnh4c: R("R_ADA") + R("R_ADD") + R("R_ASNN") - R("R_ASNS2") +
 R("R_ASPT") - R("R_CBMK") + R("R_CSND") + R("R_CYNTAH") + R("R_CYSDS") +
 R("R_CYSTL") + R("R_CYTD") + R("R_DAAD") + R("R_DADA") + R("R_DCTPD") +
 R("R_DCYTD") + R("R_DHPPDA2") + R("R_ETHAAL") + R("R_G6PDA") -
 R("R_GLNS") + R("R_GLUDy") - R("R_GLUDy_Rev") + R("R_GLUN") + R("R_GLYCL") +
 R("R_GMPR") + R("R_GUAD") + 4 R("R_HMBS") - R("R_NADS1") + R("R_NH4t")
 - R("R_NH4t_Rev") + R("R_NMNDA") + R("R_NNAM") + R("R_NTRIR2x") +
 R("R_PEAMNO") + R("R_PYAM5PO") + 2 R("R_SADH") + R("R_SERD_D") +
 R("R_SERD_L") + R("R_THRD_L") + R("R_TRPAS2") - R("R_TRPAS2_Rev") + 2
 R("R_UGLYCH")  = 0.

Mnh4e: -R("R_EX_nh4_e_") + R("R_EX_nh4_e__Rev") - R("R_NH4t") +
 R("R_NH4t_Rev")  = 0.

Mnicrntc: R("R_NAMNPP") + R("R_NMNDA") - R("R_NNAT") -
 R("R_NNDMBRT") + R("R_NNDPR")  = 0.

Mnmnb: R("Mnmnb_Transp") + R("R_EX_nmn_e_") - R("R_EX_nmn_e__Rev") 
 = 0.

Mnmnc: R("R_NADDP") - R("R_NMNAT") - R("R_NMNDA") - R("R_NMNN") +
 R("R_NMNP")  = 0.

Mnmne: -R("R_EX_nmn_e_") + R("R_EX_nmn_e__Rev") + R("R_NADDPe") -
 R("R_NMNP") - R("R_NMNt7")  = 0.

Mno2b: R("Mno2b_Transp") + R("R_EX_no2_e_") - R("R_EX_no2_e__Rev") 
 = 0.

Mno2c: R("R_NO2t2r") - R("R_NO2t2r_Rev") + R("R_NO3R1") +
 R("R_NO3R2") - R("R_NO3t7") - R("R_NTRIR2x")  = 0.

Mno2e: -R("R_EX_no2_e_") + R("R_EX_no2_e__Rev") - R("R_NO2t2r") +
 R("R_NO2t2r_Rev") + R("R_NO3t7")  = 0.

Mno3b: R("Mno3b_Transp") + R("R_EX_no3_e_") - R("R_EX_no3_e__Rev") 
 = 0.

Mno3c: -R("R_NO3R1") - R("R_NO3R2") + R("R_NO3t7")  = 0.

Mno3e: -R("R_EX_no3_e_") + R("R_EX_no3_e__Rev") - R("R_NO3t7")  = 0.

Mo2b: R("Mo2b_Transp") + R("R_EX_o2_e_") - R("R_EX_o2_e__Rev")  = 0.

Mo2c: -R("R_3HCINNMH") - R("R_3HPPPNH") - R("R_ASPO6") + R("R_CAT")
 - R("R_CINNDO") - R("R_CPPPGO") - 0.5 R("R_CYTBD") - 0.5 R("R_CYTBO3")
 - R("R_DHCINDO") - R("R_DKMPPD") - R("R_HPPPNDO") + R("R_O2t") -
 R("R_O2t_Rev") - 0.5 R("R_OMMBLHX") - 0.5 R("R_OMPHHX") - 0.5 R("R_OPHHX")
 - R("R_PDX5PO") + R("R_PDX5PO_Rev") - R("R_PEAMNO") - 1.5
 R("R_PPPGO") - R("R_PPPNDO") - R("R_PYAM5PO") + R("R_SPODM") - R("R_TAUDO")  =
 0.

Mo2e: -R("R_EX_o2_e_") + R("R_EX_o2_e__Rev") - R("R_O2t") +
 R("R_O2t_Rev")  = 0.

Mo2SuperOxidec: -2 R("R_SPODM")  = 0.

Moaac: R("R_ASPTA") - R("R_ASPTA_Rev") + R("R_CITL") - R("R_CS") -
 R("R_MCITS") + R("R_MDH") + R("R_MDH2") + R("R_MDH3") - R("R_MDH_Rev")
 + R("R_PPC") - R("R_PPCK") + R("R_TARTD")  = 0.

Mocdcab: R("Mocdcab_Transp") + R("R_EX_ocdca_e_") -
 R("R_EX_ocdca_e__Rev")  = 0.

Mocdcac: -R("R_FAO3") + R("R_OCDCAt2")  = 0.

Mocdcae: -R("R_EX_ocdca_e_") + R("R_EX_ocdca_e__Rev") -
 R("R_OCDCAt2")  = 0.

Mocdceac: -R("R_AACPS5") + 0.5 R("R_LPLIPA1") + 0.5 R("R_LPLIPA2") +
 0.5 R("R_LPLIPA3") + 0.5 R("R_PLIPA1") + 0.5 R("R_PLIPA2") + 0.5
 R("R_PLIPA3")  = 0.

Moctdpc: -R("R_DHNAOT") - R("R_HBZOPT") + R("R_OCTDPS")  = 0.

MocteACPc: R("R_AACPS5") + R("R_C181SN") - R("R_PASYN_EC")  = 0.

Mohpbc: -R("R_OHPBAT") + R("R_OHPBAT_Rev") + R("R_PERD") -
 R("R_PERD_Rev")  = 0.

Mop4enc: R("R_HKNDDH") + R("R_HKNTDH") - R("R_OP4ENH")  = 0.

Mornb: R("Mornb_Transp") + R("R_EX_orn_e_") - R("R_EX_orn_e__Rev") 
 = 0.

Mornc: R("R_ACODA") - R("R_ARGORNt7") + R("R_ARGORNt7_Rev") -
 R("R_OCBT") + R("R_OCBT_Rev") + R("R_ORNabc") - R("R_ORNDC") - R("R_ORNTA")
 - R("R_PTRCORNt7") + R("R_PTRCORNt7_Rev")  = 0.

Morne: R("R_ARGORNt7") - R("R_ARGORNt7_Rev") - R("R_EX_orn_e_") +
 R("R_EX_orn_e__Rev") - R("R_ORNabc") + R("R_PTRCORNt7") -
 R("R_PTRCORNt7_Rev")  = 0.

Morot5pc: -R("R_OMPDC") - R("R_ORPT") + R("R_ORPT_Rev")  = 0.

Morotc: R("R_DHORD2") + R("R_DHORD5") + R("R_ORPT") -
 R("R_ORPT_Rev")  = 0.

Mpacaldc: -R("R_ALDD19x") + R("R_PEAMNO")  = 0.

Mpacc: R("R_ALDD19x") - R("R_PACCOAL")  = 0.

MpaECc: 0.02 R("R_CDAPPA_EC") + 0.02 R("R_DAGK_EC") - 0.02
 R("R_DASYN_EC") + 0.02 R("R_DASYN_EC_Rev") - 0.02 R("R_PAPA_EC") + 0.02
 R("R_PASYN_EC")  = 0.

MpalmACPc: R("R_AACPS3") + R("R_C160SN") - 0.72 R("R_PASYN_EC")  =
 0.

Mpan4pc: R("R_PPCDC") - R("R_PTPATi")  = 0.

MpantRc: R("R_DPR") - R("R_PANTS")  = 0.

Mpapc: R("R_ACPS1") - R("R_BPNT") + R("R_PAPSR")  = 0.

Mpapsc: R("R_ADSK") - R("R_PAPSR")  = 0.

MpcECc: -0.02 R("R_PLIPA3")  = 0.

Mpdx5pc: -R("R_PDX5PO") + R("R_PDX5PO_Rev") + R("R_PDX5PS") -
 R("R_PDXPP") + R("R_PYDXNK")  = 0.

Mpeamnc: -R("R_PEAMNO")  = 0.

MpeECc: -0.001935 R("R_BiomassEcoli") - 0.02 R("R_PEPT_EC") + 0.02
 R("R_PEPT_EC_Rev") - 0.02 R("R_PLIPA2") + 0.02 R("R_PSD_EC")  = 0.

Mpepc: -R("R_ACGApts") - R("R_ACMANApts") - R("R_DDPA") -
 R("R_DHAPT") + R("R_ENO") - R("R_ENO_Rev") - R("R_FRUpts") - R("R_FRUpts2") -
 R("R_GALTpts") - R("R_GAMpts") - R("R_GLCpts") - R("R_KDOPS") -
 R("R_MALTpts") - R("R_MANpts") - R("R_MNLpts") - R("R_PPC") + R("R_PPCK") +
 R("R_PPS") - R("R_PSCVT") + R("R_PSCVT_Rev") - R("R_PYK") -
 R("R_SBTpts") - R("R_SUCpts") - R("R_TREpts") - R("R_UAGCVT")  = 0.

MpeptidoECc: -0.0276 R("R_BiomassEcoli") + R("R_PPTGS")  = 0.

MpgECc: -0.000464 R("R_BiomassEcoli") - 0.04 R("R_CLPNS_EC") + 0.04
 R("R_CLPNS_EC_Rev") - 0.02 R("R_LPLIPA4") - 0.02 R("R_LPLIPA5") - 0.02
 R("R_LPLIPA6") + 0.02 R("R_PGPP_EC") - 0.02 R("R_PLIPA1")  = 0.

MpgpECc: -0.02 R("R_PGPP_EC") + 0.02 R("R_PGSA_EC") - 0.02
 R("R_PGSA_EC_Rev")  = 0.

Mphaccoac: R("R_PACCOAL")  = 0.

MpheLb: R("MpheLb_Transp") + R("R_EX_phe_L_e_") -
 R("R_EX_phe_L_e__Rev")  = 0.

MpheLc: -0.176 R("R_BiomassEcoli") + R("R_PHEt2r") -
 R("R_PHEt2r_Rev") - R("R_PHETA1") + R("R_PHETA1_Rev")  = 0.

MpheLe: -R("R_EX_phe_L_e_") + R("R_EX_phe_L_e__Rev") - R("R_PHEt2r")
 + R("R_PHEt2r_Rev")  = 0.

Mphemec: R("R_FCLT") - R("R_HEMEOS")  = 0.

Mphomc: R("R_HSK") - R("R_THRS")  = 0.

Mphpyrc: R("R_PHETA1") - R("R_PHETA1_Rev") + R("R_PPNDH")  = 0.

Mphthrc: -R("R_4HTHRS") + R("R_OHPBAT") - R("R_OHPBAT_Rev") -
 R("R_PDX5PS")  = 0.

Mpib: R("Mpib_Transp") + R("R_EX_pi_e_") - R("R_EX_pi_e__Rev")  = 0.

Mpic: R("R_4HTHRS") + R("R_ACCOACr") - R("R_ACCOACr_Rev") +
 R("R_ACCOAL") + R("R_ADSS") - R("R_AGPR") + R("R_AGPR_Rev") + R("R_AIRC2") +
 R("R_ALAabc") + R("R_ALAALAr") - R("R_ALAALAr_Rev") + R("R_ARBabc") +
 R("R_ARGabc") - R("R_ASAD") + R("R_ASAD_Rev") + R("R_ASNabc") +
 R("R_ASPabc") + R("R_ASPCT") + R("R_ATPM") - R("R_ATPS4r") +
 R("R_ATPS4r_Rev") + 45.5628 R("R_BiomassEcoli") + R("R_BPNT") + R("R_CBIAT_DELETE") -
 R("R_CBIAT_DELETE_Rev") + R("R_CBL1abc") + R("R_CBLAT_DELETE") -
 R("R_CBLAT_DELETE_Rev") + R("R_CBPS") + R("R_CHLabc") + R("R_CHORS") +
 R("R_CTPS2") + R("R_CYSabc") + R("R_DAPabc") + R("R_DBTSr") -
 R("R_DBTSr_Rev") + R("R_DDPA") + R("R_DHFS") + R("R_DHQS") + R("R_DKMPPD") +
 R("R_DKMPPD2") + R("R_DNMPPA") - R("R_DURIPP") + R("R_DURIPP_Rev") +
 R("R_FBP") + R("R_FE2abc") - 2 R("R_FUCPt6_2") + R("R_G1PP") + R("R_G5SD") -
 2 R("R_G6Pt6_2") + R("R_GALabc") - R("R_GAPD") + R("R_GAPD_Rev") +
 R("R_GART") - R("R_GLCP") + R("R_GLNabc") + R("R_GLNS") + R("R_GLUabc")
 + R("R_GLUCYS") + R("R_GLYBabc") + R("R_GLYC3Pabc") - R("R_GLYC3Pt6")
 + R("R_GMHEPPA") + R("R_GSPMDS") + R("R_GTHS") + R("R_HISabc") +
 R("R_HISTP") + R("R_HYPOE") + R("R_ILEabc") + R("R_Kabc") + R("R_KDOPP") +
 R("R_KDOPS") + R("R_LEUabc") + R("R_LYSabc") + R("R_MALTabc") +
 R("R_MALTHXabc") + R("R_MALTPTabc") + R("R_MALTTRabc") + R("R_MALTTTRabc") +
 R("R_MAN1PT2") - 2 R("R_MAN6Pt6_2") + R("R_METabc") + R("R_METAT") +
 R("R_METDabc") + R("R_MGSA") + R("R_MI1PP") - R("R_MLTP1") +
 R("R_MLTP1_Rev") - R("R_MLTP2") + R("R_MLTP2_Rev") - R("R_MLTP3") +
 R("R_MLTP3_Rev") + R("R_NADPPPS") + R("R_NAMNPP") + R("R_NTD1") + R("R_NTD10") +
 R("R_NTD11") + R("R_NTD2") + R("R_NTD3") + R("R_NTD4") + R("R_NTD5") +
 R("R_NTD6") + R("R_NTD7") + R("R_NTD8") + R("R_NTD9") + R("R_OCBT") -
 R("R_OCBT_Rev") + R("R_ORNabc") + R("R_PAPA_EC") + R("R_PDX5PS") +
 R("R_PDXPP") + R("R_PGLYCP") + R("R_PGPP_EC") + 2 R("R_PIabc") +
 R("R_PIt2r") - R("R_PIt2r_Rev") + R("R_PMDPHT") + 2 R("R_PPA") + R("R_PPC") +
 R("R_PPS") + R("R_PRAGSr") - R("R_PRAGSr_Rev") + R("R_PRAIS") +
 R("R_PRASCS") - R("R_PRASCS_Rev") + R("R_PRFGS") + R("R_PROabc") +
 R("R_PSCVT") - R("R_PSCVT_Rev") + R("R_PSP_L") - R("R_PTA2") - R("R_PTAr") +
 R("R_PTAr_Rev") + R("R_PTRCabc") - R("R_PUNP1") + R("R_PUNP1_Rev") -
 R("R_PUNP2") + R("R_PUNP2_Rev") - R("R_PUNP3") + R("R_PUNP3_Rev") -
 R("R_PUNP4") + R("R_PUNP4_Rev") - R("R_PUNP5") + R("R_PUNP5_Rev") -
 R("R_PUNP6") + R("R_PUNP6_Rev") - R("R_PUNP7") + R("R_PUNP7_Rev") +
 R("R_PYDXPP") - R("R_PYNP2r") + R("R_PYNP2r_Rev") + R("R_QULNS") +
 R("R_RBFSa") + R("R_RIBabc") + R("R_RZ5PP") + R("R_SADT2") + R("R_SELNPS") +
 R("R_SPMDabc") + R("R_SUCCabc") + R("R_SUCOAS") - R("R_SUCOAS_Rev") +
 R("R_SULabc") + R("R_TAURabc") + R("R_THMabc") + R("R_THRabc") +
 R("R_THRS") - R("R_TMDPP") + R("R_TMDPP_Rev") + R("R_TRE6PP") + R("R_TSULabc")
 + R("R_UAAGDS") + R("R_UAGCVT") + R("R_UAMAGS") + R("R_UAMAS") +
 R("R_UDCPDP") + R("R_UGMDDS") + R("R_VALabc") + R("R_XYLabc")  = 0.

Mpie: -R("R_EX_pi_e_") + R("R_EX_pi_e__Rev") + 2 R("R_FUCPt6_2") + 2
 R("R_G6Pt6_2") + R("R_GLYC3Pt6") + 2 R("R_MAN6Pt6_2") - R("R_PIabc")
 - R("R_PIt2r") + R("R_PIt2r_Rev")  = 0.

Mpmcoac: -R("R_AOXSr") + R("R_AOXSr_Rev")  = 0.

MpntoRb: R("MpntoRb_Transp") + R("R_EX_pnto_R_e_") -
 R("R_EX_pnto_R_e__Rev")  = 0.

MpntoRc: R("R_PANTS") - R("R_PNTK") + R("R_PNTOt4")  = 0.

MpntoRe: -R("R_EX_pnto_R_e_") + R("R_EX_pnto_R_e__Rev") -
 R("R_PNTOt4")  = 0.

Mppac: -R("R_ACCOAL") + R("R_PPAKr") - R("R_PPAKr_Rev") +
 R("R_PPCSCT")  = 0.

Mppapc: -R("R_PPAKr") + R("R_PPAKr_Rev") + R("R_PTA2")  = 0.

Mppbngc: -4 R("R_HMBS") + R("R_PPBNGS")  = 0.

Mppcoac: R("R_ACCOAL") - R("R_MCITS") + R("R_MMCD") + R("R_OBTFL") -
 R("R_PPCSCT") - R("R_PTA2")  = 0.

Mpphnc: R("R_CHORM") - R("R_PPND") - R("R_PPNDH")  = 0.

Mppic: R("R_AACPS1") + R("R_AACPS2") + R("R_AACPS3") + R("R_AACPS4")
 + R("R_AACPS5") + R("R_ACBIPGT") + R("R_ACS") + R("R_ADNCYC") +
 R("R_ADPT") + R("R_ANPRT") + R("R_ARGSS") + R("R_ASNS1") + R("R_ASNS2") +
 R("R_ATPPRT") + 0.7302 R("R_BiomassEcoli") + R("R_CBIAT_DELETE") -
 R("R_CBIAT_DELETE_Rev") + R("R_CBLAT_DELETE") - R("R_CBLAT_DELETE_Rev") +
 R("R_DASYN_EC") - R("R_DASYN_EC_Rev") + R("R_DHBSr") -
 R("R_DHBSr_Rev") + R("R_DHNAOT") + R("R_DHPS2") + R("R_DMATT") + R("R_DNTPPA") +
 R("R_DUTPDP") + R("R_FAO1") + R("R_FAO2") + R("R_FAO3") + R("R_FMNAT") +
 R("R_G1PTT") + R("R_GALU") + R("R_GALUi") - R("R_GALU_Rev") +
 R("R_GLGC") + R("R_GLUPRT") + R("R_GLUTRS") + R("R_GMHEPAT") + R("R_GMPS2") +
 R("R_GRTT") + R("R_GTPCII2") + R("R_GUAPRT") + R("R_HBZOPT") +
 R("R_HEMEOS") + R("R_HXPRT") + R("R_KDOCT2") + R("R_MEPCT") + R("R_METAT") +
 R("R_NADS1") + R("R_NAMNPP") + R("R_NMNAT") + R("R_NNAT") +
 R("R_NNDPR") + R("R_NTPP1") + R("R_NTPP2") + R("R_NTPP3") + R("R_NTPP4") +
 R("R_NTPP5") + R("R_NTPP6") + R("R_NTPP7") + R("R_NTPP8") + 5 R("R_OCTDPS")
 - R("R_ORPT") + R("R_ORPT_Rev") + R("R_PACCOAL") + R("R_PANTS") -
 R("R_PPA") + R("R_PPNCL2") + R("R_PRATPP") + R("R_PTPATi") + R("R_SADT2")
 + R("R_SERASr") - R("R_SERASr_Rev") + R("R_SUCBZL") + R("R_THZPSN") +
 R("R_TMPPP") + R("R_UAGDP") + 8 R("R_UDCPDPS") + R("R_UPPRT") +
 R("R_XPPT")  = 0.

Mppp9c: -R("R_FCLT") + R("R_PPPGO")  = 0.

Mpppg9c: R("R_CPPPGO") - R("R_PPPGO")  = 0.

Mpppic: R("R_NTPTP1") + R("R_NTPTP2")  = 0.

Mpppnb: R("Mpppnb_Transp") + R("R_EX_pppn_e_") -
 R("R_EX_pppn_e__Rev")  = 0.

Mpppnc: -R("R_PPPNDO") + R("R_PPPNt2r") - R("R_PPPNt2r_Rev")  = 0.

Mpppne: -R("R_EX_pppn_e_") + R("R_EX_pppn_e__Rev") - R("R_PPPNt2r")
 + R("R_PPPNt2r_Rev")  = 0.

Mpramc: R("R_GLUPRT") - R("R_PRAGSr") + R("R_PRAGSr_Rev")  = 0.

Mpranc: R("R_ANPRT") - R("R_PRAIi")  = 0.

Mprbampc: -R("R_PRAMPC") + R("R_PRATPP")  = 0.

Mprbatpc: R("R_ATPPRT") - R("R_PRATPP")  = 0.

Mprfpc: R("R_PRAMPC") - R("R_PRMICIi")  = 0.

Mprlpc: -R("R_IG3PS") + R("R_PRMICIi")  = 0.

MproLb: R("MproLb_Transp") + R("R_EX_pro_L_e_") -
 R("R_EX_pro_L_e__Rev")  = 0.

MproLc: -0.21 R("R_BiomassEcoli") + R("R_P5CR") + R("R_PROabc") -
 R("R_PROD2") + R("R_PROt2r") - R("R_PROt2r_Rev") + R("R_PROt4")  = 0.

MproLe: -R("R_EX_pro_L_e_") + R("R_EX_pro_L_e__Rev") - R("R_PROabc")
 - R("R_PROt2r") + R("R_PROt2r_Rev") - R("R_PROt4")  = 0.

Mprppc: -R("R_ADPT") - R("R_ANPRT") - R("R_ATPPRT") - R("R_GLUPRT")
 - R("R_GUAPRT") - R("R_HXPRT") - R("R_NAMNPP") - R("R_NNDPR") +
 R("R_ORPT") - R("R_ORPT_Rev") + R("R_PRPPS") - R("R_PRPPS_Rev") -
 R("R_UPPRT") - R("R_XPPT")  = 0.

MpsECc: -0.000052 R("R_BiomassEcoli") - 0.02 R("R_PSD_EC") + 0.02
 R("R_PSSA_EC") - 0.02 R("R_PSSA_EC_Rev")  = 0.

MpserLc: R("R_PSERT") - R("R_PSP_L")  = 0.

Mptrcb: R("Mptrcb_Transp") + R("R_EX_ptrc_e_") -
 R("R_EX_ptrc_e__Rev")  = 0.

Mptrcc: R("R_AGMT") - 0.035 R("R_BiomassEcoli") + R("R_ORNDC") +
 R("R_PTRCabc") + R("R_PTRCORNt7") - R("R_PTRCORNt7_Rev") + R("R_PTRCt2r")
 - R("R_PTRCt2r_Rev") - R("R_PTRCTA") - R("R_SPMS")  = 0.

Mptrce: -R("R_EX_ptrc_e_") + R("R_EX_ptrc_e__Rev") - R("R_PTRCabc")
 - R("R_PTRCORNt7") + R("R_PTRCORNt7_Rev") - R("R_PTRCt2r") +
 R("R_PTRCt2r_Rev")  = 0.

Mpyam5pc: R("R_ALATA_D2") + R("R_ALATA_L2") - R("R_HYPOE") -
 R("R_PYAM5PO") + R("R_PYDAMK")  = 0.

Mpydamc: R("R_HYPOE") - R("R_PYDAMK")  = 0.

Mpydx5pc: -R("R_ALATA_D2") - R("R_ALATA_L2") + R("R_PDX5PO") -
 R("R_PDX5PO_Rev") + R("R_PYAM5PO") + R("R_PYDXK") - R("R_PYDXPP")  = 0.

Mpydxc: -R("R_PYDXK") + R("R_PYDXPP")  = 0.

Mpydxnc: R("R_PDXPP") - R("R_PYDXNK")  = 0.

Mpyrb: R("Mpyrb_Transp") + R("R_EX_pyr_e_") - R("R_EX_pyr_e__Rev") 
 = 0.

Mpyrc: R("R_ACGApts") - R("R_ACHBS") - 2 R("R_ACLS") +
 R("R_ACMANApts") + R("R_ACNML") + R("R_ADCL") + R("R_ALATA_D2") + R("R_ALATA_L") +
 R("R_ALATA_L2") - R("R_ALATA_L_Rev") + R("R_ANS") + R("R_CHRPL") +
 R("R_CYSDS") + R("R_CYSTL") + R("R_DAAD") + R("R_DDPGALA") -
 R("R_DDPGALA_Rev") + R("R_DHAPT") - R("R_DHDPS") - R("R_DXPS") + R("R_EDA") +
 R("R_FRUpts") + R("R_FRUpts2") + R("R_GALTpts") + R("R_GAMpts") +
 R("R_GLCpts") + R("R_GLCRAL") + R("R_HOPNTAL") + R("R_ICHORT") + R("R_LDH_D")
 + R("R_LDH_D2") - R("R_LDH_D_Rev") + R("R_L_LACD2") + R("R_L_LACD3") +
 R("R_MALTpts") + R("R_MANpts") + R("R_MCITL2") - R("R_MCITL2_Rev") +
 R("R_ME1") + R("R_ME2") + R("R_MNLpts") - R("R_PDH") - R("R_PFL") -
 R("R_POX") - R("R_PPS") + R("R_PYK") + R("R_PYRt2r") - R("R_PYRt2r_Rev")
 + R("R_SBTpts") + R("R_SERD_D") + R("R_SERD_L") + R("R_SHCHCS2") +
 R("R_SUCpts") + R("R_TREpts") + R("R_TRPAS2") - R("R_TRPAS2_Rev") +
 R("R_VPAMT")  = 0.

Mpyre: -R("R_EX_pyr_e_") + R("R_EX_pyr_e__Rev") - R("R_PYRt2r") +
 R("R_PYRt2r_Rev")  = 0.

Mq8c: -R("R_ASPO3") + R("R_CYTBD") + R("R_CYTBO3") - R("R_DHORD2") -
 R("R_FDH2") - R("R_G3PD5") - R("R_GLCDe") - R("R_GLYCTO2") -
 R("R_HYD1") - R("R_LDH_D2") - R("R_L_LACD2") - R("R_MDH2") - R("R_NADH5") -
 R("R_NADH6") + R("R_NO3R1") - R("R_POX") - R("R_SUCD4") +
 R("R_SUCD4_Rev")  = 0.

Mq8h2c: R("R_ASPO3") - R("R_CYTBD") - R("R_CYTBO3") + R("R_DHORD2")
 + R("R_DMQMT") + R("R_FDH2") + R("R_G3PD5") + R("R_GLCDe") +
 R("R_GLYCTO2") + R("R_HYD1") + R("R_LDH_D2") + R("R_L_LACD2") + R("R_MDH2") +
 R("R_NADH5") + R("R_NADH6") - R("R_NO3R1") + R("R_POX") + R("R_SUCD4")
 - R("R_SUCD4_Rev")  = 0.

Mqulnc: -R("R_NNDPR") + R("R_QULNS")  = 0.

Mr1pc: -R("R_PPM") + R("R_PPM_Rev") + R("R_PUNP1") -
 R("R_PUNP1_Rev") + R("R_PUNP3") - R("R_PUNP3_Rev") + R("R_PUNP5") - R("R_PUNP5_Rev")
 + R("R_PUNP7") - R("R_PUNP7_Rev") + R("R_PYNP2r") - R("R_PYNP2r_Rev") 
 = 0.

Mr5pc: R("R_AMPN") + R("R_CMPN") + R("R_NMNN") + R("R_NMNt7") +
 R("R_PPM") - R("R_PPM_Rev") - R("R_PRPPS") + R("R_PRPPS_Rev") + R("R_RBK")
 - R("R_RPI") + R("R_RPI_Rev") - R("R_TKT1") + R("R_TKT1_Rev")  = 0.

MrblLc: R("R_ARAI") - R("R_ARAI_Rev") - R("R_RBK_L1")  = 0.

Mrdmbzic: -R("R_ADOCBLS") + R("R_RZ5PP")  = 0.

Mrhcysc: R("R_AHCYSNS") - R("R_RHCCE")  = 0.

MribDb: R("MribDb_Transp") + R("R_EX_rib_D_e_") -
 R("R_EX_rib_D_e__Rev")  = 0.

MribDc: -R("R_RBK") + R("R_RIBabc")  = 0.

MribDe: -R("R_EX_rib_D_e_") + R("R_EX_rib_D_e__Rev") - R("R_RIBabc")
  = 0.

Mribflvc: -R("R_RBFK") + R("R_RBFSb")  = 0.

Mrml1pc: R("R_RMK") - R("R_RMPA") + R("R_RMPA_Rev")  = 0.

Mrmlc: R("R_RMI") - R("R_RMI_Rev") - R("R_RMK")  = 0.

Mrmnb: R("Mrmnb_Transp") + R("R_EX_rmn_e_") - R("R_EX_rmn_e__Rev") 
 = 0.

Mrmnc: -R("R_RMI") + R("R_RMI_Rev") + R("R_RMNt")  = 0.

Mrmne: -R("R_EX_rmn_e_") + R("R_EX_rmn_e__Rev") - R("R_RMNt")  = 0.

Mru5pDc: -R("R_A5PISO") + R("R_A5PISO_Rev") - R("R_DB4PS") +
 R("R_GND") - R("R_RPE") + R("R_RPE_Rev") + R("R_RPI") - R("R_RPI_Rev")  = 0.

Mru5pLc: R("R_RBK_L1") - R("R_RBP4E") + R("R_RBP4E_Rev") +
 R("R_X5PL3E")  = 0.

Ms7pc: -R("R_S7PI") - R("R_TALA") + R("R_TALA_Rev") + R("R_TKT1") -
 R("R_TKT1_Rev")  = 0.

Msbt6pc: -R("R_SBTPD") + R("R_SBTPD_Rev") + R("R_SBTpts")  = 0.

MsbtDb: R("MsbtDb_Transp") + R("R_EX_sbt_D_e_") -
 R("R_EX_sbt_D_e__Rev")  = 0.

MsbtDe: -R("R_EX_sbt_D_e_") + R("R_EX_sbt_D_e__Rev") - R("R_SBTpts")
  = 0.

Msbzcoac: -R("R_NPHS") + R("R_SUCBZL")  = 0.

Msclc: R("R_SHCHD2") - R("R_SHCHF")  = 0.

Mselnc: -R("R_SELNPS")  = 0.

Mselnpc: R("R_SELNPS")  = 0.

Mserampc: -3 R("R_ENTCS") + R("R_SERASr") - R("R_SERASr_Rev")  = 0.

MserDb: R("MserDb_Transp") + R("R_EX_ser_D_e_") -
 R("R_EX_ser_D_e__Rev")  = 0.

MserDc: R("R_DSERt2r") - R("R_DSERt2r_Rev") - R("R_SERD_D")  = 0.

MserDe: -R("R_DSERt2r") + R("R_DSERt2r_Rev") - R("R_EX_ser_D_e_") +
 R("R_EX_ser_D_e__Rev")  = 0.

MserLb: R("MserLb_Transp") + R("R_EX_ser_L_e_") -
 R("R_EX_ser_L_e__Rev")  = 0.

MserLc: -0.205 R("R_BiomassEcoli") - R("R_GHMT2") + R("R_GPDDA3") +
 R("R_PSP_L") - R("R_PSSA_EC") + R("R_PSSA_EC_Rev") - R("R_SERASr") +
 R("R_SERASr_Rev") - R("R_SERAT") + R("R_SERAT_Rev") - R("R_SERD_L") +
 R("R_SERt2r") - R("R_SERt2r_Rev") + R("R_SERt4") - R("R_TRPS1") -
 R("R_TRPS2")  = 0.

MserLe: -R("R_EX_ser_L_e_") + R("R_EX_ser_L_e__Rev") - R("R_SERt2r")
 + R("R_SERt2r_Rev") - R("R_SERt4")  = 0.

Mshemec: R("R_SHCHF")  = 0.

Mskm5pc: -R("R_PSCVT") + R("R_PSCVT_Rev") + R("R_SHKK")  = 0.

Mskmc: R("R_SHK3Dr") - R("R_SHK3Dr_Rev") - R("R_SHKK")  = 0.

Msl26dac: -R("R_SDPDS") - R("R_SDPTA") + R("R_SDPTA_Rev")  = 0.

Msl2a6oc: R("R_SDPTA") - R("R_SDPTA_Rev") + R("R_THDPS")  = 0.

Mso3c: R("R_CYANST") + R("R_PAPSR") + R("R_SULR") - R("R_SULR_Rev")
 + R("R_TAUDO")  = 0.

Mso4b: R("Mso4b_Transp") + R("R_EX_so4_e_") - R("R_EX_so4_e__Rev") 
 = 0.

Mso4c: -R("R_SADT2") + R("R_SULabc")  = 0.

Mso4e: -R("R_EX_so4_e_") + R("R_EX_so4_e__Rev") - R("R_SULabc")  =
 0.

Mspmdb: R("Mspmdb_Transp") + R("R_EX_spmd_e_") -
 R("R_EX_spmd_e__Rev")  = 0.

Mspmdc: -0.007 R("R_BiomassEcoli") + R("R_GSPMDA") - R("R_GSPMDS") +
 R("R_SPMDabc") - R("R_SPMDAT1") - R("R_SPMDAT2") + R("R_SPMS")  = 0.

Mspmde: -R("R_EX_spmd_e_") + R("R_EX_spmd_e__Rev") - R("R_SPMDabc") 
 = 0.

Mssaltppc: R("R_OXGDC2") - R("R_SHCHCS2")  = 0.

Msuc6pc: -R("R_FFSD") + R("R_SUCpts")  = 0.

Msucargc: R("R_AST") - R("R_SADH")  = 0.

Msucbzc: -R("R_SUCBZL") + R("R_SUCBZS")  = 0.

Msuccb: R("Msuccb_Transp") + R("R_EX_succ_e_") -
 R("R_EX_succ_e__Rev")  = 0.

Msuccc: R("R_ASPO5") - R("R_CITt7") + R("R_FRD2") + R("R_FRD3") +
 R("R_HKNDDH") + R("R_ICL") + R("R_MCITL2") - R("R_MCITL2_Rev") -
 R("R_PPCSCT") + R("R_SDPDS") + R("R_SGDS") + R("R_SHSL1") + R("R_SSALx") +
 R("R_SSALy") + R("R_SUCCabc") + R("R_SUCCt2_2") + R("R_SUCCt2_3") -
 R("R_SUCCt2b") - R("R_SUCD1i") - R("R_SUCFUMt") + R("R_SUCFUMt_Rev") -
 R("R_SUCOAS") + R("R_SUCOAS_Rev") - R("R_TARTRt7") + R("R_TARTRt7_Rev") +
 R("R_TAUDO")  = 0.

Msucce: R("R_CITt7") - R("R_EX_succ_e_") + R("R_EX_succ_e__Rev") -
 R("R_SUCCabc") - R("R_SUCCt2_2") - R("R_SUCCt2_3") + R("R_SUCCt2b") +
 R("R_SUCFUMt") - R("R_SUCFUMt_Rev") + R("R_TARTRt7") -
 R("R_TARTRt7_Rev")  = 0.

Msuccoac: R("R_AKGDH") - R("R_AST") - 3.e-6 R("R_BiomassEcoli") -
 R("R_HSST") - R("R_MMM2") + R("R_PPCSCT") + R("R_SUCOAS") -
 R("R_SUCOAS_Rev") - R("R_THDPS")  = 0.

Msucgluc: -R("R_SGDS") + R("R_SGSAD")  = 0.

Msucgsac: -R("R_SGSAD") + R("R_SOTA")  = 0.

Msuchmsc: R("R_HSST") - R("R_SHSL1")  = 0.

Msucornc: R("R_SADH") - R("R_SOTA")  = 0.

Msucrb: R("Msucrb_Transp") + R("R_EX_sucr_e_") -
 R("R_EX_sucr_e__Rev")  = 0.

Msucre: -R("R_EX_sucr_e_") + R("R_EX_sucr_e__Rev") - R("R_SUCpts") 
 = 0.

Msucsalc: R("R_ABTA") - R("R_SSALx") - R("R_SSALy")  = 0.

Mtag6pDc: R("R_GLTPD") - R("R_GLTPD_Rev") - R("R_PFK_2")  = 0.

MtagdpDc: R("R_PFK_2") - R("R_TGBPA") + R("R_TGBPA_Rev")  = 0.

Mtagurc: R("R_GUI2") - R("R_GUI2_Rev") + R("R_TAGURr") -
 R("R_TAGURr_Rev")  = 0.

MtartrLb: R("MtartrLb_Transp") + R("R_EX_tartr_L_e_") -
 R("R_EX_tartr_L_e__Rev")  = 0.

MtartrLc: -R("R_TARTD") + R("R_TARTRt7") - R("R_TARTRt7_Rev")  = 0.

MtartrLe: -R("R_EX_tartr_L_e_") + R("R_EX_tartr_L_e__Rev") -
 R("R_TARTRt7") + R("R_TARTRt7_Rev")  = 0.

Mtaurb: R("Mtaurb_Transp") + R("R_EX_taur_e_") -
 R("R_EX_taur_e__Rev")  = 0.

Mtaurc: -R("R_TAUDO") + R("R_TAURabc")  = 0.

Mtaure: -R("R_EX_taur_e_") + R("R_EX_taur_e__Rev") - R("R_TAURabc") 
 = 0.

Mtcyntc: R("R_CYANST")  = 0.

MtdeACPc: R("R_AACPS2") + R("R_C141SN") - 0.1 R("R_PASYN_EC")  = 0.

Mthdpc: R("R_DHDPRy") - R("R_THDPS")  = 0.

Mthfc: R("R_AICART") - R("R_AICART_Rev") + R("R_DHFR") -
 R("R_DHFR_Rev") + R("R_FTHFD") + R("R_GARFT") - R("R_GARFT_Rev") - R("R_GHMT2") -
 R("R_GLYCL") + R("R_METS") + R("R_MOHMT")  = 0.

Mthmb: R("Mthmb_Transp") + R("R_EX_thm_e_") - R("R_EX_thm_e__Rev") 
 = 0.

Mthmc: R("R_THMabc") - R("R_TMKr") + R("R_TMKr_Rev")  = 0.

Mthme: -R("R_EX_thm_e_") + R("R_EX_thm_e__Rev") - R("R_THMabc")  =
 0.

Mthmmpc: R("R_TMKr") - R("R_TMKr_Rev") - R("R_TMPKr") +
 R("R_TMPKr_Rev") + R("R_TMPPP")  = 0.

Mthmppc: -R("R_OXGDC2") + R("R_SHCHCS2") + R("R_TMPKr") -
 R("R_TMPKr_Rev")  = 0.

MthrLb: R("MthrLb_Transp") + R("R_EX_thr_L_e_") -
 R("R_EX_thr_L_e__Rev")  = 0.

MthrLc: -0.241 R("R_BiomassEcoli") + R("R_THRabc") - R("R_THRAr") +
 R("R_THRAr_Rev") - R("R_THRD") - R("R_THRD_L") + R("R_THRS") +
 R("R_THRt2r") - R("R_THRt2r_Rev") + R("R_THRt4")  = 0.

MthrLe: -R("R_EX_thr_L_e_") + R("R_EX_thr_L_e__Rev") - R("R_THRabc")
 - R("R_THRt2r") + R("R_THRt2r_Rev") - R("R_THRt4")  = 0.

Mthymc: R("R_TMDPP") - R("R_TMDPP_Rev")  = 0.

Mthymdb: R("Mthymdb_Transp") + R("R_EX_thymd_e_") -
 R("R_EX_thymd_e__Rev")  = 0.

Mthymdc: R("R_NTD5") + R("R_THMDt2") + R("R_THMDt2r") -
 R("R_THMDt2r_Rev") - R("R_TMDK1") - R("R_TMDPP") + R("R_TMDPP_Rev")  = 0.

Mthymde: -R("R_EX_thymd_e_") + R("R_EX_thymd_e__Rev") -
 R("R_THMDt2") - R("R_THMDt2r") + R("R_THMDt2r_Rev")  = 0.

Mtmab: R("Mtmab_Transp") + R("R_EX_tma_e_") - R("R_EX_tma_e__Rev") 
 = 0.

Mtmac: R("R_TMAOR1") + R("R_TMAOR2")  = 0.

Mtmae: -R("R_EX_tma_e_") + R("R_EX_tma_e__Rev") + R("R_TMAOR1e") +
 R("R_TMAOR2e")  = 0.

Mtmaob: R("Mtmaob_Transp") + R("R_EX_tmao_e_") -
 R("R_EX_tmao_e__Rev")  = 0.

Mtmaoc: -R("R_TMAOR1") - R("R_TMAOR2")  = 0.

Mtmaoe: -R("R_EX_tmao_e_") + R("R_EX_tmao_e__Rev") - R("R_TMAOR1e")
 - R("R_TMAOR2e")  = 0.

Mtrdoxc: R("R_PAPSR") + R("R_RNDR1") + R("R_RNDR2") + R("R_RNDR3") +
 R("R_RNDR4") + R("R_RNTR1") + R("R_RNTR2") + R("R_RNTR3") +
 R("R_RNTR4") - R("R_TRDR")  = 0.

Mtrdrdc: -R("R_PAPSR") - R("R_RNDR1") - R("R_RNDR2") - R("R_RNDR3")
 - R("R_RNDR4") - R("R_RNTR1") - R("R_RNTR2") - R("R_RNTR3") -
 R("R_RNTR4") + R("R_TRDR")  = 0.

Mtre6pc: -R("R_TRE6PH") - R("R_TRE6PP") + R("R_TRE6PS") +
 R("R_TREpts")  = 0.

Mtreb: R("Mtreb_Transp") + R("R_EX_tre_e_") - R("R_EX_tre_e__Rev") 
 = 0.

Mtrec: R("R_TRE6PP") - R("R_TREH")  = 0.

Mtree: -R("R_EX_tre_e_") + R("R_EX_tre_e__Rev") - R("R_TREHe") -
 R("R_TREpts")  = 0.

Mtrnagluc: R("R_GLUTRR") - R("R_GLUTRS")  = 0.

MtrpLb: R("MtrpLb_Transp") + R("R_EX_trp_L_e_") -
 R("R_EX_trp_L_e__Rev")  = 0.

MtrpLc: -0.054 R("R_BiomassEcoli") - R("R_TRPAS2") +
 R("R_TRPAS2_Rev") + R("R_TRPS1") + R("R_TRPS2") + R("R_TRPt2r") - R("R_TRPt2r_Rev") 
 = 0.

MtrpLe: -R("R_EX_trp_L_e_") + R("R_EX_trp_L_e__Rev") - R("R_TRPt2r")
 + R("R_TRPt2r_Rev")  = 0.

Mtsulb: R("Mtsulb_Transp") + R("R_EX_tsul_e_") -
 R("R_EX_tsul_e__Rev")  = 0.

Mtsulc: -R("R_CYANST") + R("R_TSULabc")  = 0.

Mtsule: -R("R_EX_tsul_e_") + R("R_EX_tsul_e__Rev") - R("R_TSULabc") 
 = 0.

Mttdcab: R("Mttdcab_Transp") + R("R_EX_ttdca_e_") -
 R("R_EX_ttdca_e__Rev")  = 0.

Mttdcac: -R("R_AACPS1") - R("R_FAO1") + 0.02 R("R_LPLIPA1") + 0.02
 R("R_LPLIPA2") + 0.02 R("R_LPLIPA3") + 0.02 R("R_PLIPA1") + 0.02
 R("R_PLIPA2") + 0.02 R("R_PLIPA3") + R("R_TTDCAt2")  = 0.

Mttdcae: -R("R_EX_ttdca_e_") + R("R_EX_ttdca_e__Rev") -
 R("R_TTDCAt2")  = 0.

Mttdceac: -R("R_AACPS2") + 0.05 R("R_LPLIPA1") + 0.05 R("R_LPLIPA2")
 + 0.05 R("R_LPLIPA3") + 0.05 R("R_PLIPA1") + 0.05 R("R_PLIPA2") +
 0.05 R("R_PLIPA3")  = 0.

MtyrLb: R("MtyrLb_Transp") + R("R_EX_tyr_L_e_") -
 R("R_EX_tyr_L_e__Rev")  = 0.

MtyrLc: -0.131 R("R_BiomassEcoli") - R("R_THZPSN") + R("R_TYRt2r") -
 R("R_TYRt2r_Rev") - R("R_TYRTA") + R("R_TYRTA_Rev")  = 0.

MtyrLe: -R("R_EX_tyr_L_e_") + R("R_EX_tyr_L_e__Rev") - R("R_TYRt2r")
 + R("R_TYRt2r_Rev")  = 0.

Mu23gac: -R("R_LPADSS") + R("R_U23GAAT") - R("R_USHD")  = 0.

Mu3agac: R("R_UAGAAT") - R("R_UAGAAT_Rev") - R("R_UHGADA")  = 0.

Mu3hgac: -R("R_U23GAAT") + R("R_UHGADA")  = 0.

Muaagmdac: -R("R_PPTGS") + R("R_UAGPT3")  = 0.

Muaccgc: R("R_UAGCVT") - R("R_UAPGR")  = 0.

Muacgamc: -R("R_ACGAMT") - R("R_UAG2Ei") - R("R_UAGAAT") +
 R("R_UAGAAT_Rev") - R("R_UAGCVT") + R("R_UAGDP") - R("R_UAGPT3")  = 0.

Muacmamc: -R("R_UACMAMO") + R("R_UAG2Ei")  = 0.

Muacmamuc: -R("R_ACMAMUT") + R("R_UACMAMO")  = 0.

Muagmdac: R("R_PAPPT3") - R("R_UAGPT3")  = 0.

Muamac: -R("R_UAMAGS") + R("R_UAMAS")  = 0.

Muamagc: -R("R_UAAGDS") + R("R_UAMAGS")  = 0.

Muamrc: -R("R_UAMAS") + R("R_UAPGR")  = 0.

Mudcpdpc: R("R_ECAP_EC") + R("R_PPTGS") - R("R_UDCPDP") +
 R("R_UDCPDPS")  = 0.

Mudcppc: -R("R_ACGAMT") - R("R_PAPPT3") + R("R_UDCPDP")  = 0.

Mudpc: R("R_ACMAMUT") + R("R_LPADSS") + 2 R("R_LPSSYN_EC") -
 R("R_NDPK2") + R("R_NDPK2_Rev") - R("R_RNDR4") + R("R_TRE6PS") +
 R("R_UAGPT3") + R("R_UMPK") - R("R_UMPK_Rev")  = 0.

Mudpgalc: R("R_UDPG4E") - R("R_UDPG4E_Rev") - R("R_UDPGALM") +
 R("R_UGLT") - R("R_UGLT_Rev")  = 0.

Mudpgalfurc: R("R_UDPGALM")  = 0.

Mudpgc: -0.003 R("R_BiomassEcoli") + R("R_GALU") + R("R_GALUi") -
 R("R_GALU_Rev") - 2 R("R_LPSSYN_EC") - R("R_TRE6PS") - R("R_UDPG4E") +
 R("R_UDPG4E_Rev") - R("R_UDPGD") - R("R_UGLT") + R("R_UGLT_Rev")  = 0.

Mudpglcurc: R("R_UDPGD")  = 0.

Mugmdac: -R("R_PAPPT3") + R("R_UGMDDS")  = 0.

Mugmdc: R("R_UAAGDS") - R("R_UGMDDS")  = 0.

Mumpc: R("R_ACGAMT") - R("R_NTD2") + R("R_NTPP8") + R("R_OMPDC") +
 R("R_PAPPT3") - R("R_UMPK") + R("R_UMPK_Rev") + R("R_UPPRT") +
 R("R_URIK2") + R("R_USHD")  = 0.

Munagac: R("R_ACGAMT") - R("R_ACMAMUT")  = 0.

Munagamuc: -R("R_AADDGT") + R("R_ACMAMUT")  = 0.

Munagamufc: R("R_AADDGT") - R("R_ECAP_EC")  = 0.

Muppg3c: -R("R_UPP3MT") + R("R_UPP3S") - R("R_UPPDC1")  = 0.

Murab: R("Murab_Transp") + R("R_EX_ura_e_") - R("R_EX_ura_e__Rev") 
 = 0.

Murac: R("R_CSND") + R("R_DURIPP") - R("R_DURIPP_Rev") +
 R("R_PYNP2r") - R("R_PYNP2r_Rev") - R("R_UPPRT") + R("R_URAt2") + R("R_URAt2r") -
 R("R_URAt2r_Rev")  = 0.

Murae: -R("R_EX_ura_e_") + R("R_EX_ura_e__Rev") - R("R_URAt2") -
 R("R_URAt2r") + R("R_URAt2r_Rev")  = 0.

Murdglycc: R("R_ALLTAH") - R("R_UGLYCH")  = 0.

Mureab: R("Mureab_Transp") + R("R_EX_urea_e_") -
 R("R_EX_urea_e__Rev")  = 0.

Mureac: R("R_AGMT") + R("R_ALLTAH") + R("R_UREAt") -
 R("R_UREAt_Rev")  = 0.

Mureae: -R("R_EX_urea_e_") + R("R_EX_urea_e__Rev") - R("R_UREAt") +
 R("R_UREAt_Rev")  = 0.

Murib: R("Murib_Transp") + R("R_EX_uri_e_") - R("R_EX_uri_e__Rev") 
 = 0.

Muric: R("R_CYTD") + R("R_NTD2") - R("R_PYNP2r") + R("R_PYNP2r_Rev")
 - R("R_URIK2") + R("R_URIt2") + R("R_URIt2r") - R("R_URIt2r_Rev")  =
 0.

Murie: -R("R_EX_uri_e_") + R("R_EX_uri_e__Rev") - R("R_URIt2") -
 R("R_URIt2r") + R("R_URIt2r_Rev")  = 0.

Mutpc: -0.136 R("R_BiomassEcoli") - R("R_CTPS2") - R("R_GALU") -
 R("R_GALUi") + R("R_GALU_Rev") + R("R_NDPK2") - R("R_NDPK2_Rev") -
 R("R_NTPP8") - R("R_RNTR4") - R("R_UAGDP")  = 0.

MvalLb: R("MvalLb_Transp") + R("R_EX_val_L_e_") -
 R("R_EX_val_L_e__Rev")  = 0.

MvalLc: -0.402 R("R_BiomassEcoli") + R("R_VALabc") + R("R_VALt2r") -
 R("R_VALt2r_Rev") - R("R_VALTA") + R("R_VALTA_Rev") + R("R_VPAMT")  =
 0.

MvalLe: -R("R_EX_val_L_e_") + R("R_EX_val_L_e__Rev") - R("R_VALabc")
 - R("R_VALt2r") + R("R_VALt2r_Rev")  = 0.

Mxanb: R("Mxanb_Transp") + R("R_EX_xan_e_") - R("R_EX_xan_e__Rev") 
 = 0.

Mxanc: R("R_GUAD") + R("R_PUNP7") - R("R_PUNP7_Rev") + R("R_XANt") +
 R("R_XANt2") - R("R_XANt_Rev") - R("R_XPPT")  = 0.

Mxane: -R("R_EX_xan_e_") + R("R_EX_xan_e__Rev") - R("R_XANt") -
 R("R_XANt2") + R("R_XANt_Rev")  = 0.

Mxmpc: -R("R_GMPS2") + R("R_IMPD") - R("R_NTD10") + R("R_XPPT")  =
 0.

Mxtsnb: R("Mxtsnb_Transp") + R("R_EX_xtsn_e_") -
 R("R_EX_xtsn_e__Rev")  = 0.

Mxtsnc: R("R_NTD10") - R("R_PUNP7") + R("R_PUNP7_Rev") +
 R("R_XTSNt2r") - R("R_XTSNt2r_Rev")  = 0.

Mxtsne: -R("R_EX_xtsn_e_") + R("R_EX_xtsn_e__Rev") - R("R_XTSNt2r")
 + R("R_XTSNt2r_Rev")  = 0.

Mxu5pDc: R("R_RBP4E") - R("R_RBP4E_Rev") + R("R_RPE") -
 R("R_RPE_Rev") - R("R_TKT1") + R("R_TKT1_Rev") - R("R_TKT2") + R("R_TKT2_Rev") +
 R("R_XYLK")  = 0.

Mxu5pLc: R("R_KG6PDC") - R("R_X5PL3E")  = 0.

MxylDb: R("MxylDb_Transp") + R("R_EX_xyl_D_e_") -
 R("R_EX_xyl_D_e__Rev")  = 0.

MxylDc: R("R_XYLabc") - R("R_XYLI1") + R("R_XYLI1_Rev") +
 R("R_XYLt2")  = 0.

MxylDe: -R("R_EX_xyl_D_e_") + R("R_EX_xyl_D_e__Rev") - R("R_XYLabc")
 - R("R_XYLt2")  = 0.

MxyluDc: R("R_XYLI1") - R("R_XYLI1_Rev") - R("R_XYLK")  = 0.

Bounds

0 <= R("R_EX_12ppd_S_e_") <= 999999
0 <= R("R_EX_12ppd_S_e__Rev") <= 999999
0 <= R("R_12PPDt") <= 999999
0 <= R("R_12PPDt_Rev") <= 999999
0 <= R("R_EX_15dap_e_") <= 999999
0 <= R("R_EX_15dap_e__Rev") <= 999999
0 <= R("R_ADSL2r") <= 999999
0 <= R("R_ADSL2r_Rev") <= 999999
0 <= R("R_DAPE") <= 999999
0 <= R("R_DAPE_Rev") <= 999999
0 <= R("R_EX_26dap_M_e_") <= 999999
0 <= R("R_EX_26dap_M_e__Rev") <= 999999
0 <= R("R_EX_2ddglcn_e_") <= 999999
0 <= R("R_EX_2ddglcn_e__Rev") <= 999999
0 <= R("R_DDPGALA") <= 999999
0 <= R("R_DDPGALA_Rev") <= 999999
0 <= R("R_PPM2") <= 999999
0 <= R("R_PPM2_Rev") <= 999999
0 <= R("R_PGM") <= 999999
0 <= R("R_PGM_Rev") <= 999999
0 <= R("R_ENO") <= 999999
0 <= R("R_ENO_Rev") <= 999999
0 <= R("R_IPPMIa") <= 999999
0 <= R("R_IPPMIa_Rev") <= 999999
0 <= R("R_DHQD") <= 999999
0 <= R("R_DHQD_Rev") <= 999999
0 <= R("R_EX_3hcinnm_e_") <= 999999
0 <= R("R_EX_3hcinnm_e__Rev") <= 999999
0 <= R("R_EX_3hpppn_e_") <= 999999
0 <= R("R_EX_3hpppn_e__Rev") <= 999999
0 <= R("R_EX_4abut_e_") <= 999999
0 <= R("R_EX_4abut_e__Rev") <= 999999
0 <= R("R_AIRC3") <= 999999
0 <= R("R_AIRC3_Rev") <= 999999
0 <= R("R_MTRI") <= 999999
0 <= R("R_MTRI_Rev") <= 999999
0 <= R("R_EX_acac_e_") <= 999999
0 <= R("R_EX_acac_e__Rev") <= 999999
0 <= R("R_EX_acald_e_") <= 999999
0 <= R("R_EX_acald_e__Rev") <= 999999
0 <= R("R_ACALDt") <= 999999
0 <= R("R_ACALDt_Rev") <= 999999
0 <= R("R_ACACT1r") <= 999999
0 <= R("R_ACACT1r_Rev") <= 999999
0 <= R("R_EX_ac_e_") <= 999999
0 <= R("R_EX_ac_e__Rev") <= 999999
0 <= R("R_EX_acgam_e_") <= 999999
0 <= R("R_EX_acgam_e__Rev") <= 999999
0 <= R("R_EX_acmana_e_") <= 999999
0 <= R("R_EX_acmana_e__Rev") <= 999999
0 <= R("R_EX_acnam_e_") <= 999999
0 <= R("R_EX_acnam_e__Rev") <= 999999
0 <= R("R_ACOATA") <= 999999
0 <= R("R_ACOATA_Rev") <= 999999
0 <= R("R_EX_ade_e_") <= 999999
0 <= R("R_EX_ade_e__Rev") <= 999999
0 <= R("R_EX_adn_e_") <= 999999
0 <= R("R_EX_adn_e__Rev") <= 999999
0 <= R("R_AICART") <= 999999
0 <= R("R_AICART_Rev") <= 999999
0 <= R("R_ACOTA") <= 999999
0 <= R("R_ACOTA_Rev") <= 999999
0 <= R("R_EX_akg_e_") <= 999999
0 <= R("R_EX_akg_e__Rev") <= 999999
0 <= R("R_EX_ala_D_e_") <= 999999
0 <= R("R_EX_ala_D_e__Rev") <= 999999
0 <= R("R_ALAR") <= 999999
0 <= R("R_ALAR_Rev") <= 999999
0 <= R("R_ALATA_L") <= 999999
0 <= R("R_ALATA_L_Rev") <= 999999
0 <= R("R_EX_ala_L_e_") <= 999999
0 <= R("R_EX_ala_L_e__Rev") <= 999999
0 <= R("R_EX_alltn_e_") <= 999999
0 <= R("R_EX_alltn_e__Rev") <= 999999
0 <= R("R_AMAOTr") <= 999999
0 <= R("R_AMAOTr_Rev") <= 999999
0 <= R("R_EX_amp_e_") <= 999999
0 <= R("R_EX_amp_e__Rev") <= 999999
0 <= R("R_ARAI") <= 999999
0 <= R("R_ARAI_Rev") <= 999999
0 <= R("R_EX_arab_L_e_") <= 999999
0 <= R("R_EX_arab_L_e__Rev") <= 999999
0 <= R("R_EX_arg_L_e_") <= 999999
0 <= R("R_EX_arg_L_e__Rev") <= 999999
0 <= R("R_ARGSL") <= 999999
0 <= R("R_ARGSL_Rev") <= 999999
0 <= R("R_EX_asn_L_e_") <= 999999
0 <= R("R_EX_asn_L_e__Rev") <= 999999
0 <= R("R_ASPTA") <= 999999
0 <= R("R_ASPTA_Rev") <= 999999
0 <= R("R_EX_asp_L_e_") <= 999999
0 <= R("R_EX_asp_L_e__Rev") <= 999999
0 <= R("R_DHBSr") <= 999999
0 <= R("R_DHBSr_Rev") <= 999999
0 <= R("R_PGK") <= 999999
0 <= R("R_PGK_Rev") <= 999999
0 <= R("R_ACKr") <= 999999
0 <= R("R_ACKr_Rev") <= 999999
0 <= R("R_ALAALAr") <= 999999
0 <= R("R_ALAALAr_Rev") <= 999999
0 <= R("R_ADK1") <= 999999
0 <= R("R_ADK1_Rev") <= 999999
0 <= R("R_ASPK") <= 999999
0 <= R("R_ASPK_Rev") <= 999999
0 <= R("R_PRASCS") <= 999999
0 <= R("R_PRASCS_Rev") <= 999999
0 <= R("R_EX_but_e_") <= 999999
0 <= R("R_EX_but_e__Rev") <= 999999
0 <= R("R_EX_cbl1_e_") <= 999999
0 <= R("R_EX_cbl1_e__Rev") <= 999999
0 <= R("R_NDPK3") <= 999999
0 <= R("R_NDPK3_Rev") <= 999999
0 <= R("R_EX_chol_e_") <= 999999
0 <= R("R_EX_chol_e__Rev") <= 999999
0 <= R("R_ACONT") <= 999999
0 <= R("R_ACONT_Rev") <= 999999
0 <= R("R_EX_cit_e_") <= 999999
0 <= R("R_EX_cit_e__Rev") <= 999999
0 <= R("R_CYTK1") <= 999999
0 <= R("R_CYTK1_Rev") <= 999999
0 <= R("R_EX_co2_e_") <= 999999
0 <= R("R_EX_co2_e__Rev") <= 999999
0 <= R("R_CO2t") <= 999999
0 <= R("R_CO2t_Rev") <= 999999
0 <= R("R_CRNBTCT") <= 999999
0 <= R("R_CRNBTCT_Rev") <= 999999
0 <= R("R_CRNCDH") <= 999999
0 <= R("R_CRNCDH_Rev") <= 999999
0 <= R("R_EX_crn_e_") <= 999999
0 <= R("R_EX_crn_e__Rev") <= 999999
0 <= R("R_EX_csn_e_") <= 999999
0 <= R("R_EX_csn_e__Rev") <= 999999
0 <= R("R_CRNCBCT") <= 999999
0 <= R("R_CRNCBCT_Rev") <= 999999
0 <= R("R_EX_cynt_e_") <= 999999
0 <= R("R_EX_cynt_e__Rev") <= 999999
0 <= R("R_EX_cys_L_e_") <= 999999
0 <= R("R_EX_cys_L_e__Rev") <= 999999
0 <= R("R_EX_cytd_e_") <= 999999
0 <= R("R_EX_cytd_e__Rev") <= 999999
0 <= R("R_EX_dad_2_e_") <= 999999
0 <= R("R_EX_dad_2_e__Rev") <= 999999
0 <= R("R_NDPK8") <= 999999
0 <= R("R_NDPK8_Rev") <= 999999
0 <= R("R_DADK") <= 999999
0 <= R("R_DADK_Rev") <= 999999
0 <= R("R_DBTSr") <= 999999
0 <= R("R_DBTSr_Rev") <= 999999
0 <= R("R_ADSL1r") <= 999999
0 <= R("R_ADSL1r_Rev") <= 999999
0 <= R("R_NDPK7") <= 999999
0 <= R("R_NDPK7_Rev") <= 999999
0 <= R("R_CYTK2") <= 999999
0 <= R("R_CYTK2_Rev") <= 999999
0 <= R("R_EX_dcyt_e_") <= 999999
0 <= R("R_EX_dcyt_e__Rev") <= 999999
0 <= R("R_NDPK5") <= 999999
0 <= R("R_NDPK5_Rev") <= 999999
0 <= R("R_DGK1") <= 999999
0 <= R("R_DGK1_Rev") <= 999999
0 <= R("R_EX_dgsn_e_") <= 999999
0 <= R("R_EX_dgsn_e__Rev") <= 999999
0 <= R("R_EX_dha_e_") <= 999999
0 <= R("R_EX_dha_e__Rev") <= 999999
0 <= R("R_DHAt") <= 999999
0 <= R("R_DHAt_Rev") <= 999999
0 <= R("R_TPI") <= 999999
0 <= R("R_TPI_Rev") <= 999999
0 <= R("R_EX_din_e_") <= 999999
0 <= R("R_EX_din_e__Rev") <= 999999
0 <= R("R_EX_dms_e_") <= 999999
0 <= R("R_EX_dms_e__Rev") <= 999999
0 <= R("R_EX_dmso_e_") <= 999999
0 <= R("R_EX_dmso_e__Rev") <= 999999
0 <= R("R_BTS2") <= 999999
0 <= R("R_BTS2_Rev") <= 999999
0 <= R("R_NDPK4") <= 999999
0 <= R("R_NDPK4_Rev") <= 999999
0 <= R("R_DTMPK") <= 999999
0 <= R("R_DTMPK_Rev") <= 999999
0 <= R("R_NDPK6") <= 999999
0 <= R("R_NDPK6_Rev") <= 999999
0 <= R("R_URIDK2r") <= 999999
0 <= R("R_URIDK2r_Rev") <= 999999
0 <= R("R_EX_duri_e_") <= 999999
0 <= R("R_EX_duri_e__Rev") <= 999999
0 <= R("R_EX_etoh_e_") <= 999999
0 <= R("R_EX_etoh_e__Rev") <= 999999
0 <= R("R_F6PA") <= 999999
0 <= R("R_F6PA_Rev") <= 999999
0 <= R("R_FCLPA") <= 999999
0 <= R("R_FCLPA_Rev") <= 999999
0 <= R("R_FBA") <= 999999
0 <= R("R_FBA_Rev") <= 999999
0 <= R("R_EX_fe2_e_") <= 999999
0 <= R("R_EX_fe2_e__Rev") <= 999999
0 <= R("R_EX_for_e_") <= 999999
0 <= R("R_EX_for_e__Rev") <= 999999
0 <= R("R_FORt") <= 999999
0 <= R("R_FORt_Rev") <= 999999
0 <= R("R_EX_fru_e_") <= 999999
0 <= R("R_EX_fru_e__Rev") <= 999999
0 <= R("R_EX_fuc1p_L_e_") <= 999999
0 <= R("R_EX_fuc1p_L_e__Rev") <= 999999
0 <= R("R_FCI") <= 999999
0 <= R("R_FCI_Rev") <= 999999
0 <= R("R_EX_fuc_L_e_") <= 999999
0 <= R("R_EX_fuc_L_e__Rev") <= 999999
0 <= R("R_EX_fum_e_") <= 999999
0 <= R("R_EX_fum_e__Rev") <= 999999
0 <= R("R_PGMT") <= 999999
0 <= R("R_PGMT_Rev") <= 999999
0 <= R("R_PGI") <= 999999
0 <= R("R_PGI_Rev") <= 999999
0 <= R("R_EX_g6p_e_") <= 999999
0 <= R("R_EX_g6p_e__Rev") <= 999999
0 <= R("R_GALKr") <= 999999
0 <= R("R_GALKr_Rev") <= 999999
0 <= R("R_EX_galct_D_e_") <= 999999
0 <= R("R_EX_galct_D_e__Rev") <= 999999
0 <= R("R_EX_galctn_D_e_") <= 999999
0 <= R("R_EX_galctn_D_e__Rev") <= 999999
0 <= R("R_EX_gal_e_") <= 999999
0 <= R("R_EX_gal_e__Rev") <= 999999
0 <= R("R_EX_galt_e_") <= 999999
0 <= R("R_EX_galt_e__Rev") <= 999999
0 <= R("R_GUI2") <= 999999
0 <= R("R_GUI2_Rev") <= 999999
0 <= R("R_EX_galur_e_") <= 999999
0 <= R("R_EX_galur_e__Rev") <= 999999
0 <= R("R_PGAMT") <= 999999
0 <= R("R_PGAMT_Rev") <= 999999
0 <= R("R_EX_gam_e_") <= 999999
0 <= R("R_EX_gam_e__Rev") <= 999999
0 <= R("R_GARFT") <= 999999
0 <= R("R_GARFT_Rev") <= 999999
0 <= R("R_EX_gbbtn_e_") <= 999999
0 <= R("R_EX_gbbtn_e__Rev") <= 999999
0 <= R("R_NDPK1") <= 999999
0 <= R("R_NDPK1_Rev") <= 999999
0 <= R("R_EX_glc_e_") <= 999999
0 <= R("R_EX_glc_e__Rev") <= 999999
0 <= R("R_EX_glcn_e_") <= 999999
0 <= R("R_EX_glcn_e__Rev") <= 999999
0 <= R("R_EX_glcr_e_") <= 999999
0 <= R("R_EX_glcr_e__Rev") <= 999999
0 <= R("R_GUI1") <= 999999
0 <= R("R_GUI1_Rev") <= 999999
0 <= R("R_EX_glcur_e_") <= 999999
0 <= R("R_EX_glcur_e__Rev") <= 999999
0 <= R("R_EX_gln_L_e_") <= 999999
0 <= R("R_EX_gln_L_e__Rev") <= 999999
0 <= R("R_GLUR") <= 999999
0 <= R("R_GLUR_Rev") <= 999999
0 <= R("R_EX_glu_L_e_") <= 999999
0 <= R("R_EX_glu_L_e__Rev") <= 999999
0 <= R("R_GLUABUTt7") <= 999999
0 <= R("R_GLUABUTt7_Rev") <= 999999
0 <= R("R_EX_glyald_e_") <= 999999
0 <= R("R_EX_glyald_e__Rev") <= 999999
0 <= R("R_GLYALDt") <= 999999
0 <= R("R_GLYALDt_Rev") <= 999999
0 <= R("R_EX_glyb_e_") <= 999999
0 <= R("R_EX_glyb_e__Rev") <= 999999
0 <= R("R_PGSA_EC") <= 999999
0 <= R("R_PGSA_EC_Rev") <= 999999
0 <= R("R_EX_glyc3p_e_") <= 999999
0 <= R("R_EX_glyc3p_e__Rev") <= 999999
0 <= R("R_GLYCt") <= 999999
0 <= R("R_GLYCt_Rev") <= 999999
0 <= R("R_EX_glyc_e_") <= 999999
0 <= R("R_EX_glyc_e__Rev") <= 999999
0 <= R("R_EX_glyclt_e_") <= 999999
0 <= R("R_EX_glyclt_e__Rev") <= 999999
0 <= R("R_EX_gly_e_") <= 999999
0 <= R("R_EX_gly_e__Rev") <= 999999
0 <= R("R_GK1") <= 999999
0 <= R("R_GK1_Rev") <= 999999
0 <= R("R_EX_gsn_e_") <= 999999
0 <= R("R_EX_gsn_e__Rev") <= 999999
0 <= R("R_ADK3") <= 999999
0 <= R("R_ADK3_Rev") <= 999999
0 <= R("R_EX_gua_e_") <= 999999
0 <= R("R_EX_gua_e__Rev") <= 999999
0 <= R("R_GUAt") <= 999999
0 <= R("R_GUAt_Rev") <= 999999
0 <= R("R_IPPMIb") <= 999999
0 <= R("R_IPPMIb_Rev") <= 999999
0 <= R("R_AHC") <= 999999
0 <= R("R_AHC_Rev") <= 999999
0 <= R("R_CBIAT_DELETE") <= 999999
0 <= R("R_CBIAT_DELETE_Rev") <= 999999
0 <= R("R_CBLAT_DELETE") <= 999999
0 <= R("R_CBLAT_DELETE_Rev") <= 999999
0 <= R("R_HCO3E") <= 999999
0 <= R("R_HCO3E_Rev") <= 999999
0 <= R("R_DHORTS") <= 999999
0 <= R("R_DHORTS_Rev") <= 999999
0 <= R("R_FUM") <= 999999
0 <= R("R_FUM_Rev") <= 999999
0 <= R("R_EX_h2o_e_") <= 999999
0 <= R("R_EX_h2o_e__Rev") <= 999999
0 <= R("R_H2Ot") <= 999999
0 <= R("R_H2Ot_Rev") <= 999999
0 <= R("R_ADMDCr") <= 999999
0 <= R("R_ADMDCr_Rev") <= 999999
0 <= R("R_ACCOACr") <= 999999
0 <= R("R_ACCOACr_Rev") <= 999999
0 <= R("R_EX_hdca_e_") <= 999999
0 <= R("R_EX_hdca_e__Rev") <= 999999
0 <= R("R_EX_h_e_") <= 999999
0 <= R("R_EX_h_e__Rev") <= 999999
0 <= R("R_DDGLCNt2r") <= 999999
0 <= R("R_DDGLCNt2r_Rev") <= 999999
0 <= R("R_HCINNMt2r") <= 999999
0 <= R("R_HCINNMt2r_Rev") <= 999999
0 <= R("R_HPPPNt2r") <= 999999
0 <= R("R_HPPPNt2r_Rev") <= 999999
0 <= R("R_ACACt2") <= 999999
0 <= R("R_ACACt2_Rev") <= 999999
0 <= R("R_ACt2r") <= 999999
0 <= R("R_ACt2r_Rev") <= 999999
0 <= R("R_ADEt2r") <= 999999
0 <= R("R_ADEt2r_Rev") <= 999999
0 <= R("R_ADNt2r") <= 999999
0 <= R("R_ADNt2r_Rev") <= 999999
0 <= R("R_AKGt2r") <= 999999
0 <= R("R_AKGt2r_Rev") <= 999999
0 <= R("R_DALAt2r") <= 999999
0 <= R("R_DALAt2r_Rev") <= 999999
0 <= R("R_ALAt2r") <= 999999
0 <= R("R_ALAt2r_Rev") <= 999999
0 <= R("R_ALLTNt2r") <= 999999
0 <= R("R_ALLTNt2r_Rev") <= 999999
0 <= R("R_ARBt2r") <= 999999
0 <= R("R_ARBt2r_Rev") <= 999999
0 <= R("R_ASNt2r") <= 999999
0 <= R("R_ASNt2r_Rev") <= 999999
0 <= R("R_BUTt2r") <= 999999
0 <= R("R_BUTt2r_Rev") <= 999999
0 <= R("R_CHLt2r") <= 999999
0 <= R("R_CHLt2r_Rev") <= 999999
0 <= R("R_CYTDt2r") <= 999999
0 <= R("R_CYTDt2r_Rev") <= 999999
0 <= R("R_ETOHt2r") <= 999999
0 <= R("R_ETOHt2r_Rev") <= 999999
0 <= R("R_FUCt") <= 999999
0 <= R("R_FUCt_Rev") <= 999999
0 <= R("R_GALCTt2r") <= 999999
0 <= R("R_GALCTt2r_Rev") <= 999999
0 <= R("R_GALCTNt2r") <= 999999
0 <= R("R_GALCTNt2r_Rev") <= 999999
0 <= R("R_GALURt2r") <= 999999
0 <= R("R_GALURt2r_Rev") <= 999999
0 <= R("R_GLCNt2r") <= 999999
0 <= R("R_GLCNt2r_Rev") <= 999999
0 <= R("R_GLCRt2r") <= 999999
0 <= R("R_GLCRt2r_Rev") <= 999999
0 <= R("R_GLCURt2r") <= 999999
0 <= R("R_GLCURt2r_Rev") <= 999999
0 <= R("R_GLUt2r") <= 999999
0 <= R("R_GLUt2r_Rev") <= 999999
0 <= R("R_GLYBt2r") <= 999999
0 <= R("R_GLYBt2r_Rev") <= 999999
0 <= R("R_GLYCLTt2r") <= 999999
0 <= R("R_GLYCLTt2r_Rev") <= 999999
0 <= R("R_GLYt2r") <= 999999
0 <= R("R_GLYt2r_Rev") <= 999999
0 <= R("R_EX_his_L_e_") <= 999999
0 <= R("R_EX_his_L_e__Rev") <= 999999
0 <= R("R_HISt2r") <= 999999
0 <= R("R_HISt2r_Rev") <= 999999
0 <= R("R_HPYRI") <= 999999
0 <= R("R_HPYRI_Rev") <= 999999
0 <= R("R_EX_hxan_e_") <= 999999
0 <= R("R_EX_hxan_e__Rev") <= 999999
0 <= R("R_HYXNt") <= 999999
0 <= R("R_HYXNt_Rev") <= 999999
0 <= R("R_EX_idon_L_e_") <= 999999
0 <= R("R_EX_idon_L_e__Rev") <= 999999
0 <= R("R_IDONt2r") <= 999999
0 <= R("R_IDONt2r_Rev") <= 999999
0 <= R("R_ILETA") <= 999999
0 <= R("R_ILETA_Rev") <= 999999
0 <= R("R_EX_ile_L_e_") <= 999999
0 <= R("R_EX_ile_L_e__Rev") <= 999999
0 <= R("R_ILEt2r") <= 999999
0 <= R("R_ILEt2r_Rev") <= 999999
0 <= R("R_IMPC") <= 999999
0 <= R("R_IMPC_Rev") <= 999999
0 <= R("R_EX_indole_e_") <= 999999
0 <= R("R_EX_indole_e__Rev") <= 999999
0 <= R("R_INDOLEt2r") <= 999999
0 <= R("R_INDOLEt2r_Rev") <= 999999
0 <= R("R_EX_ins_e_") <= 999999
0 <= R("R_EX_ins_e__Rev") <= 999999
0 <= R("R_INSt2r") <= 999999
0 <= R("R_INSt2r_Rev") <= 999999
0 <= R("R_ADK4") <= 999999
0 <= R("R_ADK4_Rev") <= 999999
0 <= R("R_EX_k_e_") <= 999999
0 <= R("R_EX_k_e__Rev") <= 999999
0 <= R("R_Kt2r") <= 999999
0 <= R("R_Kt2r_Rev") <= 999999
0 <= R("R_EX_lac_D_e_") <= 999999
0 <= R("R_EX_lac_D_e__Rev") <= 999999
0 <= R("R_D_LACt2") <= 999999
0 <= R("R_D_LACt2_Rev") <= 999999
0 <= R("R_EX_lac_L_e_") <= 999999
0 <= R("R_EX_lac_L_e__Rev") <= 999999
0 <= R("R_L_LACt2r") <= 999999
0 <= R("R_L_LACt2r_Rev") <= 999999
0 <= R("R_EX_lcts_e_") <= 999999
0 <= R("R_EX_lcts_e__Rev") <= 999999
0 <= R("R_LCTSt") <= 999999
0 <= R("R_LCTSt_Rev") <= 999999
0 <= R("R_EX_leu_L_e_") <= 999999
0 <= R("R_EX_leu_L_e__Rev") <= 999999
0 <= R("R_LEUt2r") <= 999999
0 <= R("R_LEUt2r_Rev") <= 999999
0 <= R("R_EX_lys_L_e_") <= 999999
0 <= R("R_EX_lys_L_e__Rev") <= 999999
0 <= R("R_LYSt2r") <= 999999
0 <= R("R_LYSt2r_Rev") <= 999999
0 <= R("R_MCOATA") <= 999999
0 <= R("R_MCOATA_Rev") <= 999999
0 <= R("R_EX_mal_L_e_") <= 999999
0 <= R("R_EX_mal_L_e__Rev") <= 999999
0 <= R("R_EX_malt_e_") <= 999999
0 <= R("R_EX_malt_e__Rev") <= 999999
0 <= R("R_EX_malthx_e_") <= 999999
0 <= R("R_EX_malthx_e__Rev") <= 999999
0 <= R("R_EX_maltpt_e_") <= 999999
0 <= R("R_EX_maltpt_e__Rev") <= 999999
0 <= R("R_EX_malttr_e_") <= 999999
0 <= R("R_EX_malttr_e__Rev") <= 999999
0 <= R("R_EX_maltttr_e_") <= 999999
0 <= R("R_EX_maltttr_e__Rev") <= 999999
0 <= R("R_PMANM") <= 999999
0 <= R("R_PMANM_Rev") <= 999999
0 <= R("R_MAN6PI") <= 999999
0 <= R("R_MAN6PI_Rev") <= 999999
0 <= R("R_EX_man6p_e_") <= 999999
0 <= R("R_EX_man6p_e__Rev") <= 999999
0 <= R("R_EX_man_e_") <= 999999
0 <= R("R_EX_man_e__Rev") <= 999999
0 <= R("R_EX_melib_e_") <= 999999
0 <= R("R_EX_melib_e__Rev") <= 999999
0 <= R("R_EX_met_D_e_") <= 999999
0 <= R("R_EX_met_D_e__Rev") <= 999999
0 <= R("R_MTHFC") <= 999999
0 <= R("R_MTHFC_Rev") <= 999999
0 <= R("R_EX_met_L_e_") <= 999999
0 <= R("R_EX_met_L_e__Rev") <= 999999
0 <= R("R_MCITL2") <= 999999
0 <= R("R_MCITL2_Rev") <= 999999
0 <= R("R_MME") <= 999999
0 <= R("R_MME_Rev") <= 999999
0 <= R("R_EX_mnl_e_") <= 999999
0 <= R("R_EX_mnl_e__Rev") <= 999999
0 <= R("R_NAt3_1") <= 999999
0 <= R("R_NAt3_1_Rev") <= 999999
0 <= R("R_EX_na1_e_") <= 999999
0 <= R("R_EX_na1_e__Rev") <= 999999
0 <= R("R_EX_nac_e_") <= 999999
0 <= R("R_EX_nac_e__Rev") <= 999999
0 <= R("R_DHBD") <= 999999
0 <= R("R_DHBD_Rev") <= 999999
0 <= R("R_PERD") <= 999999
0 <= R("R_PERD_Rev") <= 999999
0 <= R("R_TAGURr") <= 999999
0 <= R("R_TAGURr_Rev") <= 999999
0 <= R("R_GLTPD") <= 999999
0 <= R("R_GLTPD_Rev") <= 999999
0 <= R("R_E4PD") <= 999999
0 <= R("R_E4PD_Rev") <= 999999
0 <= R("R_LDH_D") <= 999999
0 <= R("R_LDH_D_Rev") <= 999999
0 <= R("R_LCAD") <= 999999
0 <= R("R_LCAD_Rev") <= 999999
0 <= R("R_MDH") <= 999999
0 <= R("R_MDH_Rev") <= 999999
0 <= R("R_MANAO") <= 999999
0 <= R("R_MANAO_Rev") <= 999999
0 <= R("R_M1PD") <= 999999
0 <= R("R_M1PD_Rev") <= 999999
0 <= R("R_EX_nad_e_") <= 999999
0 <= R("R_EX_nad_e__Rev") <= 999999
0 <= R("R_IDOND") <= 999999
0 <= R("R_IDOND_Rev") <= 999999
0 <= R("R_ALCD19") <= 999999
0 <= R("R_ALCD19_Rev") <= 999999
0 <= R("R_LCARS") <= 999999
0 <= R("R_LCARS_Rev") <= 999999
0 <= R("R_ADHEr") <= 999999
0 <= R("R_ADHEr_Rev") <= 999999
0 <= R("R_G6PDH2r") <= 999999
0 <= R("R_G6PDH2r_Rev") <= 999999
0 <= R("R_G3PD2") <= 999999
0 <= R("R_G3PD2_Rev") <= 999999
0 <= R("R_GLUDy") <= 999999
0 <= R("R_GLUDy_Rev") <= 999999
0 <= R("R_HSDy") <= 999999
0 <= R("R_HSDy_Rev") <= 999999
0 <= R("R_ICDHyr") <= 999999
0 <= R("R_ICDHyr_Rev") <= 999999
0 <= R("R_MTHFD") <= 999999
0 <= R("R_MTHFD_Rev") <= 999999
0 <= R("R_SULR") <= 999999
0 <= R("R_SULR_Rev") <= 999999
0 <= R("R_SHK3Dr") <= 999999
0 <= R("R_SHK3Dr_Rev") <= 999999
0 <= R("R_5DGLCNR") <= 999999
0 <= R("R_5DGLCNR_Rev") <= 999999
0 <= R("R_DHFR") <= 999999
0 <= R("R_DHFR_Rev") <= 999999
0 <= R("R_GTHOr") <= 999999
0 <= R("R_GTHOr_Rev") <= 999999
0 <= R("R_EX_nh4_e_") <= 999999
0 <= R("R_EX_nh4_e__Rev") <= 999999
0 <= R("R_NH4t") <= 999999
0 <= R("R_NH4t_Rev") <= 999999
0 <= R("R_EX_nmn_e_") <= 999999
0 <= R("R_EX_nmn_e__Rev") <= 999999
0 <= R("R_EX_no2_e_") <= 999999
0 <= R("R_EX_no2_e__Rev") <= 999999
0 <= R("R_NO2t2r") <= 999999
0 <= R("R_NO2t2r_Rev") <= 999999
0 <= R("R_EX_no3_e_") <= 999999
0 <= R("R_EX_no3_e__Rev") <= 999999
0 <= R("R_EX_o2_e_") <= 999999
0 <= R("R_EX_o2_e__Rev") <= 999999
0 <= R("R_O2t") <= 999999
0 <= R("R_O2t_Rev") <= 999999
0 <= R("R_EX_ocdca_e_") <= 999999
0 <= R("R_EX_ocdca_e__Rev") <= 999999
0 <= R("R_OHPBAT") <= 999999
0 <= R("R_OHPBAT_Rev") <= 999999
0 <= R("R_ARGORNt7") <= 999999
0 <= R("R_ARGORNt7_Rev") <= 999999
0 <= R("R_OCBT") <= 999999
0 <= R("R_OCBT_Rev") <= 999999
0 <= R("R_EX_orn_e_") <= 999999
0 <= R("R_EX_orn_e__Rev") <= 999999
0 <= R("R_DASYN_EC") <= 999999
0 <= R("R_DASYN_EC_Rev") <= 999999
0 <= R("R_PDX5PO") <= 999999
0 <= R("R_PDX5PO_Rev") <= 999999
0 <= R("R_PEPT_EC") <= 999999
0 <= R("R_PEPT_EC_Rev") <= 999999
0 <= R("R_CLPNS_EC") <= 999999
0 <= R("R_CLPNS_EC_Rev") <= 999999
0 <= R("R_PHETA1") <= 999999
0 <= R("R_PHETA1_Rev") <= 999999
0 <= R("R_EX_phe_L_e_") <= 999999
0 <= R("R_EX_phe_L_e__Rev") <= 999999
0 <= R("R_PHEt2r") <= 999999
0 <= R("R_PHEt2r_Rev") <= 999999
0 <= R("R_PTAr") <= 999999
0 <= R("R_PTAr_Rev") <= 999999
0 <= R("R_PUNP1") <= 999999
0 <= R("R_PUNP1_Rev") <= 999999
0 <= R("R_PUNP2") <= 999999
0 <= R("R_PUNP2_Rev") <= 999999
0 <= R("R_PUNP4") <= 999999
0 <= R("R_PUNP4_Rev") <= 999999
0 <= R("R_PUNP6") <= 999999
0 <= R("R_PUNP6_Rev") <= 999999
0 <= R("R_DURIPP") <= 999999
0 <= R("R_DURIPP_Rev") <= 999999
0 <= R("R_PUNP3") <= 999999
0 <= R("R_PUNP3_Rev") <= 999999
0 <= R("R_ATPS4r") <= 999999
0 <= R("R_ATPS4r_Rev") <= 999999
0 <= R("R_PUNP5") <= 999999
0 <= R("R_PUNP5_Rev") <= 999999
0 <= R("R_MLTP3") <= 999999
0 <= R("R_MLTP3_Rev") <= 999999
0 <= R("R_MLTP2") <= 999999
0 <= R("R_MLTP2_Rev") <= 999999
0 <= R("R_MLTP1") <= 999999
0 <= R("R_MLTP1_Rev") <= 999999
0 <= R("R_GAPD") <= 999999
0 <= R("R_GAPD_Rev") <= 999999
0 <= R("R_AGPR") <= 999999
0 <= R("R_AGPR_Rev") <= 999999
0 <= R("R_ASAD") <= 999999
0 <= R("R_ASAD_Rev") <= 999999
0 <= R("R_EX_pi_e_") <= 999999
0 <= R("R_EX_pi_e__Rev") <= 999999
0 <= R("R_PIt2r") <= 999999
0 <= R("R_PIt2r_Rev") <= 999999
0 <= R("R_AOXSr") <= 999999
0 <= R("R_AOXSr_Rev") <= 999999
0 <= R("R_EX_pnto_R_e_") <= 999999
0 <= R("R_EX_pnto_R_e__Rev") <= 999999
0 <= R("R_PPAKr") <= 999999
0 <= R("R_PPAKr_Rev") <= 999999
0 <= R("R_ORPT") <= 999999
0 <= R("R_ORPT_Rev") <= 999999
0 <= R("R_EX_pppn_e_") <= 999999
0 <= R("R_EX_pppn_e__Rev") <= 999999
0 <= R("R_PPPNt2r") <= 999999
0 <= R("R_PPPNt2r_Rev") <= 999999
0 <= R("R_PRAGSr") <= 999999
0 <= R("R_PRAGSr_Rev") <= 999999
0 <= R("R_EX_pro_L_e_") <= 999999
0 <= R("R_EX_pro_L_e__Rev") <= 999999
0 <= R("R_PROt2r") <= 999999
0 <= R("R_PROt2r_Rev") <= 999999
0 <= R("R_EX_ptrc_e_") <= 999999
0 <= R("R_EX_ptrc_e__Rev") <= 999999
0 <= R("R_PTRCt2r") <= 999999
0 <= R("R_PTRCt2r_Rev") <= 999999
0 <= R("R_PTRCORNt7") <= 999999
0 <= R("R_PTRCORNt7_Rev") <= 999999
0 <= R("R_EX_pyr_e_") <= 999999
0 <= R("R_EX_pyr_e__Rev") <= 999999
0 <= R("R_PYRt2r") <= 999999
0 <= R("R_PYRt2r_Rev") <= 999999
0 <= R("R_SUCD4") <= 999999
0 <= R("R_SUCD4_Rev") <= 999999
0 <= R("R_PPM") <= 999999
0 <= R("R_PPM_Rev") <= 999999
0 <= R("R_RPI") <= 999999
0 <= R("R_RPI_Rev") <= 999999
0 <= R("R_PRPPS") <= 999999
0 <= R("R_PRPPS_Rev") <= 999999
0 <= R("R_EX_rib_D_e_") <= 999999
0 <= R("R_EX_rib_D_e__Rev") <= 999999
0 <= R("R_RMPA") <= 999999
0 <= R("R_RMPA_Rev") <= 999999
0 <= R("R_RMI") <= 999999
0 <= R("R_RMI_Rev") <= 999999
0 <= R("R_EX_rmn_e_") <= 999999
0 <= R("R_EX_rmn_e__Rev") <= 999999
0 <= R("R_A5PISO") <= 999999
0 <= R("R_A5PISO_Rev") <= 999999
0 <= R("R_RPE") <= 999999
0 <= R("R_RPE_Rev") <= 999999
0 <= R("R_RBP4E") <= 999999
0 <= R("R_RBP4E_Rev") <= 999999
0 <= R("R_TALA") <= 999999
0 <= R("R_TALA_Rev") <= 999999
0 <= R("R_SBTPD") <= 999999
0 <= R("R_SBTPD_Rev") <= 999999
0 <= R("R_EX_sbt_D_e_") <= 999999
0 <= R("R_EX_sbt_D_e__Rev") <= 999999
0 <= R("R_EX_ser_D_e_") <= 999999
0 <= R("R_EX_ser_D_e__Rev") <= 999999
0 <= R("R_DSERt2r") <= 999999
0 <= R("R_DSERt2r_Rev") <= 999999
0 <= R("R_SERAT") <= 999999
0 <= R("R_SERAT_Rev") <= 999999
0 <= R("R_PSSA_EC") <= 999999
0 <= R("R_PSSA_EC_Rev") <= 999999
0 <= R("R_SERASr") <= 999999
0 <= R("R_SERASr_Rev") <= 999999
0 <= R("R_EX_ser_L_e_") <= 999999
0 <= R("R_EX_ser_L_e__Rev") <= 999999
0 <= R("R_SERt2r") <= 999999
0 <= R("R_SERt2r_Rev") <= 999999
0 <= R("R_PSCVT") <= 999999
0 <= R("R_PSCVT_Rev") <= 999999
0 <= R("R_SDPTA") <= 999999
0 <= R("R_SDPTA_Rev") <= 999999
0 <= R("R_EX_so4_e_") <= 999999
0 <= R("R_EX_so4_e__Rev") <= 999999
0 <= R("R_EX_spmd_e_") <= 999999
0 <= R("R_EX_spmd_e__Rev") <= 999999
0 <= R("R_SUCOAS") <= 999999
0 <= R("R_SUCOAS_Rev") <= 999999
0 <= R("R_SUCFUMt") <= 999999
0 <= R("R_SUCFUMt_Rev") <= 999999
0 <= R("R_EX_succ_e_") <= 999999
0 <= R("R_EX_succ_e__Rev") <= 999999
0 <= R("R_EX_sucr_e_") <= 999999
0 <= R("R_EX_sucr_e__Rev") <= 999999
0 <= R("R_TGBPA") <= 999999
0 <= R("R_TGBPA_Rev") <= 999999
0 <= R("R_EX_tartr_L_e_") <= 999999
0 <= R("R_EX_tartr_L_e__Rev") <= 999999
0 <= R("R_TARTRt7") <= 999999
0 <= R("R_TARTRt7_Rev") <= 999999
0 <= R("R_EX_taur_e_") <= 999999
0 <= R("R_EX_taur_e__Rev") <= 999999
0 <= R("R_TMKr") <= 999999
0 <= R("R_TMKr_Rev") <= 999999
0 <= R("R_EX_thm_e_") <= 999999
0 <= R("R_EX_thm_e__Rev") <= 999999
0 <= R("R_TMPKr") <= 999999
0 <= R("R_TMPKr_Rev") <= 999999
0 <= R("R_THRAr") <= 999999
0 <= R("R_THRAr_Rev") <= 999999
0 <= R("R_EX_thr_L_e_") <= 999999
0 <= R("R_EX_thr_L_e__Rev") <= 999999
0 <= R("R_THRt2r") <= 999999
0 <= R("R_THRt2r_Rev") <= 999999
0 <= R("R_TMDPP") <= 999999
0 <= R("R_TMDPP_Rev") <= 999999
0 <= R("R_EX_thymd_e_") <= 999999
0 <= R("R_EX_thymd_e__Rev") <= 999999
0 <= R("R_THMDt2r") <= 999999
0 <= R("R_THMDt2r_Rev") <= 999999
0 <= R("R_EX_tma_e_") <= 999999
0 <= R("R_EX_tma_e__Rev") <= 999999
0 <= R("R_EX_tmao_e_") <= 999999
0 <= R("R_EX_tmao_e__Rev") <= 999999
0 <= R("R_EX_tre_e_") <= 999999
0 <= R("R_EX_tre_e__Rev") <= 999999
0 <= R("R_TRPAS2") <= 999999
0 <= R("R_TRPAS2_Rev") <= 999999
0 <= R("R_EX_trp_L_e_") <= 999999
0 <= R("R_EX_trp_L_e__Rev") <= 999999
0 <= R("R_TRPt2r") <= 999999
0 <= R("R_TRPt2r_Rev") <= 999999
0 <= R("R_EX_tsul_e_") <= 999999
0 <= R("R_EX_tsul_e__Rev") <= 999999
0 <= R("R_EX_ttdca_e_") <= 999999
0 <= R("R_EX_ttdca_e__Rev") <= 999999
0 <= R("R_TYRTA") <= 999999
0 <= R("R_TYRTA_Rev") <= 999999
0 <= R("R_EX_tyr_L_e_") <= 999999
0 <= R("R_EX_tyr_L_e__Rev") <= 999999
0 <= R("R_TYRt2r") <= 999999
0 <= R("R_TYRt2r_Rev") <= 999999
0 <= R("R_UAGAAT") <= 999999
0 <= R("R_UAGAAT_Rev") <= 999999
0 <= R("R_NDPK2") <= 999999
0 <= R("R_NDPK2_Rev") <= 999999
0 <= R("R_UDPG4E") <= 999999
0 <= R("R_UDPG4E_Rev") <= 999999
0 <= R("R_UGLT") <= 999999
0 <= R("R_UGLT_Rev") <= 999999
0 <= R("R_UMPK") <= 999999
0 <= R("R_UMPK_Rev") <= 999999
0 <= R("R_EX_ura_e_") <= 999999
0 <= R("R_EX_ura_e__Rev") <= 999999
0 <= R("R_URAt2r") <= 999999
0 <= R("R_URAt2r_Rev") <= 999999
0 <= R("R_EX_urea_e_") <= 999999
0 <= R("R_EX_urea_e__Rev") <= 999999
0 <= R("R_UREAt") <= 999999
0 <= R("R_UREAt_Rev") <= 999999
0 <= R("R_PYNP2r") <= 999999
0 <= R("R_PYNP2r_Rev") <= 999999
0 <= R("R_EX_uri_e_") <= 999999
0 <= R("R_EX_uri_e__Rev") <= 999999
0 <= R("R_URIt2r") <= 999999
0 <= R("R_URIt2r_Rev") <= 999999
0 <= R("R_GALU") <= 999999
0 <= R("R_GALU_Rev") <= 999999
0 <= R("R_VALTA") <= 999999
0 <= R("R_VALTA_Rev") <= 999999
0 <= R("R_EX_val_L_e_") <= 999999
0 <= R("R_EX_val_L_e__Rev") <= 999999
0 <= R("R_VALt2r") <= 999999
0 <= R("R_VALt2r_Rev") <= 999999
0 <= R("R_EX_xan_e_") <= 999999
0 <= R("R_EX_xan_e__Rev") <= 999999
0 <= R("R_XANt") <= 999999
0 <= R("R_XANt_Rev") <= 999999
0 <= R("R_PUNP7") <= 999999
0 <= R("R_PUNP7_Rev") <= 999999
0 <= R("R_EX_xtsn_e_") <= 999999
0 <= R("R_EX_xtsn_e__Rev") <= 999999
0 <= R("R_XTSNt2r") <= 999999
0 <= R("R_XTSNt2r_Rev") <= 999999
0 <= R("R_TKT2") <= 999999
0 <= R("R_TKT2_Rev") <= 999999
0 <= R("R_TKT1") <= 999999
0 <= R("R_TKT1_Rev") <= 999999
0 <= R("R_XYLI1") <= 999999
0 <= R("R_XYLI1_Rev") <= 999999
0 <= R("R_EX_xyl_D_e_") <= 999999
0 <= R("R_EX_xyl_D_e__Rev") <= 999999
0 <= R("R_DHAD1") <= 999999
0 <= R("R_DHAD2") <= 999999
0 <= R("R_DHQS") <= 999999
0 <= R("R_EDA") <= 999999
0 <= R("R_DRPA") <= 999999
0 <= R("R_MCITD") <= 999999
0 <= R("R_MECDPS") <= 999999
0 <= R("R_SUCBZS") <= 999999
0 <= R("R_TRPS3") <= 999999
0 <= R("R_CHORS") <= 999999
0 <= R("R_ADCL") <= 999999
0 <= R("R_HOPNTAL") <= 999999
0 <= R("R_PPBNGS") <= 999999
0 <= R("R_GLCRAL") <= 999999
0 <= R("R_MDRPD") <= 999999
0 <= R("R_DHPS2") <= 999999
0 <= R("R_EDD") <= 999999
0 <= R("R_ACACCT") <= 999999
0 <= R("R_AMANAPE") <= 999999
0 <= R("R_ACNML") <= 999999
0 <= R("R_GLCS1") <= 999999
0 <= R("R_AGMHE") <= 999999
0 <= R("R_ABTA") <= 999999
0 <= R("R_ALARi") <= 999999
0 <= R("R_VPAMT") <= 999999
0 <= R("R_ALTRH") <= 999999
0 <= R("R_AMMQT8_2") <= 999999
0 <= R("R_OHPHM") <= 999999
0 <= R("R_OMBZLM") <= 999999
0 <= R("R_DMQMT") <= 999999
0 <= R("R_ACONMT") <= 999999
0 <= R("R_ASPT") <= 999999
0 <= R("R_ADNCYC") <= 999999
0 <= R("R_DAGK_EC") <= 999999
0 <= R("R_DDGLK") <= 999999
0 <= R("R_DDGALK") <= 999999
0 <= R("R_HMPK1") <= 999999
0 <= R("R_PMPK") <= 999999
0 <= R("R_CDPMEK") <= 999999
0 <= R("R_HETZK") <= 999999
0 <= R("R_MTRK") <= 999999
0 <= R("R_HPPK2") <= 999999
0 <= R("R_ACGK") <= 999999
0 <= R("R_AMANK") <= 999999
0 <= R("R_ADNK1") <= 999999
0 <= R("R_ADOCBIK") <= 999999
0 <= R("R_ADSK") <= 999999
0 <= R("R_BUTCT") <= 999999
0 <= R("R_ASPCT") <= 999999
0 <= R("R_ICHORSi") <= 999999
0 <= R("R_CHORM") <= 999999
0 <= R("R_CHRPL") <= 999999
0 <= R("R_CITL") <= 999999
0 <= R("R_ARGSS") <= 999999
0 <= R("R_GLYATi") <= 999999
0 <= R("R_OBTFL") <= 999999
0 <= R("R_ACPS1") <= 999999
0 <= R("R_ACS") <= 999999
0 <= R("R_PPNCL2") <= 999999
0 <= R("R_RBFSa") <= 999999
0 <= R("R_MGSA") <= 999999
0 <= R("R_DHNPA2") <= 999999
0 <= R("R_DHPTDCs") <= 999999
0 <= R("R_RBFSb") <= 999999
0 <= R("R_DMSOR2") <= 999999
0 <= R("R_DMSOR2e") <= 999999
0 <= R("R_DPCOAK") <= 999999
0 <= R("R_TDPADGAT") <= 999999
0 <= R("R_TDPDRE") <= 999999
0 <= R("R_TDPGDH") <= 999999
0 <= R("R_DURIK1") <= 999999
0 <= R("R_DXYLK") <= 999999
0 <= R("R_IGPDH") <= 999999
0 <= R("R_ETHAAL") <= 999999
0 <= R("R_FRUK") <= 999999
0 <= R("R_PFK") <= 999999
0 <= R("R_FCLK") <= 999999
0 <= R("R_PRAIS") <= 999999
0 <= R("R_XYLI2i") <= 999999
0 <= R("R_FRD3") <= 999999
0 <= R("R_ASPO5") <= 999999
0 <= R("R_GALCTD") <= 999999
0 <= R("R_GALCTND") <= 999999
0 <= R("R_G1PACT") <= 999999
0 <= R("R_GART") <= 999999
0 <= R("R_CRNt7") <= 999999
0 <= R("R_GDMANE") <= 999999
0 <= R("R_GMAND") <= 999999
0 <= R("R_HEX1") <= 999999
0 <= R("R_GNK") <= 999999
0 <= R("R_GLCRD") <= 999999
0 <= R("R_ADCS") <= 999999
0 <= R("R_ANS") <= 999999
0 <= R("R_GF6PTA") <= 999999
0 <= R("R_G1SATi") <= 999999
0 <= R("R_G5SADs") <= 999999
0 <= R("R_UNK3") <= 999999
0 <= R("R_PSERT") <= 999999
0 <= R("R_LEUTAi") <= 999999
0 <= R("R_ACGS") <= 999999
0 <= R("R_GLU5K") <= 999999
0 <= R("R_GLUCYS") <= 999999
0 <= R("R_DHFS") <= 999999
0 <= R("R_TDPAGTA") <= 999999
0 <= R("R_GTHS") <= 999999
0 <= R("R_G3PD7") <= 999999
0 <= R("R_GLYK") <= 999999
0 <= R("R_GLYCTO4") <= 999999
0 <= R("R_GLYCK") <= 999999
0 <= R("R_GMHEPK") <= 999999
0 <= R("R_GSNK") <= 999999
0 <= R("R_CYTDK2") <= 999999
0 <= R("R_CAT") <= 999999
0 <= R("R_FTHFD") <= 999999
0 <= R("R_MICITD") <= 999999
0 <= R("R_PGLYCP") <= 999999
0 <= R("R_PMDPHT") <= 999999
0 <= R("R_MTAN") <= 999999
0 <= R("R_RZ5PP") <= 999999
0 <= R("R_PGL") <= 999999
0 <= R("R_IPPS") <= 999999
0 <= R("R_NACODA") <= 999999
0 <= R("R_AGDC") <= 999999
0 <= R("R_ACODA") <= 999999
0 <= R("R_AGMT") <= 999999
0 <= R("R_LPLIPA3") <= 999999
0 <= R("R_LPLIPA2") <= 999999
0 <= R("R_LPLIPA1") <= 999999
0 <= R("R_AHCYSNS") <= 999999
0 <= R("R_DNTPPA") <= 999999
0 <= R("R_AMPMS") <= 999999
0 <= R("R_ALLTN") <= 999999
0 <= R("R_ALLTAH") <= 999999
0 <= R("R_NTD7") <= 999999
0 <= R("R_AMPN") <= 999999
0 <= R("R_AP4AH") <= 999999
0 <= R("R_AP5AH") <= 999999
0 <= R("R_AB6PGH") <= 999999
0 <= R("R_ASNN") <= 999999
0 <= R("R_ATPM") <= 999999
0 <= R("R_NTPP6") <= 999999
0 <= R("R_DAPabc") <= 999999
0 <= R("R_ALAabc") <= 999999
0 <= R("R_ARBabc") <= 999999
0 <= R("R_ARGabc") <= 999999
0 <= R("R_ASNabc") <= 999999
0 <= R("R_ASPabc") <= 999999
0 <= R("R_CBL1abc") <= 999999
0 <= R("R_CDAPPA_EC") <= 999999
0 <= R("R_CHLabc") <= 999999
0 <= R("R_NTD4") <= 999999
0 <= R("R_CMPN") <= 999999
0 <= R("R_NTPP4") <= 999999
0 <= R("R_CYSDS") <= 999999
0 <= R("R_CYSabc") <= 999999
0 <= R("R_CYSTL") <= 999999
0 <= R("R_NTD6") <= 999999
0 <= R("R_NTPP5") <= 999999
0 <= R("R_NTD3") <= 999999
0 <= R("R_NTPP3") <= 999999
0 <= R("R_NTD8") <= 999999
0 <= R("R_NTPP1") <= 999999
0 <= R("R_NTPTP1") <= 999999
0 <= R("R_DNMPPA") <= 999999
0 <= R("R_NTD5") <= 999999
0 <= R("R_NTPP7") <= 999999
0 <= R("R_NTD1") <= 999999
0 <= R("R_DUTPDP") <= 999999
0 <= R("R_DAAD") <= 999999
0 <= R("R_FBP") <= 999999
0 <= R("R_FE2abc") <= 999999
0 <= R("R_G1PP") <= 999999
0 <= R("R_GPDDA1") <= 999999
0 <= R("R_GPDDA2") <= 999999
0 <= R("R_GPDDA4") <= 999999
0 <= R("R_GPDDA5") <= 999999
0 <= R("R_GPDDA3") <= 999999
0 <= R("R_GALabc") <= 999999
0 <= R("R_G6PDA") <= 999999
0 <= R("R_GLUN") <= 999999
0 <= R("R_ASNS1") <= 999999
0 <= R("R_PRFGS") <= 999999
0 <= R("R_GLNabc") <= 999999
0 <= R("R_GLUabc") <= 999999
0 <= R("R_MALS") <= 999999
0 <= R("R_GLYBabc") <= 999999
0 <= R("R_GLYC3Pabc") <= 999999
0 <= R("R_GMHEPPA") <= 999999
0 <= R("R_NTD9") <= 999999
0 <= R("R_GP4GH") <= 999999
0 <= R("R_GTPCI") <= 999999
0 <= R("R_NTPP2") <= 999999
0 <= R("R_NTPTP2") <= 999999
0 <= R("R_GSPMDA") <= 999999
0 <= R("R_DKMPPD2") <= 999999
0 <= R("R_GTPCII2") <= 999999
0 <= R("R_CYSS") <= 999999
0 <= R("R_DAPDC") <= 999999
0 <= R("R_IGPS") <= 999999
0 <= R("R_MECDPDH") <= 999999
0 <= R("R_OMCDC") <= 999999
0 <= R("R_KG6PDC") <= 999999
0 <= R("R_OPHBDC") <= 999999
0 <= R("R_TMPPP") <= 999999
0 <= R("R_PPCDC") <= 999999
0 <= R("R_ARGDC") <= 999999
0 <= R("R_ASP1DC") <= 999999
0 <= R("R_MEPCT") <= 999999
0 <= R("R_FMNAT") <= 999999
0 <= R("R_FHL") <= 999999
0 <= R("R_GLGC") <= 999999
0 <= R("R_G1PTT") <= 999999
0 <= R("R_GLUDC") <= 999999
0 <= R("R_GLXCL") <= 999999
0 <= R("R_GMHEPAT") <= 999999
0 <= R("R_ACBIPGT") <= 999999
0 <= R("R_DHPPDA2") <= 999999
0 <= R("R_ADD") <= 999999
0 <= R("R_ADA") <= 999999
0 <= R("R_CSND") <= 999999
0 <= R("R_CYTD") <= 999999
0 <= R("R_DADA") <= 999999
0 <= R("R_DCTPD") <= 999999
0 <= R("R_DCYTD") <= 999999
0 <= R("R_GUAD") <= 999999
0 <= R("R_HYD3") <= 999999
0 <= R("R_AIRC2") <= 999999
0 <= R("R_CBPS") <= 999999
0 <= R("R_CYNTAH") <= 999999
0 <= R("R_METS") <= 999999
0 <= R("R_AACPS3") <= 999999
0 <= R("R_AACPS4") <= 999999
0 <= R("R_ABUTt2") <= 999999
0 <= R("R_ACNAMt2") <= 999999
0 <= R("R_ADNt2") <= 999999
0 <= R("R_ASPt2") <= 999999
0 <= R("R_CSNt2") <= 999999
0 <= R("R_CYNTt2") <= 999999
0 <= R("R_CYTDt2") <= 999999
0 <= R("R_DADNt2") <= 999999
0 <= R("R_DCYTt2") <= 999999
0 <= R("R_DGSNt2") <= 999999
0 <= R("R_DINSt2") <= 999999
0 <= R("R_DURIt2") <= 999999
0 <= R("R_GALt2") <= 999999
0 <= R("R_GLCt2") <= 999999
0 <= R("R_GSNt2") <= 999999
0 <= R("R_GUAt2") <= 999999
0 <= R("R_HDCAt2") <= 999999
0 <= R("R_ASPt2_2") <= 999999
0 <= R("R_FUMt2_2") <= 999999
0 <= R("R_ASPt2_3") <= 999999
0 <= R("R_FUMt2_3") <= 999999
0 <= R("R_HISabc") <= 999999
0 <= R("R_HISTP") <= 999999
0 <= R("R_HKNDDH") <= 999999
0 <= R("R_HKNTDH") <= 999999
0 <= R("R_UPP3S") <= 999999
0 <= R("R_HSK") <= 999999
0 <= R("R_QULNS") <= 999999
0 <= R("R_ICHORT") <= 999999
0 <= R("R_ICL") <= 999999
0 <= R("R_ILEabc") <= 999999
0 <= R("R_HSTPT") <= 999999
0 <= R("R_ADSS") <= 999999
0 <= R("R_NTD11") <= 999999
0 <= R("R_INSK") <= 999999
0 <= R("R_INSt2") <= 999999
0 <= R("R_IPDDIi") <= 999999
0 <= R("R_DMATT") <= 999999
0 <= R("R_GRTT") <= 999999
0 <= R("R_OCTDPS") <= 999999
0 <= R("R_UDCPDPS") <= 999999
0 <= R("R_EDTXS1") <= 999999
0 <= R("R_EDTXS3") <= 999999
0 <= R("R_KDOPP") <= 999999
0 <= R("R_KDOCT2") <= 999999
0 <= R("R_MOAT2") <= 999999
0 <= R("R_Kabc") <= 999999
0 <= R("R_LACZ") <= 999999
0 <= R("R_LEUabc") <= 999999
0 <= R("R_GLYOX") <= 999999
0 <= R("R_MOAT") <= 999999
0 <= R("R_TDSK") <= 999999
0 <= R("R_LYSDC") <= 999999
0 <= R("R_LYSabc") <= 999999
0 <= R("R_CADVt") <= 999999
0 <= R("R_MACPD") <= 999999
0 <= R("R_KAS14") <= 999999
0 <= R("R_KAS15") <= 999999
0 <= R("R_MALt2_2") <= 999999
0 <= R("R_MALt2_3") <= 999999
0 <= R("R_MALTabc") <= 999999
0 <= R("R_MLTG5") <= 999999
0 <= R("R_MLTG4") <= 999999
0 <= R("R_AMALT4") <= 999999
0 <= R("R_MALTHXabc") <= 999999
0 <= R("R_MLTG3") <= 999999
0 <= R("R_AMALT3") <= 999999
0 <= R("R_MALTPTabc") <= 999999
0 <= R("R_MLTG1") <= 999999
0 <= R("R_AMALT1") <= 999999
0 <= R("R_MALTTRabc") <= 999999
0 <= R("R_MLTG2") <= 999999
0 <= R("R_AMALT2") <= 999999
0 <= R("R_MALTTTRabc") <= 999999
0 <= R("R_MAN1PT2") <= 999999
0 <= R("R_MNNH") <= 999999
0 <= R("R_GALS3") <= 999999
0 <= R("R_MELIBt2") <= 999999
0 <= R("R_METDabc") <= 999999
0 <= R("R_METAT") <= 999999
0 <= R("R_METabc") <= 999999
0 <= R("R_MI1PP") <= 999999
0 <= R("R_TMDS") <= 999999
0 <= R("R_MOHMT") <= 999999
0 <= R("R_MMCD") <= 999999
0 <= R("R_DMSOR1") <= 999999
0 <= R("R_DMSOR1e") <= 999999
0 <= R("R_FRD2") <= 999999
0 <= R("R_ASPO4") <= 999999
0 <= R("R_DHORD5") <= 999999
0 <= R("R_G3PD6") <= 999999
0 <= R("R_GLYCTO3") <= 999999
0 <= R("R_HYD2") <= 999999
0 <= R("R_FDH3") <= 999999
0 <= R("R_L_LACD3") <= 999999
0 <= R("R_MDH3") <= 999999
0 <= R("R_LGTHL") <= 999999
0 <= R("R_EDTXS2") <= 999999
0 <= R("R_EDTXS4") <= 999999
0 <= R("R_NAt3_2") <= 999999
0 <= R("R_NAt3_15") <= 999999
0 <= R("R_GLUt4") <= 999999
0 <= R("R_NACUP") <= 999999
0 <= R("R_IPMD") <= 999999
0 <= R("R_PGCD") <= 999999
0 <= R("R_NADK") <= 999999
0 <= R("R_DHPPD") <= 999999
0 <= R("R_DHCIND") <= 999999
0 <= R("R_ACALDi") <= 999999
0 <= R("R_AKGDH") <= 999999
0 <= R("R_SHCHD2") <= 999999
0 <= R("R_GLYCDx") <= 999999
0 <= R("R_NADDP") <= 999999
0 <= R("R_ABUTD") <= 999999
0 <= R("R_ALDD2x") <= 999999
0 <= R("R_BETALDHx") <= 999999
0 <= R("R_FAO4") <= 999999
0 <= R("R_GCALDD") <= 999999
0 <= R("R_P5CD") <= 999999
0 <= R("R_IMPD") <= 999999
0 <= R("R_LCADi") <= 999999
0 <= R("R_ME1") <= 999999
0 <= R("R_HISTD") <= 999999
0 <= R("R_FAO2") <= 999999
0 <= R("R_NADDPe") <= 999999
0 <= R("R_DOGULNR") <= 999999
0 <= R("R_DKGLCNR2x") <= 999999
0 <= R("R_2DGLCNRx") <= 999999
0 <= R("R_2DGULRx") <= 999999
0 <= R("R_NADH9") <= 999999
0 <= R("R_TRSAR") <= 999999
0 <= R("R_BSORx") <= 999999
0 <= R("R_GLYCLTDx") <= 999999
0 <= R("R_DMPPS") <= 999999
0 <= R("R_IPDPS") <= 999999
0 <= R("R_NADH8") <= 999999
0 <= R("R_HPYRRx") <= 999999
0 <= R("R_MTHFR2") <= 999999
0 <= R("R_NADH10") <= 999999
0 <= R("R_NADH7") <= 999999
0 <= R("R_GND") <= 999999
0 <= R("R_NADPPPS") <= 999999
0 <= R("R_BETALDHy") <= 999999
0 <= R("R_ME2") <= 999999
0 <= R("R_THD2") <= 999999
0 <= R("R_DHDPRy") <= 999999
0 <= R("R_DKGLCNR1") <= 999999
0 <= R("R_DKGLCNR2y") <= 999999
0 <= R("R_KARA2i") <= 999999
0 <= R("R_2DGLCNRy") <= 999999
0 <= R("R_2DGULRy") <= 999999
0 <= R("R_DPR") <= 999999
0 <= R("R_APRAUR") <= 999999
0 <= R("R_IDOND2") <= 999999
0 <= R("R_KARA1i") <= 999999
0 <= R("R_BSORy") <= 999999
0 <= R("R_TDPDRR") <= 999999
0 <= R("R_DXPRIi") <= 999999
0 <= R("R_GOFUCR") <= 999999
0 <= R("R_GLUSy") <= 999999
0 <= R("R_G5SD") <= 999999
0 <= R("R_GLUTRR") <= 999999
0 <= R("R_GLYCLTDy") <= 999999
0 <= R("R_P5CR") <= 999999
0 <= R("R_GMPR") <= 999999
0 <= R("R_HPYRRy") <= 999999
0 <= R("R_KAS16") <= 999999
0 <= R("R_NADTRHD") <= 999999
0 <= R("R_C120SN") <= 999999
0 <= R("R_C141SN") <= 999999
0 <= R("R_C140SN") <= 999999
0 <= R("R_C161SN") <= 999999
0 <= R("R_C160SN") <= 999999
0 <= R("R_C181SN") <= 999999
0 <= R("R_NNAM") <= 999999
0 <= R("R_ASNS2") <= 999999
0 <= R("R_CBMK") <= 999999
0 <= R("R_NADS1") <= 999999
0 <= R("R_GLNS") <= 999999
0 <= R("R_NNDMBRT") <= 999999
0 <= R("R_NNAT") <= 999999
0 <= R("R_NMNDA") <= 999999
0 <= R("R_NMNN") <= 999999
0 <= R("R_NMNAT") <= 999999
0 <= R("R_NMNP") <= 999999
0 <= R("R_NMNt7") <= 999999
0 <= R("R_NTRIR2x") <= 999999
0 <= R("R_NO3R2") <= 999999
0 <= R("R_NO3t7") <= 999999
0 <= R("R_OMMBLHX") <= 999999
0 <= R("R_OMPHHX") <= 999999
0 <= R("R_OPHHX") <= 999999
0 <= R("R_ASPO6") <= 999999
0 <= R("R_DHCINDO") <= 999999
0 <= R("R_HPPPNDO") <= 999999
0 <= R("R_DKMPPD") <= 999999
0 <= R("R_CPPPGO") <= 999999
0 <= R("R_3HCINNMH") <= 999999
0 <= R("R_3HPPPNH") <= 999999
0 <= R("R_CINNDO") <= 999999
0 <= R("R_SPODM") <= 999999
0 <= R("R_PPCK") <= 999999
0 <= R("R_CS") <= 999999
0 <= R("R_FAO3") <= 999999
0 <= R("R_OCDCAt2") <= 999999
0 <= R("R_AACPS5") <= 999999
0 <= R("R_HBZOPT") <= 999999
0 <= R("R_DHNAOT") <= 999999
0 <= R("R_OP4ENH") <= 999999
0 <= R("R_ORNTA") <= 999999
0 <= R("R_ORNDC") <= 999999
0 <= R("R_ORNabc") <= 999999
0 <= R("R_OMPDC") <= 999999
0 <= R("R_ALDD19x") <= 999999
0 <= R("R_PACCOAL") <= 999999
0 <= R("R_PAPA_EC") <= 999999
0 <= R("R_PTPATi") <= 999999
0 <= R("R_PANTS") <= 999999
0 <= R("R_BPNT") <= 999999
0 <= R("R_PLIPA3") <= 999999
0 <= R("R_PDXPP") <= 999999
0 <= R("R_PEAMNO") <= 999999
0 <= R("R_PLIPA2") <= 999999
0 <= R("R_ACGApts") <= 999999
0 <= R("R_ACMANApts") <= 999999
0 <= R("R_DHAPT") <= 999999
0 <= R("R_FRUpts") <= 999999
0 <= R("R_FRUpts2") <= 999999
0 <= R("R_GALTpts") <= 999999
0 <= R("R_GAMpts") <= 999999
0 <= R("R_GLCpts") <= 999999
0 <= R("R_KDOPS") <= 999999
0 <= R("R_PPC") <= 999999
0 <= R("R_DDPA") <= 999999
0 <= R("R_PYK") <= 999999
0 <= R("R_MALTpts") <= 999999
0 <= R("R_MANpts") <= 999999
0 <= R("R_MNLpts") <= 999999
0 <= R("R_LPLIPA5") <= 999999
0 <= R("R_LPLIPA4") <= 999999
0 <= R("R_LPLIPA6") <= 999999
0 <= R("R_PLIPA1") <= 999999
0 <= R("R_PGPP_EC") <= 999999
0 <= R("R_HEMEOS") <= 999999
0 <= R("R_THRS") <= 999999
0 <= R("R_4HTHRS") <= 999999
0 <= R("R_PDX5PS") <= 999999
0 <= R("R_GLYC3Pt6") <= 999999
0 <= R("R_GLCP") <= 999999
0 <= R("R_FUCPt6_2") <= 999999
0 <= R("R_G6Pt6_2") <= 999999
0 <= R("R_MAN6Pt6_2") <= 999999
0 <= R("R_PIabc") <= 999999
0 <= R("R_PNTK") <= 999999
0 <= R("R_PNTOt4") <= 999999
0 <= R("R_ACCOAL") <= 999999
0 <= R("R_HMBS") <= 999999
0 <= R("R_MCITS") <= 999999
0 <= R("R_PTA2") <= 999999
0 <= R("R_PPNDH") <= 999999
0 <= R("R_PPND") <= 999999
0 <= R("R_PPA") <= 999999
0 <= R("R_FCLT") <= 999999
0 <= R("R_PPPGO") <= 999999
0 <= R("R_PPPNDO") <= 999999
0 <= R("R_PRAIi") <= 999999
0 <= R("R_PRAMPC") <= 999999
0 <= R("R_PRATPP") <= 999999
0 <= R("R_PRMICIi") <= 999999
0 <= R("R_IG3PS") <= 999999
0 <= R("R_PROD2") <= 999999
0 <= R("R_PROabc") <= 999999
0 <= R("R_PROt4") <= 999999
0 <= R("R_ADPT") <= 999999
0 <= R("R_ANPRT") <= 999999
0 <= R("R_ATPPRT") <= 999999
0 <= R("R_GUAPRT") <= 999999
0 <= R("R_GLUPRT") <= 999999
0 <= R("R_HXPRT") <= 999999
0 <= R("R_NAMNPP") <= 999999
0 <= R("R_PSD_EC") <= 999999
0 <= R("R_PSP_L") <= 999999
0 <= R("R_PTRCTA") <= 999999
0 <= R("R_SPMS") <= 999999
0 <= R("R_PTRCabc") <= 999999
0 <= R("R_HYPOE") <= 999999
0 <= R("R_PYAM5PO") <= 999999
0 <= R("R_PYDAMK") <= 999999
0 <= R("R_ALATA_D2") <= 999999
0 <= R("R_ALATA_L2") <= 999999
0 <= R("R_PYDXPP") <= 999999
0 <= R("R_PYDXK") <= 999999
0 <= R("R_PYDXNK") <= 999999
0 <= R("R_DHDPS") <= 999999
0 <= R("R_PFL") <= 999999
0 <= R("R_PPS") <= 999999
0 <= R("R_ACHBS") <= 999999
0 <= R("R_DXPS") <= 999999
0 <= R("R_PDH") <= 999999
0 <= R("R_ACLS") <= 999999
0 <= R("R_ASPO3") <= 999999
0 <= R("R_DHORD2") <= 999999
0 <= R("R_G3PD5") <= 999999
0 <= R("R_GLYCTO2") <= 999999
0 <= R("R_GLCDe") <= 999999
0 <= R("R_HYD1") <= 999999
0 <= R("R_FDH2") <= 999999
0 <= R("R_LDH_D2") <= 999999
0 <= R("R_L_LACD2") <= 999999
0 <= R("R_MDH2") <= 999999
0 <= R("R_NADH5") <= 999999
0 <= R("R_NADH6") <= 999999
0 <= R("R_POX") <= 999999
0 <= R("R_NO3R1") <= 999999
0 <= R("R_CYTBD") <= 999999
0 <= R("R_CYTBO3") <= 999999
0 <= R("R_NNDPR") <= 999999
0 <= R("R_RBK_L1") <= 999999
0 <= R("R_ADOCBLS") <= 999999
0 <= R("R_RHCCE") <= 999999
0 <= R("R_RBK") <= 999999
0 <= R("R_RIBabc") <= 999999
0 <= R("R_RBFK") <= 999999
0 <= R("R_RMK") <= 999999
0 <= R("R_RMNt") <= 999999
0 <= R("R_DB4PS") <= 999999
0 <= R("R_S7PI") <= 999999
0 <= R("R_SBTpts") <= 999999
0 <= R("R_NPHS") <= 999999
0 <= R("R_SHCHF") <= 999999
0 <= R("R_SELNPS") <= 999999
0 <= R("R_ENTCS") <= 999999
0 <= R("R_SERD_D") <= 999999
0 <= R("R_SERD_L") <= 999999
0 <= R("R_TRPS1") <= 999999
0 <= R("R_TRPS2") <= 999999
0 <= R("R_SERt4") <= 999999
0 <= R("R_SHKK") <= 999999
0 <= R("R_SDPDS") <= 999999
0 <= R("R_SADT2") <= 999999
0 <= R("R_SULabc") <= 999999
0 <= R("R_SPMDAT1") <= 999999
0 <= R("R_SPMDAT2") <= 999999
0 <= R("R_GSPMDS") <= 999999
0 <= R("R_SPMDabc") <= 999999
0 <= R("R_SHCHCS2") <= 999999
0 <= R("R_FFSD") <= 999999
0 <= R("R_SADH") <= 999999
0 <= R("R_SUCBZL") <= 999999
0 <= R("R_CITt7") <= 999999
0 <= R("R_SUCD1i") <= 999999
0 <= R("R_SUCCt2b") <= 999999
0 <= R("R_PPCSCT") <= 999999
0 <= R("R_SUCCabc") <= 999999
0 <= R("R_SUCCt2_2") <= 999999
0 <= R("R_SUCCt2_3") <= 999999
0 <= R("R_MMM2") <= 999999
0 <= R("R_AST") <= 999999
0 <= R("R_HSST") <= 999999
0 <= R("R_SGDS") <= 999999
0 <= R("R_SGSAD") <= 999999
0 <= R("R_SHSL1") <= 999999
0 <= R("R_SOTA") <= 999999
0 <= R("R_SUCpts") <= 999999
0 <= R("R_SSALx") <= 999999
0 <= R("R_SSALy") <= 999999
0 <= R("R_PFK_2") <= 999999
0 <= R("R_TARTD") <= 999999
0 <= R("R_TAUDO") <= 999999
0 <= R("R_TAURabc") <= 999999
0 <= R("R_PASYN_EC") <= 999999
0 <= R("R_THDPS") <= 999999
0 <= R("R_GLYCL") <= 999999
0 <= R("R_GHMT2") <= 999999
0 <= R("R_THMabc") <= 999999
0 <= R("R_OXGDC2") <= 999999
0 <= R("R_THRD_L") <= 999999
0 <= R("R_THRD") <= 999999
0 <= R("R_THRabc") <= 999999
0 <= R("R_THRt4") <= 999999
0 <= R("R_TMDK1") <= 999999
0 <= R("R_THMDt2") <= 999999
0 <= R("R_TMAOR2") <= 999999
0 <= R("R_TMAOR1") <= 999999
0 <= R("R_TMAOR2e") <= 999999
0 <= R("R_TMAOR1e") <= 999999
0 <= R("R_TRDR") <= 999999
0 <= R("R_RNDR1") <= 999999
0 <= R("R_RNTR1") <= 999999
0 <= R("R_RNDR3") <= 999999
0 <= R("R_RNTR3") <= 999999
0 <= R("R_RNDR2") <= 999999
0 <= R("R_RNTR2") <= 999999
0 <= R("R_PAPSR") <= 999999
0 <= R("R_TRE6PH") <= 999999
0 <= R("R_TRE6PP") <= 999999
0 <= R("R_TREH") <= 999999
0 <= R("R_TREHe") <= 999999
0 <= R("R_TREpts") <= 999999
0 <= R("R_GLUTRS") <= 999999
0 <= R("R_CYANST") <= 999999
0 <= R("R_TSULabc") <= 999999
0 <= R("R_AACPS1") <= 999999
0 <= R("R_FAO1") <= 999999
0 <= R("R_TTDCAt2") <= 999999
0 <= R("R_AACPS2") <= 999999
0 <= R("R_THZPSN") <= 999999
0 <= R("R_USHD") <= 999999
0 <= R("R_LPADSS") <= 999999
0 <= R("R_UHGADA") <= 999999
0 <= R("R_U23GAAT") <= 999999
0 <= R("R_PPTGS") <= 999999
0 <= R("R_UAPGR") <= 999999
0 <= R("R_UAG2Ei") <= 999999
0 <= R("R_UAGCVT") <= 999999
0 <= R("R_UACMAMO") <= 999999
0 <= R("R_UAGPT3") <= 999999
0 <= R("R_UAMAGS") <= 999999
0 <= R("R_UAAGDS") <= 999999
0 <= R("R_UAMAS") <= 999999
0 <= R("R_UDCPDP") <= 999999
0 <= R("R_ACGAMT") <= 999999
0 <= R("R_RNDR4") <= 999999
0 <= R("R_UDPGALM") <= 999999
0 <= R("R_TRE6PS") <= 999999
0 <= R("R_UDPGD") <= 999999
0 <= R("R_LPSSYN_EC") <= 999999
0 <= R("R_PAPPT3") <= 999999
0 <= R("R_UGMDDS") <= 999999
0 <= R("R_NTD2") <= 999999
0 <= R("R_ACMAMUT") <= 999999
0 <= R("R_AADDGT") <= 999999
0 <= R("R_ECAP_EC") <= 999999
0 <= R("R_UPP3MT") <= 999999
0 <= R("R_UPPDC1") <= 999999
0 <= R("R_UPPRT") <= 999999
0 <= R("R_URAt2") <= 999999
0 <= R("R_UGLYCH") <= 999999
0 <= R("R_URIK2") <= 999999
0 <= R("R_URIt2") <= 999999
0 <= R("R_NTPP8") <= 999999
0 <= R("R_CTPS2") <= 999999
0 <= R("R_UAGDP") <= 999999
0 <= R("R_GALUi") <= 999999
0 <= R("R_RNTR4") <= 999999
0 <= R("R_BiomassEcoli") <= 999999
0 <= R("R_VALabc") <= 999999
0 <= R("R_XPPT") <= 999999
0 <= R("R_XANt2") <= 999999
0 <= R("R_NTD10") <= 999999
0 <= R("R_GMPS2") <= 999999
0 <= R("R_X5PL3E") <= 999999
0 <= R("R_XYLabc") <= 999999
0 <= R("R_XYLt2") <= 999999
0 <= R("R_XYLK") <= 999999


\ Some model modifications follow

-999999 <= R("M12ppdSb_Transp") <= 0
-999999 <= R("M15dapb_Transp") <= 0
-999999 <= R("M26dapMb_Transp") <= 0
-999999 <= R("M2ddglcnb_Transp") <= 0
-999999 <= R("M3hcinnmb_Transp") <= 0
-999999 <= R("M3hpppnb_Transp") <= 0
-999999 <= R("M4abutb_Transp") <= 0
-999999 <= R("Macacb_Transp") <= 0
-999999 <= R("Macaldb_Transp") <= 0
-999999 <= R("Macb_Transp") <= 0
-999999 <= R("Macgamb_Transp") <= 0
-999999 <= R("Macmanab_Transp") <= 0
-999999 <= R("Macnamb_Transp") <= 0
-999999 <= R("Madeb_Transp") <= 0
-999999 <= R("Madnb_Transp") <= 0
-999999 <= R("Makgb_Transp") <= 0
-999999 <= R("MalaDb_Transp") <= 0
-999999 <= R("MalaLb_Transp") <= 0
-999999 <= R("Malltnb_Transp") <= 0
-999999 <= R("Mampb_Transp") <= 0
-999999 <= R("MarabLb_Transp") <= 0
-999999 <= R("MargLb_Transp") <= 0
-999999 <= R("MasnLb_Transp") <= 0
-999999 <= R("MaspLb_Transp") <= 0
-999999 <= R("Mbutb_Transp") <= 0
-999999 <= R("Mcbl1b_Transp") <= 0
-999999 <= R("Mcholb_Transp") <= 0
-999999 <= R("Mcitb_Transp") <= 0
-999999 <= R("Mco2b_Transp") <= 0
-999999 <= R("Mcrnb_Transp") <= 0
-999999 <= R("Mcsnb_Transp") <= 0
-999999 <= R("Mcyntb_Transp") <= 0
-999999 <= R("McysLb_Transp") <= 0
-999999 <= R("Mcytdb_Transp") <= 0
-999999 <= R("Mdad2b_Transp") <= 0
-999999 <= R("Mdcytb_Transp") <= 0
-999999 <= R("Mdgsnb_Transp") <= 0
-999999 <= R("Mdhab_Transp") <= 0
-999999 <= R("Mdinb_Transp") <= 0
-999999 <= R("Mdmsb_Transp") <= 0
-999999 <= R("Mdmsob_Transp") <= 0
-999999 <= R("Mdurib_Transp") <= 0
-999999 <= R("Metohb_Transp") <= 0
-999999 <= R("Mfe2b_Transp") <= 0
-999999 <= R("Mforb_Transp") <= 0
-999999 <= R("Mfrub_Transp") <= 0
-999999 <= R("Mfuc1pLb_Transp") <= 0
-999999 <= R("MfucLb_Transp") <= 0
-999999 <= R("Mfumb_Transp") <= 0
-999999 <= R("Mg6pb_Transp") <= 0
-999999 <= R("Mgalb_Transp") <= 0
-999999 <= R("MgalctDb_Transp") <= 0
-999999 <= R("MgalctnDb_Transp") <= 0
-999999 <= R("Mgaltb_Transp") <= 0
-999999 <= R("Mgalurb_Transp") <= 0
-999999 <= R("Mgamb_Transp") <= 0
-999999 <= R("Mgbbtnb_Transp") <= 0
-999999 <= R("MglcDb_Transp") <= 0
-999999 <= R("Mglcnb_Transp") <= 0
-999999 <= R("Mglcrb_Transp") <= 0
-999999 <= R("Mglcurb_Transp") <= 0
-999999 <= R("MglnLb_Transp") <= 0
-999999 <= R("MgluLb_Transp") <= 0
-999999 <= R("Mglyaldb_Transp") <= 0
-999999 <= R("Mglyb_Transp") <= 0
-999999 <= R("Mglybb_Transp") <= 0
-999999 <= R("Mglyc3pb_Transp") <= 0
-999999 <= R("Mglycb_Transp") <= 0
-999999 <= R("Mglycltb_Transp") <= 0
-999999 <= R("Mgsnb_Transp") <= 0
-999999 <= R("Mguab_Transp") <= 0
-999999 <= R("Mh2ob_Transp") <= 0
-999999 <= R("Mhb_Transp") <= 0
-999999 <= R("Mhdcab_Transp") <= 0
-999999 <= R("MhisLb_Transp") <= 0
-999999 <= R("Mhxanb_Transp") <= 0
-999999 <= R("MidonLb_Transp") <= 0
-999999 <= R("MileLb_Transp") <= 0
-999999 <= R("Mindoleb_Transp") <= 0
-999999 <= R("Minsb_Transp") <= 0
-999999 <= R("Mkb_Transp") <= 0
-999999 <= R("MlacDb_Transp") <= 0
-999999 <= R("MlacLb_Transp") <= 0
-999999 <= R("Mlctsb_Transp") <= 0
-999999 <= R("MleuLb_Transp") <= 0
-999999 <= R("MlysLb_Transp") <= 0
-999999 <= R("MmalLb_Transp") <= 0
-999999 <= R("Mmaltb_Transp") <= 0
-999999 <= R("Mmalthxb_Transp") <= 0
-999999 <= R("Mmaltptb_Transp") <= 0
-999999 <= R("Mmalttrb_Transp") <= 0
-999999 <= R("Mmaltttrb_Transp") <= 0
-999999 <= R("Mman6pb_Transp") <= 0
-999999 <= R("Mmanb_Transp") <= 0
-999999 <= R("Mmelibb_Transp") <= 0
-999999 <= R("MmetDb_Transp") <= 0
-999999 <= R("MmetLb_Transp") <= 0
-999999 <= R("Mmnlb_Transp") <= 0
-999999 <= R("Mna1b_Transp") <= 0
-999999 <= R("Mnacb_Transp") <= 0
-999999 <= R("Mnadb_Transp") <= 0
-999999 <= R("Mnh4b_Transp") <= 0
-999999 <= R("Mnmnb_Transp") <= 0
-999999 <= R("Mno2b_Transp") <= 0
-999999 <= R("Mno3b_Transp") <= 0
-999999 <= R("Mo2b_Transp") <= 0
-999999 <= R("Mocdcab_Transp") <= 0
-999999 <= R("Mornb_Transp") <= 0
-999999 <= R("MpheLb_Transp") <= 0
-999999 <= R("Mpib_Transp") <= 0
-999999 <= R("MpntoRb_Transp") <= 0
-999999 <= R("Mpppnb_Transp") <= 0
-999999 <= R("MproLb_Transp") <= 0
-999999 <= R("Mptrcb_Transp") <= 0
-999999 <= R("Mpyrb_Transp") <= 0
-999999 <= R("MribDb_Transp") <= 0
-999999 <= R("Mrmnb_Transp") <= 0
-999999 <= R("MsbtDb_Transp") <= 0
-999999 <= R("MserDb_Transp") <= 0
-999999 <= R("MserLb_Transp") <= 0
-999999 <= R("Mso4b_Transp") <= 0
-999999 <= R("Mspmdb_Transp") <= 0
-999999 <= R("Msuccb_Transp") <= 0
-999999 <= R("Msucrb_Transp") <= 0
-999999 <= R("MtartrLb_Transp") <= 0
-999999 <= R("Mtaurb_Transp") <= 0
-999999 <= R("Mthmb_Transp") <= 0
-999999 <= R("MthrLb_Transp") <= 0
-999999 <= R("Mthymdb_Transp") <= 0
-999999 <= R("Mtmab_Transp") <= 0
-999999 <= R("Mtmaob_Transp") <= 0
-999999 <= R("Mtreb_Transp") <= 0
-999999 <= R("MtrpLb_Transp") <= 0
-999999 <= R("Mtsulb_Transp") <= 0
-999999 <= R("Mttdcab_Transp") <= 0
-999999 <= R("MtyrLb_Transp") <= 0
-999999 <= R("Murab_Transp") <= 0
-999999 <= R("Mureab_Transp") <= 0
-999999 <= R("Murib_Transp") <= 0
-999999 <= R("MvalLb_Transp") <= 0
-999999 <= R("Mxanb_Transp") <= 0
-999999 <= R("Mxtsnb_Transp") <= 0
-999999 <= R("MxylDb_Transp") <= 0
