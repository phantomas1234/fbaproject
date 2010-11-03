\* Problem: Unknown *\

Maximize
 Z: + R("R_BiomassEcoli")

Subject To
 M10fthfc: + R("R_MTHFC") - R("R_GARFT") - R("R_FTHFD") - R("R_AICART")
 = 0
 M12dgrECc: + 0.02 R("R_PEPT_EC") + 0.02 R("R_PAPA_EC")
 - 0.02 R("R_DAGK_EC") = 0
 M12ppdSb: + R("R_EX_12ppd_S_e_") + R("M12ppdSb_Transp") = 0
 M12ppdSc: + R("R_LCARS") + R("R_12PPDt") = 0
 M12ppdSe: - R("R_12PPDt") - R("R_EX_12ppd_S_e_") = 0
 M13dpgc: + R("R_PGK") + R("R_GAPD") = 0
 M15dapb: + R("R_EX_15dap_e_") + R("M15dapb_Transp") = 0
 M15dapc: + R("R_LYSDC") - R("R_CADVt") = 0
 M15dape: + R("R_CADVt") - R("R_EX_15dap_e_") = 0
 M1pyr5cc: + R("R_PROD2") - R("R_P5CR") - R("R_P5CD") + R("R_G5SADs")
 = 0
 M23ddhbc: + R("R_ICHORT") - R("R_DHBD") = 0
 M23dhbac: - 3 R("R_ENTCS") + R("R_DHBSr") = 0
 M23dhbc: - R("R_DHBSr") + R("R_DHBD") = 0
 M23dhdpc: + R("R_DHDPS") - R("R_DHDPRy") = 0
 M23dhmbc: + R("R_KARA1i") - R("R_DHAD1") = 0
 M23dhmpc: + R("R_KARA2i") - R("R_DHAD2") = 0
 M23dogulnc: - R("R_DOGULNR") = 0
 M25aicsc: + R("R_PRASCS") - R("R_ADSL2r") = 0
 M25dkglcnc: - R("R_DKGLCNR2y") - R("R_DKGLCNR2x") - R("R_DKGLCNR1") = 0
 M25drappc: + R("R_GTPCII2") - R("R_DHPPDA2") = 0
 M26dapLLc: + R("R_SDPDS") - R("R_DAPE") = 0
 M26dapMb: + R("R_EX_26dap_M_e_") + R("M26dapMb_Transp") = 0
 M26dapMc: - R("R_UAAGDS") - R("R_DAPDC") + R("R_DAPabc") + R("R_DAPE")
 = 0
 M26dapMe: - R("R_DAPabc") - R("R_EX_26dap_M_e_") = 0
 M2ahbutc: + R("R_ACHBS") - R("R_KARA2i") = 0
 M2aobutc: + R("R_THRD") - R("R_GLYATi") = 0
 M2cpr5pc: + R("R_PRAIi") - R("R_IGPS") = 0
 M2dda7pc: - R("R_DHQS") + R("R_DDPA") = 0
 M2ddg6pc: + R("R_EDD") - R("R_EDA") + R("R_DDGLK") = 0
 M2ddglcnb: + R("R_EX_2ddglcn_e_") + R("M2ddglcnb_Transp") = 0
 M2ddglcnc: + R("R_MNNH") + R("R_DDGLCNt2r") + R("R_ALTRH")
 - R("R_DDGLK") = 0
 M2ddglcne: - R("R_DDGLCNt2r") - R("R_EX_2ddglcn_e_") = 0
 M2dh3dgal6pc: - R("R_DDPGALA") + R("R_DDGALK") = 0
 M2dh3dgalc: + R("R_GALCTND") - R("R_DDGALK") = 0
 M2dhglcnc: - R("R_2DGLCNRy") - R("R_2DGLCNRx") = 0
 M2dhgulnc: - R("R_2DGULRy") - R("R_2DGULRx") + R("R_DKGLCNR1") = 0
 M2dhpc: + R("R_MOHMT") - R("R_DPR") = 0
 M2dmmq8c: + R("R_TMAOR2e") + R("R_TMAOR2") - R("R_NADH9")
 - R("R_NADH8") - R("R_HYD3") - R("R_GLYCTO4") - R("R_G3PD7")
 + R("R_FRD3") + R("R_DMSOR2e") + R("R_DMSOR2") + R("R_DHNAOT")
 - R("R_AMMQT8_2") = 0
 M2dmmql8c: - R("R_TMAOR2e") - R("R_TMAOR2") + R("R_NADH9")
 + R("R_NADH8") + R("R_HYD3") + R("R_GLYCTO4") + R("R_G3PD7")
 - R("R_FRD3") - R("R_DMSOR2e") - R("R_DMSOR2") = 0
 M2dr1pc: + R("R_TMDPP") + R("R_PUNP6") + R("R_PUNP4") + R("R_PUNP2")
 - R("R_PPM2") + R("R_DURIPP") = 0
 M2dr5pc: - R("R_DRPA") + R("R_PPM2") = 0
 M2h3oppanc: - R("R_TRSAR") + R("R_HPYRI") + R("R_GLXCL")
 + R("R_GLCRAL") = 0
 M2ippmc: - R("R_IPPMIb") + R("R_IPPMIa") = 0
 M2kmbc: - R("R_UNK3") + R("R_DKMPPD2") + R("R_DKMPPD") = 0
 M2mahmpc: - R("R_TMPPP") + R("R_PMPK") = 0
 M2mcacnc: - R("R_MICITD") + R("R_MCITD") = 0
 M2mcitc: + R("R_MCITS") - R("R_MCITD") = 0
 M2me4pc: - R("R_MEPCT") + R("R_DXPRIi") = 0
 M2mecdpc: + R("R_MECDPS") - R("R_MECDPDH") = 0
 M2obutc: + R("R_THRD_L") - R("R_OBTFL") - R("R_ACHBS") = 0
 M2ohphc: + R("R_OPHHX") - R("R_OHPHM") = 0
 M2ombzlc: + R("R_OMPHHX") - R("R_OMBZLM") = 0
 M2omhmblc: + R("R_OMMBLHX") - R("R_DMQMT") = 0
 M2ommblc: - R("R_OMMBLHX") + R("R_OMBZLM") = 0
 M2omphc: - R("R_OMPHHX") + R("R_OHPHM") = 0
 M2ophc: + R("R_OPHBDC") - R("R_OPHHX") = 0
 M2p4c2mec: + R("R_CDPMEK") - R("R_MECDPS") = 0
 M2pgc: - R("R_PGM") - R("R_ENO") = 0
 M2pglycc: - R("R_PGLYCP") = 0
 M2shchcc: - R("R_SUCBZS") + R("R_SHCHCS2") = 0
 M34hppc: + R("R_TYRTA") + R("R_PPND") = 0
 M3c2hmpc: - R("R_IPMD") - R("R_IPPMIa") = 0
 M3c3hmpc: + R("R_IPPS") + R("R_IPPMIb") = 0
 M3c4mopc: - R("R_OMCDC") + R("R_IPMD") = 0
 M3dhgulnc: + R("R_DOGULNR") = 0
 M3dhgulnpc: - R("R_KG6PDC") = 0
 M3dhqc: - R("R_DHQD") + R("R_DHQS") = 0
 M3dhskc: - R("R_SHK3Dr") + R("R_DHQD") = 0
 M3hcinnmb: + R("R_EX_3hcinnm_e_") + R("M3hcinnmb_Transp") = 0
 M3hcinnmc: + R("R_HCINNMt2r") - R("R_3HCINNMH") = 0
 M3hcinnme: - R("R_HCINNMt2r") - R("R_EX_3hcinnm_e_") = 0
 M3hmrsACPc: - R("R_UAGAAT") - R("R_U23GAAT") + R("R_KAS16") = 0
 M3hpppnb: + R("R_EX_3hpppn_e_") + R("M3hpppnb_Transp") = 0
 M3hpppnc: + R("R_HPPPNt2r") - R("R_3HPPPNH") = 0
 M3hpppne: - R("R_HPPPNt2r") - R("R_EX_3hpppn_e_") = 0
 M3ig3pc: - R("R_TRPS3") - R("R_TRPS1") + R("R_IGPS") = 0
 M3mobc: - R("R_VPAMT") + R("R_VALTA") - R("R_IPPS") - R("R_MOHMT")
 + R("R_DHAD1") = 0
 M3mopc: + R("R_ILETA") + R("R_DHAD2") = 0
 M3ophbc: + R("R_HBZOPT") - R("R_OPHBDC") = 0
 M3pgc: - R("R_PGCD") + R("R_GLYCK") + R("R_PGM") - R("R_PGK") = 0
 M3phpc: - R("R_PSERT") + R("R_PGCD") = 0
 M3psmec: + R("R_PSCVT") - R("R_CHORS") = 0
 M4abutb: + R("R_EX_4abut_e_") + R("M4abutb_Transp") = 0
 M4abutc: + R("R_GLUDC") - R("R_GLUABUTt7") + R("R_ABUTt2")
 + R("R_ABUTD") - R("R_ABTA") = 0
 M4abute: + R("R_GLUABUTt7") - R("R_ABUTt2") - R("R_EX_4abut_e_") = 0
 M4abutnc: + R("R_PTRCTA") - R("R_ABUTD") = 0
 M4abzc: - R("R_DHPS2") + R("R_ADCL") = 0
 M4adchoc: + R("R_ADCS") - R("R_ADCL") = 0
 M4ahmmpc: - R("R_HMPK1") = 0
 M4ampmc: + R("R_AMPMS") + R("R_HMPK1") - R("R_PMPK") = 0
 M4c2mec: - R("R_CDPMEK") + R("R_MEPCT") = 0
 M4h2opntnc: + R("R_OP4ENH") - R("R_HOPNTAL") = 0
 M4hbac: + R("R_THZPSN") = 0
 M4hbzc: + R("R_CHRPL") - R("R_HBZOPT") = 0
 M4hthrc: + R("R_4HTHRS") = 0
 M4mhetzc: - R("R_HETZK") = 0
 M4mopc: - R("R_LEUTAi") + R("R_OMCDC") = 0
 M4mpetzc: + R("R_HETZK") + R("R_THZPSN") - R("R_TMPPP") = 0
 M4paspc: + R("R_ASPK") + R("R_ASAD") = 0
 M4perc: - R("R_PERD") + R("R_E4PD") = 0
 M4ppanc: - R("R_PPNCL2") + R("R_PNTK") = 0
 M4ppcysc: - R("R_PPCDC") + R("R_PPNCL2") = 0
 M4r5auc: + R("R_RBFSb") - R("R_RBFSa") + R("R_PMDPHT") = 0
 M5aizcc: - R("R_AIRC3") - R("R_PRASCS") = 0
 M5aopc: - 2 R("R_PPBNGS") + R("R_G1SATi") = 0
 M5aprbuc: + R("R_APRAUR") - R("R_PMDPHT") = 0
 M5apruc: - R("R_APRAUR") + R("R_DHPPDA2") = 0
 M5caizc: + R("R_AIRC2") + R("R_AIRC3") = 0
 M5dglcnc: - R("R_IDOND2") - R("R_IDOND") - R("R_5DGLCNR")
 + R("R_DKGLCNR2y") + R("R_DKGLCNR2x") = 0
 M5dh4dglcc: + R("R_GLCRD") + R("R_GALCTD") - R("R_GLCRAL") = 0
 M5mdr1pc: + R("R_MTRK") - R("R_MTRI") = 0
 M5mdru1pc: - R("R_MDRPD") + R("R_MTRI") = 0
 M5mtac: + R("R_SPMS") - R("R_MTAN") = 0
 M5mthfc: + R("R_MTHFR2") - R("R_METS") - 0.05 R("R_BiomassEcoli") = 0
 M5mtrc: + R("R_MTAN") - R("R_MTRK") = 0
 M5prdmbzc: - R("R_RZ5PP") + R("R_NNDMBRT") = 0
 M6hmhptc: - R("R_HPPK2") + R("R_DHNPA2") = 0
 M6hmhptppc: + R("R_HPPK2") - R("R_DHPS2") = 0
 M6pgcc: + R("R_PGL") + R("R_GNK") - R("R_GND") - R("R_EDD") = 0
 M6pglc: + R("R_G6PDH2r") - R("R_PGL") = 0
 M8aonnc: + R("R_AOXSr") - R("R_AMAOTr") = 0
 Maacaldc: + R("R_TAUDO") = 0
 Maacoac: + R("R_FAO4") + R("R_ACACT1r") + R("R_ACACCT") = 0
 Macacb: + R("R_EX_acac_e_") + R("Macacb_Transp") = 0
 Macacc: + R("R_ACACt2") - R("R_ACACCT") = 0
 Macace: - R("R_ACACt2") - R("R_EX_acac_e_") = 0
 MacACPc: + R("R_MACPD") - R("R_KAS14") + R("R_ACOATA") = 0
 Macaldb: + R("R_EX_acald_e_") + R("Macaldb_Transp") = 0
 Macaldc: + R("R_THRAr") + R("R_ETHAAL") - R("R_ALDD2x") + R("R_ACALDt")
 - R("R_ACALDi") + R("R_HOPNTAL") + R("R_DRPA") = 0
 Macalde: - R("R_ACALDt") - R("R_EX_acald_e_") = 0
 Macb: + R("R_EX_ac_e_") + R("Macb_Transp") = 0
 Macc: + R("R_UHGADA") + R("R_POX") + R("R_NACODA") + R("R_CYSS")
 + R("R_CITL") + R("R_BUTCT") + R("R_AGDC") + R("R_ACt2r") - R("R_ACS")
 + R("R_ACODA") - R("R_ACKr") + R("R_ALDD2x") + R("R_ACACCT") = 0
 Maccoac: - R("R_TDPADGAT") - R("R_SPMDAT2") - R("R_SPMDAT1")
 - R("R_SERAT") - R("R_PTAr") + R("R_PFL") + R("R_PDH") - R("R_MALS")
 - R("R_KAS15") - R("R_G1PACT") + 9 R("R_FAO3") + 8 R("R_FAO2")
 + 7 R("R_FAO1") - R("R_CS") - R("R_ADHEr") - R("R_ACGS")
 - R("R_ACCOACr") - R("R_BUTCT") + R("R_ACS") + R("R_ACALDi")
 - R("R_ACOATA") - 2 R("R_ACACT1r") - R("R_ACACCT") - R("R_IPPS")
 + R("R_GLYATi") - 5e-05 R("R_BiomassEcoli") = 0
 Mace: - R("R_ACt2r") - R("R_EX_ac_e_") = 0
 Macg5pc: + R("R_AGPR") + R("R_ACGK") = 0
 Macg5sac: + R("R_ACOTA") - R("R_AGPR") - R("R_NACODA") = 0
 Macgam1pc: - R("R_UAGDP") + R("R_G1PACT") = 0
 Macgam6pc: + R("R_AMANAPE") + R("R_ACGApts") - R("R_AGDC") = 0
 Macgamb: + R("R_EX_acgam_e_") + R("Macgamb_Transp") = 0
 Macgame: - R("R_EX_acgam_e_") - R("R_ACGApts") = 0
 Macgluc: - R("R_ACGK") + R("R_ACGS") = 0
 Macmanab: + R("R_EX_acmana_e_") + R("Macmanab_Transp") = 0
 Macmanac: - R("R_AMANK") + R("R_ACNML") = 0
 Macmanae: - R("R_ACMANApts") - R("R_EX_acmana_e_") = 0
 Macmanapc: + R("R_ACMANApts") + R("R_AMANK") - R("R_AMANAPE") = 0
 Macnamb: + R("R_EX_acnam_e_") + R("Macnamb_Transp") = 0
 Macnamc: + R("R_ACNAMt2") - R("R_ACNML") = 0
 Macname: - R("R_ACNAMt2") - R("R_EX_acnam_e_") = 0
 Maconmc: + R("R_ACONMT") = 0
 MaconTc: - R("R_ACONMT") = 0
 Macornc: - R("R_ACOTA") - R("R_ACODA") = 0
 MACPc: + 2 R("R_PASYN_EC") - R("R_MCOATA") + R("R_EDTXS4")
 + R("R_EDTXS3") + R("R_EDTXS2") + R("R_EDTXS1") + 7 R("R_C181SN")
 + 6 R("R_C161SN") + 6 R("R_C160SN") + 5 R("R_C141SN") + 5 R("R_C140SN")
 + 4 R("R_C120SN") + R("R_ACPS1") - R("R_AACPS5") - R("R_AACPS4")
 - R("R_AACPS3") - R("R_AACPS2") - R("R_AACPS1") + R("R_KAS14")
 - R("R_ACOATA") + R("R_UAGAAT") + R("R_U23GAAT") + R("R_KAS16") = 0
 Macserc: + R("R_SERAT") - R("R_CYSS") = 0
 MactACPc: - R("R_C181SN") - R("R_C161SN") - R("R_C160SN")
 - R("R_C141SN") - R("R_C140SN") - R("R_C120SN") + R("R_KAS15")
 + R("R_KAS14") = 0
 Mactpc: + R("R_PTAr") + R("R_ACKr") = 0
 Madeb: + R("R_EX_ade_e_") + R("Madeb_Transp") = 0
 Madec: + R("R_PUNP1") + R("R_AMPN") + R("R_AHCYSNS") - R("R_ADPT")
 + R("R_ADEt2r") - R("R_ADD") + R("R_MTAN") + R("R_PUNP2") = 0
 Madee: - R("R_ADEt2r") - R("R_EX_ade_e_") = 0
 Madnb: + R("R_EX_adn_e_") + R("Madnb_Transp") = 0
 Madnc: + R("R_NTD7") + R("R_AHC") + R("R_ADNt2r") + R("R_ADNt2")
 - R("R_ADNK1") - R("R_ADA") - R("R_PUNP1") = 0
 Madne: - R("R_ADNt2r") - R("R_ADNt2") - R("R_EX_adn_e_") = 0
 Madocbic: + R("R_CBIAT_DELETE") - R("R_ADOCBIK") = 0
 Madocbipc: - R("R_ACBIPGT") + R("R_ADOCBIK") = 0
 Madocblc: + R("R_CBLAT_DELETE") + R("R_ADOCBLS") = 0
 Madpc: + R("R_XYLK") + R("R_XYLabc") + R("R_VALabc") + R("R_URIDK2r")
 + R("R_UMPK") + R("R_UGMDDS") + R("R_UAMAS") + R("R_UAMAGS")
 + R("R_TSULabc") + R("R_TMPKr") + R("R_TMKr") + R("R_TMDK1")
 + R("R_THRabc") + R("R_THMabc") + R("R_TDSK") + R("R_TAURabc")
 + R("R_SULabc") + R("R_SUCOAS") + R("R_SUCCabc") + R("R_SPMDabc")
 + R("R_SHKK") - R("R_RNDR1") + R("R_RMK") + R("R_RIBabc")
 + R("R_RBK_L1") + R("R_RBK") + R("R_RBFK") - R("R_PYK") + R("R_PYDXNK")
 + R("R_PYDXK") + R("R_PYDAMK") + R("R_PTRCabc") + R("R_PROabc")
 + R("R_PRFGS") + R("R_PRAIS") + R("R_PRAGSr") + R("R_PPCK")
 - R("R_PPAKr") + R("R_PIabc") + R("R_PFK_2") + R("R_PFK")
 + R("R_ORNabc") + R("R_NDPK8") + R("R_NDPK7") + R("R_NDPK6")
 + R("R_NDPK5") + R("R_NDPK4") + R("R_NDPK3") + R("R_NDPK2")
 + R("R_NDPK1") + R("R_NAMNPP") + R("R_NADK") + R("R_METDabc")
 + R("R_METabc") + R("R_MALTTTRabc") + R("R_MALTTRabc")
 + R("R_MALTPTabc") + R("R_MALTHXabc") + R("R_MALTabc") + R("R_LYSabc")
 + 3 R("R_LPSSYN_EC") + R("R_LEUabc") + R("R_Kabc") + R("R_INSK")
 + R("R_ILEabc") + R("R_HSK") + R("R_HISabc") + R("R_HEX1")
 + R("R_GTHS") + R("R_GSPMDS") + R("R_GSNK") + R("R_GMHEPK")
 + R("R_GLYK") + R("R_GLYC3Pabc") + R("R_GLYBabc") + R("R_GLUCYS")
 + R("R_GLUabc") + R("R_GLU5K") + R("R_GLNS") + R("R_GLNabc")
 + R("R_GLCS1") + R("R_GK1") + R("R_GART") + R("R_GALKr")
 + R("R_GALabc") + R("R_FRUK") + R("R_FE2abc") + R("R_FCLK")
 + R("R_DXYLK") + R("R_DURIK1") + R("R_DTMPK") + R("R_DPCOAK")
 + R("R_DHFS") + R("R_DGK1") + R("R_DBTSr") + R("R_DADK") + R("R_CYTK2")
 + R("R_CYTK1") + R("R_CYSabc") + R("R_CTPS2") + R("R_CHLabc")
 + 2 R("R_CBPS") + R("R_CBMK") + R("R_CBL1abc") - R("R_ATPS4r")
 + R("R_ATPM") + R("R_ASPabc") + R("R_ASNabc") + R("R_ARGabc")
 + R("R_ARBabc") + R("R_AP5AH") + 2 R("R_AP4AH") + R("R_ALAALAr")
 + R("R_ALAabc") + R("R_ADSK") + R("R_ADK4") + R("R_ADK3")
 + 2 R("R_ADK1") + R("R_ACCOAL") + R("R_ADOCBIK") + R("R_ADNK1")
 + R("R_AMANK") + R("R_ACGK") + R("R_ACCOACr") + R("R_ACKr")
 + R("R_GNK") + R("R_MTRK") + R("R_AIRC2") + R("R_PNTK") + R("R_ASPK")
 + R("R_HETZK") + R("R_HMPK1") + R("R_GLYCK") + R("R_CDPMEK")
 + R("R_PMPK") + R("R_DDGALK") + R("R_DDGLK") + R("R_UAAGDS")
 + R("R_DAPabc") + R("R_PRASCS") + R("R_PGK") + R("R_DAGK_EC")
 + 45.5608 R("R_BiomassEcoli") = 0
 Madpglcc: + R("R_GLGC") - R("R_GLCS1") = 0
 MadphepDDc: + R("R_GMHEPAT") - R("R_AGMHE") = 0
 MadphepLDc: + R("R_AGMHE") - 3 R("R_LPSSYN_EC") = 0
 Magdpcbic: - R("R_ADOCBLS") + R("R_ACBIPGT") = 0
 Magmc: + R("R_ARGDC") - R("R_AGMT") = 0
 MagpcECc: + 0.02 R("R_PLIPA3") - 0.02 R("R_LPLIPA5")
 - 0.02 R("R_LPLIPA3") = 0
 MagpeECc: + 0.02 R("R_PLIPA2") - 0.02 R("R_LPLIPA4")
 - 0.02 R("R_LPLIPA2") = 0
 MagpgECc: + 0.02 R("R_PLIPA1") - 0.02 R("R_LPLIPA6")
 - 0.02 R("R_LPLIPA1") = 0
 Mahcysc: + 2 R("R_UPP3MT") - R("R_AHC") - R("R_AHCYSNS")
 + R("R_ACONMT") + R("R_DMQMT") + R("R_OMBZLM") + R("R_OHPHM")
 + R("R_AMMQT8_2") = 0
 Mahdtc: + R("R_GTPCI") - R("R_DNTPPA") = 0
 Maicarc: + R("R_IG3PS") + R("R_ADSL2r") - R("R_AICART") = 0
 Mairc: + R("R_PRAIS") - R("R_AIRC2") - R("R_AMPMS") = 0
 Makgb: + R("R_EX_akg_e_") + R("Makgb_Transp") = 0
 Makgc: + R("R_TDPAGTA") - R("R_SOTA") - R("R_SDPTA") - R("R_PHETA1")
 - R("R_OXGDC2") - R("R_ORNTA") + R("R_OHPBAT") + R("R_ICDHyr")
 + R("R_HSTPT") - R("R_GLUSy") + R("R_GLUDy") - R("R_ASPTA")
 - R("R_ALATA_L") + R("R_AKGt2r") - R("R_AKGDH") - R("R_ACOTA")
 - R("R_TAUDO") + R("R_LEUTAi") - R("R_PTRCTA") - R("R_ABTA")
 + R("R_PSERT") - R("R_ILETA") - R("R_VALTA") - R("R_TYRTA")
 + R("R_UNK3") = 0
 Makge: - R("R_AKGt2r") - R("R_EX_akg_e_") = 0
 Malaalac: - R("R_UGMDDS") + R("R_ALAALAr") = 0
 MalaBc: - R("R_PANTS") + R("R_ASP1DC") = 0
 MalacSc: + R("R_ACLS") - R("R_KARA1i") = 0
 MalaDb: + R("R_EX_ala_D_e_") + R("MalaDb_Transp") = 0
 MalaDc: + R("R_DALAt2r") - R("R_DAAD") - R("R_ALATA_D2") + R("R_ALARi")
 + R("R_ALAR") - 2 R("R_ALAALAr") = 0
 MalaDe: - R("R_DALAt2r") - R("R_EX_ala_D_e_") = 0
 MalaLb: + R("R_EX_ala_L_e_") + R("MalaLb_Transp") = 0
 MalaLc: + R("R_BTS2") - R("R_ALATA_L2") + R("R_ALAt2r") - R("R_ALARi")
 - R("R_ALAR") - R("R_ALATA_L") - R("R_UAMAS") + R("R_ALAabc")
 - R("R_AOXSr") + R("R_THZPSN") - R("R_VPAMT")
 - 0.488 R("R_BiomassEcoli") = 0
 MalaLe: - R("R_ALAt2r") - R("R_EX_ala_L_e_") - R("R_ALAabc") = 0
 Malltnb: + R("R_EX_alltn_e_") + R("Malltnb_Transp") = 0
 Malltnc: + R("R_ALLTNt2r") - R("R_ALLTN") = 0
 Malltne: - R("R_ALLTNt2r") - R("R_EX_alltn_e_") = 0
 Mallttc: - R("R_ALLTAH") + R("R_ALLTN") = 0
 Maltrnc: - R("R_TAGURr") - R("R_ALTRH") = 0
 Mametamc: + R("R_ADMDCr") - R("R_SPMS") = 0
 Mametc: + R("R_METAT") - R("R_ADMDCr") - 2 R("R_UPP3MT")
 - R("R_ACONMT") - R("R_AMAOTr") - R("R_DMQMT") - R("R_OMBZLM")
 - R("R_OHPHM") - R("R_AMMQT8_2") = 0
 Mamobc: + R("R_AMAOTr") = 0
 Mampb: + R("R_EX_amp_e_") + R("Mampb_Transp") = 0
 Mampc: + R("R_SUCBZL") + R("R_SELNPS") + R("R_PRPPS") + R("R_PPS")
 + R("R_PACCOAL") + R("R_NTPP6") + R("R_NADS1") + R("R_NADDP")
 + R("R_GMPS2") + R("R_GLUTRS") + R("R_BPNT") + R("R_ASNS2")
 + R("R_ASNS1") + R("R_ARGSS") + R("R_ADSL1r") + R("R_PANTS")
 - R("R_ADK4") - R("R_ADK3") - R("R_ADK1") - R("R_NTD7") + R("R_ADNK1")
 - R("R_AMPN") + R("R_ADPT") + R("R_AACPS5") + R("R_AACPS4")
 + R("R_AACPS3") + R("R_AACPS2") + R("R_AACPS1") + R("R_FAO3")
 + R("R_FAO2") + R("R_FAO1") + R("R_ACS") + R("R_HPPK2") + R("R_THZPSN")
 + 6 R("R_ENTCS") - 0.001 R("R_BiomassEcoli") = 0
 Mampe: + R("R_NADDPe") - R("R_EX_amp_e_") = 0
 Manthc: + R("R_ANS") - R("R_ANPRT") = 0
 Map4ac: - R("R_AP4AH") = 0
 Map5ac: - R("R_AP5AH") = 0
 MapgECc: + 0.02 R("R_LPLIPA6") + 0.02 R("R_LPLIPA4")
 + 0.02 R("R_LPLIPA5") = 0
 MapoACPc: - R("R_ACPS1") = 0
 Mapsc: + R("R_SADT2") - R("R_ADSK") = 0
 Mara5pc: - R("R_KDOPS") + R("R_A5PISO") = 0
 MarabLb: + R("R_EX_arab_L_e_") + R("MarabLb_Transp") = 0
 MarabLc: + R("R_ARBt2r") - R("R_ARAI") + R("R_ARBabc") = 0
 MarabLe: - R("R_ARBt2r") - R("R_EX_arab_L_e_") - R("R_ARBabc") = 0
 Marbt6pc: - R("R_AB6PGH") = 0
 MargLb: + R("R_EX_arg_L_e_") + R("MargLb_Transp") = 0
 MargLc: - R("R_AST") + R("R_ARGSL") + R("R_ARGORNt7") - R("R_ARGDC")
 + R("R_ARGabc") - 0.281 R("R_BiomassEcoli") = 0
 MargLe: - R("R_ARGORNt7") - R("R_EX_arg_L_e_") - R("R_ARGabc") = 0
 Margsucc: - R("R_ARGSL") + R("R_ARGSS") = 0
 MasnLb: + R("R_EX_asn_L_e_") + R("MasnLb_Transp") = 0
 MasnLc: + R("R_ASNt2r") - R("R_ASNN") + R("R_ASNS2") + R("R_ASNS1")
 + R("R_ASNabc") - 0.229 R("R_BiomassEcoli") = 0
 MasnLe: - R("R_ASNt2r") - R("R_EX_asn_L_e_") - R("R_ASNabc") = 0
 MaspLb: + R("R_EX_asp_L_e_") + R("MaspLb_Transp") = 0
 MaspLc: + R("R_ASPt2_3") + R("R_ASPt2_2") + R("R_ASPt2") - R("R_ASPT")
 - R("R_ASPO6") - R("R_ASPO5") - R("R_ASPO4") - R("R_ASPO3")
 - R("R_ASPCT") - R("R_ADSS") + R("R_ASNN") - R("R_ASNS2")
 - R("R_ASNS1") - R("R_ARGSS") - R("R_ASP1DC") - R("R_ASPTA")
 + R("R_ASPabc") - R("R_ASPK") - R("R_PRASCS")
 - 0.229 R("R_BiomassEcoli") = 0
 MaspLe: - R("R_ASPt2_3") - R("R_ASPt2_2") - R("R_ASPt2")
 - R("R_EX_asp_L_e_") - R("R_ASPabc") = 0
 Maspsac: + R("R_HSDy") - R("R_ASAD") - R("R_DHDPS") = 0
 Matpc: - R("R_SERASr") - R("R_RNTR1") - R("R_PTPATi") - R("R_NNAT")
 - R("R_NMNAT") - R("R_FMNAT") - R("R_ATPPRT") - R("R_ADNCYC")
 - R("R_SADT2") - R("R_SUCBZL") - R("R_SELNPS") - R("R_PRPPS")
 - R("R_PPS") - R("R_PACCOAL") - R("R_NTPP6") - R("R_NADS1")
 - R("R_GMPS2") - R("R_GLUTRS") - R("R_ASNS2") - R("R_ASNS1")
 - R("R_ARGSS") - R("R_METAT") - R("R_PANTS") - R("R_GMHEPAT")
 - R("R_GLGC") - R("R_XYLK") - R("R_XYLabc") - R("R_VALabc")
 - R("R_URIDK2r") - R("R_UMPK") - R("R_UGMDDS") - R("R_UAMAS")
 - R("R_UAMAGS") - R("R_TSULabc") - R("R_TMPKr") - R("R_TMKr")
 - R("R_TMDK1") - R("R_THRabc") - R("R_THMabc") - R("R_TDSK")
 - R("R_TAURabc") - R("R_SULabc") - R("R_SUCOAS") - R("R_SUCCabc")
 - R("R_SPMDabc") - R("R_SHKK") - R("R_RMK") - R("R_RIBabc")
 - R("R_RBK_L1") - R("R_RBK") - R("R_RBFK") + R("R_PYK") - R("R_PYDXNK")
 - R("R_PYDXK") - R("R_PYDAMK") - R("R_PTRCabc") - R("R_PROabc")
 - R("R_PRFGS") - R("R_PRAIS") - R("R_PRAGSr") - R("R_PPCK")
 + R("R_PPAKr") - R("R_PIabc") - R("R_PFK_2") - R("R_PFK")
 - R("R_ORNabc") - R("R_NDPK8") - R("R_NDPK7") - R("R_NDPK6")
 - R("R_NDPK5") - R("R_NDPK4") - R("R_NDPK3") - R("R_NDPK2")
 - R("R_NDPK1") - R("R_NAMNPP") - R("R_NADK") - R("R_METDabc")
 - R("R_METabc") - R("R_MALTTTRabc") - R("R_MALTTRabc")
 - R("R_MALTPTabc") - R("R_MALTHXabc") - R("R_MALTabc") - R("R_LYSabc")
 - R("R_LEUabc") - R("R_Kabc") - R("R_INSK") - R("R_ILEabc")
 - R("R_HSK") - R("R_HISabc") - R("R_HEX1") - R("R_GTHS")
 - R("R_GSPMDS") - R("R_GSNK") - R("R_GMHEPK") - R("R_GLYK")
 - R("R_GLYC3Pabc") - R("R_GLYBabc") - R("R_GLUCYS") - R("R_GLUabc")
 - R("R_GLU5K") - R("R_GLNS") - R("R_GLNabc") - R("R_GK1") - R("R_GART")
 - R("R_GALKr") - R("R_GALabc") - R("R_FRUK") - R("R_FE2abc")
 - R("R_FCLK") - R("R_DXYLK") - R("R_DURIK1") - R("R_DTMPK")
 - R("R_DPCOAK") - R("R_DHFS") - R("R_DGK1") - R("R_DBTSr")
 - R("R_DADK") - R("R_CYTK2") - R("R_CYTK1") - R("R_CYSabc")
 - R("R_CTPS2") - R("R_CHLabc") - 2 R("R_CBPS") - R("R_CBMK")
 - R("R_CBL1abc") + R("R_ATPS4r") - R("R_ATPM") - R("R_ASPabc")
 - R("R_ASNabc") - R("R_ARGabc") - R("R_ARBabc") + R("R_AP5AH")
 - R("R_ALAALAr") - R("R_ALAabc") - R("R_ADSK") - R("R_ADK1")
 - R("R_ACCOAL") - R("R_CBLAT_DELETE") - R("R_CBIAT_DELETE")
 - R("R_ADOCBIK") - R("R_ADNK1") - R("R_AACPS5") - R("R_AACPS4")
 - R("R_AACPS3") - R("R_AACPS2") - R("R_AACPS1") - R("R_AMANK")
 - R("R_ACGK") - R("R_FAO3") - R("R_FAO2") - R("R_FAO1")
 - R("R_ACCOACr") - R("R_ACS") - R("R_ACKr") - R("R_GNK") - R("R_HPPK2")
 - R("R_MTRK") - R("R_AIRC2") - R("R_PNTK") - R("R_ASPK") - R("R_HETZK")
 - R("R_THZPSN") - R("R_HMPK1") - R("R_GLYCK") - R("R_CDPMEK")
 - R("R_PMPK") - R("R_DDGALK") - R("R_DDGLK") - R("R_UAAGDS")
 - R("R_DAPabc") - R("R_PRASCS") - R("R_DHBSr") - R("R_PGK")
 - R("R_DAGK_EC") - 45.7318 R("R_BiomassEcoli") = 0
 Mbbtcoac: - R("R_CRNBTCT") = 0
 Mbetaldc: - R("R_BETALDHy") - R("R_BETALDHx") = 0
 Mbtcoac: + R("R_BUTCT") - R("R_FAO4") = 0
 Mbtnc: + R("R_BSORy") + R("R_BSORx") + R("R_BTS2") = 0
 Mbtnsoc: - R("R_BSORy") - R("R_BSORx") = 0
 Mbutb: + R("R_EX_but_e_") + R("Mbutb_Transp") = 0
 Mbutc: + R("R_BUTt2r") - R("R_BUTCT") = 0
 Mbute: - R("R_BUTt2r") - R("R_EX_but_e_") = 0
 Mcampc: + R("R_ADNCYC") = 0
 Mcbaspc: + R("R_DHORTS") + R("R_ASPCT") = 0
 Mcbic: - R("R_CBIAT_DELETE") = 0
 Mcbl1b: + R("R_EX_cbl1_e_") + R("Mcbl1b_Transp") = 0
 Mcbl1c: + R("R_CBL1abc") - R("R_CBLAT_DELETE") = 0
 Mcbl1e: - R("R_EX_cbl1_e_") - R("R_CBL1abc") = 0
 Mcbpc: - R("R_OCBT") - R("R_ASPCT") + R("R_CBPS") + R("R_CBMK") = 0
 Mcdpc: - R("R_RNDR3") - R("R_NDPK3") + 2 R("R_LPSSYN_EC")
 + R("R_CYTK1") = 0
 McdpdagECc: - 0.02 R("R_PSSA_EC") - 0.02 R("R_PGSA_EC")
 + 0.02 R("R_DASYN_EC") - 0.02 R("R_CDAPPA_EC") = 0
 Mcdpeac: - 2 R("R_LPSSYN_EC") + R("R_PEPT_EC") = 0
 Mcechdddc: + R("R_PPPNDO") - R("R_DHPPD") = 0
 Mcenchdddc: - R("R_DHCIND") + R("R_CINNDO") = 0
 Mcholb: + R("R_EX_chol_e_") + R("Mcholb_Transp") = 0
 Mcholc: + R("R_GPDDA1") + R("R_CHLt2r") + R("R_CHLabc") = 0
 Mchole: - R("R_CHLt2r") - R("R_EX_chol_e_") - R("R_CHLabc") = 0
 Mchorc: - R("R_ICHORSi") - R("R_CHORM") - R("R_ANS") - R("R_CHRPL")
 - R("R_ADCS") + R("R_CHORS") = 0
 Mcinnmc: - R("R_CINNDO") = 0
 Mcitb: + R("R_EX_cit_e_") + R("Mcitb_Transp") = 0
 Mcitc: + R("R_CITt7") - R("R_ACONT") + R("R_CS") - R("R_CITL") = 0
 Mcite: - R("R_CITt7") - R("R_EX_cit_e_") = 0
 McitrLc: + R("R_OCBT") - R("R_ARGSS") = 0
 Mckdoc: - R("R_MOAT2") - R("R_MOAT") + R("R_KDOCT2")
 - 3 R("R_LPSSYN_EC") = 0
 MclpnECc: + 0.02 R("R_CLPNS_EC") - 0.000129 R("R_BiomassEcoli") = 0
 Mcmpc: + R("R_NTPP4") - R("R_NTD4") + R("R_CYTDK2") - R("R_CMPN")
 + R("R_MOAT2") + R("R_MOAT") + R("R_PSSA_EC") + R("R_PGSA_EC")
 + R("R_CDAPPA_EC") + 3 R("R_LPSSYN_EC") - R("R_CYTK1") + R("R_PPNCL2")
 + R("R_MECDPS") - R("R_PEPT_EC") = 0
 Mco2b: + R("R_EX_co2_e_") + R("Mco2b_Transp") = 0
 Mco2c: + 4 R("R_UPPDC1") + R("R_UGLYCH") + R("R_SADH") + R("R_PSD_EC")
 + R("R_PPNDH") - R("R_PPC") + R("R_PDX5PS") + R("R_ORNDC")
 + R("R_OMPDC") + R("R_NNDPR") + R("R_MMCD") + R("R_ME2") + R("R_ME1")
 - R("R_HCO3E") + R("R_GLYCL") + R("R_FHL") + R("R_FDH3") + R("R_FDH2")
 + R("R_DXPS") + 2 R("R_CYNTAH") + 2 R("R_CPPPGO") + R("R_CO2t")
 + R("R_ADMDCr") + R("R_ACLS") + R("R_ASP1DC") + R("R_OXGDC2")
 + R("R_ICDHyr") + R("R_AKGDH") + R("R_ARGDC") + R("R_PPCK")
 - R("R_DBTSr") - R("R_CBMK") + 7 R("R_C181SN") + 6 R("R_C161SN")
 + 6 R("R_C160SN") + 5 R("R_C141SN") + 5 R("R_C140SN") + 4 R("R_C120SN")
 + R("R_PDH") + R("R_KAS15") + R("R_POX") + R("R_MACPD") + R("R_KAS14")
 + R("R_TAUDO") + R("R_AOXSr") + R("R_GND") + R("R_PPCDC")
 + R("R_THZPSN") + R("R_GLUDC") + R("R_KAS16") + R("R_KG6PDC")
 + R("R_OMCDC") + R("R_PPND") + R("R_OPHBDC") + R("R_GLXCL")
 + R("R_DHNAOT") + R("R_IGPS") + R("R_ACHBS") + R("R_DAPDC")
 + R("R_LYSDC") = 0
 Mco2e: - R("R_CO2t") - R("R_EX_co2_e_") = 0
 Mcoac: + R("R_THDPS") + R("R_PTA2") + R("R_NPHS") + R("R_HSST")
 + R("R_AST") - R("R_SUCBZL") - R("R_PACCOAL") - R("R_AKGDH")
 - R("R_SUCOAS") + R("R_DPCOAK") - R("R_ACCOAL") + R("R_MCOATA")
 - R("R_ACPS1") + R("R_TDPADGAT") + R("R_SPMDAT2") + R("R_SPMDAT1")
 + R("R_SERAT") + R("R_PTAr") - R("R_PFL") - R("R_PDH") + R("R_MALS")
 + R("R_KAS15") + R("R_G1PACT") - 9 R("R_FAO3") - 8 R("R_FAO2")
 - 7 R("R_FAO1") + R("R_CS") + R("R_ADHEr") + R("R_ACGS") - R("R_ACS")
 - R("R_ACALDi") + R("R_ACOATA") + R("R_ACACT1r") + R("R_AOXSr")
 + R("R_IPPS") - R("R_OBTFL") + R("R_MCITS") - R("R_GLYATi")
 - 6e-06 R("R_BiomassEcoli") = 0
 Mcpppg3c: + R("R_UPPDC1") - R("R_CPPPGO") = 0
 Mcrnb: + R("R_EX_crn_e_") + R("Mcrnb_Transp") = 0
 Mcrnc: + R("R_CRNt7") - R("R_CRNCBCT") - R("R_CRNBTCT") = 0
 Mcrncoac: - R("R_CRNCDH") + R("R_CRNCBCT") + R("R_CRNBTCT") = 0
 Mcrne: - R("R_CRNt7") - R("R_EX_crn_e_") = 0
 Mcsnb: + R("R_EX_csn_e_") + R("Mcsnb_Transp") = 0
 Mcsnc: + R("R_CSNt2") - R("R_CSND") + R("R_CMPN") = 0
 Mcsne: - R("R_CSNt2") - R("R_EX_csn_e_") = 0
 Mctbtc: + R("R_CRNCBCT") = 0
 Mctbtcoac: + R("R_CRNCDH") - R("R_CRNCBCT") = 0
 Mctpc: - R("R_RNTR3") - R("R_NTPP4") - R("R_KDOCT2") - R("R_DASYN_EC")
 + R("R_NDPK3") + R("R_CTPS2") - R("R_PPNCL2") - R("R_MEPCT")
 - 0.126 R("R_BiomassEcoli") = 0
 Mcyanc: - R("R_CYANST") = 0
 Mcyntb: + R("R_EX_cynt_e_") + R("Mcyntb_Transp") = 0
 Mcyntc: + R("R_CYNTt2") - R("R_CYNTAH") = 0
 Mcynte: - R("R_CYNTt2") - R("R_EX_cynt_e_") = 0
 McysLb: + R("R_EX_cys_L_e_") + R("McysLb_Transp") = 0
 McysLc: - R("R_SHSL1") - R("R_CYSDS") - R("R_BTS2") - R("R_GLUCYS")
 + R("R_CYSabc") + R("R_CYSS") - R("R_PPNCL2") - R("R_THZPSN")
 - 0.087 R("R_BiomassEcoli") = 0
 McysLe: - R("R_EX_cys_L_e_") - R("R_CYSabc") = 0
 McystLc: - R("R_CYSTL") + R("R_SHSL1") = 0
 Mcytdb: + R("R_EX_cytd_e_") + R("Mcytdb_Transp") = 0
 Mcytdc: + R("R_CYTDt2r") + R("R_CYTDt2") - R("R_CYTD") + R("R_NTD4")
 - R("R_CYTDK2") = 0
 Mcytde: - R("R_CYTDt2r") - R("R_CYTDt2") - R("R_EX_cytd_e_") = 0
 Mdad2b: + R("R_EX_dad_2_e_") + R("Mdad2b_Transp") = 0
 Mdad2c: + R("R_NTD6") + R("R_DADNt2") - R("R_DADA") - R("R_PUNP2") = 0
 Mdad2e: - R("R_DADNt2") - R("R_EX_dad_2_e_") = 0
 Mdadpc: + R("R_RNDR1") - R("R_NDPK8") + R("R_DADK") = 0
 Mdampc: + R("R_NTPP5") - R("R_NTD6") - R("R_DADK") = 0
 Mdannc: - R("R_DBTSr") + R("R_AMAOTr") = 0
 Mdatpc: - R("R_NTPP5") + R("R_RNTR1") + R("R_NDPK8")
 - 0.0247 R("R_BiomassEcoli") = 0
 Mdb4pc: + R("R_DB4PS") - R("R_RBFSa") = 0
 Mdcampc: + R("R_ADSS") - R("R_ADSL1r") = 0
 Mdcdpc: + R("R_RNDR3") - R("R_NDPK7") + R("R_CYTK2") = 0
 Mdcmpc: + R("R_NTPP3") - R("R_NTD3") - R("R_CYTK2") = 0
 Mdctpc: - R("R_DCTPD") - R("R_NTPP3") + R("R_RNTR3") + R("R_NDPK7")
 - 0.0254 R("R_BiomassEcoli") = 0
 Mdcytb: + R("R_EX_dcyt_e_") + R("Mdcytb_Transp") = 0
 Mdcytc: + R("R_DCYTt2") - R("R_DCYTD") + R("R_NTD3") = 0
 Mdcyte: - R("R_DCYTt2") - R("R_EX_dcyt_e_") = 0
 MddcaACPc: - R("R_EDTXS1") + R("R_C120SN") - R("R_KAS16") = 0
 Mdgdpc: + R("R_RNDR2") - R("R_NDPK5") + R("R_DGK1") = 0
 Mdgmpc: + R("R_NTPP1") - R("R_NTD8") - R("R_DGK1") = 0
 Mdgsnb: + R("R_EX_dgsn_e_") + R("Mdgsnb_Transp") = 0
 Mdgsnc: + R("R_NTPTP1") + R("R_DGSNt2") + R("R_NTD8") - R("R_PUNP4")
 = 0
 Mdgsne: - R("R_DGSNt2") - R("R_EX_dgsn_e_") = 0
 Mdgtpc: + R("R_RNTR2") - R("R_NTPTP1") - R("R_NTPP1") + R("R_NDPK5")
 - 0.0254 R("R_BiomassEcoli") = 0
 Mdhab: + R("R_EX_dha_e_") + R("Mdhab_Transp") = 0
 Mdhac: + R("R_GLYCDx") + R("R_F6PA") + R("R_DHAt") - R("R_DHAPT") = 0
 Mdhae: - R("R_DHAt") - R("R_EX_dha_e_") = 0
 Mdhapc: - R("R_TPI") + R("R_TGBPA") + R("R_RMPA") - R("R_QULNS")
 - R("R_MGSA") + R("R_G3PD6") + R("R_G3PD5") + R("R_G3PD2")
 + R("R_FCLPA") + R("R_FBA") + R("R_DHAPT") + R("R_G3PD7") = 0
 Mdhcinnmc: - R("R_DHCINDO") + R("R_DHCIND") + R("R_3HCINNMH") = 0
 Mdhfc: + R("R_TMDS") - R("R_DHFR") + R("R_DHFS") = 0
 Mdhnac: + R("R_NPHS") - R("R_DHNAOT") = 0
 Mdhnptc: + R("R_DNMPPA") - R("R_DHNPA2") = 0
 MdhorSc: - R("R_DHORD5") - R("R_DHORD2") - R("R_DHORTS") = 0
 Mdhpmpc: - R("R_DNMPPA") + R("R_DNTPPA") = 0
 Mdhpppnc: - R("R_HPPPNDO") + R("R_DHPPD") + R("R_3HPPPNH") = 0
 Mdhptc: - R("R_DHFS") + R("R_DHPS2") = 0
 Mdhptdc: + R("R_RHCCE") - R("R_DHPTDCs") = 0
 Mdinb: + R("R_EX_din_e_") + R("Mdinb_Transp") = 0
 Mdinc: + R("R_DINSt2") + R("R_DADA") - R("R_PUNP6") = 0
 Mdine: - R("R_DINSt2") - R("R_EX_din_e_") = 0
 Mdkmppc: + R("R_MDRPD") - R("R_DKMPPD2") - R("R_DKMPPD") = 0
 Mdmbzidc: - R("R_NNDMBRT") = 0
 Mdmlzc: - 2 R("R_RBFSb") + R("R_RBFSa") = 0
 Mdmppc: + R("R_IPDDIi") + R("R_DMPPS") - R("R_DMATT") = 0
 Mdmsb: + R("R_EX_dms_e_") + R("Mdmsb_Transp") = 0
 Mdmsc: + R("R_DMSOR1") + R("R_DMSOR2") = 0
 Mdmse: + R("R_DMSOR1e") - R("R_EX_dms_e_") + R("R_DMSOR2e") = 0
 Mdmsob: + R("R_EX_dmso_e_") + R("Mdmsob_Transp") = 0
 Mdmsoc: - R("R_DMSOR1") - R("R_DMSOR2") = 0
 Mdmsoe: - R("R_EX_dmso_e_") - R("R_DMSOR1e") - R("R_DMSOR2e") = 0
 Mdnadc: + R("R_NNAT") - R("R_NADS1") = 0
 Mdpcoac: + R("R_PTPATi") - R("R_DPCOAK") = 0
 Mdsclc: - R("R_SHCHD2") + R("R_UPP3MT") = 0
 Mdtbtc: - R("R_BTS2") + R("R_DBTSr") = 0
 Mdtdp4aaddgc: - R("R_AADDGT") + R("R_TDPADGAT") = 0
 Mdtdp4addgc: + R("R_TDPAGTA") - R("R_TDPADGAT") = 0
 Mdtdp4d6dgc: + R("R_TDPGDH") - R("R_TDPDRE") - R("R_TDPAGTA") = 0
 Mdtdp4d6dmc: - R("R_TDPDRR") + R("R_TDPDRE") = 0
 Mdtdpc: + R("R_AADDGT") - R("R_NDPK4") + R("R_DTMPK") = 0
 Mdtdpgluc: + R("R_G1PTT") - R("R_TDPGDH") = 0
 Mdtdprmnc: + R("R_TDPDRR") = 0
 Mdtmpc: + R("R_NTPP7") - R("R_NTD5") + R("R_TMDS") + R("R_TMDK1")
 - R("R_DTMPK") = 0
 Mdttpc: - R("R_NTPP7") - R("R_G1PTT") + R("R_NDPK4")
 - 0.0247 R("R_BiomassEcoli") = 0
 Mdudpc: + R("R_RNDR4") + R("R_URIDK2r") - R("R_NDPK6") = 0
 Mdumpc: - R("R_NTD1") + R("R_DUTPDP") - R("R_TMDS") - R("R_URIDK2r")
 + R("R_DURIK1") = 0
 Mdurib: + R("R_EX_duri_e_") + R("Mdurib_Transp") = 0
 Mduric: + R("R_DURIt2") + R("R_NTD1") + R("R_DCYTD") - R("R_DURIK1")
 - R("R_DURIPP") = 0
 Mdurie: - R("R_DURIt2") - R("R_EX_duri_e_") = 0
 Mdutpc: + R("R_RNTR4") - R("R_DUTPDP") + R("R_DCTPD") + R("R_NDPK6")
 = 0
 Mdxyl5pc: - R("R_PDX5PS") + R("R_DXPS") + R("R_DXYLK") - R("R_THZPSN")
 - R("R_DXPRIi") = 0
 Mdxylc: - R("R_DXYLK") = 0
 Me4pc: - R("R_TKT2") + R("R_TALA") - R("R_E4PD") - R("R_DDPA") = 0
 MecaECc: + R("R_ECAP_EC") = 0
 Meig3pc: - R("R_IGPDH") + R("R_IG3PS") = 0
 Menterc: + R("R_ENTCS") = 0
 Methac: + R("R_GPDDA2") - R("R_ETHAAL") = 0
 Metohb: + R("R_EX_etoh_e_") + R("Metohb_Transp") = 0
 Metohc: + R("R_ETOHt2r") + R("R_ADHEr") = 0
 Metohe: - R("R_ETOHt2r") - R("R_EX_etoh_e_") = 0
 Mf1pc: + R("R_FRUpts") - R("R_FRUK") = 0
 Mf6pc: + R("R_SBTPD") + R("R_PGI") + R("R_MAN6PI") + R("R_M1PD")
 - R("R_GF6PTA") + R("R_G6PDA") + R("R_FRUpts2") + R("R_FBP")
 + R("R_TKT2") + R("R_TALA") - R("R_F6PA") - R("R_PFK") = 0
 Mfadc: + R("R_SUCD4") - R("R_SUCD1i") + R("R_FMNAT") - R("R_DAAD")
 - 8 R("R_FAO3") - 7 R("R_FAO2") - 6 R("R_FAO1") - R("R_FAO4")
 - R("R_PROD2") - 1e-05 R("R_BiomassEcoli") = 0
 Mfadh2c: - R("R_SUCD4") + R("R_SUCD1i") + R("R_DAAD") + 8 R("R_FAO3")
 + 7 R("R_FAO2") + 6 R("R_FAO1") + R("R_FAO4") + R("R_PROD2") = 0
 Mfc1pc: - R("R_FCLPA") + R("R_FCLK") = 0
 MfclLc: + R("R_FCI") - R("R_FCLK") = 0
 Mfdpc: - R("R_FBP") - R("R_FBA") + R("R_PFK") + R("R_FRUK") = 0
 Mfe2b: + R("R_EX_fe2_e_") + R("Mfe2b_Transp") = 0
 Mfe2c: - R("R_SHCHF") - R("R_FCLT") + R("R_FE2abc") = 0
 Mfe2e: - R("R_EX_fe2_e_") - R("R_FE2abc") = 0
 Mfgamc: - R("R_PRFGS") + R("R_GART") + R("R_GARFT") = 0
 Mfmnc: - R("R_FMNAT") + R("R_RBFK") = 0
 Mforb: + R("R_EX_for_e_") + R("Mforb_Transp") = 0
 Mforc: + R("R_FORt") + R("R_DB4PS") - R("R_FHL") - R("R_FDH3")
 - R("R_FDH2") + R("R_GTPCI") - R("R_GART") + R("R_PFL")
 + 2 R("R_AMPMS") + R("R_OBTFL") + R("R_DKMPPD2") + R("R_DKMPPD")
 + R("R_GTPCII2") + R("R_FTHFD") = 0
 Mfore: - R("R_FORt") - R("R_EX_for_e_") = 0
 Mfpramc: + R("R_PRFGS") - R("R_PRAIS") = 0
 Mfpricac: + R("R_IMPC") + R("R_AICART") = 0
 Mfrdpc: - R("R_UDCPDPS") - R("R_OCTDPS") - R("R_HEMEOS") + R("R_GRTT")
 = 0
 Mfrub: + R("R_EX_fru_e_") + R("Mfrub_Transp") = 0
 Mfruc: - R("R_XYLI2i") + R("R_FFSD") = 0
 Mfrue: - R("R_EX_fru_e_") - R("R_FRUpts2") - R("R_FRUpts") = 0
 Mfruurc: + R("R_MANAO") + R("R_GUI1") = 0
 Mfuc1pLb: + R("R_EX_fuc1p_L_e_") + R("Mfuc1pLb_Transp") = 0
 Mfuc1pLc: + R("R_FUCPt6_2") = 0
 Mfuc1pLe: - R("R_FUCPt6_2") - R("R_EX_fuc1p_L_e_") = 0
 MfucLb: + R("R_EX_fuc_L_e_") + R("MfucLb_Transp") = 0
 MfucLc: + R("R_FUCt") - R("R_FCI") = 0
 MfucLe: - R("R_FUCt") - R("R_EX_fuc_L_e_") = 0
 Mfumb: + R("R_EX_fum_e_") + R("Mfumb_Transp") = 0
 Mfumc: + R("R_SUCFUMt") + R("R_HKNTDH") + R("R_FUMt2_3")
 + R("R_FUMt2_2") - R("R_FUM") - R("R_FRD2") + R("R_SUCD1i")
 + R("R_ASPT") - R("R_ASPO5") + R("R_ARGSL") + R("R_ADSL1r")
 - R("R_FRD3") + R("R_ADSL2r") = 0
 Mfume: - R("R_SUCFUMt") - R("R_FUMt2_3") - R("R_FUMt2_2")
 - R("R_EX_fum_e_") = 0
 Mg1pc: + R("R_UGLT") - R("R_PGMT") + R("R_MLTP3") + R("R_MLTP2")
 + R("R_MLTP1") + R("R_GLCP") - R("R_GALUi") - R("R_GALU") - R("R_G1PP")
 - R("R_G1PTT") - R("R_GLGC") = 0
 Mg3pc: + R("R_TKT1") + R("R_TKT2") - R("R_TALA") + R("R_TPI")
 + R("R_TGBPA") + R("R_FBA") + R("R_F6PA") - R("R_DXPS") + R("R_TRPS3")
 + R("R_TRPS1") + R("R_DRPA") + R("R_DDPGALA") + R("R_EDA")
 - R("R_GAPD") = 0
 Mg3pcc: - R("R_GPDDA1") + R("R_LPLIPA5") + R("R_LPLIPA3") = 0
 Mg3pec: - R("R_GPDDA2") + R("R_LPLIPA4") + R("R_LPLIPA2") = 0
 Mg3pgc: - R("R_GPDDA4") + R("R_LPLIPA6") + R("R_LPLIPA1") = 0
 Mg3pic: - R("R_GPDDA5") = 0
 Mg3psc: - R("R_GPDDA3") = 0
 Mg6pb: + R("R_EX_g6p_e_") + R("Mg6pb_Transp") = 0
 Mg6pc: - R("R_TRE6PS") + R("R_TRE6PH") + R("R_GLCpts") + R("R_G6Pt6_2")
 + R("R_PGMT") + R("R_FFSD") - R("R_PGI") + R("R_AB6PGH") + R("R_HEX1")
 - R("R_G6PDH2r") = 0
 Mg6pe: - R("R_G6Pt6_2") - R("R_EX_g6p_e_") = 0
 Mgal1pc: - R("R_UGLT") + R("R_GALKr") = 0
 Mgalb: + R("R_EX_gal_e_") + R("Mgalb_Transp") = 0
 Mgalc: + R("R_LACZ") + R("R_GALt2") + R("R_GALS3") - R("R_GALKr")
 + R("R_GALabc") = 0
 MgalctDb: + R("R_EX_galct_D_e_") + R("MgalctDb_Transp") = 0
 MgalctDc: + R("R_GALCTt2r") - R("R_GALCTD") = 0
 MgalctDe: - R("R_GALCTt2r") - R("R_EX_galct_D_e_") = 0
 MgalctnDb: + R("R_EX_galctn_D_e_") + R("MgalctnDb_Transp") = 0
 MgalctnDc: + R("R_GALCTNt2r") - R("R_GALCTND") = 0
 MgalctnDe: - R("R_GALCTNt2r") - R("R_EX_galctn_D_e_") = 0
 Mgale: - R("R_GALt2") - R("R_EX_gal_e_") - R("R_GALabc") = 0
 Mgalt1pc: - R("R_GLTPD") + R("R_GALTpts") = 0
 Mgaltb: + R("R_EX_galt_e_") + R("Mgaltb_Transp") = 0
 Mgalte: - R("R_EX_galt_e_") - R("R_GALTpts") = 0
 Mgalurb: + R("R_EX_galur_e_") + R("Mgalurb_Transp") = 0
 Mgalurc: - R("R_GUI2") + R("R_GALURt2r") = 0
 Mgalure: - R("R_GALURt2r") - R("R_EX_galur_e_") = 0
 Mgam1pc: - R("R_PGAMT") - R("R_G1PACT") = 0
 Mgam6pc: + R("R_GAMpts") + R("R_PGAMT") + R("R_GF6PTA") - R("R_G6PDA")
 + R("R_AGDC") = 0
 Mgamb: + R("R_EX_gam_e_") + R("Mgamb_Transp") = 0
 Mgame: - R("R_EX_gam_e_") - R("R_GAMpts") = 0
 Mgarc: + R("R_PRAGSr") - R("R_GART") - R("R_GARFT") = 0
 Mgbbtnb: + R("R_EX_gbbtn_e_") + R("Mgbbtnb_Transp") = 0
 Mgbbtnc: - R("R_CRNt7") + R("R_CRNBTCT") = 0
 Mgbbtne: - R("R_EX_gbbtn_e_") + R("R_CRNt7") = 0
 Mgcaldc: - R("R_GCALDD") + R("R_DHNPA2") = 0
 Mgdpc: + R("R_URIK2") - R("R_MAN1PT2") + 2 R("R_GP4GH") - R("R_RNDR2")
 + R("R_CYTDK2") + R("R_ADSS") + R("R_SADT2") - R("R_NDPK1")
 + R("R_GK1") + R("R_ADK3") = 0
 Mgdpddmanc: + R("R_GMAND") - R("R_GDMANE") = 0
 Mgdpfucc: + R("R_GOFUCR") = 0
 Mgdpmannc: - R("R_GMAND") + R("R_MAN1PT2") = 0
 Mgdpofucc: - R("R_GOFUCR") + R("R_GDMANE") = 0
 MglcDb: + R("R_EX_glc_e_") + R("MglcDb_Transp") = 0
 MglcDc: + 2 R("R_TREH") + R("R_MLTG5") + R("R_MLTG4") + R("R_MLTG3")
 + R("R_MLTG2") + R("R_MLTG1") + R("R_GLCt2") + R("R_AMALT4")
 + R("R_AMALT3") + R("R_AMALT2") + R("R_AMALT1") + R("R_LACZ")
 + R("R_GALS3") + R("R_TRE6PH") + R("R_G1PP") + R("R_XYLI2i")
 - R("R_HEX1") = 0
 MglcDe: + 2 R("R_TREHe") - R("R_GLCDe") - R("R_GLCt2")
 - R("R_EX_glc_e_") - R("R_GLCpts") = 0
 Mglcnb: + R("R_EX_glcn_e_") + R("Mglcnb_Transp") = 0
 Mglcnc: + R("R_GLCNt2r") - R("R_GNK") + R("R_5DGLCNR")
 + R("R_2DGLCNRy") + R("R_2DGLCNRx") = 0
 Mglcne: - R("R_GLCNt2r") - R("R_EX_glcn_e_") + R("R_GLCDe") = 0
 Mglcrb: + R("R_EX_glcr_e_") + R("Mglcrb_Transp") = 0
 Mglcrc: + R("R_GLCRt2r") - R("R_GLCRD") = 0
 Mglcre: - R("R_GLCRt2r") - R("R_EX_glcr_e_") = 0
 Mglcurb: + R("R_EX_glcur_e_") + R("Mglcurb_Transp") = 0
 Mglcurc: + R("R_GLCURt2r") - R("R_GUI1") = 0
 Mglcure: - R("R_GLCURt2r") - R("R_EX_glcur_e_") = 0
 MglnLb: + R("R_EX_gln_L_e_") + R("MglnLb_Transp") = 0
 MglnLc: - R("R_GLUPRT") - R("R_GLUN") - R("R_GF6PTA") - R("R_ANS")
 - R("R_GMPS2") - R("R_ASNS1") - R("R_GLUSy") - R("R_IG3PS")
 - R("R_PRFGS") + R("R_GLNS") + R("R_GLNabc") - R("R_CTPS2")
 - R("R_CBPS") - R("R_ADCS") - 0.25 R("R_BiomassEcoli") = 0
 MglnLe: - R("R_EX_gln_L_e_") - R("R_GLNabc") = 0
 Mglu1sac: + R("R_GLUTRR") - R("R_G1SATi") = 0
 Mglu5pc: - R("R_G5SD") + R("R_GLU5K") = 0
 Mglu5sac: + R("R_G5SD") + R("R_ORNTA") + R("R_NACODA") - R("R_G5SADs")
 = 0
 Mglucysc: - R("R_GTHS") + R("R_GLUCYS") = 0
 MgluDc: - R("R_GLUR") - R("R_UAMAGS") = 0
 MgluLb: + R("R_EX_glu_L_e_") + R("MgluLb_Transp") = 0
 MgluLc: + R("R_SGDS") + R("R_GLUt4") + R("R_GLUt2r") + R("R_GLUR")
 + R("R_GLUPRT") + R("R_GLUN") + R("R_GF6PTA") + R("R_ANS")
 + R("R_GMPS2") - R("R_GLUTRS") + R("R_ASNS1") - R("R_TDPAGTA")
 + R("R_SOTA") + R("R_SDPTA") + R("R_PHETA1") + R("R_ORNTA")
 - R("R_OHPBAT") - R("R_HSTPT") + 2 R("R_GLUSy") - R("R_GLUDy")
 + R("R_ASPTA") + R("R_ALATA_L") + R("R_IG3PS") + R("R_PRFGS")
 - R("R_GLUCYS") + R("R_GLUabc") - R("R_GLU5K") - R("R_GLNS")
 - R("R_DHFS") + R("R_CTPS2") + R("R_CBPS") + R("R_ACOTA") - R("R_ACGS")
 - R("R_LEUTAi") + R("R_ADCS") + R("R_PTRCTA") - R("R_GLUDC")
 + R("R_GLUABUTt7") + R("R_ABTA") - R("R_PSERT") + R("R_ILETA")
 + R("R_VALTA") + R("R_TYRTA") - R("R_UNK3") + R("R_P5CD")
 - 0.25 R("R_BiomassEcoli") = 0
 MgluLe: - R("R_GLUt4") - R("R_GLUt2r") - R("R_EX_glu_L_e_")
 - R("R_GLUabc") - R("R_GLUABUTt7") = 0
 Mglutrnac: - R("R_GLUTRR") + R("R_GLUTRS") = 0
 Mglxc: + R("R_ICL") + R("R_GLYCTO3") + R("R_GLYCTO2") - R("R_GLYCLTDy")
 - R("R_GLYCLTDx") + R("R_UGLYCH") - R("R_MALS") - 2 R("R_GLXCL")
 + R("R_GLYCTO4") = 0
 Mglyaldb: + R("R_EX_glyald_e_") + R("Mglyaldb_Transp") = 0
 Mglyaldc: + R("R_GLYALDt") - R("R_ALCD19") = 0
 Mglyalde: - R("R_GLYALDt") - R("R_EX_glyald_e_") = 0
 Mglyb: + R("R_EX_gly_e_") + R("Mglyb_Transp") = 0
 Mglybb: + R("R_EX_glyb_e_") + R("Mglybb_Transp") = 0
 Mglybc: + R("R_GLYBt2r") + R("R_BETALDHy") + R("R_BETALDHx")
 + R("R_GLYBabc") = 0
 Mglybe: - R("R_GLYBt2r") - R("R_EX_glyb_e_") - R("R_GLYBabc") = 0
 Mglyc: + R("R_GLYt2r") + R("R_GHMT2") - R("R_GLYCL") - R("R_PRAGSr")
 - R("R_GTHS") + R("R_THRAr") + R("R_GLYATi")
 - 0.582 R("R_BiomassEcoli") = 0
 Mglyc3pb: + R("R_EX_glyc3p_e_") + R("Mglyc3pb_Transp") = 0
 Mglyc3pc: + R("R_GLYC3Pt6") + R("R_GPDDA3") + R("R_GPDDA5")
 + R("R_GPDDA4") + R("R_GPDDA2") - R("R_G3PD6") - R("R_G3PD5")
 - R("R_G3PD2") + R("R_GPDDA1") - R("R_PGSA_EC") + R("R_GLYK")
 + R("R_GLYC3Pabc") - R("R_PASYN_EC") - R("R_G3PD7") = 0
 Mglyc3pe: - R("R_GLYC3Pt6") - R("R_EX_glyc3p_e_") - R("R_GLYC3Pabc")
 = 0
 Mglycb: + R("R_EX_glyc_e_") + R("Mglycb_Transp") = 0
 Mglycc: - R("R_GLYCt") + R("R_ALCD19") + R("R_GPDDA4") - R("R_GLYCDx")
 + R("R_CLPNS_EC") - R("R_GLYK") = 0
 Mglyce: + R("R_GLYCt") - R("R_EX_glyc_e_") = 0
 Mglycltb: + R("R_EX_glyclt_e_") + R("Mglycltb_Transp") = 0
 Mglycltc: + R("R_GLYCLTt2r") - R("R_GLYCTO3") - R("R_GLYCTO2")
 + R("R_GLYCLTDy") + R("R_GLYCLTDx") + R("R_GCALDD") + R("R_PGLYCP")
 - R("R_GLYCTO4") = 0
 Mglyclte: - R("R_GLYCLTt2r") - R("R_EX_glyclt_e_") = 0
 Mglycogenc: - R("R_GLCP") + R("R_GLCS1") - 0.154 R("R_BiomassEcoli")
 = 0
 MglycRc: + R("R_HPYRRy") + R("R_HPYRRx") - R("R_GLYCK") + R("R_TRSAR")
 = 0
 Mglye: - R("R_GLYt2r") - R("R_EX_gly_e_") = 0
 Mgmhep17bpc: - R("R_GMHEPPA") + R("R_GMHEPK") = 0
 Mgmhep1pc: + R("R_GMHEPPA") - R("R_GMHEPAT") = 0
 Mgmhep7pc: + R("R_S7PI") - R("R_GMHEPK") = 0
 Mgmpc: + R("R_NTPP2") - R("R_NTD9") + R("R_GUAPRT") - R("R_GMPR")
 + R("R_GMPS2") + R("R_GSNK") - R("R_GK1") + R("R_ADOCBLS") = 0
 Mgp4gc: - R("R_GP4GH") = 0
 Mgrdpc: - R("R_GRTT") + R("R_DMATT") = 0
 Mgsnb: + R("R_EX_gsn_e_") + R("Mgsnb_Transp") = 0
 Mgsnc: - R("R_PUNP3") + R("R_NTPTP2") + R("R_GSNt2") + R("R_NTD9")
 - R("R_GSNK") = 0
 Mgsne: - R("R_GSNt2") - R("R_EX_gsn_e_") = 0
 Mgthoxc: - R("R_GTHOr") = 0
 Mgthrdc: - R("R_LGTHL") + R("R_GSPMDA") + R("R_GLYOX") + 2 R("R_GTHOr")
 + R("R_GTHS") - R("R_GSPMDS") = 0
 Mgtpc: - R("R_NTPTP2") - R("R_NTPP2") - R("R_URIK2") - R("R_RNTR2")
 - R("R_CYTDK2") - R("R_ADSS") - R("R_SADT2") - R("R_GTPCI")
 + R("R_NDPK1") - R("R_ADK3") - R("R_ACBIPGT") - R("R_GTPCII2")
 - 0.203 R("R_BiomassEcoli") = 0
 Mgtspmdc: - R("R_GSPMDA") + R("R_GSPMDS") = 0
 Mguab: + R("R_EX_gua_e_") + R("Mguab_Transp") = 0
 Mguac: + R("R_GUAt2") + R("R_GUAt") - R("R_GUAD") + R("R_PUNP3")
 - R("R_GUAPRT") + R("R_PUNP4") = 0
 Mguae: - R("R_GUAt2") - R("R_GUAt") - R("R_EX_gua_e_") = 0
 Mh2c: - R("R_HYD2") - R("R_HYD1") + R("R_FHL") - R("R_HYD3") = 0
 Mh2mb4pc: - R("R_IPDPS") - R("R_DMPPS") + R("R_MECDPDH") = 0
 Mh2o2c: + R("R_SPODM") + R("R_PYAM5PO") + R("R_PEAMNO") + R("R_PDX5PO")
 - 2 R("R_CAT") + R("R_ASPO6") = 0
 Mh2ob: + R("R_EX_h2o_e_") + R("Mh2ob_Transp") = 0
 Mh2oc: - R("R_USHD") + R("R_UPP3S") - R("R_UDPGD") - R("R_UDCPDP")
 - R("R_UACMAMO") + R("R_TRPS2") - R("R_TRPAS2") - R("R_TRE6PP")
 + R("R_TMAOR1") - R("R_THRS") + R("R_TARTD") - 3 R("R_SULR")
 - R("R_SSALy") - R("R_SSALx") - R("R_SGSAD") - R("R_PYDXPP")
 - R("R_PSP_L") - R("R_PRATPP") - R("R_PRAMPC") + 3 R("R_PPPGO")
 - R("R_PPA") - R("R_PGPP_EC") - R("R_PDXPP") + 2 R("R_NTRIR2x")
 - R("R_NTPP8") - R("R_NTD2") - R("R_NTD11") - R("R_NTD10")
 + R("R_NO3R2") + R("R_NO3R1") - R("R_NNAM") - R("R_NMNt7")
 - R("R_NMNN") - R("R_NMNDA") - R("R_NADPPPS") - R("R_MI1PP")
 - R("R_LCADi") - R("R_LCAD") - R("R_KDOPP") - R("R_IMPD")
 - R("R_HYPOE") - R("R_HMBS") - R("R_HKNDDH") - R("R_HISTP")
 - R("R_HISTD") + R("R_H2Ot") + R("R_CYTBO3") + R("R_CYTBD")
 - R("R_ALDD19x") - R("R_PYAM5PO") - R("R_PEAMNO") + 2 R("R_CAT")
 + R("R_IPDPS") - R("R_GUAD") - R("R_GSPMDA") - R("R_GLYOX")
 - R("R_NTPTP2") - R("R_NTPP2") - R("R_NTD9") - R("R_GMHEPPA")
 + R("R_GHMT2") - R("R_SGDS") - R("R_GLUPRT") - R("R_GLUN")
 - R("R_TREH") - R("R_MLTG5") - R("R_MLTG4") - R("R_MLTG3")
 - R("R_MLTG2") - R("R_MLTG1") + R("R_GMAND") - R("R_GP4GH")
 - R("R_GCALDD") - R("R_LACZ") - R("R_GALS3") - R("R_TRE6PH")
 - R("R_GPDDA3") - R("R_GPDDA5") - R("R_GPDDA4") - R("R_G1PP")
 - R("R_HKNTDH") - R("R_FUM") - R("R_FFSD") - R("R_HEMEOS")
 - R("R_IMPC") - R("R_G6PDA") - R("R_FBP") - R("R_GPDDA2")
 + R("R_IGPDH") + R("R_RNTR4") - R("R_NTD1") - R("R_DUTPDP")
 + R("R_RNDR4") - R("R_NTPP7") - R("R_NTD5") + R("R_TDPGDH")
 + R("R_DMSOR1") + R("R_DMPPS") + R("R_DHPTDCs") - R("R_DNMPPA")
 + 2 R("R_QULNS") + R("R_RNTR2") - R("R_NTPTP1") - R("R_NTPP1")
 - R("R_NTD8") + R("R_RNDR2") - R("R_DCYTD") - R("R_DCTPD")
 - R("R_NTPP3") - R("R_NTD3") - R("R_NTPP5") - R("R_NTD6") - R("R_DADA")
 - R("R_CYTD") - R("R_CYSTL") - R("R_CYSDS") + R("R_RNTR3")
 - R("R_CSND") + R("R_CRNCDH") - R("R_THDPS") - R("R_UGLYCH")
 - 2 R("R_SADH") + R("R_PPNDH") - R("R_PPC") + 2 R("R_PDX5PS")
 - R("R_HCO3E") + 2 R("R_CPPPGO") - R("R_NTPP4") - R("R_NTD4")
 - R("R_CMPN") - R("R_GPDDA1") - R("R_CDAPPA_EC") + R("R_RNDR3")
 - R("R_DHORTS") + R("R_BSORy") + R("R_BSORx") - R("R_BETALDHy")
 - R("R_BETALDHx") + R("R_RNTR1") - R("R_ASNN") - R("R_AB6PGH")
 - R("R_KDOPS") - R("R_SADT2") - R("R_SELNPS") - R("R_PPS")
 - R("R_NTPP6") - R("R_NADDP") - R("R_GMPS2") - R("R_BPNT")
 - R("R_ASNS1") - R("R_METAT") - R("R_ALLTAH") - R("R_ALLTN")
 - R("R_DAAD") - R("R_GLUDy") - R("R_GTPCI") - R("R_DNTPPA")
 - R("R_PLIPA1") - R("R_LPLIPA1") - R("R_PLIPA2") - R("R_LPLIPA2")
 - R("R_PLIPA3") - R("R_LPLIPA3") - R("R_AGMT") - R("R_XYLabc")
 - R("R_VALabc") - R("R_TSULabc") - R("R_THRabc") - R("R_THMabc")
 - R("R_TAURabc") - R("R_SULabc") - R("R_SUCCabc") - R("R_SPMDabc")
 + R("R_RNDR1") - R("R_RIBabc") - R("R_PTRCabc") - R("R_PROabc")
 - R("R_PRFGS") - R("R_PIabc") - R("R_ORNabc") - R("R_NAMNPP")
 - R("R_METDabc") - R("R_METabc") - R("R_MALTTTRabc") - R("R_MALTTRabc")
 - R("R_MALTPTabc") - R("R_MALTHXabc") - R("R_MALTabc") - R("R_LYSabc")
 - R("R_LEUabc") - R("R_Kabc") - R("R_ILEabc") - R("R_HISabc")
 - R("R_GLYC3Pabc") - R("R_GLYBabc") - R("R_GLUabc") - R("R_GLNabc")
 - R("R_GALabc") - R("R_FE2abc") - R("R_CYSabc") - R("R_CTPS2")
 - R("R_CHLabc") - R("R_CBPS") - R("R_CBL1abc") + R("R_ATPS4r")
 - R("R_ATPM") - R("R_ASPabc") - R("R_ASNabc") - R("R_ARGabc")
 - R("R_ARBabc") - R("R_AP5AH") - R("R_AP4AH") - R("R_ALAabc")
 - R("R_CBLAT_DELETE") - R("R_CBIAT_DELETE") - R("R_NTD7") - R("R_AHC")
 - R("R_ADA") - R("R_AMPN") - R("R_AHCYSNS") - R("R_ADD")
 + 8 R("R_C181SN") + 7 R("R_C161SN") + 7 R("R_C160SN") + 6 R("R_C141SN")
 + 6 R("R_C140SN") + 5 R("R_C120SN") - R("R_MALS") - 8 R("R_FAO3")
 - 7 R("R_FAO2") - 6 R("R_FAO1") - R("R_CS") - R("R_UHGADA")
 - R("R_POX") - R("R_NACODA") - R("R_AGDC") - R("R_ACODA")
 - R("R_ALDD2x") - R("R_FAO4") - R("R_PGL") - R("R_RZ5PP") - R("R_MTAN")
 + R("R_MDRPD") + R("R_GLCRD") + R("R_GALCTD") + 2 R("R_PPBNGS")
 + 2 R("R_RBFSa") - R("R_PMDPHT") - R("R_E4PD") - R("R_4HTHRS")
 + R("R_THZPSN") - R("R_OP4ENH") - R("R_AMPMS") - R("R_ABUTD")
 + R("R_TRPS1") + R("R_3HPPPNH") + R("R_3HCINNMH") + R("R_DHQD")
 - R("R_IPPS") + R("R_SUCBZS") - R("R_PGLYCP") + R("R_ENO")
 + R("R_MECDPDH") - R("R_MCITS") - R("R_MICITD") + R("R_MCITD")
 - 3 R("R_DKMPPD2") - R("R_DKMPPD") - R("R_IPPMIb") + R("R_IPPMIa")
 + R("R_TMAOR2") + R("R_DMSOR2") - R("R_MOHMT") + R("R_GALCTND")
 + R("R_MNNH") + R("R_ALTRH") + R("R_EDD") - R("R_DDPA") + R("R_IGPS")
 - R("R_DAPabc") - R("R_SDPDS") - 3 R("R_GTPCII2") - R("R_DHPPDA2")
 + R("R_DHAD2") + R("R_DHAD1") + 2 R("R_DHDPS") - R("R_ICHORT")
 - 2 R("R_P5CD") + R("R_G5SADs") - R("R_PAPA_EC") - R("R_MTHFC")
 - R("R_FTHFD") - 45.5608 R("R_BiomassEcoli") = 0
 Mh2oe: + R("R_TMAOR1e") - R("R_H2Ot") - R("R_EX_h2o_e_") - R("R_TREHe")
 - R("R_GLCDe") + R("R_DMSOR1e") - R("R_NADDPe") + R("R_TMAOR2e")
 + R("R_DMSOR2e") = 0
 Mh2sc: - R("R_SULR") + R("R_CYSDS") - R("R_CYSS") = 0
 Mhb: + R("R_EX_h_e_") + R("Mhb_Transp") = 0
 Mhc: + R("R_XYLt2") + R("R_XTSNt2r") + R("R_XANt2") + R("R_VALt2r")
 + R("R_URIt2r") + R("R_URIt2") + R("R_URAt2r") + R("R_URAt2")
 - R("R_UAPGR") + R("R_UAGPT3") + R("R_TYRt2r") + R("R_TTDCAt2")
 + R("R_TRPt2r") - R("R_TRDR") + R("R_THRt2r") + R("R_THMDt2r")
 + R("R_THMDt2") + 2 R("R_THD2") - R("R_SUCCt2b") + 3 R("R_SUCCt2_3")
 + 2 R("R_SUCCt2_2") + R("R_SERt2r") + R("R_RMNt") + R("R_PYRt2r")
 + R("R_PTRCt2r") + R("R_PROt2r") + R("R_PPTGS") + R("R_PPPNt2r")
 + R("R_PIt2r") + R("R_PHEt2r") + 2 R("R_PAPSR") + R("R_OCDCAt2")
 + R("R_NO2t2r") + 2 R("R_NAt3_2") + 3 R("R_NAt3_15") + R("R_NAt3_1")
 - 3 R("R_NADH7") - 4.5 R("R_NADH6") - R("R_NADH5") - R("R_NADH10")
 + R("R_MTHFD") + R("R_MELIBt2") + R("R_MDH") + 3 R("R_MALt2_3")
 + 2 R("R_MALt2_2") + R("R_LYSt2r") + R("R_LPADSS") + R("R_L_LACt2r")
 + R("R_LEUt2r") + R("R_LDH_D") + R("R_LCTSt") + R("R_Kt2r")
 + R("R_INSt2r") + R("R_INSt2") + R("R_INDOLEt2r") + R("R_ILEt2r")
 + R("R_IDONt2r") + R("R_HISt2r") + R("R_HDCAt2") + R("R_DSERt2r")
 + R("R_D_LACt2") + R("R_ACMAMUT") + 2 R("R_USHD") + 3 R("R_UDPGD")
 + R("R_UDCPDP") + 3 R("R_UACMAMO") - R("R_TMAOR1") + 5 R("R_SULR")
 + 2 R("R_SSALy") + 2 R("R_SSALx") + 2 R("R_SGSAD") + R("R_PRATPP")
 + R("R_PPA") - 5 R("R_NTRIR2x") + R("R_NTPP8") - 2 R("R_NO3R2")
 - 2 R("R_NO3R1") + R("R_NMNt7") + R("R_NMNN") + 2 R("R_LCADi")
 + 2 R("R_LCAD") + R("R_IMPD") + 2 R("R_HKNDDH") + 3 R("R_HISTD")
 - 2.5 R("R_CYTBO3") - 2 R("R_CYTBD") + 2 R("R_ALDD19x")
 - 2 R("R_SPODM") - R("R_IPDPS") - 2 R("R_HYD2") - 2 R("R_HYD1")
 + R("R_GUAt2") - R("R_GUAD") + R("R_GLYOX") - R("R_GTHOr")
 + R("R_GSNt2") + R("R_NTPP2") - 2 R("R_GMPR") - R("R_HPYRRy")
 - R("R_HPYRRx") + R("R_GLYCLTt2r") + R("R_GLYt2r") + R("R_GLYBt2r")
 - R("R_ALCD19") - R("R_GLYCLTDy") - R("R_GLYCLTDx") + R("R_GLUt2r")
 - R("R_G5SD") - R("R_GLUTRR") + R("R_GLCURt2r") + R("R_GLCRt2r")
 + R("R_GLCNt2r") + R("R_GLCt2") - R("R_GOFUCR") + R("R_URIK2")
 - R("R_MAN1PT2") + 2 R("R_GP4GH") + 2 R("R_GCALDD") + R("R_GALURt2r")
 + R("R_GLTPD") + R("R_GALCTNt2r") + R("R_GALCTt2r") + R("R_GALt2")
 + R("R_TRE6PS") + R("R_GPDDA3") + R("R_GPDDA5") + R("R_GPDDA4")
 - R("R_GALUi") - R("R_GALU") + 2 R("R_HKNTDH") + 3 R("R_FUMt2_3")
 + 2 R("R_FUMt2_2") + R("R_FUCt") + R("R_MANAO") + 3 R("R_SHCHF")
 + 2 R("R_FCLT") + R("R_SBTPD") + R("R_M1PD") + R("R_ETOHt2r")
 + R("R_GPDDA2") + R("R_ECAP_EC") + R("R_DURIt2") + R("R_DUTPDP")
 + R("R_NTPP7") - R("R_G1PTT") - R("R_TDPDRR") + R("R_AADDGT")
 + R("R_SHCHD2") - R("R_DMPPS") + R("R_DINSt2") - R("R_DHFR")
 + R("R_QULNS") + R("R_G3PD2") + R("R_GLYCDx") + R("R_DGSNt2")
 + R("R_NTPP1") + R("R_DCYTt2") - R("R_DCYTD") - R("R_DCTPD")
 + R("R_NTPP3") + R("R_DB4PS") + R("R_NTPP5") + R("R_DADNt2")
 - R("R_DADA") + R("R_CYTDt2r") + R("R_CYTDt2") - R("R_CYTD")
 + R("R_SHSL1") + R("R_CYNTt2") + R("R_CYANST") + R("R_CSNt2")
 - R("R_CSND") - 4 R("R_UPPDC1") - 2 R("R_UGLYCH") - 2 R("R_SADH")
 - R("R_PSD_EC") - R("R_PPNDH") + R("R_PPC") + R("R_PDX5PS")
 - R("R_ORNDC") - R("R_OMPDC") - 2 R("R_NNDPR") - R("R_MMCD")
 + R("R_HCO3E") - R("R_FHL") - 3 R("R_FDH3") - 3 R("R_FDH2")
 - R("R_DXPS") - 3 R("R_CYNTAH") - 2 R("R_CPPPGO") + R("R_NTPP4")
 + R("R_CYTDK2") + R("R_MOAT2") + R("R_MOAT") + R("R_GPDDA1")
 + R("R_CHLt2r") + R("R_DHCIND") - R("R_CINNDO") - R("R_PPPNDO")
 + R("R_DHPPD") + R("R_PSSA_EC") + R("R_PGSA_EC") - R("R_DASYN_EC")
 + 2 R("R_CDAPPA_EC") + R("R_OCBT") + R("R_DHORTS") + R("R_BUTt2r")
 - R("R_BSORy") - R("R_BSORx") + 2 R("R_BETALDHy") + 2 R("R_BETALDHx")
 - R("R_SERASr") - R("R_PTPATi") - R("R_NNAT") - R("R_NMNAT")
 - R("R_FMNAT") + R("R_HSDy") + 3 R("R_ASPt2_3") + 2 R("R_ASPt2_2")
 + R("R_ASPt2") + R("R_ASPCT") + 2 R("R_ADSS") + R("R_ASNt2r")
 + R("R_AST") + R("R_ARBt2r") + R("R_ANS") + R("R_PRPPS") + 2 R("R_PPS")
 + R("R_NTPP6") + R("R_NADS1") + 2 R("R_NADDP") + 2 R("R_GMPS2")
 + R("R_ASNS2") + R("R_ASNS1") + R("R_ARGSS") - R("R_ADMDCr")
 + R("R_TAGURr") + R("R_ALLTNt2r") + R("R_ALLTN") + 2 R("R_BTS2")
 + R("R_ALAt2r") + R("R_DALAt2r") - R("R_ACLS") + R("R_PANTS")
 - R("R_ASP1DC") - R("R_OXGDC2") - R("R_GLUSy") + R("R_GLUDy")
 + R("R_AKGt2r") + R("R_IG3PS") + R("R_DNTPPA") + R("R_UPP3MT")
 + R("R_PLIPA1") + R("R_LPLIPA1") + R("R_PLIPA2") + R("R_LPLIPA2")
 + R("R_PLIPA3") + R("R_LPLIPA3") - R("R_ARGDC") - R("R_GMHEPAT")
 - R("R_GLGC") + R("R_XYLK") + R("R_XYLabc") + R("R_VALabc")
 + R("R_UGMDDS") + R("R_UAMAS") + R("R_UAMAGS") + R("R_TSULabc")
 + R("R_TMKr") + R("R_TMDK1") + R("R_THRabc") + R("R_THMabc")
 + R("R_TDSK") + R("R_TAURabc") + R("R_SULabc") + R("R_SUCCabc")
 + R("R_SPMDabc") + R("R_SHKK") + R("R_RMK") + R("R_RIBabc")
 + R("R_RBK_L1") + R("R_RBK") + R("R_RBFK") - R("R_PYK") + R("R_PYDXNK")
 + R("R_PYDXK") + R("R_PYDAMK") + R("R_PTRCabc") + R("R_PROabc")
 + R("R_PRFGS") + 2 R("R_PRAIS") + R("R_PRAGSr") + R("R_PIabc")
 + R("R_PFK_2") + R("R_PFK") + R("R_ORNabc") + R("R_NADK")
 + R("R_METDabc") + R("R_METabc") + R("R_MALTTTRabc") + R("R_MALTTRabc")
 + R("R_MALTPTabc") + R("R_MALTHXabc") + R("R_MALTabc") + R("R_LYSabc")
 + 10 R("R_LPSSYN_EC") + R("R_LEUabc") + R("R_Kabc") + R("R_INSK")
 + R("R_ILEabc") + R("R_HSK") + R("R_HISabc") + R("R_HEX1")
 + R("R_GTHS") + R("R_GSPMDS") + R("R_GSNK") + R("R_GMHEPK")
 + R("R_GLYK") + R("R_GLYC3Pabc") + R("R_GLYBabc") + R("R_GLUCYS")
 + R("R_GLUabc") + R("R_GLNS") + R("R_GLNabc") + R("R_GLCS1")
 + R("R_GART") + R("R_GALKr") + R("R_GALabc") + R("R_FRUK")
 + R("R_FE2abc") + R("R_FCLK") + R("R_DXYLK") + R("R_DURIK1")
 + R("R_DPCOAK") + 3 R("R_DBTSr") + R("R_CYSabc") + 2 R("R_CTPS2")
 + R("R_CHLabc") + 2 R("R_CBPS") + 2 R("R_CBMK") + R("R_CBL1abc")
 + 3 R("R_ATPS4r") + R("R_ATPM") + R("R_ASPabc") + R("R_ASNabc")
 + R("R_ARGabc") + R("R_ARBabc") + 2 R("R_AP5AH") + 2 R("R_AP4AH")
 + R("R_ALAALAr") + R("R_ALAabc") + R("R_ADSK") + R("R_ADOCBLS")
 - R("R_ACBIPGT") + R("R_ADOCBIK") + R("R_ADNt2r") + R("R_ADNt2")
 + R("R_ADNK1") - R("R_ADA") + R("R_ADEt2r") - R("R_ADD")
 - 22 R("R_C181SN") - 19 R("R_C161SN") - 20 R("R_C160SN")
 - 16 R("R_C141SN") - 17 R("R_C140SN") - 14 R("R_C120SN") + R("R_ACPS1")
 + R("R_ACNAMt2") + R("R_AMANK") - R("R_UAGDP") + R("R_AGPR")
 + R("R_TDPADGAT") + R("R_SPMDAT2") + R("R_SPMDAT1") + R("R_MALS")
 - R("R_KAS15") + R("R_G1PACT") + 8 R("R_FAO3") + 7 R("R_FAO2")
 + 6 R("R_FAO1") + R("R_CS") - 2 R("R_ADHEr") + R("R_ACGS")
 + R("R_ACCOACr") + R("R_CYSS") + R("R_ACt2r") + 2 R("R_ALDD2x")
 + R("R_ACALDi") - R("R_MACPD") - R("R_KAS14") + R("R_ACACt2")
 + R("R_FAO4") + R("R_TAUDO") - R("R_AOXSr") + R("R_G6PDH2r")
 + R("R_PGL") + R("R_GNK") + R("R_HPPK2") + R("R_NNDMBRT")
 - R("R_MTHFR2") + R("R_SPMS") + R("R_MTRK") - R("R_IDOND2")
 - R("R_IDOND") - R("R_5DGLCNR") + R("R_AIRC2") - R("R_APRAUR")
 + R("R_PPBNGS") - R("R_PPCDC") + R("R_PPNCL2") + R("R_PNTK")
 + R("R_PERD") + 2 R("R_E4PD") + R("R_ASAD") + R("R_HETZK")
 + R("R_THZPSN") + 4 R("R_AMPMS") + R("R_HMPK1") + R("R_DHPS2")
 + R("R_ADCL") - R("R_GLUDC") + R("R_ABUTt2") + 2 R("R_ABUTD")
 + R("R_PGCD") + R("R_GLYCK") + R("R_HPPPNt2r") - R("R_3HPPPNH")
 + R("R_U23GAAT") - 2 R("R_KAS16") + R("R_HCINNMt2r") - R("R_3HCINNMH")
 - R("R_SHK3Dr") - R("R_KG6PDC") - R("R_OMCDC") + R("R_IPPS")
 + R("R_IPMD") + R("R_CDPMEK") - R("R_OPHBDC") + R("R_DMQMT")
 + R("R_OMBZLM") + R("R_OHPHM") - R("R_MECDPDH") - R("R_MEPCT")
 - R("R_DXPRIi") + R("R_MCITS") - R("R_TMPPP") + 6 R("R_DKMPPD2")
 + 2 R("R_DKMPPD") - R("R_TRSAR") - R("R_GLXCL") - R("R_TMAOR2")
 - R("R_NADH9") - 3.8 R("R_NADH8") - 2 R("R_HYD3") + R("R_DHNAOT")
 + R("R_AMMQT8_2") - R("R_DPR") - R("R_2DGULRy") - R("R_2DGULRx")
 - R("R_2DGLCNRy") - R("R_2DGLCNRx") + R("R_DDGALK") + R("R_DDGLCNt2r")
 + R("R_DDGLK") - R("R_IGPS") + R("R_THRD") - R("R_ACHBS")
 + R("R_UAAGDS") - R("R_DAPDC") + R("R_DAPabc") + 2 R("R_GTPCII2")
 - R("R_DHPPDA2") - R("R_DKGLCNR2y") - R("R_DKGLCNR2x")
 - R("R_DKGLCNR1") + R("R_PRASCS") - R("R_DOGULNR") - R("R_KARA2i")
 - R("R_KARA1i") + R("R_DHDPS") - R("R_DHDPRy") + 6 R("R_ENTCS")
 + R("R_DHBD") + R("R_PROD2") - 2 R("R_P5CR") + R("R_P5CD")
 + R("R_G5SADs") - R("R_LYSDC") + R("R_CADVt") + R("R_GAPD")
 - R("R_LCARS") - R("R_PEPT_EC") + R("R_DAGK_EC") + R("R_GARFT")
 + R("R_FTHFD") + 45.56035 R("R_BiomassEcoli") = 0
 Mhco3c: + R("R_HCO3E") - R("R_CYNTAH") - R("R_CBPS") - R("R_ACCOACr")
 - R("R_AIRC2") = 0
 MhcysLc: + R("R_RHCCE") + R("R_CYSTL") + R("R_AHC") - R("R_METS") = 0
 Mhdcab: + R("R_EX_hdca_e_") + R("Mhdcab_Transp") = 0
 Mhdcac: + R("R_HDCAt2") + 0.36 R("R_PLIPA1") + 0.36 R("R_LPLIPA1")
 + 0.36 R("R_PLIPA2") + 0.36 R("R_LPLIPA2") + 0.36 R("R_PLIPA3")
 + 0.36 R("R_LPLIPA3") - R("R_AACPS3") - R("R_FAO2") = 0
 Mhdcae: - R("R_EX_hdca_e_") - R("R_HDCAt2") = 0
 Mhdceac: + 0.07 R("R_PLIPA1") + 0.07 R("R_LPLIPA1")
 + 0.07 R("R_PLIPA2") + 0.07 R("R_LPLIPA2") + 0.07 R("R_PLIPA3")
 + 0.07 R("R_LPLIPA3") - R("R_AACPS4") = 0
 MhdeACPc: - 0.14 R("R_PASYN_EC") - R("R_EDTXS3") + R("R_C161SN")
 + R("R_AACPS4") = 0
 Mhe: - R("R_XYLt2") - R("R_XTSNt2r") - R("R_XANt2") - R("R_VALt2r")
 - R("R_URIt2r") - R("R_URIt2") - R("R_URAt2r") - R("R_URAt2")
 - R("R_TYRt2r") - R("R_TTDCAt2") - R("R_TRPt2r") - R("R_THRt2r")
 - R("R_THMDt2r") - R("R_THMDt2") - 2 R("R_THD2") + R("R_SUCCt2b")
 - 3 R("R_SUCCt2_3") - 2 R("R_SUCCt2_2") - R("R_SERt2r") - R("R_RMNt")
 - R("R_PYRt2r") - R("R_PTRCt2r") - R("R_PROt2r") - R("R_PPPNt2r")
 - R("R_PIt2r") - R("R_PHEt2r") - R("R_OCDCAt2") - R("R_NO2t2r")
 - 2 R("R_NAt3_2") - 3 R("R_NAt3_15") - R("R_NAt3_1") + 2 R("R_NADH7")
 + 3.5 R("R_NADH6") - R("R_MELIBt2") - 3 R("R_MALt2_3")
 - 2 R("R_MALt2_2") - R("R_LYSt2r") - R("R_L_LACt2r") - R("R_LEUt2r")
 - R("R_LCTSt") - R("R_Kt2r") - R("R_INSt2r") - R("R_INSt2")
 - R("R_INDOLEt2r") - R("R_ILEt2r") - R("R_IDONt2r") - R("R_HISt2r")
 - R("R_HDCAt2") - R("R_DSERt2r") - R("R_D_LACt2") - R("R_EX_h_e_")
 - R("R_TMAOR1e") + 2 R("R_NO3R2") + 2 R("R_NO3R1") + 2.5 R("R_CYTBO3")
 + 2 R("R_CYTBD") + 2 R("R_HYD2") + 2 R("R_HYD1") - R("R_GUAt2")
 - R("R_GSNt2") - R("R_GLYCLTt2r") - R("R_GLYt2r") - R("R_GLYBt2r")
 - R("R_GLUt2r") - R("R_GLCURt2r") - R("R_GLCRt2r") - R("R_GLCNt2r")
 + R("R_GLCDe") - R("R_GLCt2") - R("R_GALURt2r") - R("R_GALCTNt2r")
 - R("R_GALCTt2r") - R("R_GALt2") - 3 R("R_FUMt2_3") - 2 R("R_FUMt2_2")
 - R("R_FUCt") - R("R_ETOHt2r") - R("R_DURIt2") - R("R_DINSt2")
 - R("R_DGSNt2") - R("R_DCYTt2") - R("R_DADNt2") - R("R_CYTDt2r")
 - R("R_CYTDt2") - R("R_CYNTt2") - R("R_CSNt2") + 2 R("R_FDH3")
 + 2 R("R_FDH2") - R("R_CHLt2r") - R("R_BUTt2r") - 3 R("R_ASPt2_3")
 - 2 R("R_ASPt2_2") - R("R_ASPt2") - R("R_ASNt2r") - R("R_ARBt2r")
 + 2 R("R_NADDPe") - R("R_ALLTNt2r") - R("R_ALAt2r") - R("R_DALAt2r")
 - R("R_AKGt2r") - 4 R("R_ATPS4r") - R("R_ADNt2r") - R("R_ADNt2")
 - R("R_ADEt2r") - R("R_ACNAMt2") - R("R_ACt2r") - R("R_ACACt2")
 - R("R_ABUTt2") - R("R_HPPPNt2r") - R("R_HCINNMt2r") - R("R_TMAOR2e")
 + 2.8 R("R_NADH8") + 2 R("R_HYD3") - R("R_DDGLCNt2r") - R("R_CADVt")
 = 0
 MhemeOc: + R("R_HEMEOS") = 0
 MhisLb: + R("R_EX_his_L_e_") + R("MhisLb_Transp") = 0
 MhisLc: + R("R_HISt2r") + R("R_HISTD") + R("R_HISabc")
 - 0.09 R("R_BiomassEcoli") = 0
 MhisLe: - R("R_EX_his_L_e_") - R("R_HISt2r") - R("R_HISabc") = 0
 Mhispc: - R("R_HISTP") + R("R_HSTPT") = 0
 Mhistdc: + R("R_HISTP") - R("R_HISTD") = 0
 Mhknddc: - R("R_HKNDDH") + R("R_HPPPNDO") = 0
 Mhkntdc: - R("R_HKNTDH") + R("R_DHCINDO") = 0
 Mhmbilc: - R("R_UPP3S") + R("R_HMBS") = 0
 Mhmfurnc: + R("R_DHPTDCs") = 0
 MhomLc: - R("R_HSST") - R("R_HSDy") - R("R_HSK") = 0
 Mhpyrc: - R("R_HPYRRy") - R("R_HPYRRx") - R("R_HPYRI") = 0
 Mhqnc: + R("R_AB6PGH") = 0
 Mhxanb: + R("R_EX_hxan_e_") + R("Mhxanb_Transp") = 0
 Mhxanc: + R("R_PUNP5") + R("R_HYXNt") - R("R_HXPRT") + R("R_ADD")
 + R("R_PUNP6") = 0
 Mhxane: - R("R_HYXNt") - R("R_EX_hxan_e_") = 0
 Miaspc: - R("R_QULNS") + R("R_ASPO6") + R("R_ASPO5") + R("R_ASPO4")
 + R("R_ASPO3") = 0
 Michorc: + R("R_ICHORSi") - R("R_SHCHCS2") - R("R_ICHORT") = 0
 Micitc: - R("R_ICL") + R("R_ACONT") - R("R_ICDHyr") = 0
 MidonLb: + R("R_EX_idon_L_e_") + R("MidonLb_Transp") = 0
 MidonLc: + R("R_IDONt2r") + R("R_IDOND2") + R("R_IDOND")
 + R("R_2DGULRy") + R("R_2DGULRx") = 0
 MidonLe: - R("R_EX_idon_L_e_") - R("R_IDONt2r") = 0
 Midpc: + R("R_ADK4") = 0
 MileLb: + R("R_EX_ile_L_e_") + R("MileLb_Transp") = 0
 MileLc: + R("R_ILEt2r") + R("R_ILEabc") - R("R_ILETA")
 - 0.276 R("R_BiomassEcoli") = 0
 MileLe: - R("R_EX_ile_L_e_") - R("R_ILEt2r") - R("R_ILEabc") = 0
 Mimacpc: + R("R_IGPDH") - R("R_HSTPT") = 0
 Mimpc: + R("R_HXPRT") - R("R_NTD11") - R("R_IMPD") + R("R_GMPR")
 - R("R_IMPC") - R("R_ADSS") + R("R_INSK") = 0
 Mindoleb: + R("R_EX_indole_e_") + R("Mindoleb_Transp") = 0
 Mindolec: + R("R_INDOLEt2r") - R("R_TRPS2") + R("R_TRPAS2")
 + R("R_TRPS3") = 0
 Mindolee: - R("R_EX_indole_e_") - R("R_INDOLEt2r") = 0
 Minostc: + R("R_MI1PP") + R("R_GPDDA5") = 0
 Minsb: + R("R_EX_ins_e_") + R("Minsb_Transp") = 0
 Minsc: - R("R_PUNP5") + R("R_INSt2r") + R("R_INSt2") + R("R_NTD11")
 - R("R_INSK") + R("R_ADA") = 0
 Minse: - R("R_EX_ins_e_") - R("R_INSt2r") - R("R_INSt2") = 0
 Mipdpc: + R("R_IPDPS") - 8 R("R_UDCPDPS") - 5 R("R_OCTDPS")
 - R("R_GRTT") - R("R_IPDDIi") - R("R_DMATT") = 0
 Mitpc: - R("R_ADK4") = 0
 Mkb: + R("R_EX_k_e_") + R("Mkb_Transp") = 0
 Mkc: + R("R_Kt2r") + R("R_Kabc") = 0
 Mkdo2lipid4c: + R("R_MOAT2") - R("R_EDTXS3") - R("R_EDTXS1") = 0
 Mkdo2lipid4Lc: - R("R_EDTXS2") + R("R_EDTXS1") = 0
 Mkdo2lipid4pc: - R("R_EDTXS4") + R("R_EDTXS3") = 0
 Mkdo8pc: - R("R_KDOPP") + R("R_KDOPS") = 0
 Mkdoc: + R("R_KDOPP") - R("R_KDOCT2") = 0
 Mkdolipid4c: - R("R_MOAT2") + R("R_MOAT") = 0
 Mke: - R("R_EX_k_e_") - R("R_Kt2r") - R("R_Kabc") = 0
 MlacDb: + R("R_EX_lac_D_e_") + R("MlacDb_Transp") = 0
 MlacDc: - R("R_LDH_D2") - R("R_LDH_D") + R("R_D_LACt2") + R("R_GLYOX")
 = 0
 MlacDe: - R("R_EX_lac_D_e_") - R("R_D_LACt2") = 0
 MlacLb: + R("R_EX_lac_L_e_") + R("MlacLb_Transp") = 0
 MlacLc: - R("R_L_LACD3") - R("R_L_LACD2") + R("R_L_LACt2r")
 + R("R_LCADi") + R("R_LCAD") = 0
 MlacLe: - R("R_EX_lac_L_e_") - R("R_L_LACt2r") = 0
 MlaldLc: - R("R_LCADi") - R("R_LCAD") + R("R_RMPA") + R("R_FCLPA")
 - R("R_LCARS") = 0
 Mlctsb: + R("R_EX_lcts_e_") + R("Mlctsb_Transp") = 0
 Mlctsc: + R("R_LCTSt") - R("R_LACZ") = 0
 Mlctse: - R("R_EX_lcts_e_") - R("R_LCTSt") = 0
 MleuLb: + R("R_EX_leu_L_e_") + R("MleuLb_Transp") = 0
 MleuLc: + R("R_LEUt2r") + R("R_LEUabc") + R("R_LEUTAi")
 - 0.428 R("R_BiomassEcoli") = 0
 MleuLe: - R("R_EX_leu_L_e_") - R("R_LEUt2r") - R("R_LEUabc") = 0
 MlgtSc: + R("R_LGTHL") - R("R_GLYOX") = 0
 Mlipac: - R("R_LPSSYN_EC") + R("R_EDTXS2") = 0
 Mlipacoldc: + R("R_EDTXS4") = 0
 MlipidAc: - R("R_MOAT") + R("R_TDSK") = 0
 MlipidAdsc: + R("R_LPADSS") - R("R_TDSK") = 0
 MlipidXc: - R("R_LPADSS") + R("R_USHD") = 0
 MlpsECc: + R("R_LPSSYN_EC") - 0.0084 R("R_BiomassEcoli") = 0
 MlysLb: + R("R_EX_lys_L_e_") + R("MlysLb_Transp") = 0
 MlysLc: + R("R_LYSt2r") + R("R_LYSabc") + R("R_DAPDC") - R("R_LYSDC")
 + R("R_CADVt") - 0.326 R("R_BiomassEcoli") = 0
 MlysLe: - R("R_EX_lys_L_e_") - R("R_LYSt2r") - R("R_LYSabc")
 - R("R_CADVt") = 0
 MmalACPc: + R("R_MCOATA") - 7 R("R_C181SN") - 6 R("R_C161SN")
 - 6 R("R_C160SN") - 5 R("R_C141SN") - 5 R("R_C140SN") - 4 R("R_C120SN")
 - R("R_KAS15") - R("R_MACPD") - R("R_KAS14") - R("R_KAS16") = 0
 Mmalcoac: - R("R_MCOATA") + R("R_ACCOACr") = 0
 MmalLb: + R("R_EX_mal_L_e_") + R("MmalLb_Transp") = 0
 MmalLc: - R("R_MDH3") - R("R_MDH2") - R("R_MDH") + R("R_MALt2_3")
 + R("R_MALt2_2") + R("R_FUM") - R("R_ME2") - R("R_ME1") + R("R_MALS")
 = 0
 MmalLe: - R("R_EX_mal_L_e_") - R("R_MALt2_3") - R("R_MALt2_2") = 0
 Mmalt6pc: + R("R_MALTpts") = 0
 Mmaltb: + R("R_EX_malt_e_") + R("Mmaltb_Transp") = 0
 Mmaltc: + R("R_MLTG1") - R("R_AMALT4") - R("R_AMALT3") - R("R_AMALT2")
 - R("R_AMALT1") + R("R_MALTabc") = 0
 Mmalte: - R("R_EX_malt_e_") - R("R_MALTpts") - R("R_MALTabc") = 0
 Mmalthpc: - R("R_MLTG5") + R("R_AMALT4") - R("R_MLTP3") = 0
 Mmalthxb: + R("R_EX_malthx_e_") + R("Mmalthxb_Transp") = 0
 Mmalthxc: + R("R_MLTG5") - R("R_MLTG4") - R("R_AMALT4") + R("R_AMALT3")
 + R("R_MLTP3") - R("R_MLTP2") + R("R_MALTHXabc") = 0
 Mmalthxe: - R("R_EX_malthx_e_") - R("R_MALTHXabc") = 0
 Mmaltptb: + R("R_EX_maltpt_e_") + R("Mmaltptb_Transp") = 0
 Mmaltptc: + R("R_MLTG4") - R("R_MLTG3") - R("R_AMALT3") + R("R_AMALT2")
 + R("R_MLTP2") - R("R_MLTP1") + R("R_MALTPTabc") = 0
 Mmaltpte: - R("R_EX_maltpt_e_") - R("R_MALTPTabc") = 0
 Mmalttrb: + R("R_EX_malttr_e_") + R("Mmalttrb_Transp") = 0
 Mmalttrc: + R("R_MLTG2") - R("R_MLTG1") - R("R_AMALT1")
 + R("R_MALTTRabc") = 0
 Mmalttre: - R("R_EX_malttr_e_") - R("R_MALTTRabc") = 0
 Mmaltttrb: + R("R_EX_maltttr_e_") + R("Mmaltttrb_Transp") = 0
 Mmaltttrc: + R("R_MLTG3") - R("R_MLTG2") - R("R_AMALT2")
 + R("R_AMALT1") + R("R_MLTP1") + R("R_MALTTTRabc") = 0
 Mmaltttre: - R("R_EX_maltttr_e_") - R("R_MALTTTRabc") = 0
 Mman1pc: - R("R_PMANM") - R("R_MAN1PT2") = 0
 Mman6pb: + R("R_EX_man6p_e_") + R("Mman6pb_Transp") = 0
 Mman6pc: + R("R_MANpts") + R("R_MAN6Pt6_2") + R("R_PMANM")
 - R("R_MAN6PI") = 0
 Mman6pe: - R("R_MAN6Pt6_2") - R("R_EX_man6p_e_") = 0
 Mmanac: - R("R_MANAO") - R("R_MNNH") = 0
 Mmanb: + R("R_EX_man_e_") + R("Mmanb_Transp") = 0
 Mmane: - R("R_EX_man_e_") - R("R_MANpts") = 0
 Mmelibb: + R("R_EX_melib_e_") + R("Mmelibb_Transp") = 0
 Mmelibc: + R("R_MELIBt2") - R("R_GALS3") = 0
 Mmelibe: - R("R_EX_melib_e_") - R("R_MELIBt2") = 0
 MmetDb: + R("R_EX_met_D_e_") + R("MmetDb_Transp") = 0
 MmetDc: + R("R_METDabc") = 0
 MmetDe: - R("R_EX_met_D_e_") - R("R_METDabc") = 0
 Mmethfc: + R("R_MTHFD") - R("R_MTHFC") = 0
 MmetLb: + R("R_EX_met_L_e_") + R("MmetLb_Transp") = 0
 MmetLc: - R("R_METAT") + R("R_METabc") + R("R_METS") + R("R_UNK3")
 - 0.146 R("R_BiomassEcoli") = 0
 MmetLe: - R("R_EX_met_L_e_") - R("R_METabc") = 0
 Mmi1pDc: - R("R_MI1PP") = 0
 Mmicitc: - R("R_MCITL2") + R("R_MICITD") = 0
 Mmlthfc: - R("R_MTHFD") + R("R_GHMT2") - R("R_TMDS") + R("R_GLYCL")
 - R("R_MTHFR2") - R("R_MOHMT") = 0
 MmmcoaRc: + R("R_MMM2") - R("R_MME") = 0
 MmmcoaSc: + R("R_MME") - R("R_MMCD") = 0
 Mmnl1pc: + R("R_MNLpts") - R("R_M1PD") = 0
 Mmnlb: + R("R_EX_mnl_e_") + R("Mmnlb_Transp") = 0
 Mmnle: - R("R_EX_mnl_e_") - R("R_MNLpts") = 0
 Mmql8c: + R("R_MDH3") + R("R_L_LACD3") + R("R_NADH7") + R("R_NADH10")
 - R("R_TMAOR1e") - R("R_TMAOR1") - R("R_NO3R2") + R("R_HYD2")
 + R("R_GLYCTO3") - R("R_FRD2") - R("R_DMSOR1e") - R("R_DMSOR1")
 + R("R_DHORD5") + R("R_G3PD6") + R("R_FDH3") + R("R_ASPO4") = 0
 Mmqn8c: - R("R_MDH3") - R("R_L_LACD3") - R("R_NADH7") - R("R_NADH10")
 + R("R_TMAOR1e") + R("R_TMAOR1") + R("R_NO3R2") - R("R_HYD2")
 - R("R_GLYCTO3") + R("R_FRD2") + R("R_DMSOR1e") + R("R_DMSOR1")
 - R("R_DHORD5") - R("R_G3PD6") - R("R_FDH3") - R("R_ASPO4")
 + R("R_AMMQT8_2") = 0
 Mmthgxlc: - R("R_LGTHL") + R("R_MGSA") = 0
 MmyrsACPc: - 0.04 R("R_PASYN_EC") - R("R_EDTXS4") - R("R_EDTXS2")
 + R("R_C140SN") + R("R_AACPS1") = 0
 MN1aspmdc: + R("R_SPMDAT1") = 0
 Mn8aspmdc: + R("R_SPMDAT2") = 0
 Mna1b: + R("R_EX_na1_e_") + R("Mna1b_Transp") = 0
 Mna1c: + R("R_THRt4") + R("R_SERt4") + R("R_PROt4") + R("R_PNTOt4")
 - R("R_NAt3_2") - 2 R("R_NAt3_15") - R("R_NAt3_1") + R("R_GLUt4") = 0
 Mna1e: - R("R_THRt4") - R("R_SERt4") - R("R_PROt4") - R("R_PNTOt4")
 - R("R_EX_na1_e_") + R("R_NAt3_2") + 2 R("R_NAt3_15") + R("R_NAt3_1")
 - R("R_GLUt4") = 0
 Mnacb: + R("R_EX_nac_e_") + R("Mnacb_Transp") = 0
 Mnacc: + R("R_NACUP") + R("R_NNAM") - R("R_NAMNPP") + R("R_NNDMBRT")
 = 0
 Mnace: - R("R_NACUP") - R("R_EX_nac_e_") = 0
 Mnadb: + R("R_EX_nad_e_") + R("Mnadb_Transp") = 0
 Mnadc: - R("R_NADTRHD") + R("R_THD2") + R("R_NADH7") + R("R_NADH6")
 + R("R_NADH5") + R("R_NADH10") - R("R_MDH") - R("R_LDH_D")
 - 2 R("R_UDPGD") - 2 R("R_UACMAMO") - R("R_SSALx") - R("R_SGSAD")
 + 3 R("R_NTRIR2x") + R("R_NADPPPS") - R("R_LCADi") - R("R_LCAD")
 - R("R_IMPD") - 2 R("R_HISTD") - R("R_ALDD19x") + R("R_IPDPS")
 + R("R_HPYRRx") + R("R_ALCD19") + R("R_GLYCLTDx") - R("R_GCALDD")
 - R("R_GLTPD") - R("R_MANAO") - R("R_SBTPD") - R("R_M1PD")
 - R("R_SHCHD2") + R("R_DMPPS") - R("R_GLYCDx") - R("R_PDX5PS")
 - R("R_ME1") - R("R_GLYCL") - R("R_DHCIND") + R("R_CINNDO")
 + R("R_PPPNDO") - R("R_DHPPD") + R("R_BSORx") - R("R_BETALDHx")
 + R("R_NMNAT") + R("R_NADS1") - R("R_NADDP") - R("R_TAGURr")
 - R("R_AKGDH") - R("R_NADK") - R("R_PDH") - 8 R("R_FAO3")
 - 7 R("R_FAO2") - 6 R("R_FAO1") + 2 R("R_ADHEr") - R("R_ALDD2x")
 - R("R_ACALDi") - R("R_FAO4") + R("R_MTHFR2") + R("R_IDOND")
 - R("R_PERD") - R("R_E4PD") - R("R_ABUTD") - R("R_PGCD")
 + R("R_3HPPPNH") + R("R_3HCINNMH") - R("R_IPMD") - R("R_PPND")
 + R("R_TRSAR") + R("R_NADH9") + R("R_NADH8") + R("R_2DGULRx")
 + R("R_2DGLCNRx") - R("R_THRD") + R("R_DKGLCNR2x") + R("R_DOGULNR")
 - R("R_DHBD") - R("R_P5CD") - R("R_GAPD") + R("R_LCARS")
 - 0.00215 R("R_BiomassEcoli") = 0
 Mnade: - R("R_EX_nad_e_") - R("R_NADDPe") = 0
 Mnadhc: + R("R_NADTRHD") - R("R_THD2") - R("R_NADH7") - R("R_NADH6")
 - R("R_NADH5") - R("R_NADH10") + R("R_MDH") + R("R_LDH_D")
 + 2 R("R_UDPGD") + 2 R("R_UACMAMO") + R("R_SSALx") + R("R_SGSAD")
 - 3 R("R_NTRIR2x") + R("R_LCADi") + R("R_LCAD") + R("R_IMPD")
 + 2 R("R_HISTD") + R("R_ALDD19x") - R("R_IPDPS") - R("R_HPYRRx")
 - R("R_ALCD19") - R("R_GLYCLTDx") + R("R_GCALDD") + R("R_GLTPD")
 + R("R_MANAO") + R("R_SBTPD") + R("R_M1PD") + R("R_SHCHD2")
 - R("R_DMPPS") + R("R_GLYCDx") + R("R_PDX5PS") + R("R_ME1")
 + R("R_GLYCL") + R("R_DHCIND") - R("R_CINNDO") - R("R_PPPNDO")
 + R("R_DHPPD") - R("R_BSORx") + R("R_BETALDHx") + R("R_TAGURr")
 + R("R_AKGDH") + R("R_PDH") + 8 R("R_FAO3") + 7 R("R_FAO2")
 + 6 R("R_FAO1") - 2 R("R_ADHEr") + R("R_ALDD2x") + R("R_ACALDi")
 + R("R_FAO4") - R("R_MTHFR2") - R("R_IDOND") + R("R_PERD")
 + R("R_E4PD") + R("R_ABUTD") + R("R_PGCD") - R("R_3HPPPNH")
 - R("R_3HCINNMH") + R("R_IPMD") + R("R_PPND") - R("R_TRSAR")
 - R("R_NADH9") - R("R_NADH8") - R("R_2DGULRx") - R("R_2DGLCNRx")
 + R("R_THRD") - R("R_DKGLCNR2x") - R("R_DOGULNR") + R("R_DHBD")
 + R("R_P5CD") + R("R_GAPD") - R("R_LCARS") - 5e-05 R("R_BiomassEcoli")
 = 0
 Mnadpc: + R("R_NADTRHD") + R("R_UAPGR") + R("R_TRDR") - R("R_THD2")
 - R("R_MTHFD") - 3 R("R_SULR") - R("R_SSALy") - R("R_NADPPPS")
 + R("R_GTHOr") + R("R_GMPR") + R("R_HPYRRy") + R("R_GLYCLTDy")
 + R("R_G5SD") + R("R_GLUTRR") + R("R_GOFUCR") + R("R_TDPDRR")
 + R("R_DHFR") - R("R_G3PD2") - R("R_ME2") + R("R_BSORy")
 - R("R_BETALDHy") - R("R_HSDy") - R("R_ICDHyr") + R("R_GLUSy")
 - R("R_GLUDy") + R("R_NADK") + 15 R("R_C181SN") + 13 R("R_C161SN")
 + 14 R("R_C160SN") + 11 R("R_C141SN") + 12 R("R_C140SN")
 + 10 R("R_C120SN") - R("R_AGPR") - R("R_G6PDH2r") - R("R_GND")
 + R("R_IDOND2") + R("R_5DGLCNR") + R("R_APRAUR") - R("R_ASAD")
 + R("R_KAS16") + R("R_SHK3Dr") + R("R_DXPRIi") + R("R_DPR")
 + R("R_2DGULRy") + R("R_2DGLCNRy") + R("R_DKGLCNR2y") + R("R_DKGLCNR1")
 + R("R_KARA2i") + R("R_KARA1i") + R("R_DHDPRy") + R("R_P5CR")
 - 0.00013 R("R_BiomassEcoli") = 0
 Mnadphc: - R("R_NADTRHD") - R("R_UAPGR") - R("R_TRDR") + R("R_THD2")
 + R("R_MTHFD") + 3 R("R_SULR") + R("R_SSALy") - R("R_GTHOr")
 - R("R_GMPR") - R("R_HPYRRy") - R("R_GLYCLTDy") - R("R_G5SD")
 - R("R_GLUTRR") - R("R_GOFUCR") - R("R_TDPDRR") - R("R_DHFR")
 + R("R_G3PD2") + R("R_ME2") - R("R_BSORy") + R("R_BETALDHy")
 + R("R_HSDy") + R("R_ICDHyr") - R("R_GLUSy") + R("R_GLUDy")
 - 15 R("R_C181SN") - 13 R("R_C161SN") - 14 R("R_C160SN")
 - 11 R("R_C141SN") - 12 R("R_C140SN") - 10 R("R_C120SN") + R("R_AGPR")
 + R("R_G6PDH2r") + R("R_GND") - R("R_IDOND2") - R("R_5DGLCNR")
 - R("R_APRAUR") + R("R_ASAD") - R("R_KAS16") - R("R_SHK3Dr")
 - R("R_DXPRIi") - R("R_DPR") - R("R_2DGULRy") - R("R_2DGLCNRy")
 - R("R_DKGLCNR2y") - R("R_DKGLCNR1") - R("R_KARA2i") - R("R_KARA1i")
 - R("R_DHDPRy") - R("R_P5CR") - 0.0004 R("R_BiomassEcoli") = 0
 Mncamc: - R("R_NNAM") + R("R_NMNt7") + R("R_NMNN") = 0
 Mnh4b: + R("R_EX_nh4_e_") + R("Mnh4b_Transp") = 0
 Mnh4c: + R("R_SERD_L") + R("R_SERD_D") + R("R_NH4t") + R("R_TRPAS2")
 + R("R_NTRIR2x") + R("R_NNAM") + R("R_NMNDA") + 4 R("R_HMBS")
 + R("R_PYAM5PO") + R("R_PEAMNO") + R("R_GUAD") + R("R_GMPR")
 + R("R_GLUN") + R("R_G6PDA") + R("R_DCYTD") + R("R_DCTPD")
 + R("R_DADA") + R("R_CYTD") + R("R_CYSTL") + R("R_CYSDS") + R("R_CSND")
 + 2 R("R_UGLYCH") + 2 R("R_SADH") + R("R_GLYCL") + R("R_CYNTAH")
 + R("R_ASPT") + R("R_ASNN") - R("R_NADS1") - R("R_ASNS2") + R("R_DAAD")
 + R("R_GLUDy") - R("R_GLNS") - R("R_CBMK") + R("R_ADA") + R("R_ADD")
 + R("R_ETHAAL") + R("R_THRD_L") + R("R_DHPPDA2") = 0
 Mnh4e: - R("R_NH4t") - R("R_EX_nh4_e_") = 0
 Mnicrntc: + R("R_NMNDA") + R("R_NNDPR") - R("R_NNAT") + R("R_NAMNPP")
 - R("R_NNDMBRT") = 0
 Mnmnb: + R("R_EX_nmn_e_") + R("Mnmnb_Transp") = 0
 Mnmnc: + R("R_NMNP") - R("R_NMNN") - R("R_NMNDA") - R("R_NMNAT")
 + R("R_NADDP") = 0
 Mnmne: - R("R_NMNP") - R("R_EX_nmn_e_") - R("R_NMNt7") + R("R_NADDPe")
 = 0
 Mno2b: + R("R_EX_no2_e_") + R("Mno2b_Transp") = 0
 Mno2c: - R("R_NO3t7") + R("R_NO2t2r") - R("R_NTRIR2x") + R("R_NO3R2")
 + R("R_NO3R1") = 0
 Mno2e: + R("R_NO3t7") - R("R_EX_no2_e_") - R("R_NO2t2r") = 0
 Mno3b: + R("R_EX_no3_e_") + R("Mno3b_Transp") = 0
 Mno3c: + R("R_NO3t7") - R("R_NO3R2") - R("R_NO3R1") = 0
 Mno3e: - R("R_EX_no3_e_") - R("R_NO3t7") = 0
 Mo2b: + R("R_EX_o2_e_") + R("Mo2b_Transp") = 0
 Mo2c: + R("R_O2t") - 1.5 R("R_PPPGO") - 0.5 R("R_CYTBO3")
 - 0.5 R("R_CYTBD") + R("R_SPODM") - R("R_PYAM5PO") - R("R_PEAMNO")
 - R("R_PDX5PO") + R("R_CAT") - R("R_HPPPNDO") - R("R_DHCINDO")
 - R("R_CPPPGO") - R("R_CINNDO") - R("R_PPPNDO") - R("R_ASPO6")
 - R("R_TAUDO") - R("R_3HPPPNH") - R("R_3HCINNMH") - 0.5 R("R_OMMBLHX")
 - 0.5 R("R_OMPHHX") - 0.5 R("R_OPHHX") - R("R_DKMPPD") = 0
 Mo2e: - R("R_O2t") - R("R_EX_o2_e_") = 0
 Mo2SuperOxidec: - 2 R("R_SPODM") = 0
 Moaac: + R("R_MDH3") + R("R_MDH2") + R("R_MDH") + R("R_TARTD")
 + R("R_PPC") + R("R_ASPTA") - R("R_PPCK") - R("R_CS") + R("R_CITL")
 - R("R_MCITS") = 0
 Mocdcab: + R("R_EX_ocdca_e_") + R("Mocdcab_Transp") = 0
 Mocdcac: + R("R_OCDCAt2") - R("R_FAO3") = 0
 Mocdcae: - R("R_EX_ocdca_e_") - R("R_OCDCAt2") = 0
 Mocdceac: + 0.5 R("R_PLIPA1") + 0.5 R("R_LPLIPA1") + 0.5 R("R_PLIPA2")
 + 0.5 R("R_LPLIPA2") + 0.5 R("R_PLIPA3") + 0.5 R("R_LPLIPA3")
 - R("R_AACPS5") = 0
 Moctdpc: + R("R_OCTDPS") - R("R_HBZOPT") - R("R_DHNAOT") = 0
 MocteACPc: - R("R_PASYN_EC") + R("R_C181SN") + R("R_AACPS5") = 0
 Mohpbc: - R("R_OHPBAT") + R("R_PERD") = 0
 Mop4enc: + R("R_HKNDDH") + R("R_HKNTDH") - R("R_OP4ENH") = 0
 Mornb: + R("R_EX_orn_e_") + R("Mornb_Transp") = 0
 Mornc: - R("R_PTRCORNt7") - R("R_ORNDC") - R("R_OCBT")
 - R("R_ARGORNt7") - R("R_ORNTA") + R("R_ORNabc") + R("R_ACODA") = 0
 Morne: + R("R_PTRCORNt7") - R("R_EX_orn_e_") + R("R_ARGORNt7")
 - R("R_ORNabc") = 0
 Morot5pc: - R("R_ORPT") - R("R_OMPDC") = 0
 Morotc: + R("R_ORPT") + R("R_DHORD5") + R("R_DHORD2") = 0
 Mpacaldc: - R("R_ALDD19x") + R("R_PEAMNO") = 0
 Mpacc: + R("R_ALDD19x") - R("R_PACCOAL") = 0
 MpaECc: - 0.02 R("R_DASYN_EC") + 0.02 R("R_CDAPPA_EC")
 + 0.02 R("R_PASYN_EC") - 0.02 R("R_PAPA_EC") + 0.02 R("R_DAGK_EC") = 0
 MpalmACPc: - 0.72 R("R_PASYN_EC") + R("R_C160SN") + R("R_AACPS3") = 0
 Mpan4pc: - R("R_PTPATi") + R("R_PPCDC") = 0
 MpantRc: - R("R_PANTS") + R("R_DPR") = 0
 Mpapc: + R("R_PAPSR") - R("R_BPNT") + R("R_ACPS1") = 0
 Mpapsc: - R("R_PAPSR") + R("R_ADSK") = 0
 MpcECc: - 0.02 R("R_PLIPA3") = 0
 Mpdx5pc: - R("R_PDXPP") - R("R_PDX5PO") + R("R_PDX5PS") + R("R_PYDXNK")
 = 0
 Mpeamnc: - R("R_PEAMNO") = 0
 MpeECc: + 0.02 R("R_PSD_EC") - 0.02 R("R_PLIPA2") - 0.02 R("R_PEPT_EC")
 - 0.001935 R("R_BiomassEcoli") = 0
 Mpepc: - R("R_UAGCVT") - R("R_TREpts") - R("R_SUCpts") - R("R_SBTpts")
 - R("R_MNLpts") - R("R_MANpts") - R("R_MALTpts") - R("R_GAMpts")
 - R("R_GALTpts") - R("R_GLCpts") - R("R_FRUpts2") - R("R_FRUpts")
 - R("R_DHAPT") - R("R_PPC") - R("R_KDOPS") + R("R_PPS") - R("R_PYK")
 + R("R_PPCK") - R("R_ACMANApts") - R("R_ACGApts") - R("R_PSCVT")
 + R("R_ENO") - R("R_DDPA") = 0
 MpeptidoECc: + R("R_PPTGS") - 0.0276 R("R_BiomassEcoli") = 0
 MpgECc: + 0.02 R("R_PGPP_EC") - 0.04 R("R_CLPNS_EC")
 - 0.02 R("R_PLIPA1") - 0.02 R("R_LPLIPA6") - 0.02 R("R_LPLIPA4")
 - 0.02 R("R_LPLIPA5") - 0.000464 R("R_BiomassEcoli") = 0
 MpgpECc: - 0.02 R("R_PGPP_EC") + 0.02 R("R_PGSA_EC") = 0
 Mphaccoac: + R("R_PACCOAL") = 0
 MpheLb: + R("R_EX_phe_L_e_") + R("MpheLb_Transp") = 0
 MpheLc: + R("R_PHEt2r") - R("R_PHETA1") - 0.176 R("R_BiomassEcoli") = 0
 MpheLe: - R("R_EX_phe_L_e_") - R("R_PHEt2r") = 0
 Mphemec: - R("R_HEMEOS") + R("R_FCLT") = 0
 Mphomc: - R("R_THRS") + R("R_HSK") = 0
 Mphpyrc: + R("R_PPNDH") + R("R_PHETA1") = 0
 Mphthrc: - R("R_PDX5PS") + R("R_OHPBAT") - R("R_4HTHRS") = 0
 Mpib: + R("R_EX_pi_e_") + R("Mpib_Transp") = 0
 Mpic: - R("R_PYNP2r") - R("R_PUNP7") + R("R_UAGCVT")
 - 2 R("R_MAN6Pt6_2") - R("R_PUNP5") + R("R_PIt2r") + R("R_UDCPDP")
 + R("R_TRE6PP") + R("R_THRS") + R("R_PYDXPP") + R("R_PSP_L")
 + 2 R("R_PPA") + R("R_PGPP_EC") + R("R_PDXPP") + R("R_NTD2")
 + R("R_NTD11") + R("R_NTD10") + R("R_NADPPPS") + R("R_MI1PP")
 + R("R_KDOPP") + R("R_HYPOE") + R("R_HISTP") - R("R_PUNP3")
 + R("R_NTD9") + R("R_GMHEPPA") - R("R_GLYC3Pt6") + R("R_G5SD")
 + R("R_MAN1PT2") - 2 R("R_G6Pt6_2") - R("R_MLTP3") - R("R_MLTP2")
 - R("R_MLTP1") - R("R_GLCP") + R("R_G1PP") - 2 R("R_FUCPt6_2")
 + R("R_FBP") + R("R_NTD1") + R("R_NTD5") + R("R_DNMPPA") + R("R_QULNS")
 + R("R_MGSA") + R("R_NTD8") + R("R_NTD3") + R("R_NTD6") - R("R_PTA2")
 + R("R_PPC") + R("R_PDX5PS") + R("R_NTD4") + R("R_OCBT") + R("R_ASPCT")
 + R("R_ADSS") + R("R_KDOPS") + R("R_SADT2") + R("R_SELNPS")
 + R("R_PPS") + R("R_BPNT") + R("R_METAT") + R("R_XYLabc")
 + R("R_VALabc") + R("R_UGMDDS") + R("R_UAMAS") + R("R_UAMAGS")
 + R("R_TSULabc") + R("R_THRabc") + R("R_THMabc") + R("R_TAURabc")
 + R("R_SULabc") + R("R_SUCOAS") + R("R_SUCCabc") + R("R_SPMDabc")
 + R("R_RIBabc") + R("R_PTRCabc") + R("R_PROabc") + R("R_PRFGS")
 + R("R_PRAIS") + R("R_PRAGSr") + 2 R("R_PIabc") + R("R_ORNabc")
 + R("R_NAMNPP") + R("R_METDabc") + R("R_METabc") + R("R_MALTTTRabc")
 + R("R_MALTTRabc") + R("R_MALTPTabc") + R("R_MALTHXabc")
 + R("R_MALTabc") + R("R_LYSabc") + R("R_LEUabc") + R("R_Kabc")
 + R("R_ILEabc") + R("R_HISabc") + R("R_GTHS") + R("R_GSPMDS")
 + R("R_GLYC3Pabc") + R("R_GLYBabc") + R("R_GLUCYS") + R("R_GLUabc")
 + R("R_GLNS") + R("R_GLNabc") + R("R_GART") + R("R_GALabc")
 + R("R_FE2abc") + R("R_DHFS") + R("R_DBTSr") + R("R_CYSabc")
 + R("R_CTPS2") + R("R_CHLabc") + R("R_CBPS") + R("R_CBL1abc")
 - R("R_ATPS4r") + R("R_ATPM") + R("R_ASPabc") + R("R_ASNabc")
 + R("R_ARGabc") + R("R_ARBabc") + R("R_ALAALAr") + R("R_ALAabc")
 + R("R_ACCOAL") + R("R_CBLAT_DELETE") + R("R_CBIAT_DELETE")
 + R("R_NTD7") - R("R_PUNP1") - R("R_AGPR") - R("R_PTAr")
 + R("R_ACCOACr") + R("R_RZ5PP") + R("R_AIRC2") + R("R_RBFSa")
 + R("R_PMDPHT") - R("R_ASAD") + R("R_4HTHRS") + R("R_PSCVT")
 + R("R_CHORS") + R("R_PGLYCP") + R("R_DKMPPD2") + R("R_DKMPPD")
 - R("R_TMDPP") - R("R_PUNP6") - R("R_PUNP4") - R("R_PUNP2")
 - R("R_DURIPP") + R("R_DHQS") + R("R_DDPA") + R("R_UAAGDS")
 + R("R_DAPabc") + R("R_PRASCS") - R("R_GAPD") + R("R_PAPA_EC")
 + 45.5628 R("R_BiomassEcoli") = 0
 Mpie: - R("R_EX_pi_e_") + 2 R("R_MAN6Pt6_2") - R("R_PIt2r")
 + R("R_GLYC3Pt6") + 2 R("R_G6Pt6_2") + 2 R("R_FUCPt6_2") - R("R_PIabc")
 = 0
 Mpmcoac: - R("R_AOXSr") = 0
 MpntoRb: + R("R_EX_pnto_R_e_") + R("MpntoRb_Transp") = 0
 MpntoRc: + R("R_PNTOt4") + R("R_PANTS") - R("R_PNTK") = 0
 MpntoRe: - R("R_EX_pnto_R_e_") - R("R_PNTOt4") = 0
 Mppac: + R("R_PPCSCT") + R("R_PPAKr") - R("R_ACCOAL") = 0
 Mppapc: + R("R_PTA2") - R("R_PPAKr") = 0
 Mppbngc: - 4 R("R_HMBS") + R("R_PPBNGS") = 0
 Mppcoac: - R("R_PPCSCT") - R("R_PTA2") + R("R_MMCD") + R("R_ACCOAL")
 + R("R_OBTFL") - R("R_MCITS") = 0
 Mpphnc: - R("R_PPNDH") + R("R_CHORM") - R("R_PPND") = 0
 Mppic: + R("R_XPPT") + R("R_UPPRT") - R("R_ORPT") + R("R_HXPRT")
 + R("R_PRATPP") - R("R_PPA") + R("R_NTPP8") + R("R_NTPP2")
 + R("R_GUAPRT") + R("R_GLUPRT") + R("R_GALUi") + R("R_GALU")
 + 8 R("R_UDCPDPS") + 5 R("R_OCTDPS") + R("R_HEMEOS") + R("R_GRTT")
 + R("R_DUTPDP") + R("R_NTPP7") + R("R_G1PTT") + R("R_DMATT")
 + R("R_NTPP1") + R("R_NTPP3") + R("R_NTPP5") + R("R_NNDPR")
 + R("R_NTPP4") + R("R_KDOCT2") + R("R_DASYN_EC") + R("R_SERASr")
 + R("R_PTPATi") + R("R_NNAT") + R("R_NMNAT") + R("R_FMNAT")
 + R("R_ATPPRT") + R("R_ADNCYC") + R("R_SADT2") + R("R_ANPRT")
 + R("R_SUCBZL") + R("R_PACCOAL") + R("R_NTPP6") + R("R_NADS1")
 + R("R_GMPS2") + R("R_GLUTRS") + R("R_ASNS2") + R("R_ASNS1")
 + R("R_ARGSS") + R("R_METAT") + R("R_PANTS") + R("R_DNTPPA")
 + R("R_GMHEPAT") + R("R_GLGC") + R("R_NAMNPP") + R("R_CBLAT_DELETE")
 + R("R_ACBIPGT") + R("R_CBIAT_DELETE") + R("R_ADPT") + R("R_AACPS5")
 + R("R_AACPS4") + R("R_AACPS3") + R("R_AACPS2") + R("R_AACPS1")
 + R("R_UAGDP") + R("R_FAO3") + R("R_FAO2") + R("R_FAO1") + R("R_ACS")
 + R("R_PPNCL2") + R("R_THZPSN") + R("R_DHPS2") + R("R_HBZOPT")
 + R("R_MEPCT") + R("R_TMPPP") + R("R_DHNAOT") + R("R_GTPCII2")
 + R("R_DHBSr") + 0.7302 R("R_BiomassEcoli") = 0
 Mppp9c: + R("R_PPPGO") - R("R_FCLT") = 0
 Mpppg9c: - R("R_PPPGO") + R("R_CPPPGO") = 0
 Mpppic: + R("R_NTPTP2") + R("R_NTPTP1") = 0
 Mpppnb: + R("R_EX_pppn_e_") + R("Mpppnb_Transp") = 0
 Mpppnc: + R("R_PPPNt2r") - R("R_PPPNDO") = 0
 Mpppne: - R("R_EX_pppn_e_") - R("R_PPPNt2r") = 0
 Mpramc: + R("R_GLUPRT") - R("R_PRAGSr") = 0
 Mpranc: + R("R_ANPRT") - R("R_PRAIi") = 0
 Mprbampc: + R("R_PRATPP") - R("R_PRAMPC") = 0
 Mprbatpc: - R("R_PRATPP") + R("R_ATPPRT") = 0
 Mprfpc: - R("R_PRMICIi") + R("R_PRAMPC") = 0
 Mprlpc: + R("R_PRMICIi") - R("R_IG3PS") = 0
 MproLb: + R("R_EX_pro_L_e_") + R("MproLb_Transp") = 0
 MproLc: + R("R_PROt4") + R("R_PROt2r") + R("R_PROabc") - R("R_PROD2")
 + R("R_P5CR") - 0.21 R("R_BiomassEcoli") = 0
 MproLe: - R("R_EX_pro_L_e_") - R("R_PROt4") - R("R_PROt2r")
 - R("R_PROabc") = 0
 Mprppc: - R("R_XPPT") - R("R_UPPRT") + R("R_ORPT") - R("R_HXPRT")
 - R("R_GUAPRT") - R("R_GLUPRT") - R("R_NNDPR") - R("R_ATPPRT")
 - R("R_ANPRT") + R("R_PRPPS") - R("R_NAMNPP") - R("R_ADPT") = 0
 MpsECc: - 0.02 R("R_PSD_EC") + 0.02 R("R_PSSA_EC")
 - 5.2e-05 R("R_BiomassEcoli") = 0
 MpserLc: - R("R_PSP_L") + R("R_PSERT") = 0
 Mptrcb: + R("R_EX_ptrc_e_") + R("Mptrcb_Transp") = 0
 Mptrcc: + R("R_PTRCORNt7") + R("R_PTRCt2r") + R("R_ORNDC")
 + R("R_AGMT") + R("R_PTRCabc") - R("R_SPMS") - R("R_PTRCTA")
 - 0.035 R("R_BiomassEcoli") = 0
 Mptrce: - R("R_EX_ptrc_e_") - R("R_PTRCORNt7") - R("R_PTRCt2r")
 - R("R_PTRCabc") = 0
 Mpyam5pc: - R("R_HYPOE") - R("R_PYAM5PO") + R("R_ALATA_L2")
 + R("R_ALATA_D2") + R("R_PYDAMK") = 0
 Mpydamc: + R("R_HYPOE") - R("R_PYDAMK") = 0
 Mpydx5pc: - R("R_PYDXPP") + R("R_PYAM5PO") + R("R_PDX5PO")
 - R("R_ALATA_L2") - R("R_ALATA_D2") + R("R_PYDXK") = 0
 Mpydxc: + R("R_PYDXPP") - R("R_PYDXK") = 0
 Mpydxnc: + R("R_PDXPP") - R("R_PYDXNK") = 0
 Mpyrb: + R("R_EX_pyr_e_") + R("Mpyrb_Transp") = 0
 Mpyrc: + R("R_TREpts") + R("R_SUCpts") + R("R_SBTpts") + R("R_SERD_L")
 + R("R_SERD_D") + R("R_MNLpts") + R("R_MCITL2") + R("R_MANpts")
 + R("R_MALTpts") + R("R_L_LACD3") + R("R_L_LACD2") + R("R_LDH_D2")
 + R("R_PYRt2r") + R("R_LDH_D") + R("R_TRPAS2") + R("R_GAMpts")
 + R("R_GALTpts") + R("R_GLCpts") + R("R_FRUpts2") + R("R_FRUpts")
 + R("R_DHAPT") + R("R_CYSTL") + R("R_CYSDS") + R("R_ME2") + R("R_ME1")
 - R("R_DXPS") + R("R_ANS") - R("R_PPS") + R("R_ALATA_L2") + R("R_DAAD")
 + R("R_ALATA_D2") - 2 R("R_ACLS") + R("R_ALATA_L") + R("R_PYK")
 + R("R_ACMANApts") + R("R_ACNML") + R("R_ACGApts") - R("R_PFL")
 - R("R_PDH") - R("R_POX") + R("R_CHRPL") + R("R_HOPNTAL") + R("R_ADCL")
 + R("R_VPAMT") + R("R_SHCHCS2") + R("R_GLCRAL") + R("R_DDPGALA")
 + R("R_EDA") - R("R_ACHBS") - R("R_DHDPS") + R("R_ICHORT") = 0
 Mpyre: - R("R_EX_pyr_e_") - R("R_PYRt2r") = 0
 Mq8c: - R("R_MDH2") - R("R_L_LACD2") - R("R_LDH_D2") - R("R_NADH6")
 - R("R_NADH5") + R("R_NO3R1") + R("R_CYTBO3") + R("R_CYTBD")
 - R("R_HYD1") - R("R_GLYCTO2") - R("R_GLCDe") - R("R_SUCD4")
 - R("R_DHORD2") - R("R_G3PD5") - R("R_FDH2") - R("R_ASPO3")
 - R("R_POX") = 0
 Mq8h2c: + R("R_MDH2") + R("R_L_LACD2") + R("R_LDH_D2") + R("R_NADH6")
 + R("R_NADH5") - R("R_NO3R1") - R("R_CYTBO3") - R("R_CYTBD")
 + R("R_HYD1") + R("R_GLYCTO2") + R("R_GLCDe") + R("R_SUCD4")
 + R("R_DHORD2") + R("R_G3PD5") + R("R_FDH2") + R("R_ASPO3")
 + R("R_POX") + R("R_DMQMT") = 0
 Mqulnc: + R("R_QULNS") - R("R_NNDPR") = 0
 Mr1pc: - R("R_PPM") + R("R_PYNP2r") + R("R_PUNP7") + R("R_PUNP5")
 + R("R_PUNP3") + R("R_PUNP1") = 0
 Mr5pc: - R("R_RPI") + R("R_PPM") + R("R_NMNt7") + R("R_NMNN")
 - R("R_TKT1") + R("R_CMPN") - R("R_PRPPS") + R("R_RBK") + R("R_AMPN")
 = 0
 MrblLc: + R("R_ARAI") - R("R_RBK_L1") = 0
 Mrdmbzic: - R("R_ADOCBLS") + R("R_RZ5PP") = 0
 Mrhcysc: - R("R_RHCCE") + R("R_AHCYSNS") = 0
 MribDb: + R("R_EX_rib_D_e_") + R("MribDb_Transp") = 0
 MribDc: + R("R_RIBabc") - R("R_RBK") = 0
 MribDe: - R("R_EX_rib_D_e_") - R("R_RIBabc") = 0
 Mribflvc: - R("R_RBFK") + R("R_RBFSb") = 0
 Mrml1pc: - R("R_RMPA") + R("R_RMK") = 0
 Mrmlc: + R("R_RMI") - R("R_RMK") = 0
 Mrmnb: + R("R_EX_rmn_e_") + R("Mrmnb_Transp") = 0
 Mrmnc: - R("R_RMI") + R("R_RMNt") = 0
 Mrmne: - R("R_EX_rmn_e_") - R("R_RMNt") = 0
 Mru5pDc: - R("R_RPE") + R("R_RPI") - R("R_DB4PS") - R("R_A5PISO")
 + R("R_GND") = 0
 Mru5pLc: + R("R_X5PL3E") - R("R_RBP4E") + R("R_RBK_L1") = 0
 Ms7pc: - R("R_S7PI") + R("R_TKT1") - R("R_TALA") = 0
 Msbt6pc: + R("R_SBTpts") - R("R_SBTPD") = 0
 MsbtDb: + R("R_EX_sbt_D_e_") + R("MsbtDb_Transp") = 0
 MsbtDe: - R("R_EX_sbt_D_e_") - R("R_SBTpts") = 0
 Msbzcoac: - R("R_NPHS") + R("R_SUCBZL") = 0
 Msclc: - R("R_SHCHF") + R("R_SHCHD2") = 0
 Mselnc: - R("R_SELNPS") = 0
 Mselnpc: + R("R_SELNPS") = 0
 Mserampc: + R("R_SERASr") - 3 R("R_ENTCS") = 0
 MserDb: + R("R_EX_ser_D_e_") + R("MserDb_Transp") = 0
 MserDc: - R("R_SERD_D") + R("R_DSERt2r") = 0
 MserDe: - R("R_EX_ser_D_e_") - R("R_DSERt2r") = 0
 MserLb: + R("R_EX_ser_L_e_") + R("MserLb_Transp") = 0
 MserLc: - R("R_SERD_L") + R("R_SERt4") + R("R_SERt2r") - R("R_TRPS2")
 + R("R_PSP_L") - R("R_GHMT2") + R("R_GPDDA3") - R("R_PSSA_EC")
 - R("R_SERASr") - R("R_SERAT") - R("R_TRPS1")
 - 0.205 R("R_BiomassEcoli") = 0
 MserLe: - R("R_EX_ser_L_e_") - R("R_SERt4") - R("R_SERt2r") = 0
 Mshemec: + R("R_SHCHF") = 0
 Mskm5pc: + R("R_SHKK") - R("R_PSCVT") = 0
 Mskmc: - R("R_SHKK") + R("R_SHK3Dr") = 0
 Msl26dac: - R("R_SDPTA") - R("R_SDPDS") = 0
 Msl2a6oc: + R("R_THDPS") + R("R_SDPTA") = 0
 Mso3c: + R("R_PAPSR") + R("R_SULR") + R("R_CYANST") + R("R_TAUDO") = 0
 Mso4b: + R("R_EX_so4_e_") + R("Mso4b_Transp") = 0
 Mso4c: - R("R_SADT2") + R("R_SULabc") = 0
 Mso4e: - R("R_EX_so4_e_") - R("R_SULabc") = 0
 Mspmdb: + R("R_EX_spmd_e_") + R("Mspmdb_Transp") = 0
 Mspmdc: + R("R_GSPMDA") + R("R_SPMDabc") - R("R_GSPMDS")
 - R("R_SPMDAT2") - R("R_SPMDAT1") + R("R_SPMS")
 - 0.007 R("R_BiomassEcoli") = 0
 Mspmde: - R("R_EX_spmd_e_") - R("R_SPMDabc") = 0
 Mssaltppc: + R("R_OXGDC2") - R("R_SHCHCS2") = 0
 Msuc6pc: + R("R_SUCpts") - R("R_FFSD") = 0
 Msucargc: - R("R_SADH") + R("R_AST") = 0
 Msucbzc: - R("R_SUCBZL") + R("R_SUCBZS") = 0
 Msuccb: + R("R_EX_succ_e_") + R("Msuccb_Transp") = 0
 Msuccc: - R("R_TARTRt7") - R("R_PPCSCT") + R("R_MCITL2")
 - R("R_SUCCt2b") + R("R_SUCCt2_3") + R("R_SUCCt2_2") + R("R_SSALy")
 + R("R_SSALx") + R("R_HKNDDH") + R("R_ICL") + R("R_SGDS")
 - R("R_SUCFUMt") + R("R_FRD2") - R("R_SUCD1i") + R("R_SHSL1")
 - R("R_CITt7") + R("R_ASPO5") - R("R_SUCOAS") + R("R_SUCCabc")
 + R("R_TAUDO") + R("R_FRD3") + R("R_SDPDS") = 0
 Msucce: + R("R_TARTRt7") - R("R_EX_succ_e_") + R("R_SUCCt2b")
 - R("R_SUCCt2_3") - R("R_SUCCt2_2") + R("R_SUCFUMt") + R("R_CITt7")
 - R("R_SUCCabc") = 0
 Msuccoac: + R("R_PPCSCT") - R("R_MMM2") - R("R_THDPS") - R("R_HSST")
 - R("R_AST") + R("R_AKGDH") + R("R_SUCOAS") - 3e-06 R("R_BiomassEcoli")
 = 0
 Msucgluc: + R("R_SGSAD") - R("R_SGDS") = 0
 Msucgsac: - R("R_SGSAD") + R("R_SOTA") = 0
 Msuchmsc: - R("R_SHSL1") + R("R_HSST") = 0
 Msucornc: + R("R_SADH") - R("R_SOTA") = 0
 Msucrb: + R("R_EX_sucr_e_") + R("Msucrb_Transp") = 0
 Msucre: - R("R_EX_sucr_e_") - R("R_SUCpts") = 0
 Msucsalc: - R("R_SSALy") - R("R_SSALx") + R("R_ABTA") = 0
 Mtag6pDc: + R("R_GLTPD") - R("R_PFK_2") = 0
 MtagdpDc: - R("R_TGBPA") + R("R_PFK_2") = 0
 Mtagurc: + R("R_GUI2") + R("R_TAGURr") = 0
 MtartrLb: + R("R_EX_tartr_L_e_") + R("MtartrLb_Transp") = 0
 MtartrLc: + R("R_TARTRt7") - R("R_TARTD") = 0
 MtartrLe: - R("R_EX_tartr_L_e_") - R("R_TARTRt7") = 0
 Mtaurb: + R("R_EX_taur_e_") + R("Mtaurb_Transp") = 0
 Mtaurc: + R("R_TAURabc") - R("R_TAUDO") = 0
 Mtaure: - R("R_EX_taur_e_") - R("R_TAURabc") = 0
 Mtcyntc: + R("R_CYANST") = 0
 MtdeACPc: - 0.1 R("R_PASYN_EC") + R("R_C141SN") + R("R_AACPS2") = 0
 Mthdpc: - R("R_THDPS") + R("R_DHDPRy") = 0
 Mthfc: - R("R_GHMT2") + R("R_DHFR") - R("R_GLYCL") + R("R_METS")
 + R("R_MOHMT") + R("R_GARFT") + R("R_FTHFD") + R("R_AICART") = 0
 Mthmb: + R("R_EX_thm_e_") + R("Mthmb_Transp") = 0
 Mthmc: - R("R_TMKr") + R("R_THMabc") = 0
 Mthme: - R("R_EX_thm_e_") - R("R_THMabc") = 0
 Mthmmpc: - R("R_TMPKr") + R("R_TMKr") + R("R_TMPPP") = 0
 Mthmppc: - R("R_OXGDC2") + R("R_TMPKr") + R("R_SHCHCS2") = 0
 MthrLb: + R("R_EX_thr_L_e_") + R("MthrLb_Transp") = 0
 MthrLc: + R("R_THRt4") + R("R_THRt2r") + R("R_THRS") + R("R_THRabc")
 - R("R_THRAr") - R("R_THRD_L") - R("R_THRD")
 - 0.241 R("R_BiomassEcoli") = 0
 MthrLe: - R("R_EX_thr_L_e_") - R("R_THRt4") - R("R_THRt2r")
 - R("R_THRabc") = 0
 Mthymc: + R("R_TMDPP") = 0
 Mthymdb: + R("R_EX_thymd_e_") + R("Mthymdb_Transp") = 0
 Mthymdc: + R("R_THMDt2r") + R("R_THMDt2") + R("R_NTD5") - R("R_TMDK1")
 - R("R_TMDPP") = 0
 Mthymde: - R("R_EX_thymd_e_") - R("R_THMDt2r") - R("R_THMDt2") = 0
 Mtmab: + R("R_EX_tma_e_") + R("Mtmab_Transp") = 0
 Mtmac: + R("R_TMAOR1") + R("R_TMAOR2") = 0
 Mtmae: - R("R_EX_tma_e_") + R("R_TMAOR1e") + R("R_TMAOR2e") = 0
 Mtmaob: + R("R_EX_tmao_e_") + R("Mtmaob_Transp") = 0
 Mtmaoc: - R("R_TMAOR1") - R("R_TMAOR2") = 0
 Mtmaoe: - R("R_EX_tmao_e_") - R("R_TMAOR1e") - R("R_TMAOR2e") = 0
 Mtrdoxc: - R("R_TRDR") + R("R_PAPSR") + R("R_RNTR4") + R("R_RNDR4")
 + R("R_RNTR2") + R("R_RNDR2") + R("R_RNTR3") + R("R_RNDR3")
 + R("R_RNTR1") + R("R_RNDR1") = 0
 Mtrdrdc: + R("R_TRDR") - R("R_PAPSR") - R("R_RNTR4") - R("R_RNDR4")
 - R("R_RNTR2") - R("R_RNDR2") - R("R_RNTR3") - R("R_RNDR3")
 - R("R_RNTR1") - R("R_RNDR1") = 0
 Mtre6pc: + R("R_TREpts") - R("R_TRE6PP") + R("R_TRE6PS")
 - R("R_TRE6PH") = 0
 Mtreb: + R("R_EX_tre_e_") + R("Mtreb_Transp") = 0
 Mtrec: + R("R_TRE6PP") - R("R_TREH") = 0
 Mtree: - R("R_EX_tre_e_") - R("R_TREpts") - R("R_TREHe") = 0
 Mtrnagluc: + R("R_GLUTRR") - R("R_GLUTRS") = 0
 MtrpLb: + R("R_EX_trp_L_e_") + R("MtrpLb_Transp") = 0
 MtrpLc: + R("R_TRPt2r") + R("R_TRPS2") - R("R_TRPAS2") + R("R_TRPS1")
 - 0.054 R("R_BiomassEcoli") = 0
 MtrpLe: - R("R_EX_trp_L_e_") - R("R_TRPt2r") = 0
 Mtsulb: + R("R_EX_tsul_e_") + R("Mtsulb_Transp") = 0
 Mtsulc: - R("R_CYANST") + R("R_TSULabc") = 0
 Mtsule: - R("R_EX_tsul_e_") - R("R_TSULabc") = 0
 Mttdcab: + R("R_EX_ttdca_e_") + R("Mttdcab_Transp") = 0
 Mttdcac: + R("R_TTDCAt2") + 0.02 R("R_PLIPA1") + 0.02 R("R_LPLIPA1")
 + 0.02 R("R_PLIPA2") + 0.02 R("R_LPLIPA2") + 0.02 R("R_PLIPA3")
 + 0.02 R("R_LPLIPA3") - R("R_AACPS1") - R("R_FAO1") = 0
 Mttdcae: - R("R_EX_ttdca_e_") - R("R_TTDCAt2") = 0
 Mttdceac: + 0.05 R("R_PLIPA1") + 0.05 R("R_LPLIPA1")
 + 0.05 R("R_PLIPA2") + 0.05 R("R_LPLIPA2") + 0.05 R("R_PLIPA3")
 + 0.05 R("R_LPLIPA3") - R("R_AACPS2") = 0
 MtyrLb: + R("R_EX_tyr_L_e_") + R("MtyrLb_Transp") = 0
 MtyrLc: + R("R_TYRt2r") - R("R_THZPSN") - R("R_TYRTA")
 - 0.131 R("R_BiomassEcoli") = 0
 MtyrLe: - R("R_EX_tyr_L_e_") - R("R_TYRt2r") = 0
 Mu23gac: - R("R_LPADSS") - R("R_USHD") + R("R_U23GAAT") = 0
 Mu3agac: - R("R_UHGADA") + R("R_UAGAAT") = 0
 Mu3hgac: + R("R_UHGADA") - R("R_U23GAAT") = 0
 Muaagmdac: + R("R_UAGPT3") - R("R_PPTGS") = 0
 Muaccgc: + R("R_UAGCVT") - R("R_UAPGR") = 0
 Muacgamc: - R("R_UAG2Ei") - R("R_ACGAMT") - R("R_UAGCVT")
 - R("R_UAGPT3") + R("R_UAGDP") - R("R_UAGAAT") = 0
 Muacmamc: + R("R_UAG2Ei") - R("R_UACMAMO") = 0
 Muacmamuc: - R("R_ACMAMUT") + R("R_UACMAMO") = 0
 Muagmdac: + R("R_PAPPT3") - R("R_UAGPT3") = 0
 Muamac: + R("R_UAMAS") - R("R_UAMAGS") = 0
 Muamagc: + R("R_UAMAGS") - R("R_UAAGDS") = 0
 Muamrc: + R("R_UAPGR") - R("R_UAMAS") = 0
 Mudcpdpc: + R("R_PPTGS") - R("R_UDCPDP") + R("R_UDCPDPS")
 + R("R_ECAP_EC") = 0
 Mudcppc: - R("R_PAPPT3") - R("R_ACGAMT") + R("R_UDCPDP") = 0
 Mudpc: + R("R_UAGPT3") + R("R_LPADSS") + R("R_ACMAMUT") + R("R_TRE6PS")
 - R("R_RNDR4") + R("R_UMPK") - R("R_NDPK2") + 2 R("R_LPSSYN_EC") = 0
 Mudpgalc: - R("R_UDPGALM") + R("R_UDPG4E") + R("R_UGLT") = 0
 Mudpgalfurc: + R("R_UDPGALM") = 0
 Mudpgc: - R("R_UDPG4E") - R("R_UDPGD") - R("R_TRE6PS") - R("R_UGLT")
 + R("R_GALUi") + R("R_GALU") - 2 R("R_LPSSYN_EC")
 - 0.003 R("R_BiomassEcoli") = 0
 Mudpglcurc: + R("R_UDPGD") = 0
 Mugmdac: - R("R_PAPPT3") + R("R_UGMDDS") = 0
 Mugmdc: - R("R_UGMDDS") + R("R_UAAGDS") = 0
 Mumpc: + R("R_PAPPT3") + R("R_ACGAMT") + R("R_UPPRT") + R("R_USHD")
 + R("R_NTPP8") - R("R_NTD2") + R("R_URIK2") + R("R_OMPDC")
 - R("R_UMPK") = 0
 Munagac: + R("R_ACGAMT") - R("R_ACMAMUT") = 0
 Munagamuc: + R("R_ACMAMUT") - R("R_AADDGT") = 0
 Munagamufc: - R("R_ECAP_EC") + R("R_AADDGT") = 0
 Muppg3c: + R("R_UPP3S") - R("R_UPPDC1") - R("R_UPP3MT") = 0
 Murab: + R("R_EX_ura_e_") + R("Murab_Transp") = 0
 Murac: - R("R_UPPRT") + R("R_PYNP2r") + R("R_URAt2r") + R("R_URAt2")
 + R("R_CSND") + R("R_DURIPP") = 0
 Murae: - R("R_EX_ura_e_") - R("R_URAt2r") - R("R_URAt2") = 0
 Murdglycc: - R("R_UGLYCH") + R("R_ALLTAH") = 0
 Mureab: + R("R_EX_urea_e_") + R("Mureab_Transp") = 0
 Mureac: + R("R_UREAt") + R("R_ALLTAH") + R("R_AGMT") = 0
 Mureae: - R("R_UREAt") - R("R_EX_urea_e_") = 0
 Murib: + R("R_EX_uri_e_") + R("Murib_Transp") = 0
 Muric: - R("R_PYNP2r") + R("R_URIt2r") + R("R_URIt2") + R("R_NTD2")
 - R("R_URIK2") + R("R_CYTD") = 0
 Murie: - R("R_EX_uri_e_") - R("R_URIt2r") - R("R_URIt2") = 0
 Mutpc: - R("R_NTPP8") - R("R_GALUi") - R("R_GALU") - R("R_RNTR4")
 + R("R_NDPK2") - R("R_CTPS2") - R("R_UAGDP")
 - 0.136 R("R_BiomassEcoli") = 0
 MvalLb: + R("R_EX_val_L_e_") + R("MvalLb_Transp") = 0
 MvalLc: + R("R_VALt2r") + R("R_VALabc") + R("R_VPAMT") - R("R_VALTA")
 - 0.402 R("R_BiomassEcoli") = 0
 MvalLe: - R("R_EX_val_L_e_") - R("R_VALt2r") - R("R_VALabc") = 0
 Mxanb: + R("R_EX_xan_e_") + R("Mxanb_Transp") = 0
 Mxanc: + R("R_XANt") - R("R_XPPT") + R("R_PUNP7") + R("R_XANt2")
 + R("R_GUAD") = 0
 Mxane: - R("R_XANt") - R("R_EX_xan_e_") - R("R_XANt2") = 0
 Mxmpc: + R("R_XPPT") - R("R_NTD10") + R("R_IMPD") - R("R_GMPS2") = 0
 Mxtsnb: + R("R_EX_xtsn_e_") + R("Mxtsnb_Transp") = 0
 Mxtsnc: - R("R_PUNP7") + R("R_XTSNt2r") + R("R_NTD10") = 0
 Mxtsne: - R("R_EX_xtsn_e_") - R("R_XTSNt2r") = 0
 Mxu5pDc: + R("R_RBP4E") + R("R_RPE") - R("R_TKT1") - R("R_TKT2")
 + R("R_XYLK") = 0
 Mxu5pLc: - R("R_X5PL3E") + R("R_KG6PDC") = 0
 MxylDb: + R("R_EX_xyl_D_e_") + R("MxylDb_Transp") = 0
 MxylDc: - R("R_XYLI1") + R("R_XYLt2") + R("R_XYLabc") = 0
 MxylDe: - R("R_EX_xyl_D_e_") - R("R_XYLt2") - R("R_XYLabc") = 0
 MxyluDc: + R("R_XYLI1") - R("R_XYLK") = 0

Bounds
 0 <= R("R_BiomassEcoli") <= 999999
 -999999 <= R("R_AICART") <= 999999
 0 <= R("R_FTHFD") <= 999999
 -999999 <= R("R_GARFT") <= 999999
 -999999 <= R("R_MTHFC") <= 999999
 0 <= R("R_DAGK_EC") <= 999999
 0 <= R("R_PAPA_EC") <= 999999
 -999999 <= R("R_PEPT_EC") <= 999999
 -100000 <= R("M12ppdSb_Transp") <= 100000
 -999999 <= R("R_EX_12ppd_S_e_") <= 999999
 -999999 <= R("R_12PPDt") <= 999999
 -999999 <= R("R_LCARS") <= 999999
 -999999 <= R("R_GAPD") <= 999999
 -999999 <= R("R_PGK") <= 999999
 -100000 <= R("M15dapb_Transp") <= 100000
 -999999 <= R("R_EX_15dap_e_") <= 999999
 0 <= R("R_CADVt") <= 999999
 0 <= R("R_LYSDC") <= 999999
 0 <= R("R_G5SADs") <= 999999
 0 <= R("R_P5CD") <= 999999
 0 <= R("R_P5CR") <= 999999
 0 <= R("R_PROD2") <= 999999
 -999999 <= R("R_DHBD") <= 999999
 0 <= R("R_ICHORT") <= 999999
 -999999 <= R("R_DHBSr") <= 999999
 0 <= R("R_ENTCS") <= 999999
 0 <= R("R_DHDPRy") <= 999999
 0 <= R("R_DHDPS") <= 999999
 0 <= R("R_DHAD1") <= 999999
 0 <= R("R_KARA1i") <= 999999
 0 <= R("R_DHAD2") <= 999999
 0 <= R("R_KARA2i") <= 999999
 0 <= R("R_DOGULNR") <= 999999
 -999999 <= R("R_ADSL2r") <= 999999
 -999999 <= R("R_PRASCS") <= 999999
 0 <= R("R_DKGLCNR1") <= 999999
 0 <= R("R_DKGLCNR2x") <= 999999
 0 <= R("R_DKGLCNR2y") <= 999999
 0 <= R("R_DHPPDA2") <= 999999
 0 <= R("R_GTPCII2") <= 999999
 -999999 <= R("R_DAPE") <= 999999
 0 <= R("R_SDPDS") <= 999999
 -100000 <= R("M26dapMb_Transp") <= 100000
 -999999 <= R("R_EX_26dap_M_e_") <= 999999
 0 <= R("R_DAPabc") <= 999999
 0 <= R("R_DAPDC") <= 999999
 0 <= R("R_UAAGDS") <= 999999
 0 <= R("R_ACHBS") <= 999999
 0 <= R("R_GLYATi") <= 999999
 0 <= R("R_THRD") <= 999999
 0 <= R("R_IGPS") <= 999999
 0 <= R("R_PRAIi") <= 999999
 0 <= R("R_DDPA") <= 999999
 0 <= R("R_DHQS") <= 999999
 0 <= R("R_DDGLK") <= 999999
 0 <= R("R_EDA") <= 999999
 0 <= R("R_EDD") <= 999999
 -100000 <= R("M2ddglcnb_Transp") <= 100000
 -999999 <= R("R_EX_2ddglcn_e_") <= 999999
 0 <= R("R_ALTRH") <= 999999
 -999999 <= R("R_DDGLCNt2r") <= 999999
 0 <= R("R_MNNH") <= 999999
 0 <= R("R_DDGALK") <= 999999
 -999999 <= R("R_DDPGALA") <= 999999
 0 <= R("R_GALCTND") <= 999999
 0 <= R("R_2DGLCNRx") <= 999999
 0 <= R("R_2DGLCNRy") <= 999999
 0 <= R("R_2DGULRx") <= 999999
 0 <= R("R_2DGULRy") <= 999999
 0 <= R("R_DPR") <= 999999
 0 <= R("R_MOHMT") <= 999999
 0 <= R("R_AMMQT8_2") <= 999999
 0 <= R("R_DHNAOT") <= 999999
 0 <= R("R_DMSOR2") <= 999999
 0 <= R("R_DMSOR2e") <= 999999
 0 <= R("R_FRD3") <= 999999
 0 <= R("R_G3PD7") <= 999999
 0 <= R("R_GLYCTO4") <= 999999
 0 <= R("R_HYD3") <= 999999
 0 <= R("R_NADH8") <= 999999
 0 <= R("R_NADH9") <= 999999
 0 <= R("R_TMAOR2") <= 999999
 0 <= R("R_TMAOR2e") <= 999999
 -999999 <= R("R_DURIPP") <= 999999
 -999999 <= R("R_PPM2") <= 999999
 -999999 <= R("R_PUNP2") <= 999999
 -999999 <= R("R_PUNP4") <= 999999
 -999999 <= R("R_PUNP6") <= 999999
 -999999 <= R("R_TMDPP") <= 999999
 0 <= R("R_DRPA") <= 999999
 0 <= R("R_GLCRAL") <= 999999
 0 <= R("R_GLXCL") <= 999999
 -999999 <= R("R_HPYRI") <= 999999
 0 <= R("R_TRSAR") <= 999999
 -999999 <= R("R_IPPMIa") <= 999999
 -999999 <= R("R_IPPMIb") <= 999999
 0 <= R("R_DKMPPD") <= 999999
 0 <= R("R_DKMPPD2") <= 999999
 0 <= R("R_UNK3") <= 999999
 0 <= R("R_PMPK") <= 999999
 0 <= R("R_TMPPP") <= 999999
 0 <= R("R_MCITD") <= 999999
 0 <= R("R_MICITD") <= 999999
 0 <= R("R_MCITS") <= 999999
 0 <= R("R_DXPRIi") <= 999999
 0 <= R("R_MEPCT") <= 999999
 0 <= R("R_MECDPDH") <= 999999
 0 <= R("R_MECDPS") <= 999999
 0 <= R("R_OBTFL") <= 999999
 0 <= R("R_THRD_L") <= 999999
 0 <= R("R_OHPHM") <= 999999
 0 <= R("R_OPHHX") <= 999999
 0 <= R("R_OMBZLM") <= 999999
 0 <= R("R_OMPHHX") <= 999999
 0 <= R("R_DMQMT") <= 999999
 0 <= R("R_OMMBLHX") <= 999999
 0 <= R("R_OPHBDC") <= 999999
 0 <= R("R_CDPMEK") <= 999999
 -999999 <= R("R_ENO") <= 999999
 -999999 <= R("R_PGM") <= 999999
 0 <= R("R_PGLYCP") <= 999999
 0 <= R("R_SHCHCS2") <= 999999
 0 <= R("R_SUCBZS") <= 999999
 0 <= R("R_PPND") <= 999999
 -999999 <= R("R_TYRTA") <= 999999
 0 <= R("R_IPMD") <= 999999
 0 <= R("R_IPPS") <= 999999
 0 <= R("R_OMCDC") <= 999999
 0 <= R("R_KG6PDC") <= 999999
 -999999 <= R("R_DHQD") <= 999999
 -999999 <= R("R_SHK3Dr") <= 999999
 -100000 <= R("M3hcinnmb_Transp") <= 100000
 -999999 <= R("R_EX_3hcinnm_e_") <= 999999
 0 <= R("R_3HCINNMH") <= 999999
 -999999 <= R("R_HCINNMt2r") <= 999999
 0 <= R("R_KAS16") <= 999999
 0 <= R("R_U23GAAT") <= 999999
 -999999 <= R("R_UAGAAT") <= 999999
 -100000 <= R("M3hpppnb_Transp") <= 100000
 -999999 <= R("R_EX_3hpppn_e_") <= 999999
 0 <= R("R_3HPPPNH") <= 999999
 -999999 <= R("R_HPPPNt2r") <= 999999
 0 <= R("R_TRPS1") <= 999999
 0 <= R("R_TRPS3") <= 999999
 -999999 <= R("R_VALTA") <= 999999
 0 <= R("R_VPAMT") <= 999999
 -999999 <= R("R_ILETA") <= 999999
 0 <= R("R_HBZOPT") <= 999999
 0 <= R("R_GLYCK") <= 999999
 0 <= R("R_PGCD") <= 999999
 0 <= R("R_PSERT") <= 999999
 0 <= R("R_CHORS") <= 999999
 -999999 <= R("R_PSCVT") <= 999999
 -100000 <= R("M4abutb_Transp") <= 100000
 -999999 <= R("R_EX_4abut_e_") <= 999999
 0 <= R("R_ABTA") <= 999999
 0 <= R("R_ABUTD") <= 999999
 0 <= R("R_ABUTt2") <= 999999
 -999999 <= R("R_GLUABUTt7") <= 999999
 0 <= R("R_GLUDC") <= 999999
 0 <= R("R_PTRCTA") <= 999999
 0 <= R("R_ADCL") <= 999999
 0 <= R("R_DHPS2") <= 999999
 0 <= R("R_ADCS") <= 999999
 0 <= R("R_HMPK1") <= 999999
 0 <= R("R_AMPMS") <= 999999
 0 <= R("R_HOPNTAL") <= 999999
 0 <= R("R_OP4ENH") <= 999999
 0 <= R("R_THZPSN") <= 999999
 0 <= R("R_CHRPL") <= 999999
 0 <= R("R_4HTHRS") <= 999999
 0 <= R("R_HETZK") <= 999999
 0 <= R("R_LEUTAi") <= 999999
 -999999 <= R("R_ASAD") <= 999999
 -999999 <= R("R_ASPK") <= 999999
 -999999 <= R("R_E4PD") <= 999999
 -999999 <= R("R_PERD") <= 999999
 0 <= R("R_PNTK") <= 999999
 0 <= R("R_PPNCL2") <= 999999
 0 <= R("R_PPCDC") <= 999999
 0 <= R("R_PMDPHT") <= 999999
 0 <= R("R_RBFSa") <= 999999
 0 <= R("R_RBFSb") <= 999999
 -999999 <= R("R_AIRC3") <= 999999
 0 <= R("R_G1SATi") <= 999999
 0 <= R("R_PPBNGS") <= 999999
 0 <= R("R_APRAUR") <= 999999
 0 <= R("R_AIRC2") <= 999999
 -999999 <= R("R_5DGLCNR") <= 999999
 -999999 <= R("R_IDOND") <= 999999
 0 <= R("R_IDOND2") <= 999999
 0 <= R("R_GALCTD") <= 999999
 0 <= R("R_GLCRD") <= 999999
 -999999 <= R("R_MTRI") <= 999999
 0 <= R("R_MTRK") <= 999999
 0 <= R("R_MDRPD") <= 999999
 0 <= R("R_MTAN") <= 999999
 0 <= R("R_SPMS") <= 999999
 0 <= R("R_METS") <= 999999
 0 <= R("R_MTHFR2") <= 999999
 0 <= R("R_NNDMBRT") <= 999999
 0 <= R("R_RZ5PP") <= 999999
 0 <= R("R_DHNPA2") <= 999999
 0 <= R("R_HPPK2") <= 999999
 0 <= R("R_GND") <= 999999
 0 <= R("R_GNK") <= 999999
 0 <= R("R_PGL") <= 999999
 -999999 <= R("R_G6PDH2r") <= 999999
 -999999 <= R("R_AMAOTr") <= 999999
 -999999 <= R("R_AOXSr") <= 999999
 0 <= R("R_TAUDO") <= 999999
 0 <= R("R_ACACCT") <= 999999
 -999999 <= R("R_ACACT1r") <= 999999
 0 <= R("R_FAO4") <= 999999
 -100000 <= R("Macacb_Transp") <= 100000
 -999999 <= R("R_EX_acac_e_") <= 999999
 -999999 <= R("R_ACACt2") <= 999999
 -999999 <= R("R_ACOATA") <= 999999
 0 <= R("R_KAS14") <= 999999
 0 <= R("R_MACPD") <= 999999
 -100000 <= R("Macaldb_Transp") <= 100000
 -999999 <= R("R_EX_acald_e_") <= 999999
 0 <= R("R_ACALDi") <= 999999
 -999999 <= R("R_ACALDt") <= 999999
 0 <= R("R_ALDD2x") <= 999999
 0 <= R("R_ETHAAL") <= 999999
 -999999 <= R("R_THRAr") <= 999999
 -100000 <= R("Macb_Transp") <= 100000
 -999999 <= R("R_EX_ac_e_") <= 999999
 -999999 <= R("R_ACKr") <= 999999
 0 <= R("R_ACODA") <= 999999
 0 <= R("R_ACS") <= 999999
 -999999 <= R("R_ACt2r") <= 999999
 0 <= R("R_AGDC") <= 999999
 0 <= R("R_BUTCT") <= 999999
 0 <= R("R_CITL") <= 999999
 0 <= R("R_CYSS") <= 999999
 0 <= R("R_NACODA") <= 999999
 0 <= R("R_POX") <= 999999
 0 <= R("R_UHGADA") <= 999999
 -999999 <= R("R_ACCOACr") <= 999999
 0 <= R("R_ACGS") <= 999999
 -999999 <= R("R_ADHEr") <= 999999
 0 <= R("R_CS") <= 999999
 0 <= R("R_FAO1") <= 999999
 0 <= R("R_FAO2") <= 999999
 0 <= R("R_FAO3") <= 999999
 0 <= R("R_G1PACT") <= 999999
 0 <= R("R_KAS15") <= 999999
 0 <= R("R_MALS") <= 999999
 0 <= R("R_PDH") <= 999999
 0 <= R("R_PFL") <= 999999
 -999999 <= R("R_PTAr") <= 999999
 -999999 <= R("R_SERAT") <= 999999
 0 <= R("R_SPMDAT1") <= 999999
 0 <= R("R_SPMDAT2") <= 999999
 0 <= R("R_TDPADGAT") <= 999999
 0 <= R("R_ACGK") <= 999999
 -999999 <= R("R_AGPR") <= 999999
 -999999 <= R("R_ACOTA") <= 999999
 0 <= R("R_UAGDP") <= 999999
 0 <= R("R_ACGApts") <= 999999
 0 <= R("R_AMANAPE") <= 999999
 -100000 <= R("Macgamb_Transp") <= 100000
 -999999 <= R("R_EX_acgam_e_") <= 999999
 -100000 <= R("Macmanab_Transp") <= 100000
 -999999 <= R("R_EX_acmana_e_") <= 999999
 0 <= R("R_ACNML") <= 999999
 0 <= R("R_AMANK") <= 999999
 0 <= R("R_ACMANApts") <= 999999
 -100000 <= R("Macnamb_Transp") <= 100000
 -999999 <= R("R_EX_acnam_e_") <= 999999
 0 <= R("R_ACNAMt2") <= 999999
 0 <= R("R_ACONMT") <= 999999
 0 <= R("R_AACPS1") <= 999999
 0 <= R("R_AACPS2") <= 999999
 0 <= R("R_AACPS3") <= 999999
 0 <= R("R_AACPS4") <= 999999
 0 <= R("R_AACPS5") <= 999999
 0 <= R("R_ACPS1") <= 999999
 0 <= R("R_C120SN") <= 999999
 0 <= R("R_C140SN") <= 999999
 0 <= R("R_C141SN") <= 999999
 0 <= R("R_C160SN") <= 999999
 0 <= R("R_C161SN") <= 999999
 0 <= R("R_C181SN") <= 999999
 0 <= R("R_EDTXS1") <= 999999
 0 <= R("R_EDTXS2") <= 999999
 0 <= R("R_EDTXS3") <= 999999
 0 <= R("R_EDTXS4") <= 999999
 -999999 <= R("R_MCOATA") <= 999999
 0 <= R("R_PASYN_EC") <= 999999
 -100000 <= R("Madeb_Transp") <= 100000
 -999999 <= R("R_EX_ade_e_") <= 999999
 0 <= R("R_ADD") <= 999999
 -999999 <= R("R_ADEt2r") <= 999999
 0 <= R("R_ADPT") <= 999999
 0 <= R("R_AHCYSNS") <= 999999
 0 <= R("R_AMPN") <= 999999
 -999999 <= R("R_PUNP1") <= 999999
 -100000 <= R("Madnb_Transp") <= 100000
 -999999 <= R("R_EX_adn_e_") <= 999999
 0 <= R("R_ADA") <= 999999
 0 <= R("R_ADNK1") <= 999999
 0 <= R("R_ADNt2") <= 999999
 -999999 <= R("R_ADNt2r") <= 999999
 -999999 <= R("R_AHC") <= 999999
 0 <= R("R_NTD7") <= 999999
 0 <= R("R_ADOCBIK") <= 999999
 -999999 <= R("R_CBIAT_DELETE") <= 999999
 0 <= R("R_ACBIPGT") <= 999999
 0 <= R("R_ADOCBLS") <= 999999
 -999999 <= R("R_CBLAT_DELETE") <= 999999
 0 <= R("R_ACCOAL") <= 999999
 -999999 <= R("R_ADK1") <= 999999
 -999999 <= R("R_ADK3") <= 999999
 -999999 <= R("R_ADK4") <= 999999
 0 <= R("R_ADSK") <= 999999
 0 <= R("R_ALAabc") <= 999999
 -999999 <= R("R_ALAALAr") <= 999999
 0 <= R("R_AP4AH") <= 999999
 0 <= R("R_AP5AH") <= 999999
 0 <= R("R_ARBabc") <= 999999
 0 <= R("R_ARGabc") <= 999999
 0 <= R("R_ASNabc") <= 999999
 0 <= R("R_ASPabc") <= 999999
 0 <= R("R_ATPM") <= 999999
 -999999 <= R("R_ATPS4r") <= 999999
 0 <= R("R_CBL1abc") <= 999999
 0 <= R("R_CBMK") <= 999999
 0 <= R("R_CBPS") <= 999999
 0 <= R("R_CHLabc") <= 999999
 0 <= R("R_CTPS2") <= 999999
 0 <= R("R_CYSabc") <= 999999
 -999999 <= R("R_CYTK1") <= 999999
 -999999 <= R("R_CYTK2") <= 999999
 -999999 <= R("R_DADK") <= 999999
 -999999 <= R("R_DBTSr") <= 999999
 -999999 <= R("R_DGK1") <= 999999
 0 <= R("R_DHFS") <= 999999
 0 <= R("R_DPCOAK") <= 999999
 -999999 <= R("R_DTMPK") <= 999999
 0 <= R("R_DURIK1") <= 999999
 0 <= R("R_DXYLK") <= 999999
 0 <= R("R_FCLK") <= 999999
 0 <= R("R_FE2abc") <= 999999
 0 <= R("R_FRUK") <= 999999
 0 <= R("R_GALabc") <= 999999
 -999999 <= R("R_GALKr") <= 999999
 0 <= R("R_GART") <= 999999
 -999999 <= R("R_GK1") <= 999999
 0 <= R("R_GLCS1") <= 999999
 0 <= R("R_GLNabc") <= 999999
 0 <= R("R_GLNS") <= 999999
 0 <= R("R_GLU5K") <= 999999
 0 <= R("R_GLUabc") <= 999999
 0 <= R("R_GLUCYS") <= 999999
 0 <= R("R_GLYBabc") <= 999999
 0 <= R("R_GLYC3Pabc") <= 999999
 0 <= R("R_GLYK") <= 999999
 0 <= R("R_GMHEPK") <= 999999
 0 <= R("R_GSNK") <= 999999
 0 <= R("R_GSPMDS") <= 999999
 0 <= R("R_GTHS") <= 999999
 0 <= R("R_HEX1") <= 999999
 0 <= R("R_HISabc") <= 999999
 0 <= R("R_HSK") <= 999999
 0 <= R("R_ILEabc") <= 999999
 0 <= R("R_INSK") <= 999999
 0 <= R("R_Kabc") <= 999999
 0 <= R("R_LEUabc") <= 999999
 0 <= R("R_LPSSYN_EC") <= 999999
 0 <= R("R_LYSabc") <= 999999
 0 <= R("R_MALTabc") <= 999999
 0 <= R("R_MALTHXabc") <= 999999
 0 <= R("R_MALTPTabc") <= 999999
 0 <= R("R_MALTTRabc") <= 999999
 0 <= R("R_MALTTTRabc") <= 999999
 0 <= R("R_METabc") <= 999999
 0 <= R("R_METDabc") <= 999999
 0 <= R("R_NADK") <= 999999
 0 <= R("R_NAMNPP") <= 999999
 -999999 <= R("R_NDPK1") <= 999999
 -999999 <= R("R_NDPK2") <= 999999
 -999999 <= R("R_NDPK3") <= 999999
 -999999 <= R("R_NDPK4") <= 999999
 -999999 <= R("R_NDPK5") <= 999999
 -999999 <= R("R_NDPK6") <= 999999
 -999999 <= R("R_NDPK7") <= 999999
 -999999 <= R("R_NDPK8") <= 999999
 0 <= R("R_ORNabc") <= 999999
 0 <= R("R_PFK") <= 999999
 0 <= R("R_PFK_2") <= 999999
 0 <= R("R_PIabc") <= 999999
 -999999 <= R("R_PPAKr") <= 999999
 0 <= R("R_PPCK") <= 999999
 -999999 <= R("R_PRAGSr") <= 999999
 0 <= R("R_PRAIS") <= 999999
 0 <= R("R_PRFGS") <= 999999
 0 <= R("R_PROabc") <= 999999
 0 <= R("R_PTRCabc") <= 999999
 0 <= R("R_PYDAMK") <= 999999
 0 <= R("R_PYDXK") <= 999999
 0 <= R("R_PYDXNK") <= 999999
 0 <= R("R_PYK") <= 999999
 0 <= R("R_RBFK") <= 999999
 0 <= R("R_RBK") <= 999999
 0 <= R("R_RBK_L1") <= 999999
 0 <= R("R_RIBabc") <= 999999
 0 <= R("R_RMK") <= 999999
 0 <= R("R_RNDR1") <= 999999
 0 <= R("R_SHKK") <= 999999
 0 <= R("R_SPMDabc") <= 999999
 0 <= R("R_SUCCabc") <= 999999
 -999999 <= R("R_SUCOAS") <= 999999
 0 <= R("R_SULabc") <= 999999
 0 <= R("R_TAURabc") <= 999999
 0 <= R("R_TDSK") <= 999999
 0 <= R("R_THMabc") <= 999999
 0 <= R("R_THRabc") <= 999999
 0 <= R("R_TMDK1") <= 999999
 -999999 <= R("R_TMKr") <= 999999
 -999999 <= R("R_TMPKr") <= 999999
 0 <= R("R_TSULabc") <= 999999
 0 <= R("R_UAMAGS") <= 999999
 0 <= R("R_UAMAS") <= 999999
 0 <= R("R_UGMDDS") <= 999999
 -999999 <= R("R_UMPK") <= 999999
 -999999 <= R("R_URIDK2r") <= 999999
 0 <= R("R_VALabc") <= 999999
 0 <= R("R_XYLabc") <= 999999
 0 <= R("R_XYLK") <= 999999
 0 <= R("R_GLGC") <= 999999
 0 <= R("R_AGMHE") <= 999999
 0 <= R("R_GMHEPAT") <= 999999
 0 <= R("R_AGMT") <= 999999
 0 <= R("R_ARGDC") <= 999999
 0 <= R("R_LPLIPA3") <= 999999
 0 <= R("R_LPLIPA5") <= 999999
 0 <= R("R_PLIPA3") <= 999999
 0 <= R("R_LPLIPA2") <= 999999
 0 <= R("R_LPLIPA4") <= 999999
 0 <= R("R_PLIPA2") <= 999999
 0 <= R("R_LPLIPA1") <= 999999
 0 <= R("R_LPLIPA6") <= 999999
 0 <= R("R_PLIPA1") <= 999999
 0 <= R("R_UPP3MT") <= 999999
 0 <= R("R_DNTPPA") <= 999999
 0 <= R("R_GTPCI") <= 999999
 0 <= R("R_IG3PS") <= 999999
 -100000 <= R("Makgb_Transp") <= 100000
 -999999 <= R("R_EX_akg_e_") <= 999999
 0 <= R("R_AKGDH") <= 999999
 -999999 <= R("R_AKGt2r") <= 999999
 -999999 <= R("R_ALATA_L") <= 999999
 -999999 <= R("R_ASPTA") <= 999999
 -999999 <= R("R_GLUDy") <= 999999
 0 <= R("R_GLUSy") <= 999999
 0 <= R("R_HSTPT") <= 999999
 -999999 <= R("R_ICDHyr") <= 999999
 -999999 <= R("R_OHPBAT") <= 999999
 0 <= R("R_ORNTA") <= 999999
 0 <= R("R_OXGDC2") <= 999999
 -999999 <= R("R_PHETA1") <= 999999
 -999999 <= R("R_SDPTA") <= 999999
 0 <= R("R_SOTA") <= 999999
 0 <= R("R_TDPAGTA") <= 999999
 0 <= R("R_ASP1DC") <= 999999
 0 <= R("R_PANTS") <= 999999
 0 <= R("R_ACLS") <= 999999
 -100000 <= R("MalaDb_Transp") <= 100000
 -999999 <= R("R_EX_ala_D_e_") <= 999999
 -999999 <= R("R_ALAR") <= 999999
 0 <= R("R_ALARi") <= 999999
 0 <= R("R_ALATA_D2") <= 999999
 0 <= R("R_DAAD") <= 999999
 -999999 <= R("R_DALAt2r") <= 999999
 -100000 <= R("MalaLb_Transp") <= 100000
 -999999 <= R("R_EX_ala_L_e_") <= 999999
 -999999 <= R("R_ALAt2r") <= 999999
 0 <= R("R_ALATA_L2") <= 999999
 -999999 <= R("R_BTS2") <= 999999
 -100000 <= R("Malltnb_Transp") <= 100000
 -999999 <= R("R_EX_alltn_e_") <= 999999
 0 <= R("R_ALLTN") <= 999999
 -999999 <= R("R_ALLTNt2r") <= 999999
 0 <= R("R_ALLTAH") <= 999999
 -999999 <= R("R_TAGURr") <= 999999
 -999999 <= R("R_ADMDCr") <= 999999
 0 <= R("R_METAT") <= 999999
 -100000 <= R("Mampb_Transp") <= 100000
 -999999 <= R("R_EX_amp_e_") <= 999999
 -999999 <= R("R_ADSL1r") <= 999999
 0 <= R("R_ARGSS") <= 999999
 0 <= R("R_ASNS1") <= 999999
 0 <= R("R_ASNS2") <= 999999
 0 <= R("R_BPNT") <= 999999
 0 <= R("R_GLUTRS") <= 999999
 0 <= R("R_GMPS2") <= 999999
 0 <= R("R_NADDP") <= 999999
 0 <= R("R_NADS1") <= 999999
 0 <= R("R_NTPP6") <= 999999
 0 <= R("R_PACCOAL") <= 999999
 0 <= R("R_PPS") <= 999999
 -999999 <= R("R_PRPPS") <= 999999
 0 <= R("R_SELNPS") <= 999999
 0 <= R("R_SUCBZL") <= 999999
 0 <= R("R_NADDPe") <= 999999
 0 <= R("R_ANPRT") <= 999999
 0 <= R("R_ANS") <= 999999
 0 <= R("R_SADT2") <= 999999
 -999999 <= R("R_A5PISO") <= 999999
 0 <= R("R_KDOPS") <= 999999
 -100000 <= R("MarabLb_Transp") <= 100000
 -999999 <= R("R_EX_arab_L_e_") <= 999999
 -999999 <= R("R_ARAI") <= 999999
 -999999 <= R("R_ARBt2r") <= 999999
 0 <= R("R_AB6PGH") <= 999999
 -100000 <= R("MargLb_Transp") <= 100000
 -999999 <= R("R_EX_arg_L_e_") <= 999999
 -999999 <= R("R_ARGORNt7") <= 999999
 -999999 <= R("R_ARGSL") <= 999999
 0 <= R("R_AST") <= 999999
 -100000 <= R("MasnLb_Transp") <= 100000
 -999999 <= R("R_EX_asn_L_e_") <= 999999
 0 <= R("R_ASNN") <= 999999
 -999999 <= R("R_ASNt2r") <= 999999
 -100000 <= R("MaspLb_Transp") <= 100000
 -999999 <= R("R_EX_asp_L_e_") <= 999999
 0 <= R("R_ADSS") <= 999999
 0 <= R("R_ASPCT") <= 999999
 0 <= R("R_ASPO3") <= 999999
 0 <= R("R_ASPO4") <= 999999
 0 <= R("R_ASPO5") <= 999999
 0 <= R("R_ASPO6") <= 999999
 0 <= R("R_ASPT") <= 999999
 0 <= R("R_ASPt2") <= 999999
 0 <= R("R_ASPt2_2") <= 999999
 0 <= R("R_ASPt2_3") <= 999999
 -999999 <= R("R_HSDy") <= 999999
 0 <= R("R_ADNCYC") <= 999999
 0 <= R("R_ATPPRT") <= 999999
 0 <= R("R_FMNAT") <= 999999
 0 <= R("R_NMNAT") <= 999999
 0 <= R("R_NNAT") <= 999999
 0 <= R("R_PTPATi") <= 999999
 0 <= R("R_RNTR1") <= 999999
 -999999 <= R("R_SERASr") <= 999999
 -999999 <= R("R_CRNBTCT") <= 999999
 0 <= R("R_BETALDHx") <= 999999
 0 <= R("R_BETALDHy") <= 999999
 0 <= R("R_BSORx") <= 999999
 0 <= R("R_BSORy") <= 999999
 -100000 <= R("Mbutb_Transp") <= 100000
 -999999 <= R("R_EX_but_e_") <= 999999
 -999999 <= R("R_BUTt2r") <= 999999
 -999999 <= R("R_DHORTS") <= 999999
 -100000 <= R("Mcbl1b_Transp") <= 100000
 -999999 <= R("R_EX_cbl1_e_") <= 999999
 -999999 <= R("R_OCBT") <= 999999
 0 <= R("R_RNDR3") <= 999999
 0 <= R("R_CDAPPA_EC") <= 999999
 -999999 <= R("R_DASYN_EC") <= 999999
 -999999 <= R("R_PGSA_EC") <= 999999
 -999999 <= R("R_PSSA_EC") <= 999999
 0 <= R("R_DHPPD") <= 999999
 0 <= R("R_PPPNDO") <= 999999
 0 <= R("R_CINNDO") <= 999999
 0 <= R("R_DHCIND") <= 999999
 -100000 <= R("Mcholb_Transp") <= 100000
 -999999 <= R("R_EX_chol_e_") <= 999999
 -999999 <= R("R_CHLt2r") <= 999999
 0 <= R("R_GPDDA1") <= 999999
 0 <= R("R_CHORM") <= 999999
 0 <= R("R_ICHORSi") <= 999999
 -100000 <= R("Mcitb_Transp") <= 100000
 -999999 <= R("R_EX_cit_e_") <= 999999
 -999999 <= R("R_ACONT") <= 999999
 0 <= R("R_CITt7") <= 999999
 0 <= R("R_KDOCT2") <= 999999
 0 <= R("R_MOAT") <= 999999
 0 <= R("R_MOAT2") <= 999999
 -999999 <= R("R_CLPNS_EC") <= 999999
 0 <= R("R_CMPN") <= 999999
 0 <= R("R_CYTDK2") <= 999999
 0 <= R("R_NTD4") <= 999999
 0 <= R("R_NTPP4") <= 999999
 -100000 <= R("Mco2b_Transp") <= 100000
 -999999 <= R("R_EX_co2_e_") <= 999999
 -999999 <= R("R_CO2t") <= 999999
 0 <= R("R_CPPPGO") <= 999999
 0 <= R("R_CYNTAH") <= 999999
 0 <= R("R_DXPS") <= 999999
 0 <= R("R_FDH2") <= 999999
 0 <= R("R_FDH3") <= 999999
 0 <= R("R_FHL") <= 999999
 0 <= R("R_GLYCL") <= 999999
 -999999 <= R("R_HCO3E") <= 999999
 0 <= R("R_ME1") <= 999999
 0 <= R("R_ME2") <= 999999
 0 <= R("R_MMCD") <= 999999
 0 <= R("R_NNDPR") <= 999999
 0 <= R("R_OMPDC") <= 999999
 0 <= R("R_ORNDC") <= 999999
 0 <= R("R_PDX5PS") <= 999999
 0 <= R("R_PPC") <= 999999
 0 <= R("R_PPNDH") <= 999999
 0 <= R("R_PSD_EC") <= 999999
 0 <= R("R_SADH") <= 999999
 0 <= R("R_UGLYCH") <= 999999
 0 <= R("R_UPPDC1") <= 999999
 0 <= R("R_HSST") <= 999999
 0 <= R("R_NPHS") <= 999999
 0 <= R("R_PTA2") <= 999999
 0 <= R("R_THDPS") <= 999999
 -100000 <= R("Mcrnb_Transp") <= 100000
 -999999 <= R("R_EX_crn_e_") <= 999999
 -999999 <= R("R_CRNCBCT") <= 999999
 0 <= R("R_CRNt7") <= 999999
 -999999 <= R("R_CRNCDH") <= 999999
 -100000 <= R("Mcsnb_Transp") <= 100000
 -999999 <= R("R_EX_csn_e_") <= 999999
 0 <= R("R_CSND") <= 999999
 0 <= R("R_CSNt2") <= 999999
 0 <= R("R_RNTR3") <= 999999
 0 <= R("R_CYANST") <= 999999
 -100000 <= R("Mcyntb_Transp") <= 100000
 -999999 <= R("R_EX_cynt_e_") <= 999999
 0 <= R("R_CYNTt2") <= 999999
 -100000 <= R("McysLb_Transp") <= 100000
 -999999 <= R("R_EX_cys_L_e_") <= 999999
 0 <= R("R_CYSDS") <= 999999
 0 <= R("R_SHSL1") <= 999999
 0 <= R("R_CYSTL") <= 999999
 -100000 <= R("Mcytdb_Transp") <= 100000
 -999999 <= R("R_EX_cytd_e_") <= 999999
 0 <= R("R_CYTD") <= 999999
 0 <= R("R_CYTDt2") <= 999999
 -999999 <= R("R_CYTDt2r") <= 999999
 -100000 <= R("Mdad2b_Transp") <= 100000
 -999999 <= R("R_EX_dad_2_e_") <= 999999
 0 <= R("R_DADA") <= 999999
 0 <= R("R_DADNt2") <= 999999
 0 <= R("R_NTD6") <= 999999
 0 <= R("R_NTPP5") <= 999999
 0 <= R("R_DB4PS") <= 999999
 0 <= R("R_NTD3") <= 999999
 0 <= R("R_NTPP3") <= 999999
 0 <= R("R_DCTPD") <= 999999
 -100000 <= R("Mdcytb_Transp") <= 100000
 -999999 <= R("R_EX_dcyt_e_") <= 999999
 0 <= R("R_DCYTD") <= 999999
 0 <= R("R_DCYTt2") <= 999999
 0 <= R("R_RNDR2") <= 999999
 0 <= R("R_NTD8") <= 999999
 0 <= R("R_NTPP1") <= 999999
 -100000 <= R("Mdgsnb_Transp") <= 100000
 -999999 <= R("R_EX_dgsn_e_") <= 999999
 0 <= R("R_DGSNt2") <= 999999
 0 <= R("R_NTPTP1") <= 999999
 0 <= R("R_RNTR2") <= 999999
 -100000 <= R("Mdhab_Transp") <= 100000
 -999999 <= R("R_EX_dha_e_") <= 999999
 0 <= R("R_DHAPT") <= 999999
 -999999 <= R("R_DHAt") <= 999999
 -999999 <= R("R_F6PA") <= 999999
 0 <= R("R_GLYCDx") <= 999999
 -999999 <= R("R_FBA") <= 999999
 -999999 <= R("R_FCLPA") <= 999999
 -999999 <= R("R_G3PD2") <= 999999
 0 <= R("R_G3PD5") <= 999999
 0 <= R("R_G3PD6") <= 999999
 0 <= R("R_MGSA") <= 999999
 0 <= R("R_QULNS") <= 999999
 -999999 <= R("R_RMPA") <= 999999
 -999999 <= R("R_TGBPA") <= 999999
 -999999 <= R("R_TPI") <= 999999
 0 <= R("R_DHCINDO") <= 999999
 -999999 <= R("R_DHFR") <= 999999
 0 <= R("R_TMDS") <= 999999
 0 <= R("R_DNMPPA") <= 999999
 0 <= R("R_DHORD2") <= 999999
 0 <= R("R_DHORD5") <= 999999
 0 <= R("R_HPPPNDO") <= 999999
 0 <= R("R_DHPTDCs") <= 999999
 0 <= R("R_RHCCE") <= 999999
 -100000 <= R("Mdinb_Transp") <= 100000
 -999999 <= R("R_EX_din_e_") <= 999999
 0 <= R("R_DINSt2") <= 999999
 0 <= R("R_DMATT") <= 999999
 0 <= R("R_DMPPS") <= 999999
 0 <= R("R_IPDDIi") <= 999999
 -100000 <= R("Mdmsb_Transp") <= 100000
 -999999 <= R("R_EX_dms_e_") <= 999999
 0 <= R("R_DMSOR1") <= 999999
 0 <= R("R_DMSOR1e") <= 999999
 -100000 <= R("Mdmsob_Transp") <= 100000
 -999999 <= R("R_EX_dmso_e_") <= 999999
 0 <= R("R_SHCHD2") <= 999999
 0 <= R("R_AADDGT") <= 999999
 0 <= R("R_TDPDRE") <= 999999
 0 <= R("R_TDPGDH") <= 999999
 0 <= R("R_TDPDRR") <= 999999
 0 <= R("R_G1PTT") <= 999999
 0 <= R("R_NTD5") <= 999999
 0 <= R("R_NTPP7") <= 999999
 0 <= R("R_RNDR4") <= 999999
 0 <= R("R_DUTPDP") <= 999999
 0 <= R("R_NTD1") <= 999999
 -100000 <= R("Mdurib_Transp") <= 100000
 -999999 <= R("R_EX_duri_e_") <= 999999
 0 <= R("R_DURIt2") <= 999999
 0 <= R("R_RNTR4") <= 999999
 -999999 <= R("R_TALA") <= 999999
 -999999 <= R("R_TKT2") <= 999999
 0 <= R("R_ECAP_EC") <= 999999
 0 <= R("R_IGPDH") <= 999999
 0 <= R("R_GPDDA2") <= 999999
 -100000 <= R("Metohb_Transp") <= 100000
 -999999 <= R("R_EX_etoh_e_") <= 999999
 -999999 <= R("R_ETOHt2r") <= 999999
 0 <= R("R_FRUpts") <= 999999
 0 <= R("R_FBP") <= 999999
 0 <= R("R_FRUpts2") <= 999999
 0 <= R("R_G6PDA") <= 999999
 0 <= R("R_GF6PTA") <= 999999
 -999999 <= R("R_M1PD") <= 999999
 -999999 <= R("R_MAN6PI") <= 999999
 -999999 <= R("R_PGI") <= 999999
 -999999 <= R("R_SBTPD") <= 999999
 0 <= R("R_SUCD1i") <= 999999
 -999999 <= R("R_SUCD4") <= 999999
 -999999 <= R("R_FCI") <= 999999
 -100000 <= R("Mfe2b_Transp") <= 100000
 -999999 <= R("R_EX_fe2_e_") <= 999999
 0 <= R("R_FCLT") <= 999999
 0 <= R("R_SHCHF") <= 999999
 -100000 <= R("Mforb_Transp") <= 100000
 -999999 <= R("R_EX_for_e_") <= 999999
 -999999 <= R("R_FORt") <= 999999
 -999999 <= R("R_IMPC") <= 999999
 0 <= R("R_GRTT") <= 999999
 0 <= R("R_HEMEOS") <= 999999
 0 <= R("R_OCTDPS") <= 999999
 0 <= R("R_UDCPDPS") <= 999999
 -100000 <= R("Mfrub_Transp") <= 100000
 -999999 <= R("R_EX_fru_e_") <= 999999
 0 <= R("R_FFSD") <= 999999
 0 <= R("R_XYLI2i") <= 999999
 -999999 <= R("R_GUI1") <= 999999
 -999999 <= R("R_MANAO") <= 999999
 -100000 <= R("Mfuc1pLb_Transp") <= 100000
 -999999 <= R("R_EX_fuc1p_L_e_") <= 999999
 0 <= R("R_FUCPt6_2") <= 999999
 -100000 <= R("MfucLb_Transp") <= 100000
 -999999 <= R("R_EX_fuc_L_e_") <= 999999
 -999999 <= R("R_FUCt") <= 999999
 -100000 <= R("Mfumb_Transp") <= 100000
 -999999 <= R("R_EX_fum_e_") <= 999999
 0 <= R("R_FRD2") <= 999999
 -999999 <= R("R_FUM") <= 999999
 0 <= R("R_FUMt2_2") <= 999999
 0 <= R("R_FUMt2_3") <= 999999
 0 <= R("R_HKNTDH") <= 999999
 -999999 <= R("R_SUCFUMt") <= 999999
 0 <= R("R_G1PP") <= 999999
 -999999 <= R("R_GALU") <= 999999
 0 <= R("R_GALUi") <= 999999
 0 <= R("R_GLCP") <= 999999
 -999999 <= R("R_MLTP1") <= 999999
 -999999 <= R("R_MLTP2") <= 999999
 -999999 <= R("R_MLTP3") <= 999999
 -999999 <= R("R_PGMT") <= 999999
 -999999 <= R("R_UGLT") <= 999999
 -999999 <= R("R_TKT1") <= 999999
 0 <= R("R_GPDDA4") <= 999999
 0 <= R("R_GPDDA5") <= 999999
 0 <= R("R_GPDDA3") <= 999999
 -100000 <= R("Mg6pb_Transp") <= 100000
 -999999 <= R("R_EX_g6p_e_") <= 999999
 0 <= R("R_G6Pt6_2") <= 999999
 0 <= R("R_GLCpts") <= 999999
 0 <= R("R_TRE6PH") <= 999999
 0 <= R("R_TRE6PS") <= 999999
 -100000 <= R("Mgalb_Transp") <= 100000
 -999999 <= R("R_EX_gal_e_") <= 999999
 0 <= R("R_GALS3") <= 999999
 0 <= R("R_GALt2") <= 999999
 0 <= R("R_LACZ") <= 999999
 -100000 <= R("MgalctDb_Transp") <= 100000
 -999999 <= R("R_EX_galct_D_e_") <= 999999
 -999999 <= R("R_GALCTt2r") <= 999999
 -100000 <= R("MgalctnDb_Transp") <= 100000
 -999999 <= R("R_EX_galctn_D_e_") <= 999999
 -999999 <= R("R_GALCTNt2r") <= 999999
 0 <= R("R_GALTpts") <= 999999
 -999999 <= R("R_GLTPD") <= 999999
 -100000 <= R("Mgaltb_Transp") <= 100000
 -999999 <= R("R_EX_galt_e_") <= 999999
 -100000 <= R("Mgalurb_Transp") <= 100000
 -999999 <= R("R_EX_galur_e_") <= 999999
 -999999 <= R("R_GALURt2r") <= 999999
 -999999 <= R("R_GUI2") <= 999999
 -999999 <= R("R_PGAMT") <= 999999
 0 <= R("R_GAMpts") <= 999999
 -100000 <= R("Mgamb_Transp") <= 100000
 -999999 <= R("R_EX_gam_e_") <= 999999
 -100000 <= R("Mgbbtnb_Transp") <= 100000
 -999999 <= R("R_EX_gbbtn_e_") <= 999999
 0 <= R("R_GCALDD") <= 999999
 0 <= R("R_GP4GH") <= 999999
 0 <= R("R_MAN1PT2") <= 999999
 0 <= R("R_URIK2") <= 999999
 0 <= R("R_GDMANE") <= 999999
 0 <= R("R_GMAND") <= 999999
 0 <= R("R_GOFUCR") <= 999999
 -100000 <= R("MglcDb_Transp") <= 100000
 -999999 <= R("R_EX_glc_e_") <= 999999
 0 <= R("R_AMALT1") <= 999999
 0 <= R("R_AMALT2") <= 999999
 0 <= R("R_AMALT3") <= 999999
 0 <= R("R_AMALT4") <= 999999
 0 <= R("R_GLCt2") <= 999999
 0 <= R("R_MLTG1") <= 999999
 0 <= R("R_MLTG2") <= 999999
 0 <= R("R_MLTG3") <= 999999
 0 <= R("R_MLTG4") <= 999999
 0 <= R("R_MLTG5") <= 999999
 0 <= R("R_TREH") <= 999999
 0 <= R("R_GLCDe") <= 999999
 0 <= R("R_TREHe") <= 999999
 -100000 <= R("Mglcnb_Transp") <= 100000
 -999999 <= R("R_EX_glcn_e_") <= 999999
 -999999 <= R("R_GLCNt2r") <= 999999
 -100000 <= R("Mglcrb_Transp") <= 100000
 -999999 <= R("R_EX_glcr_e_") <= 999999
 -999999 <= R("R_GLCRt2r") <= 999999
 -100000 <= R("Mglcurb_Transp") <= 100000
 -999999 <= R("R_EX_glcur_e_") <= 999999
 -999999 <= R("R_GLCURt2r") <= 999999
 -100000 <= R("MglnLb_Transp") <= 100000
 -999999 <= R("R_EX_gln_L_e_") <= 999999
 0 <= R("R_GLUN") <= 999999
 0 <= R("R_GLUPRT") <= 999999
 0 <= R("R_GLUTRR") <= 999999
 0 <= R("R_G5SD") <= 999999
 -999999 <= R("R_GLUR") <= 999999
 -100000 <= R("MgluLb_Transp") <= 100000
 -999999 <= R("R_EX_glu_L_e_") <= 999999
 -999999 <= R("R_GLUt2r") <= 999999
 0 <= R("R_GLUt4") <= 999999
 0 <= R("R_SGDS") <= 999999
 0 <= R("R_GLYCLTDx") <= 999999
 0 <= R("R_GLYCLTDy") <= 999999
 0 <= R("R_GLYCTO2") <= 999999
 0 <= R("R_GLYCTO3") <= 999999
 0 <= R("R_ICL") <= 999999
 -100000 <= R("Mglyaldb_Transp") <= 100000
 -999999 <= R("R_EX_glyald_e_") <= 999999
 -999999 <= R("R_ALCD19") <= 999999
 -999999 <= R("R_GLYALDt") <= 999999
 -100000 <= R("Mglyb_Transp") <= 100000
 -999999 <= R("R_EX_gly_e_") <= 999999
 -100000 <= R("Mglybb_Transp") <= 100000
 -999999 <= R("R_EX_glyb_e_") <= 999999
 -999999 <= R("R_GLYBt2r") <= 999999
 0 <= R("R_GHMT2") <= 999999
 -999999 <= R("R_GLYt2r") <= 999999
 -100000 <= R("Mglyc3pb_Transp") <= 100000
 -999999 <= R("R_EX_glyc3p_e_") <= 999999
 0 <= R("R_GLYC3Pt6") <= 999999
 -100000 <= R("Mglycb_Transp") <= 100000
 -999999 <= R("R_EX_glyc_e_") <= 999999
 -999999 <= R("R_GLYCt") <= 999999
 -100000 <= R("Mglycltb_Transp") <= 100000
 -999999 <= R("R_EX_glyclt_e_") <= 999999
 -999999 <= R("R_GLYCLTt2r") <= 999999
 0 <= R("R_HPYRRx") <= 999999
 0 <= R("R_HPYRRy") <= 999999
 0 <= R("R_GMHEPPA") <= 999999
 0 <= R("R_S7PI") <= 999999
 0 <= R("R_GMPR") <= 999999
 0 <= R("R_GUAPRT") <= 999999
 0 <= R("R_NTD9") <= 999999
 0 <= R("R_NTPP2") <= 999999
 -100000 <= R("Mgsnb_Transp") <= 100000
 -999999 <= R("R_EX_gsn_e_") <= 999999
 0 <= R("R_GSNt2") <= 999999
 0 <= R("R_NTPTP2") <= 999999
 -999999 <= R("R_PUNP3") <= 999999
 -999999 <= R("R_GTHOr") <= 999999
 0 <= R("R_GLYOX") <= 999999
 0 <= R("R_GSPMDA") <= 999999
 0 <= R("R_LGTHL") <= 999999
 -100000 <= R("Mguab_Transp") <= 100000
 -999999 <= R("R_EX_gua_e_") <= 999999
 0 <= R("R_GUAD") <= 999999
 -999999 <= R("R_GUAt") <= 999999
 0 <= R("R_GUAt2") <= 999999
 0 <= R("R_HYD1") <= 999999
 0 <= R("R_HYD2") <= 999999
 0 <= R("R_IPDPS") <= 999999
 0 <= R("R_CAT") <= 999999
 -999999 <= R("R_PDX5PO") <= 999999
 0 <= R("R_PEAMNO") <= 999999
 0 <= R("R_PYAM5PO") <= 999999
 0 <= R("R_SPODM") <= 999999
 -100000 <= R("Mh2ob_Transp") <= 100000
 -999999 <= R("R_EX_h2o_e_") <= 999999
 0 <= R("R_ALDD19x") <= 999999
 0 <= R("R_CYTBD") <= 999999
 0 <= R("R_CYTBO3") <= 999999
 -999999 <= R("R_H2Ot") <= 999999
 0 <= R("R_HISTD") <= 999999
 0 <= R("R_HISTP") <= 999999
 0 <= R("R_HKNDDH") <= 999999
 0 <= R("R_HMBS") <= 999999
 0 <= R("R_HYPOE") <= 999999
 0 <= R("R_IMPD") <= 999999
 0 <= R("R_KDOPP") <= 999999
 -999999 <= R("R_LCAD") <= 999999
 0 <= R("R_LCADi") <= 999999
 0 <= R("R_MI1PP") <= 999999
 0 <= R("R_NADPPPS") <= 999999
 0 <= R("R_NMNDA") <= 999999
 0 <= R("R_NMNN") <= 999999
 0 <= R("R_NMNt7") <= 999999
 0 <= R("R_NNAM") <= 999999
 0 <= R("R_NO3R1") <= 999999
 0 <= R("R_NO3R2") <= 999999
 0 <= R("R_NTD10") <= 999999
 0 <= R("R_NTD11") <= 999999
 0 <= R("R_NTD2") <= 999999
 0 <= R("R_NTPP8") <= 999999
 0 <= R("R_NTRIR2x") <= 999999
 0 <= R("R_PDXPP") <= 999999
 0 <= R("R_PGPP_EC") <= 999999
 0 <= R("R_PPA") <= 999999
 0 <= R("R_PPPGO") <= 999999
 0 <= R("R_PRAMPC") <= 999999
 0 <= R("R_PRATPP") <= 999999
 0 <= R("R_PSP_L") <= 999999
 0 <= R("R_PYDXPP") <= 999999
 0 <= R("R_SGSAD") <= 999999
 0 <= R("R_SSALx") <= 999999
 0 <= R("R_SSALy") <= 999999
 -999999 <= R("R_SULR") <= 999999
 0 <= R("R_TARTD") <= 999999
 0 <= R("R_THRS") <= 999999
 0 <= R("R_TMAOR1") <= 999999
 0 <= R("R_TRE6PP") <= 999999
 -999999 <= R("R_TRPAS2") <= 999999
 0 <= R("R_TRPS2") <= 999999
 0 <= R("R_UACMAMO") <= 999999
 0 <= R("R_UDCPDP") <= 999999
 0 <= R("R_UDPGD") <= 999999
 0 <= R("R_UPP3S") <= 999999
 0 <= R("R_USHD") <= 999999
 0 <= R("R_TMAOR1e") <= 999999
 -100000 <= R("Mhb_Transp") <= 100000
 -999999 <= R("R_EX_h_e_") <= 999999
 0 <= R("R_ACMAMUT") <= 999999
 -999999 <= R("R_D_LACt2") <= 999999
 -999999 <= R("R_DSERt2r") <= 999999
 0 <= R("R_HDCAt2") <= 999999
 -999999 <= R("R_HISt2r") <= 999999
 -999999 <= R("R_IDONt2r") <= 999999
 -999999 <= R("R_ILEt2r") <= 999999
 -999999 <= R("R_INDOLEt2r") <= 999999
 0 <= R("R_INSt2") <= 999999
 -999999 <= R("R_INSt2r") <= 999999
 -999999 <= R("R_Kt2r") <= 999999
 -999999 <= R("R_LCTSt") <= 999999
 -999999 <= R("R_LDH_D") <= 999999
 -999999 <= R("R_LEUt2r") <= 999999
 -999999 <= R("R_L_LACt2r") <= 999999
 0 <= R("R_LPADSS") <= 999999
 -999999 <= R("R_LYSt2r") <= 999999
 0 <= R("R_MALt2_2") <= 999999
 0 <= R("R_MALt2_3") <= 999999
 -999999 <= R("R_MDH") <= 999999
 0 <= R("R_MELIBt2") <= 999999
 -999999 <= R("R_MTHFD") <= 999999
 0 <= R("R_NADH10") <= 999999
 0 <= R("R_NADH5") <= 999999
 0 <= R("R_NADH6") <= 999999
 0 <= R("R_NADH7") <= 999999
 -999999 <= R("R_NAt3_1") <= 999999
 0 <= R("R_NAt3_15") <= 999999
 0 <= R("R_NAt3_2") <= 999999
 -999999 <= R("R_NO2t2r") <= 999999
 0 <= R("R_OCDCAt2") <= 999999
 0 <= R("R_PAPSR") <= 999999
 -999999 <= R("R_PHEt2r") <= 999999
 -999999 <= R("R_PIt2r") <= 999999
 -999999 <= R("R_PPPNt2r") <= 999999
 0 <= R("R_PPTGS") <= 999999
 -999999 <= R("R_PROt2r") <= 999999
 -999999 <= R("R_PTRCt2r") <= 999999
 -999999 <= R("R_PYRt2r") <= 999999
 0 <= R("R_RMNt") <= 999999
 -999999 <= R("R_SERt2r") <= 999999
 0 <= R("R_SUCCt2_2") <= 999999
 0 <= R("R_SUCCt2_3") <= 999999
 0 <= R("R_SUCCt2b") <= 999999
 0 <= R("R_THD2") <= 999999
 0 <= R("R_THMDt2") <= 999999
 -999999 <= R("R_THMDt2r") <= 999999
 -999999 <= R("R_THRt2r") <= 999999
 0 <= R("R_TRDR") <= 999999
 -999999 <= R("R_TRPt2r") <= 999999
 0 <= R("R_TTDCAt2") <= 999999
 -999999 <= R("R_TYRt2r") <= 999999
 0 <= R("R_UAGPT3") <= 999999
 0 <= R("R_UAPGR") <= 999999
 0 <= R("R_URAt2") <= 999999
 -999999 <= R("R_URAt2r") <= 999999
 0 <= R("R_URIt2") <= 999999
 -999999 <= R("R_URIt2r") <= 999999
 -999999 <= R("R_VALt2r") <= 999999
 0 <= R("R_XANt2") <= 999999
 -999999 <= R("R_XTSNt2r") <= 999999
 0 <= R("R_XYLt2") <= 999999
 -100000 <= R("Mhdcab_Transp") <= 100000
 -999999 <= R("R_EX_hdca_e_") <= 999999
 -100000 <= R("MhisLb_Transp") <= 100000
 -999999 <= R("R_EX_his_L_e_") <= 999999
 -100000 <= R("Mhxanb_Transp") <= 100000
 -999999 <= R("R_EX_hxan_e_") <= 999999
 0 <= R("R_HXPRT") <= 999999
 -999999 <= R("R_HYXNt") <= 999999
 -999999 <= R("R_PUNP5") <= 999999
 -100000 <= R("MidonLb_Transp") <= 100000
 -999999 <= R("R_EX_idon_L_e_") <= 999999
 -100000 <= R("MileLb_Transp") <= 100000
 -999999 <= R("R_EX_ile_L_e_") <= 999999
 -100000 <= R("Mindoleb_Transp") <= 100000
 -999999 <= R("R_EX_indole_e_") <= 999999
 -100000 <= R("Minsb_Transp") <= 100000
 -999999 <= R("R_EX_ins_e_") <= 999999
 -100000 <= R("Mkb_Transp") <= 100000
 -999999 <= R("R_EX_k_e_") <= 999999
 -100000 <= R("MlacDb_Transp") <= 100000
 -999999 <= R("R_EX_lac_D_e_") <= 999999
 0 <= R("R_LDH_D2") <= 999999
 -100000 <= R("MlacLb_Transp") <= 100000
 -999999 <= R("R_EX_lac_L_e_") <= 999999
 0 <= R("R_L_LACD2") <= 999999
 0 <= R("R_L_LACD3") <= 999999
 -100000 <= R("Mlctsb_Transp") <= 100000
 -999999 <= R("R_EX_lcts_e_") <= 999999
 -100000 <= R("MleuLb_Transp") <= 100000
 -999999 <= R("R_EX_leu_L_e_") <= 999999
 -100000 <= R("MlysLb_Transp") <= 100000
 -999999 <= R("R_EX_lys_L_e_") <= 999999
 -100000 <= R("MmalLb_Transp") <= 100000
 -999999 <= R("R_EX_mal_L_e_") <= 999999
 0 <= R("R_MDH2") <= 999999
 0 <= R("R_MDH3") <= 999999
 0 <= R("R_MALTpts") <= 999999
 -100000 <= R("Mmaltb_Transp") <= 100000
 -999999 <= R("R_EX_malt_e_") <= 999999
 -100000 <= R("Mmalthxb_Transp") <= 100000
 -999999 <= R("R_EX_malthx_e_") <= 999999
 -100000 <= R("Mmaltptb_Transp") <= 100000
 -999999 <= R("R_EX_maltpt_e_") <= 999999
 -100000 <= R("Mmalttrb_Transp") <= 100000
 -999999 <= R("R_EX_malttr_e_") <= 999999
 -100000 <= R("Mmaltttrb_Transp") <= 100000
 -999999 <= R("R_EX_maltttr_e_") <= 999999
 -999999 <= R("R_PMANM") <= 999999
 -100000 <= R("Mman6pb_Transp") <= 100000
 -999999 <= R("R_EX_man6p_e_") <= 999999
 0 <= R("R_MAN6Pt6_2") <= 999999
 0 <= R("R_MANpts") <= 999999
 -100000 <= R("Mmanb_Transp") <= 100000
 -999999 <= R("R_EX_man_e_") <= 999999
 -100000 <= R("Mmelibb_Transp") <= 100000
 -999999 <= R("R_EX_melib_e_") <= 999999
 -100000 <= R("MmetDb_Transp") <= 100000
 -999999 <= R("R_EX_met_D_e_") <= 999999
 -100000 <= R("MmetLb_Transp") <= 100000
 -999999 <= R("R_EX_met_L_e_") <= 999999
 -999999 <= R("R_MCITL2") <= 999999
 -999999 <= R("R_MME") <= 999999
 0 <= R("R_MMM2") <= 999999
 0 <= R("R_MNLpts") <= 999999
 -100000 <= R("Mmnlb_Transp") <= 100000
 -999999 <= R("R_EX_mnl_e_") <= 999999
 -100000 <= R("Mna1b_Transp") <= 100000
 -999999 <= R("R_EX_na1_e_") <= 999999
 0 <= R("R_PNTOt4") <= 999999
 0 <= R("R_PROt4") <= 999999
 0 <= R("R_SERt4") <= 999999
 0 <= R("R_THRt4") <= 999999
 -100000 <= R("Mnacb_Transp") <= 100000
 -999999 <= R("R_EX_nac_e_") <= 999999
 0 <= R("R_NACUP") <= 999999
 -100000 <= R("Mnadb_Transp") <= 100000
 -999999 <= R("R_EX_nad_e_") <= 999999
 0 <= R("R_NADTRHD") <= 999999
 -100000 <= R("Mnh4b_Transp") <= 100000
 -999999 <= R("R_EX_nh4_e_") <= 999999
 -999999 <= R("R_NH4t") <= 999999
 0 <= R("R_SERD_D") <= 999999
 0 <= R("R_SERD_L") <= 999999
 -100000 <= R("Mnmnb_Transp") <= 100000
 -999999 <= R("R_EX_nmn_e_") <= 999999
 0 <= R("R_NMNP") <= 999999
 -100000 <= R("Mno2b_Transp") <= 100000
 -999999 <= R("R_EX_no2_e_") <= 999999
 0 <= R("R_NO3t7") <= 999999
 -100000 <= R("Mno3b_Transp") <= 100000
 -999999 <= R("R_EX_no3_e_") <= 999999
 -100000 <= R("Mo2b_Transp") <= 100000
 -999999 <= R("R_EX_o2_e_") <= 999999
 -999999 <= R("R_O2t") <= 999999
 -100000 <= R("Mocdcab_Transp") <= 100000
 -999999 <= R("R_EX_ocdca_e_") <= 999999
 -100000 <= R("Mornb_Transp") <= 100000
 -999999 <= R("R_EX_orn_e_") <= 999999
 -999999 <= R("R_PTRCORNt7") <= 999999
 -999999 <= R("R_ORPT") <= 999999
 0 <= R("R_SBTpts") <= 999999
 0 <= R("R_SUCpts") <= 999999
 0 <= R("R_TREpts") <= 999999
 0 <= R("R_UAGCVT") <= 999999
 -100000 <= R("MpheLb_Transp") <= 100000
 -999999 <= R("R_EX_phe_L_e_") <= 999999
 -100000 <= R("Mpib_Transp") <= 100000
 -999999 <= R("R_EX_pi_e_") <= 999999
 -999999 <= R("R_PUNP7") <= 999999
 -999999 <= R("R_PYNP2r") <= 999999
 -100000 <= R("MpntoRb_Transp") <= 100000
 -999999 <= R("R_EX_pnto_R_e_") <= 999999
 0 <= R("R_PPCSCT") <= 999999
 0 <= R("R_UPPRT") <= 999999
 0 <= R("R_XPPT") <= 999999
 -100000 <= R("Mpppnb_Transp") <= 100000
 -999999 <= R("R_EX_pppn_e_") <= 999999
 0 <= R("R_PRMICIi") <= 999999
 -100000 <= R("MproLb_Transp") <= 100000
 -999999 <= R("R_EX_pro_L_e_") <= 999999
 -100000 <= R("Mptrcb_Transp") <= 100000
 -999999 <= R("R_EX_ptrc_e_") <= 999999
 -100000 <= R("Mpyrb_Transp") <= 100000
 -999999 <= R("R_EX_pyr_e_") <= 999999
 -999999 <= R("R_PPM") <= 999999
 -999999 <= R("R_RPI") <= 999999
 -100000 <= R("MribDb_Transp") <= 100000
 -999999 <= R("R_EX_rib_D_e_") <= 999999
 -999999 <= R("R_RMI") <= 999999
 -100000 <= R("Mrmnb_Transp") <= 100000
 -999999 <= R("R_EX_rmn_e_") <= 999999
 -999999 <= R("R_RPE") <= 999999
 -999999 <= R("R_RBP4E") <= 999999
 0 <= R("R_X5PL3E") <= 999999
 -100000 <= R("MsbtDb_Transp") <= 100000
 -999999 <= R("R_EX_sbt_D_e_") <= 999999
 -100000 <= R("MserDb_Transp") <= 100000
 -999999 <= R("R_EX_ser_D_e_") <= 999999
 -100000 <= R("MserLb_Transp") <= 100000
 -999999 <= R("R_EX_ser_L_e_") <= 999999
 -100000 <= R("Mso4b_Transp") <= 100000
 -999999 <= R("R_EX_so4_e_") <= 999999
 -100000 <= R("Mspmdb_Transp") <= 100000
 -999999 <= R("R_EX_spmd_e_") <= 999999
 -100000 <= R("Msuccb_Transp") <= 100000
 -999999 <= R("R_EX_succ_e_") <= 999999
 -999999 <= R("R_TARTRt7") <= 999999
 -100000 <= R("Msucrb_Transp") <= 100000
 -999999 <= R("R_EX_sucr_e_") <= 999999
 -100000 <= R("MtartrLb_Transp") <= 100000
 -999999 <= R("R_EX_tartr_L_e_") <= 999999
 -100000 <= R("Mtaurb_Transp") <= 100000
 -999999 <= R("R_EX_taur_e_") <= 999999
 -100000 <= R("Mthmb_Transp") <= 100000
 -999999 <= R("R_EX_thm_e_") <= 999999
 -100000 <= R("MthrLb_Transp") <= 100000
 -999999 <= R("R_EX_thr_L_e_") <= 999999
 -100000 <= R("Mthymdb_Transp") <= 100000
 -999999 <= R("R_EX_thymd_e_") <= 999999
 -100000 <= R("Mtmab_Transp") <= 100000
 -999999 <= R("R_EX_tma_e_") <= 999999
 -100000 <= R("Mtmaob_Transp") <= 100000
 -999999 <= R("R_EX_tmao_e_") <= 999999
 -100000 <= R("Mtreb_Transp") <= 100000
 -999999 <= R("R_EX_tre_e_") <= 999999
 -100000 <= R("MtrpLb_Transp") <= 100000
 -999999 <= R("R_EX_trp_L_e_") <= 999999
 -100000 <= R("Mtsulb_Transp") <= 100000
 -999999 <= R("R_EX_tsul_e_") <= 999999
 -100000 <= R("Mttdcab_Transp") <= 100000
 -999999 <= R("R_EX_ttdca_e_") <= 999999
 -100000 <= R("MtyrLb_Transp") <= 100000
 -999999 <= R("R_EX_tyr_L_e_") <= 999999
 0 <= R("R_ACGAMT") <= 999999
 0 <= R("R_UAG2Ei") <= 999999
 0 <= R("R_PAPPT3") <= 999999
 -999999 <= R("R_UDPG4E") <= 999999
 0 <= R("R_UDPGALM") <= 999999
 -100000 <= R("Murab_Transp") <= 100000
 -999999 <= R("R_EX_ura_e_") <= 999999
 -100000 <= R("Mureab_Transp") <= 100000
 -999999 <= R("R_EX_urea_e_") <= 999999
 -999999 <= R("R_UREAt") <= 999999
 -100000 <= R("Murib_Transp") <= 100000
 -999999 <= R("R_EX_uri_e_") <= 999999
 -100000 <= R("MvalLb_Transp") <= 100000
 -999999 <= R("R_EX_val_L_e_") <= 999999
 -100000 <= R("Mxanb_Transp") <= 100000
 -999999 <= R("R_EX_xan_e_") <= 999999
 -999999 <= R("R_XANt") <= 999999
 -100000 <= R("Mxtsnb_Transp") <= 100000
 -999999 <= R("R_EX_xtsn_e_") <= 999999
 -100000 <= R("MxylDb_Transp") <= 100000
 -999999 <= R("R_EX_xyl_D_e_") <= 999999
 -999999 <= R("R_XYLI1") <= 999999

End
