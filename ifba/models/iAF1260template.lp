Maximize
Z: R("R_Ec_biomass_iAF1260_core_59p81M")

Subject To
M10fthfc: -R("R_AICART") + R("R_AICART_Rev") - 0.000223
 R("R_Ec_biomass_iAF1260_core_59p81M") - R("R_FMETTRS") - R("R_FTHFD") -
 R("R_GARFT") + R("R_GARFT_Rev") + R("R_MTHFC") - R("R_MTHFC_Rev") -
 R("R_ULA4NFT")  = 0.

M12dgr120c: R("R_12DGR120tipp") - R("R_DAGK120") + R("R_PAPA120")  =
 0.

M12dgr120p: -R("R_12DGR120tipp") + R("R_PAPA120pp")  = 0.

M12dgr140c: R("R_12DGR140tipp") - R("R_DAGK140") + R("R_PAPA140")  =
 0.

M12dgr140p: -R("R_12DGR140tipp") + R("R_PAPA140pp")  = 0.

M12dgr141c: R("R_12DGR141tipp") - R("R_DAGK141") + R("R_PAPA141")  =
 0.

M12dgr141p: -R("R_12DGR141tipp") + R("R_PAPA141pp")  = 0.

M12dgr160c: R("R_12DGR160tipp") - R("R_DAGK160") + R("R_PAPA160")  =
 0.

M12dgr160p: -R("R_12DGR160tipp") + R("R_PAPA160pp")  = 0.

M12dgr161c: R("R_12DGR161tipp") - R("R_DAGK161") + R("R_PAPA161")  =
 0.

M12dgr161p: -R("R_12DGR161tipp") + R("R_PAPA161pp") +
 R("R_PETNT161pp")  = 0.

M12dgr180c: R("R_12DGR180tipp") - R("R_DAGK180") + R("R_PAPA180")  =
 0.

M12dgr180p: -R("R_12DGR180tipp") + R("R_PAPA180pp")  = 0.

M12dgr181c: R("R_12DGR181tipp") - R("R_DAGK181") + R("R_PAPA181")  =
 0.

M12dgr181p: -R("R_12DGR181tipp") + R("R_PAPA181pp") +
 R("R_PETNT181pp")  = 0.

M12ppdRb: R("M12ppdRb_Transp") + R("R_EX_12ppd_R_e_") -
 R("R_EX_12ppd_R_e__Rev")  = 0.

M12ppdRc: R("R_12PPDRtpp") - R("R_12PPDRtpp_Rev") + R("R_ALR4x") +
 R("R_LCARR") - R("R_LCARR_Rev")  = 0.

M12ppdRe: -R("R_12PPDRtex") + R("R_12PPDRtex_Rev") -
 R("R_EX_12ppd_R_e_") + R("R_EX_12ppd_R_e__Rev")  = 0.

M12ppdRp: R("R_12PPDRtex") - R("R_12PPDRtex_Rev") - R("R_12PPDRtpp")
 + R("R_12PPDRtpp_Rev")  = 0.

M12ppdSb: R("M12ppdSb_Transp") + R("R_EX_12ppd_S_e_") -
 R("R_EX_12ppd_S_e__Rev")  = 0.

M12ppdSc: R("R_12PPDStpp") - R("R_12PPDStpp_Rev") + R("R_LCARS") -
 R("R_LCARS_Rev")  = 0.

M12ppdSe: -R("R_12PPDStex") + R("R_12PPDStex_Rev") -
 R("R_EX_12ppd_S_e_") + R("R_EX_12ppd_S_e__Rev")  = 0.

M12ppdSp: R("R_12PPDStex") - R("R_12PPDStex_Rev") - R("R_12PPDStpp")
 + R("R_12PPDStpp_Rev")  = 0.

M13dpgc: R("R_GAPD") - R("R_GAPD_Rev") + R("R_PGK") - R("R_PGK_Rev")
  = 0.

M14glucanb: R("M14glucanb_Transp") + R("R_EX_14glucan_e_") -
 R("R_EX_14glucan_e__Rev")  = 0.

M14glucanc: R("R_14GLUCANabcpp") - R("R_AAMYL")  = 0.

M14glucane: -R("R_14GLUCANtexi") - R("R_EX_14glucan_e_") +
 R("R_EX_14glucan_e__Rev")  = 0.

M14glucanp: -R("R_14GLUCANabcpp") + R("R_14GLUCANtexi") -
 R("R_AAMYLpp")  = 0.

M15dapb: R("M15dapb_Transp") + R("R_EX_15dap_e_") -
 R("R_EX_15dap_e__Rev")  = 0.

M15dapc: -R("R_CADVtpp") + R("R_LYSDC")  = 0.

M15dape: -R("R_DAPtex") + R("R_DAPtex_Rev") - R("R_EX_15dap_e_") +
 R("R_EX_15dap_e__Rev")  = 0.

M15dapp: R("R_CADVtpp") + R("R_DAPtex") - R("R_DAPtex_Rev")  = 0.

M1agpe120p: -R("R_LPLIPAL1E120pp") + R("R_PLIPA2E120pp")  = 0.

M1agpe140p: -R("R_LPLIPAL1E140pp") + R("R_PLIPA2E140pp")  = 0.

M1agpe141p: -R("R_LPLIPAL1E141pp") + R("R_PLIPA2E141pp")  = 0.

M1agpe160p: -R("R_LPLIPAL1E160pp") + R("R_PLIPA2E160pp")  = 0.

M1agpe161p: -R("R_LPLIPAL1E161pp") + R("R_PLIPA2E161pp")  = 0.

M1agpe180p: -R("R_LPLIPAL1E180pp") + R("R_PLIPA2E180pp")  = 0.

M1agpe181p: -R("R_LPLIPAL1E181pp") + R("R_PLIPA2E181pp")  = 0.

M1agpg120p: -R("R_LPLIPAL1G120pp") + R("R_PLIPA2G120pp")  = 0.

M1agpg140p: -R("R_LPLIPAL1G140pp") + R("R_PLIPA2G140pp")  = 0.

M1agpg141p: -R("R_LPLIPAL1G141pp") + R("R_PLIPA2G141pp")  = 0.

M1agpg160p: -R("R_LPLIPAL1G160pp") + R("R_PLIPA2G160pp")  = 0.

M1agpg161p: -R("R_LPLIPAL1G161pp") + R("R_PLIPA2G161pp")  = 0.

M1agpg180p: -R("R_LPLIPAL1G180pp") + R("R_PLIPA2G180pp")  = 0.

M1agpg181p: -R("R_LPLIPAL1G181pp") + R("R_PLIPA2G181pp")  = 0.

M1ddecg3pc: -R("R_AGPAT120") + R("R_G3PAT120")  = 0.

M1ddecg3pp: -R("R_LPLIPAL1A120pp") + R("R_PLIPA2A120pp")  = 0.

M1hdec9eg3pc: -R("R_AGPAT161") + R("R_G3PAT161")  = 0.

M1hdec9eg3pp: -R("R_LPLIPAL1A161pp") + R("R_PLIPA2A161pp")  = 0.

M1hdecg3pc: -R("R_AGPAT160") + R("R_G3PAT160")  = 0.

M1hdecg3pp: -R("R_LPLIPAL1A160pp") + R("R_PLIPA2A160pp")  = 0.

M1odec11eg3pc: -R("R_AGPAT181") + R("R_G3PAT181")  = 0.

M1odec11eg3pp: -R("R_LPLIPAL1A181pp") + R("R_PLIPA2A181pp")  = 0.

M1odecg3pc: -R("R_AGPAT180") + R("R_G3PAT180")  = 0.

M1odecg3pp: -R("R_LPLIPAL1A180pp") + R("R_PLIPA2A180pp")  = 0.

M1pyr5cc: R("R_G5SADs") - R("R_P5CD") - R("R_P5CR") + R("R_PROD2") 
 = 0.

M1tdec7eg3pc: -R("R_AGPAT141") + R("R_G3PAT141")  = 0.

M1tdec7eg3pp: -R("R_LPLIPAL1A141pp") + R("R_PLIPA2A141pp")  = 0.

M1tdecg3pc: -R("R_AGPAT140") + R("R_G3PAT140")  = 0.

M1tdecg3pp: -R("R_LPLIPAL1A140pp") + R("R_PLIPA2A140pp")  = 0.

M23campb: R("M23campb_Transp") + R("R_EX_23camp_e_") -
 R("R_EX_23camp_e__Rev")  = 0.

M23campe: -R("R_23CAMPtex") + R("R_23CAMPtex_Rev") -
 R("R_EX_23camp_e_") + R("R_EX_23camp_e__Rev")  = 0.

M23campp: R("R_23CAMPtex") - R("R_23CAMPtex_Rev") - R("R_23PDE7pp") 
 = 0.

M23ccmpb: R("M23ccmpb_Transp") + R("R_EX_23ccmp_e_") -
 R("R_EX_23ccmp_e__Rev")  = 0.

M23ccmpe: -R("R_23CCMPtex") + R("R_23CCMPtex_Rev") -
 R("R_EX_23ccmp_e_") + R("R_EX_23ccmp_e__Rev")  = 0.

M23ccmpp: R("R_23CCMPtex") - R("R_23CCMPtex_Rev") - R("R_23PDE4pp") 
 = 0.

M23cgmpb: R("M23cgmpb_Transp") + R("R_EX_23cgmp_e_") -
 R("R_EX_23cgmp_e__Rev")  = 0.

M23cgmpe: -R("R_23CGMPtex") + R("R_23CGMPtex_Rev") -
 R("R_EX_23cgmp_e_") + R("R_EX_23cgmp_e__Rev")  = 0.

M23cgmpp: R("R_23CGMPtex") - R("R_23CGMPtex_Rev") - R("R_23PDE9pp") 
 = 0.

M23cumpb: R("M23cumpb_Transp") + R("R_EX_23cump_e_") -
 R("R_EX_23cump_e__Rev")  = 0.

M23cumpe: -R("R_23CUMPtex") + R("R_23CUMPtex_Rev") -
 R("R_EX_23cump_e_") + R("R_EX_23cump_e__Rev")  = 0.

M23cumpp: R("R_23CUMPtex") - R("R_23CUMPtex_Rev") - R("R_23PDE2pp") 
 = 0.

M23dappab: R("M23dappab_Transp") + R("R_EX_23dappa_e_") -
 R("R_EX_23dappa_e__Rev")  = 0.

M23dappac: R("R_23DAPPAt2pp") - R("R_DAPAL")  = 0.

M23dappae: -R("R_23DAPPAtex") + R("R_23DAPPAtex_Rev") -
 R("R_EX_23dappa_e_") + R("R_EX_23dappa_e__Rev")  = 0.

M23dappap: -R("R_23DAPPAt2pp") + R("R_23DAPPAtex") -
 R("R_23DAPPAtex_Rev")  = 0.

M23ddhbc: -R("R_DHBD") + R("R_DHBD_Rev") + R("R_ICHORT")  = 0.

M23dhbac: R("R_DHBS") - 3 R("R_ENTCS")  = 0.

M23dhbc: R("R_DHBD") - R("R_DHBD_Rev") - R("R_DHBS") + R("R_DHBSH") 
 = 0.

M23dhbzsc: -R("R_DHBSH") + 3 R("R_ENTERES") + 3 R("R_ENTERES2") +
 R("R_FE3DHBZR")  = 0.

M23dhdpc: -R("R_DHDPRy") + R("R_DHDPS")  = 0.

M23dhmbc: -R("R_DHAD1") - R("R_KARA1") + R("R_KARA1_Rev")  = 0.

M23dhmpc: -R("R_DHAD2") + R("R_KARA2") - R("R_KARA2_Rev")  = 0.

M23dogulnc: -R("R_DOGULNR")  = 0.

M25aicsc: -R("R_ADSL2r") + R("R_ADSL2r_Rev") + R("R_PRASCSi")  = 0.

M25dkglcnc: -R("R_DKGLCNR1") - R("R_DKGLCNR2x") - R("R_DKGLCNR2y") 
 = 0.

M25drappc: -R("R_DHPPDA2") + R("R_GTPCII2")  = 0.

M26dapLLc: -R("R_DAPE") + R("R_DAPE_Rev") + R("R_SDPDS")  = 0.

M26dapMb: R("M26dapMb_Transp") + R("R_EX_26dap_M_e_") -
 R("R_EX_26dap_M_e__Rev")  = 0.

M26dapMc: R("R_DAPabcpp") - R("R_DAPDC") + R("R_DAPE") -
 R("R_DAPE_Rev") + R("R_LADGMDH") - R("R_UAAGDS")  = 0.

M26dapMe: -R("R_26DAHtex") + R("R_26DAHtex_Rev") -
 R("R_EX_26dap_M_e_") + R("R_EX_26dap_M_e__Rev")  = 0.

M26dapMp: R("R_26DAHtex") - R("R_26DAHtex_Rev") - R("R_DAPabcpp")  =
 0.

M2agpe120c: R("R_2AGPE120tipp") - R("R_2AGPEAT120") -
 R("R_LPLIPAL2ATE120") - R("R_LPLIPAL2E120")  = 0.

M2agpe120p: -R("R_2AGPE120tipp") + R("R_PLIPA1E120pp")  = 0.

M2agpe140c: R("R_2AGPE140tipp") - R("R_2AGPEAT140") -
 R("R_LPLIPAL2ATE140") - R("R_LPLIPAL2E140")  = 0.

M2agpe140p: -R("R_2AGPE140tipp") + R("R_PLIPA1E140pp")  = 0.

M2agpe141c: R("R_2AGPE141tipp") - R("R_2AGPEAT141") -
 R("R_LPLIPAL2ATE141") - R("R_LPLIPAL2E141")  = 0.

M2agpe141p: -R("R_2AGPE141tipp") + R("R_PLIPA1E141pp")  = 0.

M2agpe160c: R("R_2AGPE160tipp") - R("R_2AGPEAT160") -
 R("R_LPLIPAL2ATE160") - R("R_LPLIPAL2E160")  = 0.

M2agpe160p: -R("R_2AGPE160tipp") + R("R_ALPATE160pp") +
 R("R_PLIPA1E160pp")  = 0.

M2agpe161c: R("R_2AGPE161tipp") - R("R_2AGPEAT161") -
 R("R_LPLIPAL2ATE161") - R("R_LPLIPAL2E161")  = 0.

M2agpe161p: -R("R_2AGPE161tipp") + R("R_PLIPA1E161pp")  = 0.

M2agpe180c: R("R_2AGPE180tipp") - R("R_2AGPEAT180") -
 R("R_LPLIPAL2ATE180") - R("R_LPLIPAL2E180")  = 0.

M2agpe180p: -R("R_2AGPE180tipp") + R("R_PLIPA1E180pp")  = 0.

M2agpe181c: R("R_2AGPE181tipp") - R("R_2AGPEAT181") -
 R("R_LPLIPAL2ATE181") - R("R_LPLIPAL2E181")  = 0.

M2agpe181p: -R("R_2AGPE181tipp") + R("R_PLIPA1E181pp")  = 0.

M2agpg120c: R("R_2AGPG120tipp") - R("R_2AGPGAT120") -
 R("R_LPLIPAL2ATG120") - R("R_LPLIPAL2G120")  = 0.

M2agpg120p: -R("R_2AGPG120tipp") + R("R_PLIPA1G120pp")  = 0.

M2agpg140c: R("R_2AGPG140tipp") - R("R_2AGPGAT140") -
 R("R_LPLIPAL2ATG140") - R("R_LPLIPAL2G140")  = 0.

M2agpg140p: -R("R_2AGPG140tipp") + R("R_PLIPA1G140pp")  = 0.

M2agpg141c: R("R_2AGPG141tipp") - R("R_2AGPGAT141") -
 R("R_LPLIPAL2ATG141") - R("R_LPLIPAL2G141")  = 0.

M2agpg141p: -R("R_2AGPG141tipp") + R("R_PLIPA1G141pp")  = 0.

M2agpg160c: R("R_2AGPG160tipp") - R("R_2AGPGAT160") -
 R("R_LPLIPAL2ATG160") - R("R_LPLIPAL2G160")  = 0.

M2agpg160p: -R("R_2AGPG160tipp") + R("R_ALPATG160pp") +
 R("R_PLIPA1G160pp")  = 0.

M2agpg161c: R("R_2AGPG161tipp") - R("R_2AGPGAT161") -
 R("R_LPLIPAL2ATG161") - R("R_LPLIPAL2G161")  = 0.

M2agpg161p: -R("R_2AGPG161tipp") + R("R_PLIPA1G161pp")  = 0.

M2agpg180c: R("R_2AGPG180tipp") - R("R_2AGPGAT180") -
 R("R_LPLIPAL2ATG180") - R("R_LPLIPAL2G180")  = 0.

M2agpg180p: -R("R_2AGPG180tipp") + R("R_PLIPA1G180pp")  = 0.

M2agpg181c: R("R_2AGPG181tipp") - R("R_2AGPGAT181") -
 R("R_LPLIPAL2ATG181") - R("R_LPLIPAL2G181")  = 0.

M2agpg181p: -R("R_2AGPG181tipp") + R("R_PLIPA1G181pp")  = 0.

M2ahbutc: R("R_ACHBS") - R("R_KARA2") + R("R_KARA2_Rev")  = 0.

M2amsac: R("R_DSERDHr") - R("R_DSERDHr_Rev") + R("R_LSERDHr") -
 R("R_LSERDHr_Rev")  = 0.

M2aobutc: -R("R_AOBUTDs") + R("R_ATHRDHr") - R("R_ATHRDHr_Rev") +
 R("R_GLYAT") - R("R_GLYAT_Rev") + R("R_THRD")  = 0.

M2cpr5pc: -R("R_IGPS") + R("R_PRAIi")  = 0.

M2dda7pc: R("R_DDPA") - R("R_DHQS")  = 0.

M2ddecg3pc: R("R_2AGPA120tipp") - R("R_LPLIPAL2A120")  = 0.

M2ddecg3pp: -R("R_2AGPA120tipp") + R("R_PLIPA1A120pp")  = 0.

M2ddg6pc: R("R_DDGLK") - R("R_EDA") + R("R_EDD")  = 0.

M2ddglcnb: R("M2ddglcnb_Transp") + R("R_EX_2ddglcn_e_") -
 R("R_EX_2ddglcn_e__Rev")  = 0.

M2ddglcnc: R("R_ALTRH") + R("R_DDGLCNt2rpp") -
 R("R_DDGLCNt2rpp_Rev") - R("R_DDGLK") + R("R_MNNH")  = 0.

M2ddglcne: -R("R_DDGLCNtex") + R("R_DDGLCNtex_Rev") -
 R("R_EX_2ddglcn_e_") + R("R_EX_2ddglcn_e__Rev")  = 0.

M2ddglcnp: -R("R_DDGLCNt2rpp") + R("R_DDGLCNt2rpp_Rev") +
 R("R_DDGLCNtex") - R("R_DDGLCNtex_Rev")  = 0.

M2dh3dgal6pc: R("R_DDGALK") - R("R_DDPGALA") + R("R_DDPGALA_Rev")  =
 0.

M2dh3dgalc: -R("R_DDGALK") + R("R_GALCTND")  = 0.

M2dhglcnc: -R("R_2DGLCNRx") - R("R_2DGLCNRy")  = 0.

M2dhgulnc: -R("R_2DGULRx") - R("R_2DGULRy") + R("R_DKGLCNR1")  = 0.

M2dhpc: -R("R_DPR") + R("R_MOHMT")  = 0.

M2dmmq8c: R("R_DMSOR2") + R("R_DMSOR2pp") + R("R_FRD3") -
 R("R_G3PD7") - R("R_GLYCTO4") - R("R_HYD3pp") - R("R_NADH18pp") - R("R_NADH9") -
 R("R_NADPHQR4") + R("R_TMAOR2") + R("R_TMAOR2pp")  = 0.

M2dmmql8c: -R("R_AMMQLT8") + R("R_DHNAOT4") - R("R_DMSOR2") -
 R("R_DMSOR2pp") - R("R_FRD3") + R("R_G3PD7") + R("R_GLYCTO4") +
 R("R_HYD3pp") + R("R_NADH18pp") + R("R_NADH9") + R("R_NADPHQR4") - R("R_TMAOR2") -
 R("R_TMAOR2pp")  = 0.

M2dr1pc: R("R_DURIPP") - R("R_DURIPP_Rev") - R("R_PPM2") +
 R("R_PPM2_Rev") + R("R_PUNP2") - R("R_PUNP2_Rev") + R("R_PUNP4") -
 R("R_PUNP4_Rev") + R("R_PUNP6") - R("R_PUNP6_Rev") + R("R_TMDPP") -
 R("R_TMDPP_Rev")  = 0.

M2dr5pc: -R("R_DRPA") + R("R_PPM2") - R("R_PPM2_Rev")  = 0.

M2h3oppanc: R("R_GLCRAL") + R("R_GLXCL") + R("R_HPYRI") -
 R("R_HPYRI_Rev") - R("R_TRSARr") + R("R_TRSARr_Rev")  = 0.

M2hdec9eg3pc: R("R_2AGPA161tipp") - R("R_LPLIPAL2A161")  = 0.

M2hdec9eg3pp: -R("R_2AGPA161tipp") + R("R_PLIPA1A161pp")  = 0.

M2hdecg3pc: R("R_2AGPA160tipp") - R("R_LPLIPAL2A160")  = 0.

M2hdecg3pp: -R("R_2AGPA160tipp") + R("R_PLIPA1A160pp")  = 0.

M2ippmc: R("R_IPPMIa") - R("R_IPPMIa_Rev") - R("R_IPPMIb") +
 R("R_IPPMIb_Rev")  = 0.

M2kmbc: R("R_DKMPPD") + R("R_DKMPPD3") - R("R_UNK3")  = 0.

M2mahmpc: -R("R_2MAHMP") + R("R_PMPK") - R("R_TMPPP")  = 0.

M2mcacnc: R("R_MCITD") - R("R_MICITD")  = 0.

M2mcitc: -R("R_MCITD") + R("R_MCITS")  = 0.

M2me4pc: R("R_DXPRIi") - R("R_MEPCT")  = 0.

M2mecdpc: -R("R_MECDPDH2") + R("R_MECDPS")  = 0.

M2obutc: -R("R_ACHBS") - R("R_OBTFL") + R("R_THRD_L")  = 0.

M2odec11eg3pc: R("R_2AGPA181tipp") - R("R_LPLIPAL2A181")  = 0.

M2odec11eg3pp: -R("R_2AGPA181tipp") + R("R_PLIPA1A181pp")  = 0.

M2odecg3pc: R("R_2AGPA180tipp") - R("R_LPLIPAL2A180")  = 0.

M2odecg3pp: -R("R_2AGPA180tipp") + R("R_PLIPA1A180pp")  = 0.

M2ohphc: -0.000223 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_OHPHM") + R("R_OPHHX") + R("R_OPHHX3")  = 0.

M2ombzlc: -R("R_OMBZLM") + R("R_OMPHHX") + R("R_OMPHHX3")  = 0.

M2omhmblc: -R("R_DMQMT") + R("R_OMMBLHX") + R("R_OMMBLHX3")  = 0.

M2ommblc: R("R_OMBZLM") - R("R_OMMBLHX") - R("R_OMMBLHX3")  = 0.

M2omphc: R("R_OHPHM") - R("R_OMPHHX") - R("R_OMPHHX3")  = 0.

M2ophc: R("R_OPHBDC") - R("R_OPHHX") - R("R_OPHHX3")  = 0.

M2p4c2mec: R("R_CDPMEK") - R("R_MECDPS")  = 0.

M2pgc: -R("R_ENO") + R("R_ENO_Rev") + R("R_GLYCK2") - R("R_PGM") +
 R("R_PGM_Rev")  = 0.

M2pglycc: -R("R_PGLYCP")  = 0.

M2shchcc: R("R_SHCHCS2") - R("R_SUCBZS")  = 0.

M2tdec7eg3pc: R("R_2AGPA141tipp") - R("R_LPLIPAL2A141")  = 0.

M2tdec7eg3pp: -R("R_2AGPA141tipp") + R("R_PLIPA1A141pp")  = 0.

M2tdecg3pc: R("R_2AGPA140tipp") - R("R_LPLIPAL2A140")  = 0.

M2tdecg3pp: -R("R_2AGPA140tipp") + R("R_PLIPA1A140pp")  = 0.

M2tpr3dpcoac: R("R_TPRDCOAS")  = 0.

M34dhpacb: R("M34dhpacb_Transp") + R("R_EX_34dhpac_e_") -
 R("R_EX_34dhpac_e__Rev")  = 0.

M34dhpace: -R("R_34dhpactex") + R("R_34dhpactex_Rev") -
 R("R_EX_34dhpac_e_") + R("R_EX_34dhpac_e__Rev")  = 0.

M34dhpacp: R("R_34dhpactex") - R("R_34dhpactex_Rev") +
 R("R_42A12BOOXpp")  = 0.

M34hppc: R("R_PPND") + R("R_TYRTA") - R("R_TYRTA_Rev")  = 0.

M3ampb: R("M3ampb_Transp") + R("R_EX_3amp_e_") -
 R("R_EX_3amp_e__Rev")  = 0.

M3ampe: -R("R_3AMPtex") + R("R_3AMPtex_Rev") - R("R_EX_3amp_e_") +
 R("R_EX_3amp_e__Rev")  = 0.

M3ampp: R("R_23PDE7pp") + R("R_3AMPtex") - R("R_3AMPtex_Rev") -
 R("R_3NTD7pp")  = 0.

M3c2hmpc: -R("R_IPMD") - R("R_IPPMIa") + R("R_IPPMIa_Rev")  = 0.

M3c3hmpc: R("R_IPPMIb") - R("R_IPPMIb_Rev") + R("R_IPPS")  = 0.

M3c4mopc: R("R_IPMD") - R("R_OMCDC")  = 0.

M3cmpb: R("M3cmpb_Transp") + R("R_EX_3cmp_e_") -
 R("R_EX_3cmp_e__Rev")  = 0.

M3cmpe: -R("R_3CMPtex") + R("R_3CMPtex_Rev") - R("R_EX_3cmp_e_") +
 R("R_EX_3cmp_e__Rev")  = 0.

M3cmpp: R("R_23PDE4pp") + R("R_3CMPtex") - R("R_3CMPtex_Rev") -
 R("R_3NTD4pp")  = 0.

M3dhgulnc: -R("R_3KGK") + R("R_DOGULNR")  = 0.

M3dhgulnpc: R("R_3KGK") + R("R_ASCBPL") - R("R_KG6PDC")  = 0.

M3dhqc: R("R_DHQS") - R("R_DHQTi")  = 0.

M3dhskc: R("R_DHQTi") - R("R_SHK3Dr") + R("R_SHK3Dr_Rev")  = 0.

M3gmpb: R("M3gmpb_Transp") + R("R_EX_3gmp_e_") -
 R("R_EX_3gmp_e__Rev")  = 0.

M3gmpe: -R("R_3GMPtex") + R("R_3GMPtex_Rev") - R("R_EX_3gmp_e_") +
 R("R_EX_3gmp_e__Rev")  = 0.

M3gmpp: R("R_23PDE9pp") + R("R_3GMPtex") - R("R_3GMPtex_Rev") -
 R("R_3NTD9pp")  = 0.

M3haACPc: -R("R_3HAD40") + R("R_3OAR40") - R("R_3OAR40_Rev")  = 0.

M3hbcoac: -R("R_ECOAH1") + R("R_ECOAH1_Rev") - R("R_HACD1i")  = 0.

M3hcddec5eACPc: -R("R_3HAD121") + R("R_3OAR121")  = 0.

M3hcinnmb: R("M3hcinnmb_Transp") + R("R_EX_3hcinnm_e_") -
 R("R_EX_3hcinnm_e__Rev")  = 0.

M3hcinnmc: -R("R_3HCINNMH") + R("R_HCINNMt2rpp") -
 R("R_HCINNMt2rpp_Rev")  = 0.

M3hcinnme: -R("R_EX_3hcinnm_e_") + R("R_EX_3hcinnm_e__Rev") -
 R("R_HCINNMtex") + R("R_HCINNMtex_Rev")  = 0.

M3hcinnmp: -R("R_HCINNMt2rpp") + R("R_HCINNMt2rpp_Rev") +
 R("R_HCINNMtex") - R("R_HCINNMtex_Rev")  = 0.

M3hcmrs7eACPc: -R("R_3HAD141") + R("R_3OAR141")  = 0.

M3hcpalm9eACPc: -R("R_3HAD161") + R("R_3OAR161")  = 0.

M3hcvac11eACPc: -R("R_3HAD181") + R("R_3OAR181")  = 0.

M3hdcoac: -R("R_ECOAH4") + R("R_ECOAH4_Rev") - R("R_HACD4i")  = 0.

M3hddcoac: -R("R_ECOAH5") + R("R_ECOAH5_Rev") - R("R_HACD5i")  = 0.

M3hddecACPc: -R("R_3HAD120") + R("R_3OAR120") - R("R_3OAR120_Rev") 
 = 0.

M3hdecACPc: -R("R_3HAD100") + R("R_3OAR100") - R("R_3OAR100_Rev")  =
 0.

M3hhcoac: -R("R_ECOAH2") + R("R_ECOAH2_Rev") - R("R_HACD2i")  = 0.

M3hhdcoac: -R("R_ECOAH7") + R("R_ECOAH7_Rev") - R("R_HACD7i")  = 0.

M3hhexACPc: -R("R_3HAD60") + R("R_3OAR60") - R("R_3OAR60_Rev")  = 0.

M3hmrsACPc: -R("R_3HAD140") + R("R_3OAR140") - R("R_3OAR140_Rev") -
 R("R_U23GAAT") - R("R_UAGAAT") + R("R_UAGAAT_Rev")  = 0.

M3hocoac: -R("R_ECOAH3") + R("R_ECOAH3_Rev") - R("R_HACD3i")  = 0.

M3hoctaACPc: -R("R_3HAD180") + R("R_3OAR180") - R("R_3OAR180_Rev") 
 = 0.

M3hoctACPc: -R("R_3HAD80") + R("R_3OAR80") - R("R_3OAR80_Rev")  = 0.

M3hodcoac: -R("R_ECOAH8") + R("R_ECOAH8_Rev") - R("R_HACD8i")  = 0.

M3hpalmACPc: -R("R_3HAD160") + R("R_3OAR160") - R("R_3OAR160_Rev") 
 = 0.

M3hpppnb: R("M3hpppnb_Transp") + R("R_EX_3hpppn_e_") -
 R("R_EX_3hpppn_e__Rev")  = 0.

M3hpppnc: -R("R_3HPPPNH") + R("R_HPPPNt2rpp") -
 R("R_HPPPNt2rpp_Rev")  = 0.

M3hpppne: -R("R_EX_3hpppn_e_") + R("R_EX_3hpppn_e__Rev") -
 R("R_HPPPNtex") + R("R_HPPPNtex_Rev")  = 0.

M3hpppnp: -R("R_HPPPNt2rpp") + R("R_HPPPNt2rpp_Rev") +
 R("R_HPPPNtex") - R("R_HPPPNtex_Rev")  = 0.

M3htdcoac: -R("R_ECOAH6") + R("R_ECOAH6_Rev") - R("R_HACD6i")  = 0.

M3ig3pc: R("R_IGPS") - R("R_TRPS1") - R("R_TRPS3")  = 0.

M3mobc: R("R_DHAD1") - R("R_IPPS") - R("R_MOHMT") + R("R_VALTA") -
 R("R_VALTA_Rev") - R("R_VPAMT")  = 0.

M3mopc: R("R_DHAD2") + R("R_ILETA") - R("R_ILETA_Rev")  = 0.

M3ocddec5eACPc: -R("R_3OAR121") + R("R_3OAS121")  = 0.

M3ocmrs7eACPc: -R("R_3OAR141") + R("R_3OAS141")  = 0.

M3ocpalm9eACPc: -R("R_3OAR161") + R("R_3OAS161")  = 0.

M3ocvac11eACPc: -R("R_3OAR181") + R("R_3OAS181")  = 0.

M3odcoac: R("R_HACD4i") - R("R_KAT4")  = 0.

M3oddcoac: R("R_HACD5i") - R("R_KAT5")  = 0.

M3oddecACPc: -R("R_3OAR120") + R("R_3OAR120_Rev") + R("R_3OAS120") 
 = 0.

M3odecACPc: -R("R_3OAR100") + R("R_3OAR100_Rev") + R("R_3OAS100")  =
 0.

M3ohcoac: R("R_HACD2i") - R("R_KAT2")  = 0.

M3ohdcoac: R("R_HACD7i") - R("R_KAT7")  = 0.

M3ohexACPc: -R("R_3OAR60") + R("R_3OAR60_Rev") + R("R_3OAS60")  = 0.

M3omrsACPc: -R("R_3OAR140") + R("R_3OAR140_Rev") + R("R_3OAS140")  =
 0.

M3oocoac: R("R_HACD3i") - R("R_KAT3")  = 0.

M3ooctACPc: -R("R_3OAR80") + R("R_3OAR80_Rev") + R("R_3OAS80")  = 0.

M3ooctdACPc: -R("R_3OAR180") + R("R_3OAR180_Rev") + R("R_3OAS180") 
 = 0.

M3oodcoac: R("R_HACD8i") - R("R_KAT8")  = 0.

M3opalmACPc: -R("R_3OAR160") + R("R_3OAR160_Rev") + R("R_3OAS160") 
 = 0.

M3ophbc: R("R_HBZOPT") - R("R_OPHBDC")  = 0.

M3otdcoac: R("R_HACD6i") - R("R_KAT6")  = 0.

M3pgc: R("R_GLYCK") - R("R_PGCD") - R("R_PGK") + R("R_PGK_Rev") +
 R("R_PGM") - R("R_PGM_Rev")  = 0.

M3phpc: R("R_PGCD") - R("R_PSERT")  = 0.

M3psmec: -R("R_CHORS") + R("R_PSCVT") - R("R_PSCVT_Rev")  = 0.

M3salac: -R("R_CYSSADS")  = 0.

M3umpb: R("M3umpb_Transp") + R("R_EX_3ump_e_") -
 R("R_EX_3ump_e__Rev")  = 0.

M3umpe: -R("R_3UMPtex") + R("R_3UMPtex_Rev") - R("R_EX_3ump_e_") +
 R("R_EX_3ump_e__Rev")  = 0.

M3umpp: R("R_23PDE2pp") - R("R_3NTD2pp") + R("R_3UMPtex") -
 R("R_3UMPtex_Rev")  = 0.

M4abutb: R("M4abutb_Transp") + R("R_EX_4abut_e_") -
 R("R_EX_4abut_e__Rev")  = 0.

M4abutc: -R("R_ABTA") + R("R_ABUTD") + R("R_ABUTt2pp") +
 R("R_GGGABAH") - R("R_GLUABUTt7pp") + R("R_GLUABUTt7pp_Rev") + R("R_GLUDC")  =
 0.

M4abute: -R("R_ABUTtex") + R("R_ABUTtex_Rev") - R("R_EX_4abut_e_") +
 R("R_EX_4abut_e__Rev")  = 0.

M4abutnc: -R("R_ABUTD") + R("R_PTRCTA")  = 0.

M4abutp: -R("R_ABUTt2pp") + R("R_ABUTtex") - R("R_ABUTtex_Rev") +
 R("R_GLUABUTt7pp") - R("R_GLUABUTt7pp_Rev")  = 0.

M4abzc: R("R_ADCL") - R("R_DHPS2")  = 0.

M4adchoc: -R("R_ADCL") + R("R_ADCS")  = 0.

M4ahmmpc: -R("R_HMPK1")  = 0.

M4ampmc: R("R_2MAHMP") + R("R_AMPMS2") + R("R_HMPK1") - R("R_PMPK") 
 = 0.

M4c2mec: -R("R_CDPMEK") + R("R_MEPCT")  = 0.

M4h2opntnc: -R("R_HOPNTAL") + R("R_OP4ENH")  = 0.

M4hbab: R("M4hbab_Transp") + R("R_DM_4HBA")  = 0.

M4hbac: -R("R_DM_4HBA") + R("R_THZPSN")  = 0.

M4hbzc: R("R_CHRPL") - R("R_HBZOPT")  = 0.

M4hoxpacdb: R("M4hoxpacdb_Transp") + R("R_EX_4hoxpacd_e_") -
 R("R_EX_4hoxpacd_e__Rev")  = 0.

M4hoxpacde: -R("R_4HOXPACDtex") + R("R_4HOXPACDtex_Rev") -
 R("R_EX_4hoxpacd_e_") + R("R_EX_4hoxpacd_e__Rev")  = 0.

M4hoxpacdp: R("R_4HOXPACDtex") - R("R_4HOXPACDtex_Rev") +
 R("R_TYROXDApp")  = 0.

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

M5aizcc: -R("R_AIRC3") + R("R_AIRC3_Rev") - R("R_PRASCSi")  = 0.

M5aopc: R("R_G1SAT") - R("R_G1SAT_Rev") - 2 R("R_PPBNGS")  = 0.

M5aprbuc: R("R_APRAUR") - R("R_PMDPHT")  = 0.

M5apruc: -R("R_APRAUR") + R("R_DHPPDA2")  = 0.

M5caizc: R("R_AIRC2") + R("R_AIRC3") - R("R_AIRC3_Rev")  = 0.

M5dglcnb: R("M5dglcnb_Transp") + R("R_EX_5dglcn_e_") -
 R("R_EX_5dglcn_e__Rev")  = 0.

M5dglcnc: -R("R_5DGLCNR") + R("R_5DGLCNR_Rev") + R("R_5DGLCNt2rpp")
 - R("R_5DGLCNt2rpp_Rev") + R("R_DKGLCNR2x") + R("R_DKGLCNR2y") -
 R("R_IDOND") - R("R_IDOND2") + R("R_IDOND_Rev")  = 0.

M5dglcne: -R("R_5DGLCNtex") + R("R_5DGLCNtex_Rev") -
 R("R_EX_5dglcn_e_") + R("R_EX_5dglcn_e__Rev")  = 0.

M5dglcnp: -R("R_5DGLCNt2rpp") + R("R_5DGLCNt2rpp_Rev") +
 R("R_5DGLCNtex") - R("R_5DGLCNtex_Rev")  = 0.

M5dh4dglcc: R("R_GALCTD") - R("R_GLCRAL") + R("R_GLCRD")  = 0.

M5dribb: R("M5dribb_Transp") + R("R_DM_5DRIB")  = 0.

M5dribc: R("R_5DOAN") - R("R_DM_5DRIB")  = 0.

M5mdr1pc: -R("R_MTRI") + R("R_MTRI_Rev") + R("R_MTRK")  = 0.

M5mdru1pc: -R("R_MDRPD") + R("R_MTRI") - R("R_MTRI_Rev")  = 0.

M5mtac: -R("R_MTAN") + R("R_SPMS")  = 0.

M5mthfc: -R("R_METS") + R("R_MTHFR2")  = 0.

M5mtrc: R("R_MTAN") - R("R_MTRK")  = 0.

M5prdmbzc: R("R_NNDMBRT") - R("R_RZ5PP")  = 0.

M6hmhptc: R("R_DHNPA2") - R("R_HPPK2")  = 0.

M6hmhptppc: -R("R_DHPS2") + R("R_HPPK2")  = 0.

M6pgcc: -R("R_EDD") - R("R_GND") + R("R_GNK") + R("R_PGL")  = 0.

M6pglc: R("R_G6PDH2r") - R("R_G6PDH2r_Rev") - R("R_PGL")  = 0.

M8aonnc: -R("R_AMAOTr") + R("R_AMAOTr_Rev") + R("R_AOXSr") -
 R("R_AOXSr_Rev")  = 0.

Maacaldb: R("Maacaldb_Transp") + R("R_DM_AACALD")  = 0.

Maacaldc: -R("R_DM_AACALD") + R("R_TAUDO")  = 0.

Maacoac: R("R_ACACCT") + R("R_ACACT1r") - R("R_ACACT1r_Rev") +
 R("R_HACD1i") - R("R_KAT1")  = 0.

Maactc: -R("R_AACTOOR") + R("R_AOBUTDs")  = 0.

Macacb: R("Macacb_Transp") + R("R_EX_acac_e_") -
 R("R_EX_acac_e__Rev")  = 0.

Macacc: -R("R_ACACCT") + R("R_ACACt2pp") - R("R_ACACt2pp_Rev")  = 0.

Macace: -R("R_ACACtex") + R("R_ACACtex_Rev") - R("R_EX_acac_e_") +
 R("R_EX_acac_e__Rev")  = 0.

Macacp: -R("R_ACACt2pp") + R("R_ACACt2pp_Rev") + R("R_ACACtex") -
 R("R_ACACtex_Rev")  = 0.

MacACPc: R("R_ACOATA") - R("R_ACOATA_Rev") - R("R_KAS14") +
 R("R_MACPD")  = 0.

Macaldb: R("Macaldb_Transp") + R("R_EX_acald_e_") -
 R("R_EX_acald_e__Rev")  = 0.

Macaldc: -R("R_ACALD") + R("R_ACALD_Rev") + R("R_ACALDtpp") -
 R("R_ACALDtpp_Rev") + R("R_ALCD2x") - R("R_ALCD2x_Rev") - R("R_ALDD2x") -
 R("R_ALDD2y") + R("R_DRPA") + R("R_ETHAAL") + R("R_FDMO3") +
 R("R_HOPNTAL") + R("R_THRA2i") + R("R_THRAi")  = 0.

Macalde: -R("R_ACALDtex") + R("R_ACALDtex_Rev") - R("R_EX_acald_e_")
 + R("R_EX_acald_e__Rev")  = 0.

Macaldp: R("R_ACALDtex") - R("R_ACALDtex_Rev") - R("R_ACALDtpp") +
 R("R_ACALDtpp_Rev")  = 0.

Macanthc: R("R_ACANTHAT")  = 0.

Macb: R("Macb_Transp") + R("R_EX_ac_e_") - R("R_EX_ac_e__Rev")  = 0.

Macc: R("R_ACACCT") - R("R_ACKr") + R("R_ACKr_Rev") + R("R_ACODA") -
 R("R_ACS") + R("R_ACt2rpp") - R("R_ACt2rpp_Rev") + R("R_ACt4pp") +
 R("R_AGDC") + R("R_ALDD2x") + R("R_ALDD2y") + R("R_BUTCT") + R("R_CITL")
 + R("R_CYSS") + R("R_HXCT") + R("R_NACODA") + R("R_POX") +
 R("R_UHGADA")  = 0.

Maccoac: -R("R_ACACCT") - 2 R("R_ACACT1r") + 2 R("R_ACACT1r_Rev") +
 R("R_ACALD") - R("R_ACALD_Rev") - R("R_ACANTHAT") - R("R_ACCOAC") -
 R("R_ACGS") - R("R_ACOATA") + R("R_ACOATA_Rev") + R("R_ACS") -
 R("R_BUTCT") - R("R_CS") - R("R_G1PACT") - R("R_GLCATr") + R("R_GLCATr_Rev") -
 R("R_GLYAT") + R("R_GLYAT_Rev") - R("R_HXCT") - R("R_IPPS") -
 R("R_KAS15") + 2 R("R_KAT1") + R("R_KAT2") + R("R_KAT3") + R("R_KAT4") +
 R("R_KAT5") + R("R_KAT6") + R("R_KAT7") + R("R_KAT8") - R("R_MALS") -
 R("R_MALTATr") + R("R_MALTATr_Rev") - R("R_O16AT") + R("R_PDH") +
 R("R_PFL") - R("R_PTAr") + R("R_PTAr_Rev") - R("R_SERAT") + R("R_SERAT_Rev") -
 R("R_SPMDAT1") - R("R_SPMDAT2") - R("R_TDPADGAT")  = 0.

Mace: -R("R_ACtex") + R("R_ACtex_Rev") - R("R_EX_ac_e_") +
 R("R_EX_ac_e__Rev")  = 0.

Macetolc: R("R_ALR2") - R("R_ALR4x")  = 0.

Macg5pc: R("R_ACGK") + R("R_AGPR") - R("R_AGPR_Rev")  = 0.

Macg5sac: R("R_ACOTA") - R("R_ACOTA_Rev") - R("R_AGPR") +
 R("R_AGPR_Rev") - R("R_NACODA")  = 0.

Macgal1pb: R("Macgal1pb_Transp") + R("R_EX_acgal1p_e_") -
 R("R_EX_acgal1p_e__Rev")  = 0.

Macgal1pe: -R("R_ACGAL1Ptex") + R("R_ACGAL1Ptex_Rev") -
 R("R_EX_acgal1p_e_") + R("R_EX_acgal1p_e__Rev")  = 0.

Macgal1pp: -R("R_ACGAL1PPpp") + R("R_ACGAL1Ptex") -
 R("R_ACGAL1Ptex_Rev") + R("R_UACGALPpp")  = 0.

Macgalb: R("Macgalb_Transp") + R("R_EX_acgal_e_") -
 R("R_EX_acgal_e__Rev")  = 0.

Macgale: -R("R_ACGALtex") + R("R_ACGALtex_Rev") - R("R_EX_acgal_e_")
 + R("R_EX_acgal_e__Rev")  = 0.

Macgalp: R("R_ACGAL1PPpp") + R("R_ACGALtex") - R("R_ACGALtex_Rev") 
 = 0.

Macgam1pb: R("Macgam1pb_Transp") + R("R_EX_acgam1p_e_") -
 R("R_EX_acgam1p_e__Rev")  = 0.

Macgam1pc: R("R_G1PACT") - R("R_UAGDP")  = 0.

Macgam1pe: -R("R_ACGAM1Ptex") + R("R_ACGAM1Ptex_Rev") -
 R("R_EX_acgam1p_e_") + R("R_EX_acgam1p_e__Rev")  = 0.

Macgam1pp: -R("R_ACGAM1PPpp") + R("R_ACGAM1Ptex") -
 R("R_ACGAM1Ptex_Rev") + R("R_UACGAMPpp")  = 0.

Macgam6pc: R("R_ACGAMK") + R("R_ACGAptspp") + R("R_ACM6PH") -
 R("R_AGDC") + R("R_AMANAPEr") - R("R_AMANAPEr_Rev")  = 0.

Macgamb: R("Macgamb_Transp") + R("R_EX_acgam_e_") -
 R("R_EX_acgam_e__Rev")  = 0.

Macgamc: -R("R_ACGAMK") + R("R_AGM3PH") + R("R_AGM4PH") +
 R("R_AGMH")  = 0.

Macgame: -R("R_ACGAtex") + R("R_ACGAtex_Rev") - R("R_EX_acgam_e_") +
 R("R_EX_acgam_e__Rev")  = 0.

Macgamp: R("R_ACGAM1PPpp") - R("R_ACGAptspp") + R("R_ACGAtex") -
 R("R_ACGAtex_Rev")  = 0.

MacglcDc: R("R_GLCATr") - R("R_GLCATr_Rev")  = 0.

Macgluc: -R("R_ACGK") + R("R_ACGS")  = 0.

Macmaltc: R("R_MALTATr") - R("R_MALTATr_Rev")  = 0.

Macmanab: R("Macmanab_Transp") + R("R_EX_acmana_e_") -
 R("R_EX_acmana_e__Rev")  = 0.

Macmanac: R("R_ACNML") - R("R_AMANK")  = 0.

Macmanae: -R("R_ACMANAtex") + R("R_ACMANAtex_Rev") -
 R("R_EX_acmana_e_") + R("R_EX_acmana_e__Rev")  = 0.

Macmanap: -R("R_ACMANAptspp") + R("R_ACMANAtex") -
 R("R_ACMANAtex_Rev")  = 0.

Macmanapc: R("R_ACMANAptspp") - R("R_AMANAPEr") +
 R("R_AMANAPEr_Rev") + R("R_AMANK")  = 0.

Macmum6pc: -R("R_ACM6PH") + R("R_ACMUMptspp") + R("R_ANHMK")  = 0.

Macmumb: R("Macmumb_Transp") + R("R_EX_acmum_e_") -
 R("R_EX_acmum_e__Rev")  = 0.

Macmume: -R("R_ACMUMtex") - R("R_EX_acmum_e_") +
 R("R_EX_acmum_e__Rev")  = 0.

Macmump: -R("R_ACMUMptspp") + R("R_ACMUMtex")  = 0.

Macnamb: R("Macnamb_Transp") + R("R_EX_acnam_e_") -
 R("R_EX_acnam_e__Rev")  = 0.

Macnamc: R("R_ACNAMt2pp") - R("R_ACNML")  = 0.

Macname: -R("R_ACNAMtex") + R("R_ACNAMtex_Rev") - R("R_EX_acnam_e_")
 + R("R_EX_acnam_e__Rev")  = 0.

Macnamp: -R("R_ACNAMt2pp") + R("R_ACNAMtex") - R("R_ACNAMtex_Rev") 
 = 0.

Macolipab: R("Macolipab_Transp") + R("R_EX_acolipa_e_") -
 R("R_EX_acolipa_e__Rev")  = 0.

Macolipae: R("R_ACOLIPAtex") - R("R_EX_acolipa_e_") +
 R("R_EX_acolipa_e__Rev")  = 0.

Macolipap: -R("R_ACOLIPAtex") + R("R_LA4NTpp")  = 0.

MaconCc: R("R_ACONIs") - R("R_ACONIs_Rev") + R("R_ACONTa") -
 R("R_ACONTa_Rev") - R("R_ACONTb") + R("R_ACONTb_Rev")  = 0.

Maconmc: R("R_ACONMT")  = 0.

MaconTc: -R("R_ACONIs") + R("R_ACONIs_Rev") - R("R_ACONMT")  = 0.

Macornc: -R("R_ACODA") - R("R_ACOTA") + R("R_ACOTA_Rev")  = 0.

Macp: -R("R_ACt2rpp") + R("R_ACt2rpp_Rev") - R("R_ACt4pp") +
 R("R_ACtex") - R("R_ACtex_Rev")  = 0.

MACPc: R("R_3OAS100") + R("R_3OAS120") + R("R_3OAS121") +
 R("R_3OAS140") + R("R_3OAS141") + R("R_3OAS160") + R("R_3OAS161") +
 R("R_3OAS180") + R("R_3OAS181") + R("R_3OAS60") + R("R_3OAS80") - R("R_AACPS1") -
 R("R_AACPS2") - R("R_AACPS3") - R("R_AACPS4") - R("R_AACPS5") -
 R("R_AACPS6") - R("R_AACPS7") - R("R_AACPS8") - R("R_AACPS9") -
 R("R_ACOATA") + R("R_ACOATA_Rev") + R("R_ACPS1") + R("R_AGPAT120") +
 R("R_AGPAT140") + R("R_AGPAT141") + R("R_AGPAT160") + R("R_AGPAT161") +
 R("R_AGPAT180") + R("R_AGPAT181") + R("R_EDTXS1") + R("R_EDTXS2") + R("R_EDTXS3")
 + R("R_EDTXS4") + R("R_FA100ACPHi") + R("R_FA120ACPHi") +
 R("R_FA140ACPHi") + R("R_FA141ACPHi") + R("R_FA160ACPHi") + R("R_FA161ACPHi") +
 R("R_FA80ACPHi") + R("R_G3PAT120") + R("R_G3PAT140") + R("R_G3PAT141")
 + R("R_G3PAT160") + R("R_G3PAT161") + R("R_G3PAT180") +
 R("R_G3PAT181") + R("R_KAS14") - R("R_MCOATA") + R("R_MCOATA_Rev") + R("R_U23GAAT")
 + R("R_UAGAAT") - R("R_UAGAAT_Rev")  = 0.

Macserb: R("Macserb_Transp") + R("R_EX_acser_e_") -
 R("R_EX_acser_e__Rev")  = 0.

Macserc: -R("R_ACSERtpp") - R("R_CYSS") + R("R_SERAT") -
 R("R_SERAT_Rev")  = 0.

Macsere: -R("R_ACSERtex") + R("R_ACSERtex_Rev") - R("R_EX_acser_e_")
 + R("R_EX_acser_e__Rev")  = 0.

Macserp: R("R_ACSERtex") - R("R_ACSERtex_Rev") + R("R_ACSERtpp")  =
 0.

MactACPc: -R("R_3OAR40") + R("R_3OAR40_Rev") + R("R_KAS14") +
 R("R_KAS15")  = 0.

Mactpc: R("R_ACKr") - R("R_ACKr_Rev") + R("R_PTAr") -
 R("R_PTAr_Rev")  = 0.

Madeb: R("Madeb_Transp") + R("R_EX_ade_e_") - R("R_EX_ade_e__Rev") 
 = 0.

Madec: R("R_5DOAN") - R("R_ADD") + R("R_ADEt2rpp") -
 R("R_ADEt2rpp_Rev") + R("R_ADNUC") - R("R_ADPT") + R("R_AHCYSNS") + R("R_AMPN") +
 R("R_MTAN") + R("R_PUNP1") - R("R_PUNP1_Rev") + R("R_PUNP2") -
 R("R_PUNP2_Rev") + R("R_TPRDCOAS")  = 0.

Madee: -R("R_ADEtex") + R("R_ADEtex_Rev") - R("R_EX_ade_e_") +
 R("R_EX_ade_e__Rev")  = 0.

Madep: -R("R_ADEt2rpp") + R("R_ADEt2rpp_Rev") + R("R_ADEtex") -
 R("R_ADEtex_Rev")  = 0.

Madnb: R("Madnb_Transp") + R("R_EX_adn_e_") - R("R_EX_adn_e__Rev") 
 = 0.

Madnc: -R("R_ADA") - R("R_ADNK1") + R("R_ADNt2pp") + R("R_ADNt2rpp")
 - R("R_ADNt2rpp_Rev") - R("R_ADNUC") + R("R_NTD7") - R("R_PUNP1") +
 R("R_PUNP1_Rev")  = 0.

Madne: -R("R_ADNtex") + R("R_ADNtex_Rev") - R("R_EX_adn_e_") +
 R("R_EX_adn_e__Rev")  = 0.

Madnp: R("R_3NTD7pp") - R("R_ADNt2pp") - R("R_ADNt2rpp") +
 R("R_ADNt2rpp_Rev") + R("R_ADNtex") - R("R_ADNtex_Rev") + R("R_NTD7pp")  = 0.

Madocbic: -R("R_ADOCBIK") + R("R_CBIAT") - R("R_CBIAT_Rev")  = 0.

Madocbipc: -R("R_ACBIPGT") + R("R_ADOCBIK")  = 0.

Madocblb: R("Madocblb_Transp") + R("R_EX_adocbl_e_") -
 R("R_EX_adocbl_e__Rev")  = 0.

Madocblc: R("R_ADOCBLabcpp") + R("R_ADOCBLS") + R("R_CBLAT") -
 R("R_CBLAT_Rev")  = 0.

Madocble: -R("R_ADOCBLtonex") - R("R_EX_adocbl_e_") +
 R("R_EX_adocbl_e__Rev")  = 0.

Madocblp: -R("R_ADOCBLabcpp") + R("R_ADOCBLtonex")  = 0.

Madpc: R("R_14GLUCANabcpp") + R("R_3KGK") + R("R_3PEPTabcpp") +
 R("R_4PEPTabcpp") + R("R_ACCOAC") + R("R_ACCOAL") + R("R_ACGAMK") +
 R("R_ACGK") + R("R_ACKr") - R("R_ACKr_Rev") + 2 R("R_ADK1") - 2
 R("R_ADK1_Rev") + R("R_ADK3") - R("R_ADK3_Rev") + R("R_ADK4") - R("R_ADK4_Rev") +
 R("R_ADNK1") + R("R_ADOCBIK") + R("R_ADOCBLabcpp") + R("R_ADSK") +
 R("R_AIRC2") + R("R_ALAabcpp") + R("R_ALAALAabcpp") + R("R_ALAALAr") -
 R("R_ALAALAr_Rev") + R("R_ALLabcpp") + R("R_ALLK") + R("R_AMANK") +
 R("R_ANHMK") + 2 R("R_AP4AH") + R("R_AP5AH") + R("R_ARBabcpp") +
 R("R_ARBTNabcpp") + R("R_ARGabcpp") + R("R_ASNabcpp") + R("R_ASO3t8pp") +
 R("R_ASPabcpp") + R("R_ASPK") - R("R_ASPK_Rev") + R("R_ATPM") -
 R("R_ATPS4rpp") + R("R_ATPS4rpp_Rev") + R("R_BUTSO3abcpp") + R("R_CBIuabcpp") +
 R("R_CBL1abcpp") + R("R_CBMKr") - R("R_CBMKr_Rev") + 2 R("R_CBPS") +
 R("R_CD2abcpp") + R("R_CDPMEK") + R("R_CGLYabcpp") + R("R_CHLabcpp") +
 R("R_COBALT2abcpp") + R("R_COLIPAabcpp") + R("R_CPGNabcpp") +
 R("R_CRNabcpp") + R("R_CRNCAL2") + R("R_CRNDabcpp") + R("R_CRNDCAL2") +
 R("R_CTBTabcpp") + R("R_CTBTCAL2") + R("R_CTPS2") + R("R_CU1abcpp") +
 R("R_CU2abcpp") + R("R_CYSabc2pp") + R("R_CYSabcpp") + R("R_CYSDabcpp") +
 R("R_CYTK1") - R("R_CYTK1_Rev") + R("R_CYTK2") - R("R_CYTK2_Rev") +
 R("R_DADK") - R("R_DADK_Rev") + R("R_DAGK120") + R("R_DAGK140") +
 R("R_DAGK141") + R("R_DAGK160") + R("R_DAGK161") + R("R_DAGK180") +
 R("R_DAGK181") + R("R_DAPabcpp") + R("R_DBTS") + R("R_DDGALK") + R("R_DDGLK") +
 R("R_DGK1") - R("R_DGK1_Rev") + R("R_DHFS") + R("R_DPCOAK") +
 R("R_DTMPK") - R("R_DTMPK_Rev") + R("R_DURIK1") + R("R_DXYLK") + 59.81
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_ETHSO3abcpp") + R("R_FCLK") +
 R("R_FE2abcpp") + R("R_FE3abcpp") + R("R_FE3DCITabcpp") +
 R("R_FE3DHBZSabcpp") + R("R_FE3HOXabcpp") + R("R_FECRMabcpp") + R("R_FEENTERabcpp") +
 R("R_FEOXAMabcpp") + R("R_FRUK") + R("R_FRULYSK") + R("R_G3PCabcpp")
 + R("R_G3PEabcpp") + R("R_G3PGabcpp") + R("R_G3PIabcpp") +
 R("R_G3PSabcpp") + R("R_GALabcpp") + R("R_GALKr") - R("R_GALKr_Rev") +
 R("R_GART") + R("R_GGPTRCS") + R("R_GK1") - R("R_GK1_Rev") + R("R_GLCabcpp") +
 R("R_GLCS1") + R("R_GLNabcpp") + R("R_GLNS") + R("R_GLU5K") +
 R("R_GLUabcpp") + R("R_GLUCYS") + R("R_GLYBabcpp") + R("R_GLYC3Pabcpp") +
 R("R_GLYCK") + R("R_GLYCK2") + R("R_GLYK") + R("R_GMHEPK") + R("R_GNK") +
 R("R_GSNK") + R("R_GSPMDS") + R("R_GTHRDabc2pp") + R("R_GTHRDabcpp") +
 R("R_GTHS") + R("R_HEPK1") + R("R_HEPK2") + R("R_HEPT1") +
 R("R_HEPT2") + R("R_HEPT3") + R("R_HEPT4") + R("R_HETZK") + R("R_HEX1") +
 R("R_HEX4") + R("R_HEX7") + R("R_HG2abcpp") + R("R_HISabcpp") + R("R_HMPK1")
 + R("R_HSK") + R("R_ILEabcpp") + R("R_INSK") + R("R_ISETACabcpp") +
 R("R_K2L4Aabcpp") + R("R_Kabcpp") + R("R_LEUabcpp") + R("R_LIPAabcpp") +
 R("R_LIPACabcpp") + R("R_LYSabcpp") + R("R_MALTabcpp") +
 R("R_MALTHXabcpp") + R("R_MALTPTabcpp") + R("R_MALTTRabcpp") + R("R_MALTTTRabcpp")
 + R("R_METabcpp") + R("R_METDabcpp") + R("R_METSOX1abcpp") +
 R("R_METSOX2abcpp") + R("R_MG2uabcpp") + R("R_MOBDabcpp") + R("R_MSO3abcpp") +
 R("R_MTRK") + R("R_NADK") + R("R_NAMNPP") + R("R_NDPK1") -
 R("R_NDPK1_Rev") + R("R_NDPK2") - R("R_NDPK2_Rev") + R("R_NDPK3") -
 R("R_NDPK3_Rev") + R("R_NDPK4") - R("R_NDPK4_Rev") + R("R_NDPK5") -
 R("R_NDPK5_Rev") + R("R_NDPK6") - R("R_NDPK6_Rev") + R("R_NDPK7") - R("R_NDPK7_Rev")
 + R("R_NDPK8") - R("R_NDPK8_Rev") + R("R_NI2abcpp") +
 R("R_NI2uabcpp") + R("R_NTP1") + 2 R("R_OMMBLHX3") + 2 R("R_OMPHHX3") + 2
 R("R_OPHHX3") + R("R_ORNabcpp") + R("R_PA120abcpp") + R("R_PA140abcpp") +
 R("R_PA141abcpp") + R("R_PA160abcpp") + R("R_PA161abcpp") + R("R_PA180abcpp")
 + R("R_PA181abcpp") + R("R_PE120abcpp") + R("R_PE140abcpp") +
 R("R_PE141abcpp") + R("R_PE160abcpp") + R("R_PE161abcpp") + R("R_PE180abcpp")
 + R("R_PE181abcpp") + R("R_PFK") + R("R_PFK_2") + R("R_PG120abcpp") +
 R("R_PG140abcpp") + R("R_PG141abcpp") + R("R_PG160abcpp") +
 R("R_PG161abcpp") + R("R_PG180abcpp") + R("R_PG181abcpp") + R("R_PGK") -
 R("R_PGK_Rev") + R("R_PGP120abcpp") + R("R_PGP140abcpp") +
 R("R_PGP141abcpp") + R("R_PGP160abcpp") + R("R_PGP161abcpp") + R("R_PGP180abcpp") +
 R("R_PGP181abcpp") + R("R_PHEMEabcpp") + R("R_PIuabcpp") + R("R_PMPK") +
 R("R_PNTK") - R("R_PPAKr") + R("R_PPAKr_Rev") + R("R_PPCK") +
 R("R_PPK2r") - R("R_PPK2r_Rev") + R("R_PPKr") - R("R_PPKr_Rev") +
 R("R_PRAGSr") - R("R_PRAGSr_Rev") + R("R_PRAIS") + R("R_PRASCSi") + R("R_PRFGS") +
 R("R_PROabcpp") + R("R_PROGLYabcpp") + R("R_PTRCabcpp") +
 R("R_PYDAMK") + R("R_PYDXK") + R("R_PYDXNK") - R("R_PYK") + R("R_R15BPK") +
 R("R_R1PK") + R("R_RBFK") + R("R_RBK") + R("R_RBK_L1") + R("R_RIBabcpp") +
 R("R_RMK") - R("R_RNDR1") - R("R_RNDR1b") + R("R_SHKK") +
 R("R_SPMDabcpp") + R("R_SUCOAS") - R("R_SUCOAS_Rev") + R("R_SULabcpp") +
 R("R_SULFACabcpp") + R("R_TAURabcpp") + R("R_TDSK") + R("R_THMabcpp") +
 R("R_THRabcpp") + R("R_TMDK1") + R("R_TMK") + R("R_TMPK") + R("R_TSULabcpp")
 + R("R_TUNGSabcpp") + R("R_UAAGDS") + R("R_UAMAGS") + R("R_UAMAS") +
 R("R_UGMDDS") + R("R_UM3PL") + R("R_UM4PL") + R("R_UMPK") -
 R("R_UMPK_Rev") + R("R_URIDK2r") - R("R_URIDK2r_Rev") + R("R_VALabcpp") +
 R("R_XYLabcpp") + R("R_XYLK") + R("R_XYLK2") + R("R_ZN2abcpp") +
 R("R_ZNabcpp")  = 0.

Madpglcc: -R("R_GLCS1") + R("R_GLGC")  = 0.

MadphepDDc: -R("R_AGMHE") + R("R_GMHEPAT")  = 0.

MadphepLDc: R("R_AGMHE") - R("R_HEPT1") - R("R_HEPT2") -
 R("R_HEPT3") - R("R_HEPT4")  = 0.

Madpribc: -R("R_ADPRDP") + R("R_NADN")  = 0.

Magb: R("Magb_Transp") + R("R_EX_ag_e_") - R("R_EX_ag_e__Rev")  = 0.

Magc: -R("R_AGt3")  = 0.

Magdpcbic: R("R_ACBIPGT") - R("R_ADOCBLS")  = 0.

Mage: R("R_AGt3") - R("R_EX_ag_e_") + R("R_EX_ag_e__Rev")  = 0.

Magmb: R("Magmb_Transp") + R("R_EX_agm_e_") - R("R_EX_agm_e__Rev") 
 = 0.

Magmc: -R("R_AGMT") - R("R_ARGAGMt7pp") + R("R_ARGAGMt7pp_Rev") +
 R("R_ARGDC")  = 0.

Magme: -R("R_AGMtex") + R("R_AGMtex_Rev") - R("R_EX_agm_e_") +
 R("R_EX_agm_e__Rev")  = 0.

Magmp: R("R_AGMtex") - R("R_AGMtex_Rev") + R("R_ARGAGMt7pp") -
 R("R_ARGAGMt7pp_Rev") + R("R_ARGDCpp")  = 0.

Mahcysc: R("R_ACONMT") - R("R_AHCYSNS") + R("R_AMMQLT8") + 2
 R("R_CFAS160E") + 2 R("R_CFAS160G") + 2 R("R_CFAS180E") + 2 R("R_CFAS180G") +
 R("R_DMQMT") + R("R_HCYSMT") + R("R_OHPHM") + R("R_OMBZLM") + 2
 R("R_UPP3MT")  = 0.

Mahdtc: -R("R_DHPTPE") + R("R_DHPTPE_Rev") - R("R_DNTPPA") +
 R("R_GTPCI")  = 0.

Maicarc: R("R_ADSL2r") - R("R_ADSL2r_Rev") - R("R_AICART") +
 R("R_AICART_Rev") + R("R_IG3PS")  = 0.

Mairc: -R("R_AIRC2") - R("R_AMPMS2") + R("R_PRAIS")  = 0.

Makgb: R("Makgb_Transp") + R("R_EX_akg_e_") - R("R_EX_akg_e__Rev") 
 = 0.

Makgc: -R("R_ABTA") - R("R_ACOTA") + R("R_ACOTA_Rev") - R("R_AKGDH")
 + R("R_AKGt2rpp") - R("R_AKGt2rpp_Rev") - R("R_ALATA_L") +
 R("R_ALATA_L_Rev") - R("R_ASPTA") + R("R_ASPTA_Rev") + R("R_GLUDy") -
 R("R_GLUDy_Rev") - R("R_GLUSy") + R("R_HSTPT") + R("R_ICDHyr") -
 R("R_ICDHyr_Rev") - R("R_ILETA") + R("R_ILETA_Rev") + R("R_LEUTAi") + R("R_OHPBAT") -
 R("R_OHPBAT_Rev") - R("R_OXGDC2") - R("R_PHETA1") + R("R_PHETA1_Rev")
 + R("R_PSERT") - R("R_PTRCTA") - R("R_SDPTA") + R("R_SDPTA_Rev") -
 R("R_SOTA") - R("R_TAUDO") + R("R_TDPAGTA") - R("R_TYRTA") +
 R("R_TYRTA_Rev") + R("R_UDPKAAT") - R("R_UDPKAAT_Rev") + R("R_UNK3") -
 R("R_VALTA") + R("R_VALTA_Rev")  = 0.

Makge: -R("R_AKGtex") + R("R_AKGtex_Rev") - R("R_EX_akg_e_") +
 R("R_EX_akg_e__Rev")  = 0.

Makgp: -R("R_AKGt2rpp") + R("R_AKGt2rpp_Rev") + R("R_AKGtex") -
 R("R_AKGtex_Rev")  = 0.

Malaalab: R("Malaalab_Transp") + R("R_EX_alaala_e_") -
 R("R_EX_alaala_e__Rev")  = 0.

Malaalac: R("R_ALAALAabcpp") - R("R_ALAALAD") + R("R_ALAALAr") -
 R("R_ALAALAr_Rev") - R("R_UGMDDS")  = 0.

Malaalae: -R("R_ALAALAtex") + R("R_ALAALAtex_Rev") -
 R("R_EX_alaala_e_") + R("R_EX_alaala_e__Rev")  = 0.

Malaalap: -R("R_ALAALAabcpp") + R("R_ALAALAtex") -
 R("R_ALAALAtex_Rev") + R("R_MCTP1Bpp") + R("R_MLDCP1App") + R("R_MLDCP2App") +
 R("R_MLDCP3App")  = 0.

MalaBb: R("MalaBb_Transp") + R("R_EX_ala_B_e_") -
 R("R_EX_ala_B_e__Rev")  = 0.

MalaBc: R("R_ASP1DC") + R("R_BALAt2pp") - R("R_PANTS")  = 0.

MalaBe: -R("R_BALAtex") + R("R_BALAtex_Rev") - R("R_EX_ala_B_e_") +
 R("R_EX_ala_B_e__Rev")  = 0.

MalaBp: -R("R_BALAt2pp") + R("R_BALAtex") - R("R_BALAtex_Rev")  = 0.

MalacSc: R("R_ACLS") + R("R_KARA1") - R("R_KARA1_Rev")  = 0.

MalaDb: R("MalaDb_Transp") + R("R_EX_ala_D_e_") -
 R("R_EX_ala_D_e__Rev")  = 0.

MalaDc: R("R_4PCP") + R("R_AGM4PCP") + 2 R("R_ALAALAD") - 2
 R("R_ALAALAr") + 2 R("R_ALAALAr_Rev") + R("R_ALAR") - R("R_ALAR_Rev") -
 R("R_ALATA_D2") + R("R_AM4PCP") - R("R_DAAD") + R("R_DALAt2pp") +
 R("R_UM4PCP")  = 0.

MalaDe: -R("R_DALAtex") + R("R_DALAtex_Rev") - R("R_EX_ala_D_e_") +
 R("R_EX_ala_D_e__Rev")  = 0.

MalaDp: R("R_4PCPpp") + R("R_AGM4PCPpp") - R("R_DALAt2pp") +
 R("R_DALAtex") - R("R_DALAtex_Rev") + R("R_MCTP1App") + 2 R("R_MCTP2App") +
 R("R_MDDCP1pp") + R("R_MDDCP2pp") + R("R_MDDCP3pp") + R("R_MDDCP4pp") +
 R("R_MDDCP5pp") + R("R_MLDCP1Bpp") + R("R_MLDCP2Bpp")  = 0.

MalaLb: R("MalaLb_Transp") + R("R_EX_ala_L_e_") -
 R("R_EX_ala_L_e__Rev")  = 0.

MalaLc: R("R_ALAabcpp") - R("R_ALAR") + R("R_ALAR_Rev") +
 R("R_ALAt2pp") + R("R_ALAt4pp") - R("R_ALATA_L") - R("R_ALATA_L2") +
 R("R_ALATA_L_Rev") - R("R_ALATRS") - R("R_AOXSr") + R("R_AOXSr_Rev") +
 R("R_CYSSADS") - 0.5137 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_LALGP") +
 R("R_THZPSN") - R("R_UAMAS") - R("R_VPAMT")  = 0.

MalaLe: -R("R_ALAtex") + R("R_ALAtex_Rev") - R("R_EX_ala_L_e_") +
 R("R_EX_ala_L_e__Rev")  = 0.

MalaLp: -R("R_ALAabcpp") - R("R_ALAt2pp") - R("R_ALAt4pp") +
 R("R_ALAtex") - R("R_ALAtex_Rev")  = 0.

Malatrnac: R("R_ALATRS")  = 0.

Mall6pc: R("R_ALLK") - R("R_ALLPI") + R("R_ALLPI_Rev")  = 0.

MallDb: R("MallDb_Transp") + R("R_EX_all_D_e_") -
 R("R_EX_all_D_e__Rev")  = 0.

MallDc: R("R_ALLabcpp") - R("R_ALLK")  = 0.

MallDe: -R("R_ALLtex") + R("R_ALLtex_Rev") - R("R_EX_all_D_e_") +
 R("R_EX_all_D_e__Rev")  = 0.

MallDp: -R("R_ALLabcpp") + R("R_ALLtex") - R("R_ALLtex_Rev")  = 0.

Malltnb: R("Malltnb_Transp") + R("R_EX_alltn_e_") -
 R("R_EX_alltn_e__Rev")  = 0.

Malltnc: -R("R_ALLTN") + R("R_ALLTNt2rpp") - R("R_ALLTNt2rpp_Rev") +
 R("R_URIC")  = 0.

Malltne: -R("R_ALLTNtex") + R("R_ALLTNtex_Rev") - R("R_EX_alltn_e_")
 + R("R_EX_alltn_e__Rev")  = 0.

Malltnp: -R("R_ALLTNt2rpp") + R("R_ALLTNt2rpp_Rev") +
 R("R_ALLTNtex") - R("R_ALLTNtex_Rev")  = 0.

Mallttc: -R("R_ALLTAMH") + R("R_ALLTN")  = 0.

Mallul6pc: R("R_ALLPI") - R("R_ALLPI_Rev") - R("R_ALLULPE") +
 R("R_ALLULPE_Rev")  = 0.

Malppp: -R("R_ALPATE160pp") - R("R_ALPATG160pp")  = 0.

Maltrnc: -R("R_ALTRH") - R("R_TAGURr") + R("R_TAGURr_Rev")  = 0.

Mametamc: R("R_ADMDC") - R("R_SPMS")  = 0.

Mametc: -R("R_ACONMT") - R("R_ADMDC") - R("R_AMAOTr") +
 R("R_AMAOTr_Rev") - R("R_AMMQLT8") - R("R_BTS4") - 2 R("R_CFAS160E") - 2
 R("R_CFAS160G") - 2 R("R_CFAS180E") - 2 R("R_CFAS180G") - 2 R("R_CPPPGO2") -
 R("R_DMQMT") - 0.000223 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_HCYSMT") + R("R_METAT") - R("R_OHPHM") - R("R_OMBZLM") - 2 R("R_UPP3MT") 
 = 0.

Mamobc: R("R_AMAOTr") - R("R_AMAOTr_Rev")  = 0.

Mampb: R("Mampb_Transp") + R("R_EX_amp_e_") - R("R_EX_amp_e__Rev") 
 = 0.

Mampc: R("R_2AGPEAT120") + R("R_2AGPEAT140") + R("R_2AGPEAT141") +
 R("R_2AGPEAT160") + R("R_2AGPEAT161") + R("R_2AGPEAT180") +
 R("R_2AGPEAT181") + R("R_2AGPGAT120") + R("R_2AGPGAT140") + R("R_2AGPGAT141") +
 R("R_2AGPGAT160") + R("R_2AGPGAT161") + R("R_2AGPGAT180") +
 R("R_2AGPGAT181") + R("R_AACPS1") + R("R_AACPS2") + R("R_AACPS3") + R("R_AACPS4")
 + R("R_AACPS5") + R("R_AACPS6") + R("R_AACPS7") + R("R_AACPS8") +
 R("R_AACPS9") + R("R_ACS") - R("R_ADK1") + R("R_ADK1_Rev") - R("R_ADK3")
 + R("R_ADK3_Rev") - R("R_ADK4") + R("R_ADK4_Rev") + R("R_ADNK1") +
 R("R_ADPRDP") + R("R_ADPT") + R("R_ADSL1r") - R("R_ADSL1r_Rev") +
 R("R_ALATRS") - R("R_AMPN") + R("R_ARGSS") + R("R_ARGTRS") + R("R_ASNS1") +
 R("R_ASNS2") + R("R_ASNTRS") + R("R_ASPTRS") + R("R_BPNT") +
 R("R_CYSTRS") + 6 R("R_ENTCS") + R("R_FACOAL100t2pp") + R("R_FACOAL120t2pp") +
 R("R_FACOAL140t2pp") + R("R_FACOAL141t2pp") + R("R_FACOAL160t2pp") +
 R("R_FACOAL161t2pp") + R("R_FACOAL180t2pp") + R("R_FACOAL181t2pp") +
 R("R_FACOAL60t2pp") + R("R_FACOAL80t2pp") + R("R_GDPDPK") + R("R_GLNTRS")
 + R("R_GLUTRS") + R("R_GLYTRS") + R("R_GMPS2") + R("R_GTPDPK") +
 R("R_HISTRS") + R("R_HPPK2") + R("R_ILETRS") + R("R_LEUTRS") +
 R("R_LYSTRS") + R("R_METTRS") + R("R_NADDP") + R("R_NADS1") - R("R_NTD7") +
 R("R_NTPP6") + R("R_PACCOAL") + R("R_PANTS") + R("R_PDE1") + R("R_PHETRS")
 + R("R_PPS") + R("R_PROTRS") + R("R_PRPPS") - R("R_PRPPS_Rev") +
 R("R_SELNPS") + R("R_SERTRS") + R("R_SERTRS2") + R("R_SUCBZL") +
 R("R_THRTRS") + R("R_THZPSN") + R("R_TRPTRS") + R("R_TYRTRS") + R("R_VALTRS")  =
 0.

Mampe: -R("R_AMPtex") + R("R_AMPtex_Rev") - R("R_EX_amp_e_") +
 R("R_EX_amp_e__Rev")  = 0.

Mampp: R("R_AMPtex") - R("R_AMPtex_Rev") - R("R_NTD7pp")  = 0.

Manhgm3pc: -R("R_AGM3PA") - R("R_AGM3PH") + R("R_AGM3Pt2pp") +
 R("R_AGM4PCP")  = 0.

Manhgm3pp: -R("R_AGM3PApp") - R("R_AGM3Pt2pp") + R("R_AGM4PCPpp") +
 R("R_MLTGY2pp") + 2 R("R_MLTGY3pp")  = 0.

Manhgm4pc: -R("R_AGM4PA") - R("R_AGM4PCP") - R("R_AGM4PH") +
 R("R_AGM4Pt2pp")  = 0.

Manhgm4pp: -R("R_AGM4PApp") - R("R_AGM4PCPpp") - R("R_AGM4Pt2pp") +
 2 R("R_MLTGY1pp") + R("R_MLTGY2pp") + R("R_MLTGY4pp")  = 0.

Manhgmb: R("Manhgmb_Transp") + R("R_EX_anhgm_e_") -
 R("R_EX_anhgm_e__Rev")  = 0.

Manhgmc: R("R_AGM3PA") + R("R_AGM4PA") - R("R_AGMH") +
 R("R_AGMt2pp")  = 0.

Manhgme: -R("R_ANHGMtex") + R("R_ANHGMtex_Rev") - R("R_EX_anhgm_e_")
 + R("R_EX_anhgm_e__Rev")  = 0.

Manhgmp: R("R_AGM3PApp") + R("R_AGM4PApp") - R("R_AGMt2pp") +
 R("R_ANHGMtex") - R("R_ANHGMtex_Rev")  = 0.

Manhm3pc: R("R_AGM3PH") - R("R_AM3PA") + R("R_AM4PCP")  = 0.

Manhm4pc: R("R_AGM4PH") - R("R_AM4PA") - R("R_AM4PCP")  = 0.

Manhmc: R("R_AGMH") + R("R_AM3PA") + R("R_AM4PA") - R("R_ANHMK")  =
 0.

Manthc: -R("R_ACANTHAT") - R("R_ANPRT") + R("R_ANS")  = 0.

Map4ac: -R("R_AP4AH")  = 0.

Map5ac: -R("R_AP5AH")  = 0.

Mapg120c: R("R_LPLIPAL2ATE120") + R("R_LPLIPAL2ATG120")  = 0.

Mapg140c: R("R_LPLIPAL2ATE140") + R("R_LPLIPAL2ATG140")  = 0.

Mapg141c: R("R_LPLIPAL2ATE141") + R("R_LPLIPAL2ATG141")  = 0.

Mapg160c: R("R_LPLIPAL2ATE160") + R("R_LPLIPAL2ATG160")  = 0.

Mapg161c: R("R_LPLIPAL2ATE161") + R("R_LPLIPAL2ATG161")  = 0.

Mapg180c: R("R_LPLIPAL2ATE180") + R("R_LPLIPAL2ATG180")  = 0.

Mapg181c: R("R_LPLIPAL2ATE181") + R("R_LPLIPAL2ATG181")  = 0.

MapoACPc: -R("R_ACPS1")  = 0.

Mapsc: -R("R_ADSK") + R("R_SADT2")  = 0.

Mara5pc: R("R_A5PISO") - R("R_A5PISO_Rev") - R("R_KDOPS")  = 0.

MarabLb: R("MarabLb_Transp") + R("R_EX_arab_L_e_") -
 R("R_EX_arab_L_e__Rev")  = 0.

MarabLc: -R("R_ARAI") + R("R_ARAI_Rev") + R("R_ARBabcpp") +
 R("R_ARBt2rpp") - R("R_ARBt2rpp_Rev") - R("R_ARBt3ipp")  = 0.

MarabLe: -R("R_ARBtex") + R("R_ARBtex_Rev") - R("R_EX_arab_L_e_") +
 R("R_EX_arab_L_e__Rev")  = 0.

MarabLp: -R("R_ARBabcpp") - R("R_ARBt2rpp") + R("R_ARBt2rpp_Rev") +
 R("R_ARBt3ipp") + R("R_ARBtex") - R("R_ARBtex_Rev")  = 0.

Maragundc: R("R_O16AT") - R("R_O16GLCT1")  = 0.

Marbt6pc: -R("R_AB6PGH")  = 0.

Marbtnb: R("Marbtnb_Transp") + R("R_EX_arbtn_e_") -
 R("R_EX_arbtn_e__Rev")  = 0.

Marbtnc: 2 R("R_ARBTNR1") + 2 R("R_ARBTNR2") + 2 R("R_ARBTNR3") -
 R("R_ARBTNtpp")  = 0.

Marbtne: -R("R_ARBTNexs") + R("R_ARBTNtex") - R("R_EX_arbtn_e_") +
 R("R_EX_arbtn_e__Rev")  = 0.

Marbtnfe3b: R("Marbtnfe3b_Transp") + R("R_EX_arbtn_fe3_e_") -
 R("R_EX_arbtn_fe3_e__Rev")  = 0.

Marbtnfe3c: R("R_ARBTNabcpp") - 2 R("R_ARBTNR1") - 2 R("R_ARBTNR2")
 - 2 R("R_ARBTNR3")  = 0.

Marbtnfe3e: R("R_ARBTNexs") - R("R_ARBTNtonex") -
 R("R_EX_arbtn_fe3_e_") + R("R_EX_arbtn_fe3_e__Rev")  = 0.

Marbtnfe3p: -R("R_ARBTNabcpp") + R("R_ARBTNtonex")  = 0.

Marbtnp: -R("R_ARBTNtex") + R("R_ARBTNtpp")  = 0.

MargLb: R("MargLb_Transp") + R("R_EX_arg_L_e_") -
 R("R_EX_arg_L_e__Rev")  = 0.

MargLc: R("R_ARGabcpp") + R("R_ARGAGMt7pp") - R("R_ARGAGMt7pp_Rev")
 - R("R_ARGDC") + R("R_ARGORNt7pp") - R("R_ARGORNt7pp_Rev") +
 R("R_ARGSL") - R("R_ARGSL_Rev") - R("R_ARGt3pp") - R("R_ARGTRS") - R("R_AST") -
 0.2958 R("R_Ec_biomass_iAF1260_core_59p81M")  = 0.

MargLe: -R("R_ARGtex") + R("R_ARGtex_Rev") - R("R_EX_arg_L_e_") +
 R("R_EX_arg_L_e__Rev")  = 0.

MargLp: -R("R_ARGabcpp") - R("R_ARGAGMt7pp") + R("R_ARGAGMt7pp_Rev")
 - R("R_ARGDCpp") - R("R_ARGORNt7pp") + R("R_ARGORNt7pp_Rev") +
 R("R_ARGt3pp") + R("R_ARGtex") - R("R_ARGtex_Rev")  = 0.

Margsucc: -R("R_ARGSL") + R("R_ARGSL_Rev") + R("R_ARGSS")  = 0.

Margtrnac: R("R_ARGTRS")  = 0.

Mascb6pc: -R("R_ASCBPL") + R("R_ASCBptspp")  = 0.

MascbLb: R("MascbLb_Transp") + R("R_EX_ascb_L_e_") -
 R("R_EX_ascb_L_e__Rev")  = 0.

MascbLe: -R("R_ASCBtex") + R("R_ASCBtex_Rev") - R("R_EX_ascb_L_e_")
 + R("R_EX_ascb_L_e__Rev")  = 0.

MascbLp: -R("R_ASCBptspp") + R("R_ASCBtex") - R("R_ASCBtex_Rev")  =
 0.

MasnLb: R("MasnLb_Transp") + R("R_EX_asn_L_e_") -
 R("R_EX_asn_L_e__Rev")  = 0.

MasnLc: R("R_ASNabcpp") - R("R_ASNN") + R("R_ASNS1") + R("R_ASNS2")
 + R("R_ASNt2rpp") - R("R_ASNt2rpp_Rev") - R("R_ASNTRS") - 0.2411
 R("R_Ec_biomass_iAF1260_core_59p81M")  = 0.

MasnLe: -R("R_ASNtex") + R("R_ASNtex_Rev") - R("R_EX_asn_L_e_") +
 R("R_EX_asn_L_e__Rev")  = 0.

MasnLp: -R("R_ASNabcpp") - R("R_ASNNpp") - R("R_ASNt2rpp") +
 R("R_ASNt2rpp_Rev") + R("R_ASNtex") - R("R_ASNtex_Rev")  = 0.

Masntrnac: R("R_ASNTRS")  = 0.

Maso3b: R("Maso3b_Transp") + R("R_EX_aso3_e_") -
 R("R_EX_aso3_e__Rev")  = 0.

Maso3c: -R("R_ASO3t8pp") + R("R_ASR")  = 0.

Maso3e: -R("R_ASO3tex") + R("R_ASO3tex_Rev") - R("R_EX_aso3_e_") +
 R("R_EX_aso3_e__Rev")  = 0.

Maso3p: R("R_ASO3t8pp") + R("R_ASO3tex") - R("R_ASO3tex_Rev")  = 0.

Maso4c: -R("R_ASR")  = 0.

MaspLb: R("MaspLb_Transp") + R("R_EX_asp_L_e_") -
 R("R_EX_asp_L_e__Rev")  = 0.

MaspLc: -R("R_ADSS") - R("R_ARGSS") + R("R_ASNN") - R("R_ASNS1") -
 R("R_ASNS2") - R("R_ASP1DC") + R("R_ASPabcpp") - R("R_ASPCT") -
 R("R_ASPK") + R("R_ASPK_Rev") - R("R_ASPO3") - R("R_ASPO4") - R("R_ASPO5") -
 R("R_ASPO6") - R("R_ASPT") + R("R_ASPt2_2pp") + R("R_ASPt2_3pp") +
 R("R_ASPt2pp") - R("R_ASPTA") + R("R_ASPTA_Rev") - R("R_ASPTRS") - 0.2411
 R("R_Ec_biomass_iAF1260_core_59p81M") - R("R_PRASCSi")  = 0.

MaspLe: -R("R_ASPtex") + R("R_ASPtex_Rev") - R("R_EX_asp_L_e_") +
 R("R_EX_asp_L_e__Rev")  = 0.

MaspLp: R("R_ASNNpp") - R("R_ASPabcpp") - R("R_ASPt2_2pp") -
 R("R_ASPt2_3pp") - R("R_ASPt2pp") + R("R_ASPtex") - R("R_ASPtex_Rev")  = 0.

Maspsac: -R("R_ASAD") + R("R_ASAD_Rev") - R("R_DHDPS") + R("R_HSDy")
 - R("R_HSDy_Rev")  = 0.

Masptrnac: R("R_ASPTRS")  = 0.

MathrLc: -R("R_ATHRDHr") + R("R_ATHRDHr_Rev") - R("R_THRA2i")  = 0.

Matpc: -R("R_14GLUCANabcpp") - R("R_2AGPEAT120") - R("R_2AGPEAT140")
 - R("R_2AGPEAT141") - R("R_2AGPEAT160") - R("R_2AGPEAT161") -
 R("R_2AGPEAT180") - R("R_2AGPEAT181") - R("R_2AGPGAT120") - R("R_2AGPGAT140")
 - R("R_2AGPGAT141") - R("R_2AGPGAT160") - R("R_2AGPGAT161") -
 R("R_2AGPGAT180") - R("R_2AGPGAT181") - R("R_3KGK") - R("R_3PEPTabcpp") -
 R("R_4PEPTabcpp") - R("R_AACPS1") - R("R_AACPS2") - R("R_AACPS3") -
 R("R_AACPS4") - R("R_AACPS5") - R("R_AACPS6") - R("R_AACPS7") -
 R("R_AACPS8") - R("R_AACPS9") - R("R_ACCOAC") - R("R_ACCOAL") - R("R_ACGAMK") -
 R("R_ACGK") - R("R_ACKr") + R("R_ACKr_Rev") - R("R_ACS") - R("R_ADK1")
 + R("R_ADK1_Rev") - R("R_ADNCYC") - R("R_ADNK1") - R("R_ADOCBIK") -
 R("R_ADOCBLabcpp") - R("R_ADSK") - R("R_AIRC2") - R("R_ALAabcpp") -
 R("R_ALAALAabcpp") - R("R_ALAALAr") + R("R_ALAALAr_Rev") - R("R_ALATRS") -
 R("R_ALLabcpp") - R("R_ALLK") - R("R_AMANK") - R("R_ANHMK") +
 R("R_AP5AH") - R("R_ARBabcpp") - R("R_ARBTNabcpp") - R("R_ARGabcpp") -
 R("R_ARGSS") - R("R_ARGTRS") - R("R_ASNabcpp") - R("R_ASNS1") - R("R_ASNS2")
 - R("R_ASNTRS") - R("R_ASO3t8pp") - R("R_ASPabcpp") - R("R_ASPK") +
 R("R_ASPK_Rev") - R("R_ASPTRS") - R("R_ATPHs") - R("R_ATPM") -
 R("R_ATPPRT") + R("R_ATPS4rpp") - R("R_ATPS4rpp_Rev") - R("R_BUTSO3abcpp") -
 R("R_CBIAT") + R("R_CBIAT_Rev") - R("R_CBIuabcpp") - R("R_CBL1abcpp") -
 R("R_CBLAT") + R("R_CBLAT_Rev") - R("R_CBMKr") + R("R_CBMKr_Rev") - 2
 R("R_CBPS") - R("R_CD2abcpp") - R("R_CDPMEK") - R("R_CGLYabcpp") -
 R("R_CHLabcpp") - R("R_COBALT2abcpp") - R("R_COLIPAabcpp") -
 R("R_CPGNabcpp") - R("R_CRNabcpp") - R("R_CRNCAL2") - R("R_CRNDabcpp") -
 R("R_CRNDCAL2") - R("R_CTBTabcpp") - R("R_CTBTCAL2") - R("R_CTPS2") -
 R("R_CU1abcpp") - R("R_CU2abcpp") - R("R_CYSabc2pp") - R("R_CYSabcpp") -
 R("R_CYSDabcpp") - R("R_CYSTRS") - R("R_CYTK1") + R("R_CYTK1_Rev") -
 R("R_CYTK2") + R("R_CYTK2_Rev") - R("R_DADK") + R("R_DADK_Rev") -
 R("R_DAGK120") - R("R_DAGK140") - R("R_DAGK141") - R("R_DAGK160") -
 R("R_DAGK161") - R("R_DAGK180") - R("R_DAGK181") - R("R_DAPabcpp") - R("R_DBTS") -
 R("R_DDGALK") - R("R_DDGLK") - R("R_DGK1") + R("R_DGK1_Rev") -
 R("R_DHBS") - R("R_DHFS") - R("R_DPCOAK") - R("R_DTMPK") + R("R_DTMPK_Rev")
 - R("R_DURIK1") - R("R_DXYLK") - 59.984
 R("R_Ec_biomass_iAF1260_core_59p81M") - R("R_ETHSO3abcpp") - R("R_FACOAL100t2pp") -
 R("R_FACOAL120t2pp") - R("R_FACOAL140t2pp") - R("R_FACOAL141t2pp") -
 R("R_FACOAL160t2pp") - R("R_FACOAL161t2pp") - R("R_FACOAL180t2pp") -
 R("R_FACOAL181t2pp") - R("R_FACOAL60t2pp") - R("R_FACOAL80t2pp") - R("R_FCLK") -
 R("R_FE2abcpp") - R("R_FE3abcpp") - R("R_FE3DCITabcpp") -
 R("R_FE3DHBZSabcpp") - R("R_FE3HOXabcpp") - R("R_FECRMabcpp") - R("R_FEENTERabcpp") -
 R("R_FEOXAMabcpp") - R("R_FMNAT") - R("R_FRUK") - R("R_FRULYSK") -
 R("R_G3PCabcpp") - R("R_G3PEabcpp") - R("R_G3PGabcpp") - R("R_G3PIabcpp") -
 R("R_G3PSabcpp") - R("R_GALabcpp") - R("R_GALKr") + R("R_GALKr_Rev") -
 R("R_GART") - R("R_GDPDPK") - R("R_GGPTRCS") - R("R_GK1") +
 R("R_GK1_Rev") - R("R_GLCabcpp") - R("R_GLGC") - R("R_GLNabcpp") - R("R_GLNS")
 - R("R_GLNTRS") - R("R_GLU5K") - R("R_GLUabcpp") - R("R_GLUCYS") -
 R("R_GLUTRS") - R("R_GLYBabcpp") - R("R_GLYC3Pabcpp") - R("R_GLYCK") -
 R("R_GLYCK2") - R("R_GLYK") - R("R_GLYTRS") - R("R_GMHEPAT") -
 R("R_GMHEPK") - R("R_GMPS2") - R("R_GNK") - R("R_GSNK") - R("R_GSPMDS") -
 R("R_GTHRDabc2pp") - R("R_GTHRDabcpp") - R("R_GTHS") - R("R_GTPDPK") -
 R("R_HEPK1") - R("R_HEPK2") - R("R_HETZK") - R("R_HEX1") - R("R_HEX4") -
 R("R_HEX7") - R("R_HG2abcpp") - R("R_HISabcpp") - R("R_HISTRS") -
 R("R_HMPK1") - R("R_HPPK2") - R("R_HSK") - R("R_ILEabcpp") - R("R_ILETRS")
 - R("R_INSK") - R("R_ISETACabcpp") - R("R_K2L4Aabcpp") - R("R_Kabcpp")
 - R("R_LEUabcpp") - R("R_LEUTRS") - R("R_LIPAabcpp") -
 R("R_LIPACabcpp") - R("R_LYSabcpp") - R("R_LYSTRS") - R("R_MALTabcpp") -
 R("R_MALTHXabcpp") - R("R_MALTPTabcpp") - R("R_MALTTRabcpp") -
 R("R_MALTTTRabcpp") - R("R_METabcpp") - R("R_METAT") - R("R_METDabcpp") -
 R("R_METSOX1abcpp") - R("R_METSOX2abcpp") - R("R_METTRS") - R("R_MG2uabcpp") -
 R("R_MOBDabcpp") - R("R_MSO3abcpp") - R("R_MTRK") - R("R_NADK") -
 R("R_NADS1") - R("R_NAMNPP") - R("R_NDPK1") + R("R_NDPK1_Rev") - R("R_NDPK2") +
 R("R_NDPK2_Rev") - R("R_NDPK3") + R("R_NDPK3_Rev") - R("R_NDPK4") +
 R("R_NDPK4_Rev") - R("R_NDPK5") + R("R_NDPK5_Rev") - R("R_NDPK6") +
 R("R_NDPK6_Rev") - R("R_NDPK7") + R("R_NDPK7_Rev") - R("R_NDPK8") +
 R("R_NDPK8_Rev") - R("R_NI2abcpp") - R("R_NI2uabcpp") - R("R_NMNAT") -
 R("R_NNATr") + R("R_NNATr_Rev") - R("R_NTP1") - R("R_NTPP6") - 2
 R("R_OMMBLHX3") - 2 R("R_OMPHHX3") - 2 R("R_OPHHX3") - R("R_ORNabcpp") -
 R("R_PA120abcpp") - R("R_PA140abcpp") - R("R_PA141abcpp") -
 R("R_PA160abcpp") - R("R_PA161abcpp") - R("R_PA180abcpp") - R("R_PA181abcpp") -
 R("R_PACCOAL") - R("R_PANTS") - R("R_PE120abcpp") - R("R_PE140abcpp") -
 R("R_PE141abcpp") - R("R_PE160abcpp") - R("R_PE161abcpp") -
 R("R_PE180abcpp") - R("R_PE181abcpp") - R("R_PFK") - R("R_PFK_2") -
 R("R_PG120abcpp") - R("R_PG140abcpp") - R("R_PG141abcpp") - R("R_PG160abcpp") -
 R("R_PG161abcpp") - R("R_PG180abcpp") - R("R_PG181abcpp") - R("R_PGK") +
 R("R_PGK_Rev") - R("R_PGP120abcpp") - R("R_PGP140abcpp") -
 R("R_PGP141abcpp") - R("R_PGP160abcpp") - R("R_PGP161abcpp") - R("R_PGP180abcpp") -
 R("R_PGP181abcpp") - R("R_PHEMEabcpp") - R("R_PHETRS") -
 R("R_PIuabcpp") - R("R_PMPK") - R("R_PNTK") + R("R_PPAKr") - R("R_PPAKr_Rev") -
 R("R_PPCK") - R("R_PPK2r") + R("R_PPK2r_Rev") - R("R_PPKr") +
 R("R_PPKr_Rev") - R("R_PPS") - R("R_PRAGSr") + R("R_PRAGSr_Rev") - R("R_PRAIS") -
 R("R_PRASCSi") - R("R_PRFGS") - R("R_PROabcpp") - R("R_PROGLYabcpp")
 - R("R_PROTRS") - R("R_PRPPS") + R("R_PRPPS_Rev") - R("R_PTPATi") -
 R("R_PTRCabcpp") - R("R_PYDAMK") - R("R_PYDXK") - R("R_PYDXNK") +
 R("R_PYK") - R("R_R15BPK") - R("R_R1PK") - R("R_RBFK") - R("R_RBK") -
 R("R_RBK_L1") - R("R_RIBabcpp") - R("R_RMK") - R("R_RNTR1c") - R("R_SADT2")
 - R("R_SELNPS") - R("R_SERASr") + R("R_SERASr_Rev") - R("R_SERTRS") -
 R("R_SERTRS2") - R("R_SHKK") - R("R_SPMDabcpp") - R("R_SUCBZL") -
 R("R_SUCOAS") + R("R_SUCOAS_Rev") - R("R_SULabcpp") - R("R_SULFACabcpp") -
 R("R_TAURabcpp") - R("R_TDSK") - R("R_THMabcpp") - R("R_THRabcpp") -
 R("R_THRTRS") - R("R_THZPSN") - R("R_TMDK1") - R("R_TMK") -
 R("R_TMPK") - R("R_TPRDCOAS") - R("R_TRPTRS") - R("R_TSULabcpp") -
 R("R_TUNGSabcpp") - R("R_TYRTRS") - R("R_UAAGDS") - R("R_UAMAGS") - R("R_UAMAS") -
 R("R_UGMDDS") - R("R_UM3PL") - R("R_UM4PL") - R("R_UMPK") +
 R("R_UMPK_Rev") - R("R_URIDK2r") + R("R_URIDK2r_Rev") - R("R_VALabcpp") -
 R("R_VALTRS") - R("R_XYLabcpp") - R("R_XYLK") - R("R_XYLK2") -
 R("R_ZN2abcpp") - R("R_ZNabcpp")  = 0.

Mb2coac: R("R_ACOAD1f") - R("R_ACOAD1f_Rev") + R("R_ECOAH1") -
 R("R_ECOAH1_Rev")  = 0.

Mbbtcoac: -R("R_CRNBTCT") + R("R_CRNBTCT_Rev")  = 0.

Mbetaldc: -R("R_BETALDHx") - R("R_BETALDHy")  = 0.

Mbglycogenc: R("R_GLBRAN2") - R("R_GLCP2") - R("R_GLDBRAN2")  = 0.

Mbtalc: -R("R_ALDD4") + R("R_FDMO4")  = 0.

Mbtcoac: -R("R_ACOAD1f") + R("R_ACOAD1f_Rev") + R("R_BUTCT") +
 R("R_KAT2")  = 0.

Mbtnc: R("R_BSORx") + R("R_BSORy") + R("R_BTS4")  = 0.

Mbtnsoc: -R("R_BSORx") - R("R_BSORy")  = 0.

Mbut2eACPc: R("R_3HAD40") - R("R_EAR40x") - R("R_EAR40y")  = 0.

MbutACPc: -R("R_3OAS60") + R("R_EAR40x") + R("R_EAR40y")  = 0.

Mbutb: R("Mbutb_Transp") + R("R_EX_but_e_") - R("R_EX_but_e__Rev") 
 = 0.

Mbutc: R("R_ALDD4") - R("R_BUTCT") + R("R_BUTt2rpp") -
 R("R_BUTt2rpp_Rev")  = 0.

Mbute: -R("R_BUTtex") + R("R_BUTtex_Rev") - R("R_EX_but_e_") +
 R("R_EX_but_e__Rev")  = 0.

Mbutp: -R("R_BUTt2rpp") + R("R_BUTt2rpp_Rev") + R("R_BUTtex") -
 R("R_BUTtex_Rev")  = 0.

Mbutso3b: R("Mbutso3b_Transp") + R("R_EX_butso3_e_") -
 R("R_EX_butso3_e__Rev")  = 0.

Mbutso3c: R("R_BUTSO3abcpp") - R("R_FDMO4")  = 0.

Mbutso3e: -R("R_BUTSO3tex") + R("R_BUTSO3tex_Rev") -
 R("R_EX_butso3_e_") + R("R_EX_butso3_e__Rev")  = 0.

Mbutso3p: -R("R_BUTSO3abcpp") + R("R_BUTSO3tex") -
 R("R_BUTSO3tex_Rev")  = 0.

Mca2b: R("Mca2b_Transp") + R("R_EX_ca2_e_") - R("R_EX_ca2_e__Rev") 
 = 0.

Mca2c: -R("R_CA2t3pp") - R("R_CAt6pp") + R("R_CAt6pp_Rev") -
 0.004737 R("R_Ec_biomass_iAF1260_core_59p81M")  = 0.

Mca2e: -R("R_CA2tex") + R("R_CA2tex_Rev") - R("R_EX_ca2_e_") +
 R("R_EX_ca2_e__Rev")  = 0.

Mca2p: R("R_CA2t3pp") + R("R_CA2tex") - R("R_CA2tex_Rev") +
 R("R_CAt6pp") - R("R_CAt6pp_Rev")  = 0.

Mcampc: R("R_ADNCYC") - R("R_PDE1")  = 0.

Mcbaspc: R("R_ASPCT") + R("R_DHORTS") - R("R_DHORTS_Rev")  = 0.

Mcbib: R("Mcbib_Transp") + R("R_EX_cbi_e_") - R("R_EX_cbi_e__Rev") 
 = 0.

Mcbic: -R("R_CBIAT") + R("R_CBIAT_Rev") + R("R_CBIuabcpp")  = 0.

Mcbie: -R("R_CBItonex") - R("R_EX_cbi_e_") + R("R_EX_cbi_e__Rev")  =
 0.

Mcbip: R("R_CBItonex") - R("R_CBIuabcpp")  = 0.

Mcbl1b: R("Mcbl1b_Transp") + R("R_EX_cbl1_e_") -
 R("R_EX_cbl1_e__Rev")  = 0.

Mcbl1c: R("R_CBL1abcpp") - R("R_CBLAT") + R("R_CBLAT_Rev")  = 0.

Mcbl1e: -R("R_CBL1tonex") - R("R_EX_cbl1_e_") +
 R("R_EX_cbl1_e__Rev")  = 0.

Mcbl1p: -R("R_CBL1abcpp") + R("R_CBL1tonex")  = 0.

Mcbpc: -R("R_ASPCT") + R("R_CBMKr") - R("R_CBMKr_Rev") + R("R_CBPS")
 - R("R_OCBT") + R("R_OCBT_Rev") + R("R_OXAMTC")  = 0.

Mcd2b: R("Mcd2b_Transp") + R("R_EX_cd2_e_") - R("R_EX_cd2_e__Rev") 
 = 0.

Mcd2c: -R("R_CD2abcpp") - R("R_CD2t3pp") + R("R_CD2tpp")  = 0.

Mcd2e: -R("R_CD2tex") + R("R_CD2tex_Rev") - R("R_EX_cd2_e_") +
 R("R_EX_cd2_e__Rev")  = 0.

Mcd2p: R("R_CD2abcpp") + R("R_CD2t3pp") + R("R_CD2tex") -
 R("R_CD2tex_Rev") - R("R_CD2tpp")  = 0.

Mcddec5eACPc: -R("R_3OAS141") + R("R_EAR121x") + R("R_EAR121y")  =
 0.

Mcdec3eACPc: -R("R_3OAS121") + R("R_T2DECAI") - R("R_T2DECAI_Rev") 
 = 0.

Mcdpc: R("R_CYTK1") - R("R_CYTK1_Rev") - R("R_NDPK3") +
 R("R_NDPK3_Rev") + R("R_NTP5") - R("R_RNDR3") - R("R_RNDR3b")  = 0.

Mcdpdddecgc: -R("R_CDAPPA120") + R("R_DASYN120") - R("R_PGSA120") -
 R("R_PSSA120")  = 0.

Mcdpdhdec9egc: -R("R_CDAPPA161") + R("R_DASYN161") - R("R_PGSA161")
 - R("R_PSSA161")  = 0.

Mcdpdhdecgc: -R("R_CDAPPA160") + R("R_DASYN160") - R("R_PGSA160") -
 R("R_PSSA160")  = 0.

Mcdpdodec11egc: -R("R_CDAPPA181") + R("R_DASYN181") - R("R_PGSA181")
 - R("R_PSSA181")  = 0.

Mcdpdodecgc: -R("R_CDAPPA180") + R("R_DASYN180") - R("R_PGSA180") -
 R("R_PSSA180")  = 0.

Mcdpdtdec7egc: -R("R_CDAPPA141") + R("R_DASYN141") - R("R_PGSA141")
 - R("R_PSSA141")  = 0.

Mcdpdtdecgc: -R("R_CDAPPA140") + R("R_DASYN140") - R("R_PGSA140") -
 R("R_PSSA140")  = 0.

Mcechdddc: -R("R_DHPPD") + R("R_PPPNDO")  = 0.

Mcenchdddc: R("R_CINNDO") - R("R_DHCIND")  = 0.

Mcglyb: R("Mcglyb_Transp") + R("R_EX_cgly_e_") -
 R("R_EX_cgly_e__Rev")  = 0.

Mcglyc: -R("R_AMPTASECG") + R("R_CGLYabcpp")  = 0.

Mcglye: -R("R_CGLYtex") + R("R_CGLYtex_Rev") - R("R_EX_cgly_e_") +
 R("R_EX_cgly_e__Rev")  = 0.

Mcglyp: -R("R_CGLYabcpp") + R("R_CGLYtex") - R("R_CGLYtex_Rev") +
 R("R_GTHRDHpp")  = 0.

Mcholb: R("Mcholb_Transp") + R("R_EX_chol_e_") -
 R("R_EX_chol_e__Rev")  = 0.

Mcholc: R("R_CHLabcpp") + R("R_CHLt2pp") + R("R_GPDDA1")  = 0.

Mchole: -R("R_CHLtex") + R("R_CHLtex_Rev") - R("R_EX_chol_e_") +
 R("R_EX_chol_e__Rev")  = 0.

Mcholp: -R("R_CHLabcpp") - R("R_CHLt2pp") + R("R_CHLtex") -
 R("R_CHLtex_Rev") + R("R_GPDDA1pp")  = 0.

Mchorc: -R("R_ADCS") - R("R_ANS") - R("R_CHORM") + R("R_CHORS") -
 R("R_CHRPL") - R("R_ICHORS") - R("R_ICHORSi") + R("R_ICHORS_Rev")  = 0.

Mcinnmc: -R("R_CINNDO")  = 0.

Mcitb: R("Mcitb_Transp") + R("R_EX_cit_e_") - R("R_EX_cit_e__Rev") 
 = 0.

Mcitc: -R("R_ACONTa") + R("R_ACONTa_Rev") - R("R_CITL") +
 R("R_CITt7pp") + R("R_CS") + 2 R("R_FE3DCITabcpp")  = 0.

Mcite: -R("R_CITtex") + R("R_CITtex_Rev") - R("R_EX_cit_e_") +
 R("R_EX_cit_e__Rev")  = 0.

Mcitp: -R("R_CITt7pp") + R("R_CITtex") - R("R_CITtex_Rev")  = 0.

McitrLc: -R("R_ARGSS") + R("R_OCBT") - R("R_OCBT_Rev")  = 0.

Mckdoc: R("R_KDOCT2") - R("R_MOAT") - R("R_MOAT2") - R("R_MOAT3C") 
 = 0.

Mclb: R("Mclb_Transp") + R("R_EX_cl_e_") - R("R_EX_cl_e__Rev")  = 0.

Mclc: 2 R("R_CLt3_2pp") - 0.004737
 R("R_Ec_biomass_iAF1260_core_59p81M")  = 0.

Mcle: -R("R_CLtex") + R("R_CLtex_Rev") - R("R_EX_cl_e_") +
 R("R_EX_cl_e__Rev")  = 0.

Mclp: -2 R("R_CLt3_2pp") + R("R_CLtex") - R("R_CLtex_Rev")  = 0.

Mclpn120p: -R("R_CLPNH120pp") + R("R_CLPNS120pp") -
 R("R_CLPNS120pp_Rev")  = 0.

Mclpn140p: -R("R_CLPNH140pp") + R("R_CLPNS140pp") -
 R("R_CLPNS140pp_Rev")  = 0.

Mclpn141p: -R("R_CLPNH141pp") + R("R_CLPNS141pp") -
 R("R_CLPNS141pp_Rev")  = 0.

Mclpn160p: -R("R_CLPNH160pp") + R("R_CLPNS160pp") -
 R("R_CLPNS160pp_Rev")  = 0.

Mclpn161p: -R("R_CLPNH161pp") + R("R_CLPNS161pp") -
 R("R_CLPNS161pp_Rev")  = 0.

Mclpn180p: -R("R_CLPNH180pp") + R("R_CLPNS180pp") -
 R("R_CLPNS180pp_Rev")  = 0.

Mclpn181p: -R("R_CLPNH181pp") + R("R_CLPNS181pp") -
 R("R_CLPNS181pp_Rev")  = 0.

Mcmpb: R("Mcmpb_Transp") + R("R_EX_cmp_e_") - R("R_EX_cmp_e__Rev") 
 = 0.

Mcmpc: R("R_CDAPPA120") + R("R_CDAPPA140") + R("R_CDAPPA141") +
 R("R_CDAPPA160") + R("R_CDAPPA161") + R("R_CDAPPA180") + R("R_CDAPPA181")
 - R("R_CMPN") + R("R_CYTDK2") - R("R_CYTK1") + R("R_CYTK1_Rev") +
 R("R_MECDPS") + R("R_MOAT") + R("R_MOAT2") + R("R_MOAT3C") - R("R_NTD4") +
 R("R_NTPP4") + R("R_PGSA120") + R("R_PGSA140") + R("R_PGSA141") +
 R("R_PGSA160") + R("R_PGSA161") + R("R_PGSA180") + R("R_PGSA181") +
 R("R_PPNCL2") + R("R_PSSA120") + R("R_PSSA140") + R("R_PSSA141") +
 R("R_PSSA160") + R("R_PSSA161") + R("R_PSSA180") + R("R_PSSA181")  = 0.

Mcmpe: -R("R_CMPtex") + R("R_CMPtex_Rev") - R("R_EX_cmp_e_") +
 R("R_EX_cmp_e__Rev")  = 0.

Mcmpp: R("R_CMPtex") - R("R_CMPtex_Rev") - R("R_NTD4pp")  = 0.

Mco2b: R("Mco2b_Transp") + R("R_EX_co2_e_") - R("R_EX_co2_e__Rev") 
 = 0.

Mco2c: R("R_3OAS100") + R("R_3OAS120") + R("R_3OAS121") +
 R("R_3OAS140") + R("R_3OAS141") + R("R_3OAS160") + R("R_3OAS161") +
 R("R_3OAS180") + R("R_3OAS181") + R("R_3OAS60") + R("R_3OAS80") + R("R_ACHBS") +
 R("R_ACLS") + R("R_ADMDC") + R("R_AKGDH") + R("R_ALLTAMH") +
 R("R_AOBUTDs") + R("R_AOXSr") - R("R_AOXSr_Rev") + R("R_ARGDC") + R("R_ASP1DC")
 - R("R_CBMKr") + R("R_CBMKr_Rev") + R("R_CO2tpp") - R("R_CO2tpp_Rev")
 + 2 R("R_CPPPGO") + 2 R("R_CPPPGO2") + 2 R("R_CYNTAH") + R("R_DAPDC")
 - R("R_DBTS") + R("R_DHNAOT4") + R("R_DXPS") + R("R_FDH4pp") +
 R("R_FDH5pp") + R("R_FHL") + R("R_GLUDC") + R("R_GLXCL") + R("R_GLYCL") +
 R("R_GND") - R("R_HCO3E") + R("R_HCO3E_Rev") + R("R_ICDHyr") -
 R("R_ICDHyr_Rev") + R("R_IGPS") + R("R_KAS14") + R("R_KAS15") + R("R_KG6PDC") +
 R("R_LYSDC") + R("R_MACPD") + R("R_MALDDH") + R("R_ME1") + R("R_ME2")
 + R("R_MMCD") + R("R_NNDPR") + R("R_OAADC") + R("R_OMCDC") +
 R("R_OMPDC") + R("R_OPHBDC") + R("R_ORNDC") + R("R_OXGDC2") + R("R_PDH") +
 R("R_PDX5PS") + R("R_POX") - R("R_PPC") + R("R_PPCDC") + R("R_PPCK") +
 R("R_PPND") + R("R_PPNDH") + R("R_PSD120") + R("R_PSD140") +
 R("R_PSD141") + R("R_PSD160") + R("R_PSD161") + R("R_PSD180") + R("R_PSD181") +
 R("R_SADH") + R("R_TAUDO") + R("R_THZPSN") + R("R_UDPGDC") +
 R("R_UGLYCH") + 4 R("R_UPPDC1") + R("R_URIC")  = 0.

Mco2e: -R("R_CO2tex") + R("R_CO2tex_Rev") - R("R_EX_co2_e_") +
 R("R_EX_co2_e__Rev")  = 0.

Mco2p: R("R_ARGDCpp") + R("R_CO2tex") - R("R_CO2tex_Rev") -
 R("R_CO2tpp") + R("R_CO2tpp_Rev")  = 0.

Mcoac: R("R_ACACT1r") - R("R_ACACT1r_Rev") - R("R_ACALD") +
 R("R_ACALD_Rev") + R("R_ACANTHAT") - R("R_ACCOAL") + R("R_ACGS") +
 R("R_ACOATA") - R("R_ACOATA_Rev") - R("R_ACPS1") - R("R_ACS") - R("R_AKGDH") +
 R("R_AOXSr") - R("R_AOXSr_Rev") + R("R_AST") - R("R_CRNCAL2") -
 R("R_CRNDCAL2") + R("R_CS") - R("R_CTBTCAL2") + R("R_DPCOAK") - 0.000576
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_FACOAE100") + R("R_FACOAE120")
 + R("R_FACOAE140") + R("R_FACOAE141") + R("R_FACOAE160") +
 R("R_FACOAE161") + R("R_FACOAE180") + R("R_FACOAE181") + R("R_FACOAE60") +
 R("R_FACOAE80") - R("R_FACOAL100t2pp") - R("R_FACOAL120t2pp") -
 R("R_FACOAL140t2pp") - R("R_FACOAL141t2pp") - R("R_FACOAL160t2pp") -
 R("R_FACOAL161t2pp") - R("R_FACOAL180t2pp") - R("R_FACOAL181t2pp") -
 R("R_FACOAL60t2pp") - R("R_FACOAL80t2pp") + R("R_G1PACT") + R("R_GLCATr") -
 R("R_GLCATr_Rev") + R("R_GLYAT") - R("R_GLYAT_Rev") + R("R_HSST") +
 R("R_IPPS") + R("R_KAS15") - R("R_KAT1") - R("R_KAT2") - R("R_KAT3") -
 R("R_KAT4") - R("R_KAT5") - R("R_KAT6") - R("R_KAT7") - R("R_KAT8") +
 R("R_MALS") + R("R_MALTATr") - R("R_MALTATr_Rev") + R("R_MCITS") +
 R("R_MCOATA") - R("R_MCOATA_Rev") + R("R_NPHS") + R("R_O16AT") - R("R_OBTFL") -
 R("R_PACCOAL") - R("R_PDH") - R("R_PFL") + R("R_PTA2") + R("R_PTAr") -
 R("R_PTAr_Rev") + R("R_SERAT") - R("R_SERAT_Rev") + R("R_SPMDAT1") +
 R("R_SPMDAT2") - R("R_SUCBZL") - R("R_SUCOAS") + R("R_SUCOAS_Rev") +
 R("R_TDPADGAT") + R("R_THDPS")  = 0.

Mcobalt2b: R("Mcobalt2b_Transp") + R("R_EX_cobalt2_e_") -
 R("R_EX_cobalt2_e__Rev")  = 0.

Mcobalt2c: -R("R_COBALT2abcpp") - R("R_COBALT2t3pp") +
 R("R_COBALT2tpp") - 0.003158 R("R_Ec_biomass_iAF1260_core_59p81M")  = 0.

Mcobalt2e: -R("R_COBALT2tex") + R("R_COBALT2tex_Rev") -
 R("R_EX_cobalt2_e_") + R("R_EX_cobalt2_e__Rev")  = 0.

Mcobalt2p: R("R_COBALT2abcpp") + R("R_COBALT2t3pp") +
 R("R_COBALT2tex") - R("R_COBALT2tex_Rev") - R("R_COBALT2tpp")  = 0.

Mcolipab: R("Mcolipab_Transp") + R("R_EX_colipa_e_") -
 R("R_EX_colipa_e__Rev")  = 0.

Mcolipac: -R("R_COLIPAabcpp") + R("R_HEPT4")  = 0.

Mcolipae: R("R_COLIPAtex") - R("R_EX_colipa_e_") +
 R("R_EX_colipa_e__Rev") - R("R_LIPAHT2ex")  = 0.

Mcolipap: R("R_COLIPAabcpp") - R("R_COLIPAtex") - R("R_ECA4OALpp") -
 R("R_LA4NTpp") - R("R_O16A4Lpp")  = 0.

Mcpe160c: R("R_CFAS160E")  = 0.

Mcpe180c: R("R_CFAS180E")  = 0.

Mcpg160c: R("R_CFAS160G")  = 0.

Mcpg180c: R("R_CFAS180G")  = 0.

Mcpgnb: R("Mcpgnb_Transp") + R("R_EX_cpgn_e_") -
 R("R_EX_cpgn_e__Rev")  = 0.

Mcpgnc: R("R_CPGNabcpp") - 2 R("R_CPGNR1") - 2 R("R_CPGNR2") - 2
 R("R_CPGNR3")  = 0.

Mcpgne: R("R_CPGNexs") - R("R_CPGNtonex") - R("R_EX_cpgn_e_") +
 R("R_EX_cpgn_e__Rev")  = 0.

Mcpgnp: -R("R_CPGNabcpp") + R("R_CPGNtonex")  = 0.

Mcpgnunb: R("Mcpgnunb_Transp") + R("R_EX_cpgn_un_e_") -
 R("R_EX_cpgn_un_e__Rev")  = 0.

Mcpgnunc: 2 R("R_CPGNR1") + 2 R("R_CPGNR2") + 2 R("R_CPGNR3") -
 R("R_CPGNUtpp")  = 0.

Mcpgnune: -R("R_CPGNexs") + R("R_CPGNUtex") - R("R_EX_cpgn_un_e_") +
 R("R_EX_cpgn_un_e__Rev")  = 0.

Mcpgnunp: -R("R_CPGNUtex") + R("R_CPGNUtpp")  = 0.

Mcpppg3c: -R("R_CPPPGO") - R("R_CPPPGO2") + R("R_UPPDC1")  = 0.

Mcrnb: R("Mcrnb_Transp") + R("R_EX_crn_e_") - R("R_EX_crn_e__Rev") 
 = 0.

Mcrnc: R("R_CRNabcpp") - R("R_CRNBTCT") + R("R_CRNBTCT_Rev") -
 R("R_CRNCAL2") - R("R_CRNCBCT") + R("R_CRNCBCT_Rev") + R("R_CRNt2rpp") -
 R("R_CRNt2rpp_Rev") + R("R_CRNt7pp") + R("R_CRNt8pp")  = 0.

Mcrncoac: R("R_CRNBTCT") - R("R_CRNBTCT_Rev") + R("R_CRNCAL2") -
 R("R_CRNCAR") + R("R_CRNCAR_Rev") + R("R_CRNCBCT") - R("R_CRNCBCT_Rev") -
 R("R_CRNCDH") + R("R_CRNCDH_Rev")  = 0.

McrnDc: R("R_CRNDabcpp") - R("R_CRNDCAL2") + R("R_CRNDt2rpp") -
 R("R_CRNDt2rpp_Rev") - R("R_CRNt8pp")  = 0.

McrnDcoac: R("R_CRNCAR") - R("R_CRNCAR_Rev") + R("R_CRNDCAL2")  = 0.

McrnDp: -R("R_CRNDabcpp") - R("R_CRNDt2rpp") + R("R_CRNDt2rpp_Rev")
 + R("R_CRNt8pp")  = 0.

Mcrne: -R("R_CRNtex") + R("R_CRNtex_Rev") - R("R_EX_crn_e_") +
 R("R_EX_crn_e__Rev")  = 0.

Mcrnp: -R("R_CRNabcpp") - R("R_CRNt2rpp") + R("R_CRNt2rpp_Rev") -
 R("R_CRNt7pp") - R("R_CRNt8pp") + R("R_CRNtex") - R("R_CRNtex_Rev")  =
 0.

Mcsnb: R("Mcsnb_Transp") + R("R_EX_csn_e_") - R("R_EX_csn_e__Rev") 
 = 0.

Mcsnc: R("R_CMPN") - R("R_CSND") + R("R_CSNt2pp") + R("R_CYTDH")  =
 0.

Mcsne: -R("R_CSNtex") + R("R_CSNtex_Rev") - R("R_EX_csn_e_") +
 R("R_EX_csn_e__Rev")  = 0.

Mcsnp: -R("R_CSNt2pp") + R("R_CSNtex") - R("R_CSNtex_Rev")  = 0.

Mctbtc: R("R_CRNCBCT") - R("R_CRNCBCT_Rev") + R("R_CTBTabcpp") -
 R("R_CTBTCAL2") + R("R_CTBTt2rpp") - R("R_CTBTt2rpp_Rev")  = 0.

Mctbtcoac: -R("R_CRNCBCT") + R("R_CRNCBCT_Rev") + R("R_CRNCDH") -
 R("R_CRNCDH_Rev") + R("R_CTBTCAL2")  = 0.

Mctbtp: -R("R_CTBTabcpp") - R("R_CTBTt2rpp") + R("R_CTBTt2rpp_Rev") 
 = 0.

Mctpc: R("R_CTPS2") - R("R_DASYN120") - R("R_DASYN140") -
 R("R_DASYN141") - R("R_DASYN160") - R("R_DASYN161") - R("R_DASYN180") -
 R("R_DASYN181") - 0.1335 R("R_Ec_biomass_iAF1260_core_59p81M") - R("R_KDOCT2")
 - R("R_MEPCT") + R("R_NDPK3") - R("R_NDPK3_Rev") - R("R_NTP5") -
 R("R_NTPP4") - R("R_PPNCL2") - R("R_RNTR3c")  = 0.

Mcu2b: R("Mcu2b_Transp") + R("R_EX_cu2_e_") - R("R_EX_cu2_e__Rev") 
 = 0.

Mcu2c: -R("R_CU2abcpp") + R("R_CU2tpp") - 0.003158
 R("R_Ec_biomass_iAF1260_core_59p81M")  = 0.

Mcu2e: -R("R_CU2tex") + R("R_CU2tex_Rev") - R("R_EX_cu2_e_") +
 R("R_EX_cu2_e__Rev")  = 0.

Mcu2p: 4 R("R_CU1Opp") + R("R_CU2abcpp") + R("R_CU2tex") -
 R("R_CU2tex_Rev") - R("R_CU2tpp")  = 0.

Mcub: R("Mcub_Transp") + R("R_EX_cu_e_") - R("R_EX_cu_e__Rev")  = 0.

Mcuc: -R("R_CU1abcpp") - R("R_CUt3")  = 0.

Mcue: R("R_CUt3") - R("R_CUtex") + R("R_CUtex_Rev") -
 R("R_EX_cu_e_") + R("R_EX_cu_e__Rev")  = 0.

Mcup: R("R_CU1abcpp") - 4 R("R_CU1Opp") + R("R_CUtex") -
 R("R_CUtex_Rev")  = 0.

Mcyanb: R("Mcyanb_Transp") + R("R_EX_cyan_e_") -
 R("R_EX_cyan_e__Rev")  = 0.

Mcyanc: -R("R_CYANST") - R("R_MCPST")  = 0.

Mcyane: -R("R_CYANtex") + R("R_CYANtex_Rev") - R("R_EX_cyan_e_") +
 R("R_EX_cyan_e__Rev")  = 0.

Mcyanp: -R("R_CYANSTpp") + R("R_CYANtex") - R("R_CYANtex_Rev")  = 0.

Mcyntb: R("Mcyntb_Transp") + R("R_EX_cynt_e_") -
 R("R_EX_cynt_e__Rev")  = 0.

Mcyntc: -R("R_CYNTAH") + R("R_CYNTt2pp")  = 0.

Mcynte: -R("R_CYNTtex") + R("R_CYNTtex_Rev") - R("R_EX_cynt_e_") +
 R("R_EX_cynt_e__Rev")  = 0.

Mcyntp: -R("R_CYNTt2pp") + R("R_CYNTtex") - R("R_CYNTtex_Rev")  = 0.

McysDb: R("McysDb_Transp") + R("R_EX_cys_D_e_") -
 R("R_EX_cys_D_e__Rev")  = 0.

McysDc: R("R_CYSDabcpp") - R("R_CYSDDS")  = 0.

McysDe: -R("R_CYSDtex") + R("R_CYSDtex_Rev") - R("R_EX_cys_D_e_") +
 R("R_EX_cys_D_e__Rev")  = 0.

McysDp: -R("R_CYSDabcpp") + R("R_CYSDtex") - R("R_CYSDtex_Rev")  =
 0.

McysLb: R("McysLb_Transp") + R("R_EX_cys_L_e_") -
 R("R_EX_cys_L_e__Rev")  = 0.

McysLc: R("R_AMPTASECG") - R("R_CYSabc2pp") + R("R_CYSabcpp") -
 R("R_CYSDS") + R("R_CYSS") - R("R_CYStpp") - R("R_CYSTRS") - 0.09158
 R("R_Ec_biomass_iAF1260_core_59p81M") - R("R_GLUCYS") - R("R_PPNCL2") -
 R("R_SHSL1") - R("R_THZPSN")  = 0.

McysLe: -R("R_CYStex") + R("R_CYStex_Rev") - R("R_EX_cys_L_e_") +
 R("R_EX_cys_L_e__Rev")  = 0.

McysLp: R("R_CYSabc2pp") - R("R_CYSabcpp") + R("R_CYStex") -
 R("R_CYStex_Rev") + R("R_CYStpp")  = 0.

McystLc: -R("R_CYSTL") + R("R_SHSL1")  = 0.

Mcystrnac: R("R_CYSTRS")  = 0.

Mcytdb: R("Mcytdb_Transp") + R("R_EX_cytd_e_") -
 R("R_EX_cytd_e__Rev")  = 0.

Mcytdc: -R("R_CYTD") - R("R_CYTDH") - R("R_CYTDK2") +
 R("R_CYTDt2pp") + R("R_CYTDt2rpp") - R("R_CYTDt2rpp_Rev") + R("R_NTD4")  = 0.

Mcytde: -R("R_CYTDtex") + R("R_CYTDtex_Rev") - R("R_EX_cytd_e_") +
 R("R_EX_cytd_e__Rev")  = 0.

Mcytdp: R("R_3NTD4pp") - R("R_CYTDt2pp") - R("R_CYTDt2rpp") +
 R("R_CYTDt2rpp_Rev") + R("R_CYTDtex") - R("R_CYTDtex_Rev") + R("R_NTD4pp") 
 = 0.

Mdad2b: R("Mdad2b_Transp") + R("R_EX_dad_2_e_") -
 R("R_EX_dad_2_e__Rev")  = 0.

Mdad2c: -R("R_DADA") + R("R_DADNt2pp") + R("R_NTD6") - R("R_PUNP2")
 + R("R_PUNP2_Rev")  = 0.

Mdad2e: -R("R_DADNtex") + R("R_DADNtex_Rev") - R("R_EX_dad_2_e_") +
 R("R_EX_dad_2_e__Rev")  = 0.

Mdad2p: -R("R_DADNt2pp") + R("R_DADNtex") - R("R_DADNtex_Rev") +
 R("R_NTD6pp")  = 0.

Mdad5c: -R("R_5DOAN") + R("R_BTS4") + 2 R("R_CPPPGO2")  = 0.

Mdadpc: R("R_DADK") - R("R_DADK_Rev") - R("R_NDPK8") +
 R("R_NDPK8_Rev") + R("R_RNDR1") + R("R_RNDR1b")  = 0.

Mdampb: R("Mdampb_Transp") + R("R_EX_damp_e_") -
 R("R_EX_damp_e__Rev")  = 0.

Mdampc: -R("R_DADK") + R("R_DADK_Rev") - R("R_NTD6") + R("R_NTPP5") 
 = 0.

Mdampe: -R("R_DAMPtex") + R("R_DAMPtex_Rev") - R("R_EX_damp_e_") +
 R("R_EX_damp_e__Rev")  = 0.

Mdampp: R("R_DAMPtex") - R("R_DAMPtex_Rev") - R("R_NTD6pp")  = 0.

Mdannc: R("R_AMAOTr") - R("R_AMAOTr_Rev") - R("R_DBTS")  = 0.

Mdatpc: -R("R_DATPHs") - 0.02617
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_NDPK8") - R("R_NDPK8_Rev") - R("R_NTPP5") + R("R_RNTR1c")  =
 0.

Mdb4pc: R("R_DB4PS") - R("R_RBFSa")  = 0.

Mdc2coac: R("R_ACOAD4f") - R("R_ACOAD4f_Rev") + R("R_ECOAH4") -
 R("R_ECOAH4_Rev")  = 0.

MdcaACPc: -R("R_3OAS120") + R("R_AACPS8") + R("R_EAR100x") +
 R("R_EAR100y") - R("R_FA100ACPHi")  = 0.

Mdcab: R("Mdcab_Transp") + R("R_EX_dca_e_") - R("R_EX_dca_e__Rev") 
 = 0.

Mdcac: -R("R_AACPS8") + R("R_FA100ACPHi") + R("R_FACOAE100")  = 0.

Mdcacoac: -R("R_ACOAD4f") + R("R_ACOAD4f_Rev") - R("R_FACOAE100") +
 R("R_FACOAL100t2pp") + R("R_KAT5")  = 0.

Mdcae: -R("R_DCAtex") + R("R_DCAtex_Rev") - R("R_EX_dca_e_") +
 R("R_EX_dca_e__Rev")  = 0.

Mdcampc: -R("R_ADSL1r") + R("R_ADSL1r_Rev") + R("R_ADSS")  = 0.

Mdcap: R("R_DCAtex") - R("R_DCAtex_Rev") - R("R_FACOAL100t2pp")  =
 0.

Mdcdpc: R("R_CYTK2") - R("R_CYTK2_Rev") - R("R_NDPK7") +
 R("R_NDPK7_Rev") + R("R_RNDR3") + R("R_RNDR3b")  = 0.

Mdcmpb: R("Mdcmpb_Transp") + R("R_EX_dcmp_e_") -
 R("R_EX_dcmp_e__Rev")  = 0.

Mdcmpc: -R("R_CYTK2") + R("R_CYTK2_Rev") - R("R_NTD3") +
 R("R_NTPP3")  = 0.

Mdcmpe: -R("R_DCMPtex") + R("R_DCMPtex_Rev") - R("R_EX_dcmp_e_") +
 R("R_EX_dcmp_e__Rev")  = 0.

Mdcmpp: R("R_DCMPtex") - R("R_DCMPtex_Rev") - R("R_NTD3pp")  = 0.

Mdctpc: -R("R_DCTPD") - 0.02702
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_NDPK7") - R("R_NDPK7_Rev") - R("R_NTPP3") + R("R_RNTR3c")  =
 0.

Mdcytb: R("Mdcytb_Transp") + R("R_EX_dcyt_e_") -
 R("R_EX_dcyt_e__Rev")  = 0.

Mdcytc: -R("R_DCYTD") + R("R_DCYTt2pp") + R("R_NTD3")  = 0.

Mdcyte: -R("R_DCYTtex") + R("R_DCYTtex_Rev") - R("R_EX_dcyt_e_") +
 R("R_EX_dcyt_e__Rev")  = 0.

Mdcytp: -R("R_DCYTt2pp") + R("R_DCYTtex") - R("R_DCYTtex_Rev") +
 R("R_NTD3pp")  = 0.

Mdd2coac: R("R_ACOAD5f") - R("R_ACOAD5f_Rev") + R("R_ECOAH5") -
 R("R_ECOAH5_Rev")  = 0.

MddcaACPc: -R("R_3OAS140") + R("R_AACPS7") - R("R_AGPAT120") +
 R("R_EAR120x") + R("R_EAR120y") - R("R_EDTXS1") - R("R_FA120ACPHi") -
 R("R_G3PAT120")  = 0.

Mddcab: R("Mddcab_Transp") + R("R_EX_ddca_e_") -
 R("R_EX_ddca_e__Rev")  = 0.

Mddcac: -R("R_2AGPEAT120") - R("R_2AGPGAT120") - R("R_AACPS7") +
 R("R_FA120ACPHi") + R("R_FACOAE120") + R("R_LPLIPAL2A120") +
 R("R_LPLIPAL2E120") + R("R_LPLIPAL2G120")  = 0.

Mddcacoac: -R("R_ACOAD5f") + R("R_ACOAD5f_Rev") - R("R_FACOAE120") +
 R("R_FACOAL120t2pp") + R("R_KAT6")  = 0.

Mddcae: -R("R_DDCAtexi") - R("R_EX_ddca_e_") + R("R_EX_ddca_e__Rev")
  = 0.

Mddcap: R("R_DDCAtexi") - R("R_FACOAL120t2pp") +
 R("R_LPLIPAL1A120pp") + R("R_LPLIPAL1E120pp") + R("R_LPLIPAL1G120pp") +
 R("R_PLIPA1A120pp") + R("R_PLIPA1E120pp") + R("R_PLIPA1G120pp") + R("R_PLIPA2A120pp") +
 R("R_PLIPA2E120pp") + R("R_PLIPA2G120pp")  = 0.

Mdgdpc: R("R_DGK1") - R("R_DGK1_Rev") - R("R_NDPK5") +
 R("R_NDPK5_Rev") + R("R_RNDR2") + R("R_RNDR2b")  = 0.

Mdgmpb: R("Mdgmpb_Transp") + R("R_EX_dgmp_e_") -
 R("R_EX_dgmp_e__Rev")  = 0.

Mdgmpc: -R("R_DGK1") + R("R_DGK1_Rev") - R("R_NTD8") + R("R_NTPP1") 
 = 0.

Mdgmpe: -R("R_DGMPtex") + R("R_DGMPtex_Rev") - R("R_EX_dgmp_e_") +
 R("R_EX_dgmp_e__Rev")  = 0.

Mdgmpp: R("R_DGMPtex") - R("R_DGMPtex_Rev") - R("R_NTD8pp")  = 0.

Mdgsnb: R("Mdgsnb_Transp") + R("R_EX_dgsn_e_") -
 R("R_EX_dgsn_e__Rev")  = 0.

Mdgsnc: R("R_DGSNt2pp") + R("R_NTD8") + R("R_NTPTP1") - R("R_PUNP4")
 + R("R_PUNP4_Rev")  = 0.

Mdgsne: -R("R_DGSNtex") + R("R_DGSNtex_Rev") - R("R_EX_dgsn_e_") +
 R("R_EX_dgsn_e__Rev")  = 0.

Mdgsnp: -R("R_DGSNt2pp") + R("R_DGSNtex") - R("R_DGSNtex_Rev") +
 R("R_NTD8pp")  = 0.

Mdgtpc: -0.02702 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_NDPK5") - R("R_NDPK5_Rev") - R("R_NTPP1") - R("R_NTPTP1") + R("R_RNTR2c")  =
 0.

Mdhab: R("Mdhab_Transp") + R("R_EX_dha_e_") - R("R_EX_dha_e__Rev") 
 = 0.

Mdhac: -R("R_DHAPT") + R("R_DHAtpp") - R("R_DHAtpp_Rev") +
 R("R_F6PA") - R("R_F6PA_Rev") + R("R_GLYCDx")  = 0.

Mdhae: -R("R_DHAtex") + R("R_DHAtex_Rev") - R("R_EX_dha_e_") +
 R("R_EX_dha_e__Rev")  = 0.

Mdhap: R("R_DHAtex") - R("R_DHAtex_Rev") - R("R_DHAtpp") +
 R("R_DHAtpp_Rev")  = 0.

Mdhapc: R("R_DHAPT") + R("R_FBA") - R("R_FBA_Rev") + R("R_FCLPA") -
 R("R_FCLPA_Rev") + R("R_G3PD2") - R("R_G3PD2_Rev") + R("R_G3PD5") +
 R("R_G3PD6") + R("R_G3PD7") - R("R_MGSA") - R("R_QULNS") + R("R_RMPA") -
 R("R_RMPA_Rev") + R("R_TGBPA") - R("R_TGBPA_Rev") - R("R_TPI") +
 R("R_TPI_Rev")  = 0.

Mdhcinnmc: R("R_3HCINNMH") + R("R_DHCIND") - R("R_DHCINDO")  = 0.

Mdhfc: -R("R_DHFR") + R("R_DHFR_Rev") + R("R_DHFS") + R("R_TMDS")  =
 0.

Mdhmptpc: R("R_DHPTPE") - R("R_DHPTPE_Rev")  = 0.

Mdhnac: -R("R_DHNAOT4") + R("R_NPHS")  = 0.

Mdhnptc: -R("R_DHNPA2") + R("R_DNMPPA")  = 0.

MdhorSc: -R("R_DHORD2") - R("R_DHORD5") - R("R_DHORTS") +
 R("R_DHORTS_Rev")  = 0.

Mdhpmpc: -R("R_DNMPPA") + R("R_DNTPPA")  = 0.

Mdhpppnc: R("R_3HPPPNH") + R("R_DHPPD") - R("R_HPPPNDO")  = 0.

Mdhptc: -R("R_DHFS") + R("R_DHPS2")  = 0.

Mdhptdc: -R("R_DHPTDCs") + R("R_RHCCE")  = 0.

Mdidpc: R("R_NTP11")  = 0.

Mdimpb: R("Mdimpb_Transp") + R("R_EX_dimp_e_") -
 R("R_EX_dimp_e__Rev")  = 0.

Mdimpc: -R("R_NTD12") + R("R_NTPP10")  = 0.

Mdimpe: -R("R_DIMPtex") + R("R_DIMPtex_Rev") - R("R_EX_dimp_e_") +
 R("R_EX_dimp_e__Rev")  = 0.

Mdimpp: R("R_DIMPtex") - R("R_DIMPtex_Rev") - R("R_NTD12pp")  = 0.

Mdinb: R("Mdinb_Transp") + R("R_EX_din_e_") - R("R_EX_din_e__Rev") 
 = 0.

Mdinc: R("R_DADA") + R("R_DINSt2pp") + R("R_NTD12") - R("R_PUNP6") +
 R("R_PUNP6_Rev")  = 0.

Mdine: -R("R_DINStex") + R("R_DINStex_Rev") - R("R_EX_din_e_") +
 R("R_EX_din_e__Rev")  = 0.

Mdinp: -R("R_DINSt2pp") + R("R_DINStex") - R("R_DINStex_Rev") +
 R("R_NTD12pp")  = 0.

Mditpc: R("R_DATPHs") - R("R_NTP11") - R("R_NTPP10")  = 0.

Mdkmppc: -R("R_DKMPPD") - R("R_DKMPPD3") + R("R_MDRPD")  = 0.

Mdmbzidc: -R("R_NNDMBRT")  = 0.

Mdmlzc: R("R_RBFSa") - 2 R("R_RBFSb")  = 0.

Mdmppc: -R("R_DMATT") + R("R_DMPPS") + R("R_IPDDI") -
 R("R_IPDDI_Rev")  = 0.

Mdmsb: R("Mdmsb_Transp") + R("R_EX_dms_e_") - R("R_EX_dms_e__Rev") 
 = 0.

Mdmsc: R("R_DMSOR1") + R("R_DMSOR2")  = 0.

Mdmse: -R("R_DMStex") + R("R_DMStex_Rev") - R("R_EX_dms_e_") +
 R("R_EX_dms_e__Rev")  = 0.

Mdmsob: R("Mdmsob_Transp") + R("R_EX_dmso_e_") -
 R("R_EX_dmso_e__Rev")  = 0.

Mdmsoc: -R("R_DMSOR1") - R("R_DMSOR2")  = 0.

Mdmsoe: -R("R_DMSOtex") + R("R_DMSOtex_Rev") - R("R_EX_dmso_e_") +
 R("R_EX_dmso_e__Rev")  = 0.

Mdmsop: -R("R_DMSOR1pp") - R("R_DMSOR2pp") + R("R_DMSOtex") -
 R("R_DMSOtex_Rev")  = 0.

Mdmsp: R("R_DMSOR1pp") + R("R_DMSOR2pp") + R("R_DMStex") -
 R("R_DMStex_Rev")  = 0.

Mdnadc: -R("R_NADS1") + R("R_NNATr") - R("R_NNATr_Rev")  = 0.

Mdopab: R("Mdopab_Transp") + R("R_EX_dopa_e_") -
 R("R_EX_dopa_e__Rev")  = 0.

Mdopae: -R("R_DOPAtex") + R("R_DOPAtex_Rev") - R("R_EX_dopa_e_") +
 R("R_EX_dopa_e__Rev")  = 0.

Mdopap: -R("R_42A12BOOXpp") + R("R_DOPAtex") - R("R_DOPAtex_Rev")  =
 0.

Mdpcoac: -R("R_DPCOAK") + R("R_PTPATi") - R("R_TPRDCOAS")  = 0.

Mdsbaoxp: R("R_DSBAO1") + R("R_DSBAO2")  = 0.

Mdsbardp: -R("R_DSBAO1") - R("R_DSBAO2")  = 0.

Mdsbcoxp: -R("R_DSBCGT") - R("R_TDSR1")  = 0.

Mdsbcrdp: R("R_DSBCGT") + R("R_TDSR1")  = 0.

Mdsbdoxc: -R("R_DSBDR") + R("R_TDSR1") + R("R_TDSR2")  = 0.

Mdsbdrdc: R("R_DSBDR") - R("R_TDSR1") - R("R_TDSR2")  = 0.

Mdsbgoxp: -R("R_DSBGGT") - R("R_TDSR2")  = 0.

Mdsbgrdp: R("R_DSBGGT") + R("R_TDSR2")  = 0.

Mdsclc: -R("R_SHCHD2") + R("R_UPP3MT")  = 0.

Mdtbtc: -R("R_BTS4") + R("R_DBTS")  = 0.

Mdtdp4aaddgc: -R("R_AADDGT") + R("R_TDPADGAT")  = 0.

Mdtdp4addgc: -R("R_TDPADGAT") + R("R_TDPAGTA")  = 0.

Mdtdp4d6dgc: -R("R_TDPAGTA") - R("R_TDPDRE") + R("R_TDPGDH")  = 0.

Mdtdp4d6dmc: R("R_TDPDRE") - R("R_TDPDRR")  = 0.

Mdtdpc: R("R_AADDGT") + R("R_DTMPK") - R("R_DTMPK_Rev") -
 R("R_NDPK4") + R("R_NDPK4_Rev") + R("R_RHAT1")  = 0.

Mdtdpgluc: R("R_G1PTT") - R("R_TDPGDH")  = 0.

Mdtdprmnc: -R("R_RHAT1") + R("R_TDPDRR")  = 0.

Mdtmpb: R("Mdtmpb_Transp") + R("R_EX_dtmp_e_") -
 R("R_EX_dtmp_e__Rev")  = 0.

Mdtmpc: -R("R_DTMPK") + R("R_DTMPK_Rev") - R("R_NTD5") +
 R("R_NTPP7") + R("R_TMDK1") + R("R_TMDS")  = 0.

Mdtmpe: -R("R_DTMPtex") + R("R_DTMPtex_Rev") - R("R_EX_dtmp_e_") +
 R("R_EX_dtmp_e__Rev")  = 0.

Mdtmpp: R("R_DTMPtex") - R("R_DTMPtex_Rev") - R("R_NTD5pp")  = 0.

Mdttpc: -0.02617 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_G1PTT") + R("R_NDPK4") - R("R_NDPK4_Rev") - R("R_NTPP7")  = 0.

Mdudpc: -R("R_NDPK6") + R("R_NDPK6_Rev") + R("R_RNDR4") +
 R("R_RNDR4b") + R("R_URIDK2r") - R("R_URIDK2r_Rev")  = 0.

Mdumpb: R("Mdumpb_Transp") + R("R_EX_dump_e_") -
 R("R_EX_dump_e__Rev")  = 0.

Mdumpc: R("R_DURIK1") + R("R_DUTPDP") - R("R_NTD1") - R("R_TMDS") -
 R("R_URIDK2r") + R("R_URIDK2r_Rev")  = 0.

Mdumpe: -R("R_DUMPtex") + R("R_DUMPtex_Rev") - R("R_EX_dump_e_") +
 R("R_EX_dump_e__Rev")  = 0.

Mdumpp: R("R_DUMPtex") - R("R_DUMPtex_Rev") - R("R_NTD1pp")  = 0.

Mdurib: R("Mdurib_Transp") + R("R_EX_duri_e_") -
 R("R_EX_duri_e__Rev")  = 0.

Mduric: R("R_DCYTD") - R("R_DURIK1") - R("R_DURIPP") +
 R("R_DURIPP_Rev") + R("R_DURIt2pp") + R("R_NTD1")  = 0.

Mdurie: -R("R_DURItex") + R("R_DURItex_Rev") - R("R_EX_duri_e_") +
 R("R_EX_duri_e__Rev")  = 0.

Mdurip: -R("R_DURIt2pp") + R("R_DURItex") - R("R_DURItex_Rev") +
 R("R_NTD1pp")  = 0.

Mdutpc: R("R_DCTPD") - R("R_DUTPDP") + R("R_NDPK6") -
 R("R_NDPK6_Rev") + R("R_RNTR4c")  = 0.

Mdxyl5pc: -R("R_DXPRIi") + R("R_DXPS") + R("R_DXYLK") -
 R("R_PDX5PS") - R("R_THZPSN")  = 0.

Mdxylc: -R("R_DXYLK")  = 0.

Me4pc: -R("R_DDPA") - R("R_E4PD") + R("R_E4PD_Rev") + R("R_TALA") -
 R("R_TALA_Rev") - R("R_TKT2") + R("R_TKT2_Rev")  = 0.

Meca2undp: R("R_ECAP1pp") - R("R_ECAP2pp")  = 0.

Meca3undp: R("R_ECAP2pp") - R("R_ECAP3pp")  = 0.

Meca4colipab: R("Meca4colipab_Transp") + R("R_EX_eca4colipa_e_") -
 R("R_EX_eca4colipa_e__Rev")  = 0.

Meca4colipae: R("R_ECA4COLIPAtex") - R("R_EX_eca4colipa_e_") +
 R("R_EX_eca4colipa_e__Rev")  = 0.

Meca4colipap: -R("R_ECA4COLIPAtex") + R("R_ECA4OALpp")  = 0.

Meca4undp: -R("R_ECA4OALpp") + R("R_ECAP3pp")  = 0.

Meig3pc: R("R_IG3PS") - R("R_IGPDH")  = 0.

Menlipab: R("Menlipab_Transp") + R("R_EX_enlipa_e_") -
 R("R_EX_enlipa_e__Rev")  = 0.

Menlipae: R("R_ENLIPAtex") - R("R_EX_enlipa_e_") +
 R("R_EX_enlipa_e__Rev")  = 0.

Menlipap: -R("R_ENLIPAtex") + R("R_PETNT161pp") + R("R_PETNT181pp") 
 = 0.

Menterb: R("Menterb_Transp") + R("R_EX_enter_e_") -
 R("R_EX_enter_e__Rev")  = 0.

Menterc: R("R_ENTCS") - R("R_ENTERES") + 2 R("R_FEENTERR1") + 2
 R("R_FEENTERR2") + 2 R("R_FEENTERR3") - R("R_FEENTERtpp")  = 0.

Mentere: -R("R_EX_enter_e_") + R("R_EX_enter_e__Rev") -
 R("R_FEENTERexs") + R("R_FEENTERtex")  = 0.

Menterp: -R("R_FEENTERtex") + R("R_FEENTERtpp")  = 0.

Methab: R("Methab_Transp") + R("R_EX_etha_e_") -
 R("R_EX_etha_e__Rev")  = 0.

Methac: -R("R_ETHAAL") + R("R_ETHAt2pp") + R("R_GPDDA2")  = 0.

Methae: -R("R_ETHAtex") + R("R_ETHAtex_Rev") - R("R_EX_etha_e_") +
 R("R_EX_etha_e__Rev")  = 0.

Methap: -R("R_ETHAt2pp") + R("R_ETHAtex") - R("R_ETHAtex_Rev") +
 R("R_GPDDA2pp")  = 0.

Methso3b: R("Methso3b_Transp") + R("R_EX_ethso3_e_") -
 R("R_EX_ethso3_e__Rev")  = 0.

Methso3c: R("R_ETHSO3abcpp") - R("R_FDMO3")  = 0.

Methso3e: -R("R_ETHSO3tex") + R("R_ETHSO3tex_Rev") -
 R("R_EX_ethso3_e_") + R("R_EX_ethso3_e__Rev")  = 0.

Methso3p: -R("R_ETHSO3abcpp") + R("R_ETHSO3tex") -
 R("R_ETHSO3tex_Rev")  = 0.

Metohb: R("Metohb_Transp") + R("R_EX_etoh_e_") -
 R("R_EX_etoh_e__Rev")  = 0.

Metohc: -R("R_ALCD2x") + R("R_ALCD2x_Rev") + R("R_ETOHt2rpp") -
 R("R_ETOHt2rpp_Rev")  = 0.

Metohe: -R("R_ETOHtex") + R("R_ETOHtex_Rev") - R("R_EX_etoh_e_") +
 R("R_EX_etoh_e__Rev")  = 0.

Metohp: -R("R_ETOHt2rpp") + R("R_ETOHt2rpp_Rev") + R("R_ETOHtex") -
 R("R_ETOHtex_Rev")  = 0.

Mf1pc: -R("R_FRUK") + R("R_FRUptspp")  = 0.

Mf6pb: R("Mf6pb_Transp") + R("R_EX_f6p_e_") - R("R_EX_f6p_e__Rev") 
 = 0.

Mf6pc: R("R_ALLULPE") - R("R_ALLULPE_Rev") - R("R_F6PA") +
 R("R_F6PA_Rev") - R("R_F6PP") + R("R_F6Pt6_2pp") + R("R_FBP") +
 R("R_FRUpts2pp") + R("R_G6PDA") - R("R_GF6PTA") + R("R_HEX7") + R("R_M1PD") -
 R("R_M1PD_Rev") + R("R_MAN6PI") - R("R_MAN6PI_Rev") - R("R_PFK") + R("R_PGI")
 - R("R_PGI_Rev") + R("R_SBTPD") - R("R_SBTPD_Rev") + R("R_TALA") -
 R("R_TALA_Rev") + R("R_TKT2") - R("R_TKT2_Rev")  = 0.

Mf6pe: -R("R_EX_f6p_e_") + R("R_EX_f6p_e__Rev") - R("R_F6Ptex") +
 R("R_F6Ptex_Rev")  = 0.

Mf6pp: -R("R_F6Pt6_2pp") + R("R_F6Ptex") - R("R_F6Ptex_Rev")  = 0.

Mfadc: -R("R_ACOAD1f") + R("R_ACOAD1f_Rev") - R("R_ACOAD2f") +
 R("R_ACOAD2f_Rev") - R("R_ACOAD3f") + R("R_ACOAD3f_Rev") - R("R_ACOAD4f") +
 R("R_ACOAD4f_Rev") - R("R_ACOAD5f") + R("R_ACOAD5f_Rev") -
 R("R_ACOAD6f") + R("R_ACOAD6f_Rev") - R("R_ACOAD7f") + R("R_ACOAD7f_Rev") -
 R("R_ACOAD8f") + R("R_ACOAD8f_Rev") + R("R_ARBTNR1") + R("R_CPGNR1") -
 R("R_DAAD") - 0.000223 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_FADRx") - R("R_FADRx2") + R("R_FE3HOXR1") + R("R_FE3Ri") + R("R_FECRMR1") +
 R("R_FEENTERR1") + R("R_FEOXAMR1") + R("R_FMNAT") - R("R_PROD2")  =
 0.

Mfadh2c: R("R_ACOAD1f") - R("R_ACOAD1f_Rev") + R("R_ACOAD2f") -
 R("R_ACOAD2f_Rev") + R("R_ACOAD3f") - R("R_ACOAD3f_Rev") + R("R_ACOAD4f")
 - R("R_ACOAD4f_Rev") + R("R_ACOAD5f") - R("R_ACOAD5f_Rev") +
 R("R_ACOAD6f") - R("R_ACOAD6f_Rev") + R("R_ACOAD7f") - R("R_ACOAD7f_Rev") +
 R("R_ACOAD8f") - R("R_ACOAD8f_Rev") - R("R_ARBTNR1") - R("R_CPGNR1") +
 R("R_DAAD") + R("R_FADRx") + R("R_FADRx2") - R("R_FE3HOXR1") -
 R("R_FE3Ri") - R("R_FECRMR1") - R("R_FEENTERR1") - R("R_FEOXAMR1") +
 R("R_PROD2")  = 0.

Mfaldb: R("Mfaldb_Transp") + R("R_EX_fald_e_") -
 R("R_EX_fald_e__Rev")  = 0.

Mfaldc: R("R_FALDtpp") - R("R_FALDtpp_Rev") - R("R_FALGTHLs") +
 R("R_FALGTHLs_Rev") + R("R_FDMO2") + R("R_MTRPOX") + R("R_SARCOX")  = 0.

Mfalde: -R("R_EX_fald_e_") + R("R_EX_fald_e__Rev") - R("R_FALDtex")
 + R("R_FALDtex_Rev")  = 0.

Mfaldp: R("R_FALDtex") - R("R_FALDtex_Rev") - R("R_FALDtpp") +
 R("R_FALDtpp_Rev")  = 0.

Mfc1pc: R("R_FCLK") - R("R_FCLPA") + R("R_FCLPA_Rev")  = 0.

MfclLc: R("R_FCI") - R("R_FCI_Rev") - R("R_FCLK")  = 0.

Mfdpc: -R("R_FBA") + R("R_FBA_Rev") - R("R_FBP") + R("R_FRUK") +
 R("R_PFK")  = 0.

Mfe2b: R("Mfe2b_Transp") + R("R_EX_fe2_e_") - R("R_EX_fe2_e__Rev") 
 = 0.

Mfe2c: 2 R("R_ARBTNR1") + 2 R("R_ARBTNR2") + 2 R("R_ARBTNR3") + 2
 R("R_CPGNR1") + 2 R("R_CPGNR2") + 2 R("R_CPGNR3") - 0.007106
 R("R_Ec_biomass_iAF1260_core_59p81M") - R("R_FCLT") + R("R_FE2abcpp") +
 R("R_FE2t2pp") - R("R_FE2t3pp") + R("R_FE2tpp") + 2 R("R_FE3HOXR1") + 2
 R("R_FE3HOXR2") + 2 R("R_FE3HOXR3") + 2 R("R_FE3Ri") + 2 R("R_FECRMR1") + 2
 R("R_FECRMR2") + 2 R("R_FECRMR3") + 2 R("R_FEENTERR1") + 2
 R("R_FEENTERR2") + 2 R("R_FEENTERR3") + 2 R("R_FEOXAMR1") + 2 R("R_FEOXAMR2") + 2
 R("R_FEOXAMR3") - R("R_SHCHF")  = 0.

Mfe2e: -R("R_EX_fe2_e_") + R("R_EX_fe2_e__Rev") - R("R_FE2tex") +
 R("R_FE2tex_Rev")  = 0.

Mfe2p: -R("R_FE2abcpp") - R("R_FE2t2pp") + R("R_FE2t3pp") +
 R("R_FE2tex") - R("R_FE2tex_Rev") - R("R_FE2tpp") - 4 R("R_FEROpp")  = 0.

Mfe3b: R("Mfe3b_Transp") + R("R_EX_fe3_e_") - R("R_EX_fe3_e__Rev") 
 = 0.

Mfe3c: -0.007106 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_ENTERES2") + R("R_FE3abcpp") + R("R_FE3DCITabcpp") + R("R_FE3DHBZR") - 2
 R("R_FE3Ri")  = 0.

Mfe3dcitb: R("Mfe3dcitb_Transp") + R("R_EX_fe3dcit_e_") -
 R("R_EX_fe3dcit_e__Rev")  = 0.

Mfe3dcite: -R("R_EX_fe3dcit_e_") + R("R_EX_fe3dcit_e__Rev") -
 R("R_FE3DCITtonex")  = 0.

Mfe3dcitp: -R("R_FE3DCITabcpp") + R("R_FE3DCITtonex")  = 0.

Mfe3dhbzsb: R("Mfe3dhbzsb_Transp") + R("R_EX_fe3dhbzs_e_") -
 R("R_EX_fe3dhbzs_e__Rev")  = 0.

Mfe3dhbzsc: -R("R_FE3DHBZR") + R("R_FE3DHBZSabcpp")  = 0.

Mfe3dhbzse: -R("R_EX_fe3dhbzs_e_") + R("R_EX_fe3dhbzs_e__Rev") -
 R("R_FE3DHBZStonex")  = 0.

Mfe3dhbzsp: -R("R_FE3DHBZSabcpp") + R("R_FE3DHBZStonex")  = 0.

Mfe3e: -R("R_ARBTNexs") - R("R_CPGNexs") - R("R_EX_fe3_e_") +
 R("R_EX_fe3_e__Rev") - R("R_FE3HOXexs") - R("R_FE3tex") + R("R_FE3tex_Rev")
 - R("R_FECRMexs") - R("R_FEENTERexs") - R("R_FEOXAMexs")  = 0.

Mfe3hoxb: R("Mfe3hoxb_Transp") + R("R_EX_fe3hox_e_") -
 R("R_EX_fe3hox_e__Rev")  = 0.

Mfe3hoxc: R("R_FE3HOXabcpp") - 2 R("R_FE3HOXR1") - 2 R("R_FE3HOXR2")
 - 2 R("R_FE3HOXR3")  = 0.

Mfe3hoxe: -R("R_EX_fe3hox_e_") + R("R_EX_fe3hox_e__Rev") +
 R("R_FE3HOXexs") - R("R_FE3HOXtonex")  = 0.

Mfe3hoxp: -R("R_FE3HOXabcpp") + R("R_FE3HOXtonex")  = 0.

Mfe3hoxunb: R("Mfe3hoxunb_Transp") + R("R_EX_fe3hox_un_e_") -
 R("R_EX_fe3hox_un_e__Rev")  = 0.

Mfe3hoxunc: 2 R("R_FE3HOXR1") + 2 R("R_FE3HOXR2") + 2
 R("R_FE3HOXR3") - R("R_FE3HOXUtex")  = 0.

Mfe3hoxune: -R("R_EX_fe3hox_un_e_") + R("R_EX_fe3hox_un_e__Rev") -
 R("R_FE3HOXexs") + R("R_FE3HOXUtpp")  = 0.

Mfe3hoxunp: R("R_FE3HOXUtex") - R("R_FE3HOXUtpp")  = 0.

Mfe3p: -R("R_FE3abcpp") + R("R_FE3tex") - R("R_FE3tex_Rev") + 4
 R("R_FEROpp")  = 0.

Mfecrmb: R("Mfecrmb_Transp") + R("R_EX_fecrm_e_") -
 R("R_EX_fecrm_e__Rev")  = 0.

Mfecrmc: R("R_FECRMabcpp") - 2 R("R_FECRMR1") - 2 R("R_FECRMR2") - 2
 R("R_FECRMR3")  = 0.

Mfecrme: -R("R_EX_fecrm_e_") + R("R_EX_fecrm_e__Rev") +
 R("R_FECRMexs") - R("R_FECRMtonex")  = 0.

Mfecrmp: -R("R_FECRMabcpp") + R("R_FECRMtonex")  = 0.

Mfecrmunb: R("Mfecrmunb_Transp") + R("R_EX_fecrm_un_e_") -
 R("R_EX_fecrm_un_e__Rev")  = 0.

Mfecrmunc: 2 R("R_FECRMR1") + 2 R("R_FECRMR2") + 2 R("R_FECRMR3") -
 R("R_FECRMUtpp")  = 0.

Mfecrmune: -R("R_EX_fecrm_un_e_") + R("R_EX_fecrm_un_e__Rev") -
 R("R_FECRMexs") + R("R_FECRMUtex")  = 0.

Mfecrmunp: -R("R_FECRMUtex") + R("R_FECRMUtpp")  = 0.

Mfeenterb: R("Mfeenterb_Transp") + R("R_EX_feenter_e_") -
 R("R_EX_feenter_e__Rev")  = 0.

Mfeenterc: -R("R_ENTERES2") + R("R_FEENTERabcpp") - 2
 R("R_FEENTERR1") - 2 R("R_FEENTERR2") - 2 R("R_FEENTERR3")  = 0.

Mfeentere: -R("R_EX_feenter_e_") + R("R_EX_feenter_e__Rev") +
 R("R_FEENTERexs") - R("R_FEENTERtonex")  = 0.

Mfeenterp: -R("R_FEENTERabcpp") + R("R_FEENTERtonex")  = 0.

Mfeoxamb: R("Mfeoxamb_Transp") + R("R_EX_feoxam_e_") -
 R("R_EX_feoxam_e__Rev")  = 0.

Mfeoxamc: R("R_FEOXAMabcpp") - 2 R("R_FEOXAMR1") - 2 R("R_FEOXAMR2")
 - 2 R("R_FEOXAMR3")  = 0.

Mfeoxame: -R("R_EX_feoxam_e_") + R("R_EX_feoxam_e__Rev") +
 R("R_FEOXAMexs") - R("R_FEOXAMtonex")  = 0.

Mfeoxamp: -R("R_FEOXAMabcpp") + R("R_FEOXAMtonex")  = 0.

Mfeoxamunb: R("Mfeoxamunb_Transp") + R("R_EX_feoxam_un_e_") -
 R("R_EX_feoxam_un_e__Rev")  = 0.

Mfeoxamunc: 2 R("R_FEOXAMR1") + 2 R("R_FEOXAMR2") + 2
 R("R_FEOXAMR3") - R("R_FEOXAMUtpp")  = 0.

Mfeoxamune: -R("R_EX_feoxam_un_e_") + R("R_EX_feoxam_un_e__Rev") -
 R("R_FEOXAMexs") + R("R_FEOXAMUtex")  = 0.

Mfeoxamunp: -R("R_FEOXAMUtex") + R("R_FEOXAMUtpp")  = 0.

Mfgamc: R("R_GARFT") - R("R_GARFT_Rev") + R("R_GART") - R("R_PRFGS")
  = 0.

Mfldoxc: -R("R_FLDR") + R("R_RNTR1c") + R("R_RNTR2c") +
 R("R_RNTR3c") + R("R_RNTR4c")  = 0.

Mfldrdc: R("R_FLDR") - R("R_RNTR1c") - R("R_RNTR2c") - R("R_RNTR3c")
 - R("R_RNTR4c")  = 0.

Mfmettrnac: R("R_FMETTRS")  = 0.

Mfmnc: R("R_ARBTNR2") + R("R_CPGNR2") + R("R_FDMO") + R("R_FDMO2") +
 R("R_FDMO3") + R("R_FDMO4") + R("R_FDMO6") + R("R_FE3HOXR2") +
 R("R_FECRMR2") + R("R_FEENTERR2") + R("R_FEOXAMR2") - R("R_FMNAT") -
 R("R_FMNRx") - R("R_FMNRx2") + R("R_RBFK")  = 0.

Mfmnh2c: -R("R_ARBTNR2") - R("R_CPGNR2") - R("R_FDMO") -
 R("R_FDMO2") - R("R_FDMO3") - R("R_FDMO4") - R("R_FDMO6") - R("R_FE3HOXR2") -
 R("R_FECRMR2") - R("R_FEENTERR2") - R("R_FEOXAMR2") + R("R_FMNRx") +
 R("R_FMNRx2")  = 0.

Mforb: R("Mforb_Transp") + R("R_EX_for_e_") - R("R_EX_for_e__Rev") 
 = 0.

Mforc: 2 R("R_AMPMS2") + R("R_DB4PS") + R("R_DKMPPD") +
 R("R_DKMPPD3") - R("R_FHL") + R("R_FORt2pp") - R("R_FORtppi") + R("R_FTHFD") -
 R("R_GART") + R("R_GTPCI") + R("R_GTPCII2") + R("R_OBTFL") + R("R_PFL") +
 R("R_SFGTHi") + R("R_UPLA4FNF")  = 0.

Mfore: -R("R_EX_for_e_") + R("R_EX_for_e__Rev") - R("R_FORtex") +
 R("R_FORtex_Rev")  = 0.

Mforp: -R("R_FDH4pp") - R("R_FDH5pp") - R("R_FORt2pp") +
 R("R_FORtex") - R("R_FORtex_Rev") + R("R_FORtppi")  = 0.

Mfpramc: -R("R_PRAIS") + R("R_PRFGS")  = 0.

Mfpricac: R("R_AICART") - R("R_AICART_Rev") + R("R_IMPC") -
 R("R_IMPC_Rev")  = 0.

Mfrdpc: R("R_GRTT") - R("R_HEMEOS") - R("R_OCTDPS") - R("R_UDCPDPS")
  = 0.

Mfrub: R("Mfrub_Transp") + R("R_EX_fru_e_") - R("R_EX_fru_e__Rev") 
 = 0.

Mfruc: R("R_F6PP") + R("R_FFSD") - R("R_HEX7") + R("R_XYLI2") -
 R("R_XYLI2_Rev")  = 0.

Mfrue: -R("R_EX_fru_e_") + R("R_EX_fru_e__Rev") - R("R_FRUtex") +
 R("R_FRUtex_Rev")  = 0.

Mfrulysb: R("Mfrulysb_Transp") + R("R_EX_frulys_e_") -
 R("R_EX_frulys_e__Rev")  = 0.

Mfrulysc: R("R_FRULYSE") - R("R_FRULYSE_Rev") - R("R_FRULYSK") +
 R("R_FRULYSt2pp")  = 0.

Mfrulyse: -R("R_EX_frulys_e_") + R("R_EX_frulys_e__Rev") -
 R("R_FRULYStex") + R("R_FRULYStex_Rev")  = 0.

Mfrulysp: -R("R_FRULYSt2pp") + R("R_FRULYStex") -
 R("R_FRULYStex_Rev")  = 0.

Mfrulyspc: -R("R_FRULYSDG") + R("R_FRULYSDG_Rev") + R("R_FRULYSK") 
 = 0.

Mfrup: -R("R_FRUpts2pp") - R("R_FRUptspp") + R("R_FRUtex") -
 R("R_FRUtex_Rev")  = 0.

Mfruurb: R("Mfruurb_Transp") + R("R_EX_fruur_e_") -
 R("R_EX_fruur_e__Rev")  = 0.

Mfruurc: R("R_FRUURt2rpp") - R("R_FRUURt2rpp_Rev") + R("R_GUI1") -
 R("R_GUI1_Rev") + R("R_MANAO") - R("R_MANAO_Rev")  = 0.

Mfruure: -R("R_EX_fruur_e_") + R("R_EX_fruur_e__Rev") -
 R("R_FRUURtex") + R("R_FRUURtex_Rev")  = 0.

Mfruurp: -R("R_FRUURt2rpp") + R("R_FRUURt2rpp_Rev") +
 R("R_FRUURtex") - R("R_FRUURtex_Rev")  = 0.

MfucLb: R("MfucLb_Transp") + R("R_EX_fuc_L_e_") -
 R("R_EX_fuc_L_e__Rev")  = 0.

MfucLc: -R("R_FCI") + R("R_FCI_Rev") + R("R_FUCtpp") -
 R("R_FUCtpp_Rev")  = 0.

MfucLe: -R("R_EX_fuc_L_e_") + R("R_EX_fuc_L_e__Rev") - R("R_FUCtex")
 + R("R_FUCtex_Rev")  = 0.

MfucLp: R("R_FUCtex") - R("R_FUCtex_Rev") - R("R_FUCtpp") +
 R("R_FUCtpp_Rev")  = 0.

Mfumb: R("Mfumb_Transp") + R("R_EX_fum_e_") - R("R_EX_fum_e__Rev") 
 = 0.

Mfumc: R("R_ADSL1r") - R("R_ADSL1r_Rev") + R("R_ADSL2r") -
 R("R_ADSL2r_Rev") + R("R_ARGSL") - R("R_ARGSL_Rev") - R("R_ASPO5") +
 R("R_ASPT") - R("R_FRD2") - R("R_FRD3") - R("R_FUM") + R("R_FUM_Rev") +
 R("R_FUMt2_2pp") + R("R_FUMt2_3pp") + R("R_HKNTDH") - 3 R("R_PPPGO3") +
 R("R_SUCDi") + R("R_SUCFUMtpp") - R("R_SUCFUMtpp_Rev")  = 0.

Mfume: -R("R_EX_fum_e_") + R("R_EX_fum_e__Rev") - R("R_FUMtex") +
 R("R_FUMtex_Rev")  = 0.

Mfump: -R("R_FUMt2_2pp") - R("R_FUMt2_3pp") + R("R_FUMtex") -
 R("R_FUMtex_Rev") - R("R_SUCFUMtpp") + R("R_SUCFUMtpp_Rev")  = 0.

Mg1pb: R("Mg1pb_Transp") + R("R_EX_g1p_e_") - R("R_EX_g1p_e__Rev") 
 = 0.

Mg1pc: -R("R_G1PTT") - R("R_GALUi") + R("R_GLCP") + R("R_GLCP2") -
 R("R_GLGC") + R("R_MLTP1") - R("R_MLTP1_Rev") + R("R_MLTP2") -
 R("R_MLTP2_Rev") + R("R_MLTP3") - R("R_MLTP3_Rev") - R("R_PGMT") +
 R("R_PGMT_Rev") + R("R_UGLT") - R("R_UGLT_Rev")  = 0.

Mg1pe: -R("R_EX_g1p_e_") + R("R_EX_g1p_e__Rev") - R("R_G1Ptex") +
 R("R_G1Ptex_Rev")  = 0.

Mg1pp: -R("R_G1PPpp") + R("R_G1Ptex") - R("R_G1Ptex_Rev") +
 R("R_UDPGPpp")  = 0.

Mg3pc: R("R_DDPGALA") - R("R_DDPGALA_Rev") + R("R_DRPA") -
 R("R_DXPS") + R("R_EDA") + R("R_F6PA") - R("R_F6PA_Rev") + R("R_FBA") -
 R("R_FBA_Rev") - R("R_GAPD") + R("R_GAPD_Rev") - R("R_TALA") +
 R("R_TALA_Rev") + R("R_TGBPA") - R("R_TGBPA_Rev") + R("R_TKT1") - R("R_TKT1_Rev") +
 R("R_TKT2") - R("R_TKT2_Rev") + R("R_TPI") - R("R_TPI_Rev") +
 R("R_TRPS1") + R("R_TRPS3")  = 0.

Mg3pcb: R("Mg3pcb_Transp") + R("R_EX_g3pc_e_") -
 R("R_EX_g3pc_e__Rev")  = 0.

Mg3pcc: R("R_G3PCabcpp") - R("R_GPDDA1")  = 0.

Mg3pce: -R("R_EX_g3pc_e_") + R("R_EX_g3pc_e__Rev") - R("R_G3PCtex")
 + R("R_G3PCtex_Rev")  = 0.

Mg3pcp: -R("R_G3PCabcpp") + R("R_G3PCtex") - R("R_G3PCtex_Rev") -
 R("R_GPDDA1pp")  = 0.

Mg3peb: R("Mg3peb_Transp") + R("R_EX_g3pe_e_") -
 R("R_EX_g3pe_e__Rev")  = 0.

Mg3pec: R("R_G3PEabcpp") - R("R_GPDDA2") + R("R_LPLIPAL2ATE120") +
 R("R_LPLIPAL2ATE140") + R("R_LPLIPAL2ATE141") + R("R_LPLIPAL2ATE160") +
 R("R_LPLIPAL2ATE161") + R("R_LPLIPAL2ATE180") + R("R_LPLIPAL2ATE181")
 + R("R_LPLIPAL2E120") + R("R_LPLIPAL2E140") + R("R_LPLIPAL2E141") +
 R("R_LPLIPAL2E160") + R("R_LPLIPAL2E161") + R("R_LPLIPAL2E180") +
 R("R_LPLIPAL2E181")  = 0.

Mg3pee: -R("R_EX_g3pe_e_") + R("R_EX_g3pe_e__Rev") - R("R_G3PEtex")
 + R("R_G3PEtex_Rev")  = 0.

Mg3pep: -R("R_G3PEabcpp") + R("R_G3PEtex") - R("R_G3PEtex_Rev") -
 R("R_GPDDA2pp") + R("R_LPLIPAL1E120pp") + R("R_LPLIPAL1E140pp") +
 R("R_LPLIPAL1E141pp") + R("R_LPLIPAL1E160pp") + R("R_LPLIPAL1E161pp") +
 R("R_LPLIPAL1E180pp") + R("R_LPLIPAL1E181pp")  = 0.

Mg3pgb: R("Mg3pgb_Transp") + R("R_EX_g3pg_e_") -
 R("R_EX_g3pg_e__Rev")  = 0.

Mg3pgc: R("R_G3PGabcpp") - R("R_GPDDA4") + R("R_LPLIPAL2ATG120") +
 R("R_LPLIPAL2ATG140") + R("R_LPLIPAL2ATG141") + R("R_LPLIPAL2ATG160") +
 R("R_LPLIPAL2ATG161") + R("R_LPLIPAL2ATG180") + R("R_LPLIPAL2ATG181")
 + R("R_LPLIPAL2G120") + R("R_LPLIPAL2G140") + R("R_LPLIPAL2G141") +
 R("R_LPLIPAL2G160") + R("R_LPLIPAL2G161") + R("R_LPLIPAL2G180") +
 R("R_LPLIPAL2G181")  = 0.

Mg3pge: -R("R_EX_g3pg_e_") + R("R_EX_g3pg_e__Rev") - R("R_G3PGtex")
 + R("R_G3PGtex_Rev")  = 0.

Mg3pgp: -R("R_G3PGabcpp") + R("R_G3PGtex") - R("R_G3PGtex_Rev") -
 R("R_GPDDA4pp") + R("R_LPLIPAL1G120pp") + R("R_LPLIPAL1G140pp") +
 R("R_LPLIPAL1G141pp") + R("R_LPLIPAL1G160pp") + R("R_LPLIPAL1G161pp") +
 R("R_LPLIPAL1G180pp") + R("R_LPLIPAL1G181pp")  = 0.

Mg3pib: R("Mg3pib_Transp") + R("R_EX_g3pi_e_") -
 R("R_EX_g3pi_e__Rev")  = 0.

Mg3pic: R("R_G3PIabcpp") - R("R_GPDDA5")  = 0.

Mg3pie: -R("R_EX_g3pi_e_") + R("R_EX_g3pi_e__Rev") - R("R_G3PItex")
 + R("R_G3PItex_Rev")  = 0.

Mg3pip: -R("R_G3PIabcpp") + R("R_G3PItex") - R("R_G3PItex_Rev") -
 R("R_GPDDA5pp")  = 0.

Mg3psb: R("Mg3psb_Transp") + R("R_EX_g3ps_e_") -
 R("R_EX_g3ps_e__Rev")  = 0.

Mg3psc: R("R_G3PSabcpp") - R("R_GPDDA3")  = 0.

Mg3pse: -R("R_EX_g3ps_e_") + R("R_EX_g3ps_e__Rev") - R("R_G3PStex")
 + R("R_G3PStex_Rev")  = 0.

Mg3psp: -R("R_G3PSabcpp") + R("R_G3PStex") - R("R_G3PStex_Rev") -
 R("R_GPDDA3pp")  = 0.

Mg6pb: R("Mg6pb_Transp") + R("R_EX_g6p_e_") - R("R_EX_g6p_e__Rev") 
 = 0.

Mg6pc: R("R_AB6PGH") + R("R_FFSD") + R("R_FRULYSDG") -
 R("R_FRULYSDG_Rev") - R("R_G6PDH2r") + R("R_G6PDH2r_Rev") - R("R_G6PP") +
 R("R_G6Pt6_2pp") + R("R_GLCptspp") + R("R_HEX1") - R("R_PGI") + R("R_PGI_Rev")
 + R("R_PGMT") - R("R_PGMT_Rev") + R("R_TRE6PH") - R("R_TRE6PS")  = 0.

Mg6pe: -R("R_EX_g6p_e_") + R("R_EX_g6p_e__Rev") - R("R_G6Ptex") +
 R("R_G6Ptex_Rev")  = 0.

Mg6pp: -R("R_G6Pt6_2pp") + R("R_G6Ptex") - R("R_G6Ptex_Rev")  = 0.

Mgagicolipac: R("R_GALT1") - R("R_GLCTR2")  = 0.

Mgal1pb: R("Mgal1pb_Transp") + R("R_EX_gal1p_e_") -
 R("R_EX_gal1p_e__Rev")  = 0.

Mgal1pc: R("R_GALKr") - R("R_GALKr_Rev") - R("R_UGLT") +
 R("R_UGLT_Rev")  = 0.

Mgal1pe: -R("R_EX_gal1p_e_") + R("R_EX_gal1p_e__Rev") -
 R("R_GAL1Ptex") + R("R_GAL1Ptex_Rev")  = 0.

Mgal1pp: -R("R_GAL1PPpp") + R("R_GAL1Ptex") - R("R_GAL1Ptex_Rev") +
 R("R_UDPGALPpp")  = 0.

Mgalb: R("Mgalb_Transp") + R("R_EX_gal_e_") - R("R_EX_gal_e__Rev") 
 = 0.

MgalbDb: R("MgalbDb_Transp") + R("R_EX_gal_bD_e_") -
 R("R_EX_gal_bD_e__Rev")  = 0.

MgalbDe: -R("R_EX_gal_bD_e_") + R("R_EX_gal_bD_e__Rev") -
 R("R_GALBDtex") + R("R_GALBDtex_Rev")  = 0.

MgalbDp: R("R_GALBDtex") - R("R_GALBDtex_Rev") - R("R_GALM2pp")  =
 0.

Mgalc: R("R_GALabcpp") - R("R_GALKr") + R("R_GALKr_Rev") +
 R("R_GALS3") + R("R_GALt2pp") + R("R_LACZ")  = 0.

MgalctDb: R("MgalctDb_Transp") + R("R_EX_galct_D_e_") -
 R("R_EX_galct_D_e__Rev")  = 0.

MgalctDc: -R("R_GALCTD") + R("R_GALCTt2rpp") - R("R_GALCTt2rpp_Rev")
  = 0.

MgalctDe: -R("R_EX_galct_D_e_") + R("R_EX_galct_D_e__Rev") -
 R("R_GALCTtex") + R("R_GALCTtex_Rev")  = 0.

MgalctDp: -R("R_GALCTt2rpp") + R("R_GALCTt2rpp_Rev") +
 R("R_GALCTtex") - R("R_GALCTtex_Rev")  = 0.

MgalctnDb: R("MgalctnDb_Transp") + R("R_EX_galctn_D_e_") -
 R("R_EX_galctn_D_e__Rev")  = 0.

MgalctnDc: -R("R_GALCTND") + R("R_GALCTNt2pp")  = 0.

MgalctnDe: -R("R_EX_galctn_D_e_") + R("R_EX_galctn_D_e__Rev") -
 R("R_GALCTNtex") + R("R_GALCTNtex_Rev")  = 0.

MgalctnDp: -R("R_GALCTNt2pp") + R("R_GALCTNtex") -
 R("R_GALCTNtex_Rev")  = 0.

MgalctnLb: R("MgalctnLb_Transp") + R("R_EX_galctn_L_e_") -
 R("R_EX_galctn_L_e__Rev")  = 0.

MgalctnLc: -R("R_GALCTLO") + R("R_GALCTNLt2pp")  = 0.

MgalctnLe: -R("R_EX_galctn_L_e_") + R("R_EX_galctn_L_e__Rev") -
 R("R_GALCTNLtex") + R("R_GALCTNLtex_Rev")  = 0.

MgalctnLp: -R("R_GALCTNLt2pp") + R("R_GALCTNLtex") -
 R("R_GALCTNLtex_Rev")  = 0.

Mgale: -R("R_EX_gal_e_") + R("R_EX_gal_e__Rev") - R("R_GALtex") +
 R("R_GALtex_Rev")  = 0.

Mgalp: R("R_GAL1PPpp") - R("R_GALabcpp") + R("R_GALM2pp") -
 R("R_GALt2pp") + R("R_GALtex") - R("R_GALtex_Rev") + R("R_LACZpp")  = 0.

Mgalt1pc: R("R_GALTptspp") - R("R_GLTPD") + R("R_GLTPD_Rev")  = 0.

Mgaltb: R("Mgaltb_Transp") + R("R_EX_galt_e_") -
 R("R_EX_galt_e__Rev")  = 0.

Mgalte: -R("R_EX_galt_e_") + R("R_EX_galt_e__Rev") - R("R_GALTtex")
 + R("R_GALTtex_Rev")  = 0.

Mgaltp: -R("R_GALTptspp") + R("R_GALTtex") - R("R_GALTtex_Rev")  =
 0.

Mgalurb: R("Mgalurb_Transp") + R("R_EX_galur_e_") -
 R("R_EX_galur_e__Rev")  = 0.

Mgalurc: R("R_GALURt2rpp") - R("R_GALURt2rpp_Rev") - R("R_GUI2") +
 R("R_GUI2_Rev")  = 0.

Mgalure: -R("R_EX_galur_e_") + R("R_EX_galur_e__Rev") -
 R("R_GALURtex") + R("R_GALURtex_Rev")  = 0.

Mgalurp: -R("R_GALURt2rpp") + R("R_GALURt2rpp_Rev") +
 R("R_GALURtex") - R("R_GALURtex_Rev")  = 0.

Mgam1pc: -R("R_G1PACT") - R("R_PGAMT") + R("R_PGAMT_Rev")  = 0.

Mgam6pb: R("Mgam6pb_Transp") + R("R_EX_gam6p_e_") -
 R("R_EX_gam6p_e__Rev")  = 0.

Mgam6pc: R("R_AGDC") - R("R_G6PDA") + R("R_GAM6Pt6_2pp") +
 R("R_GAMptspp") + R("R_GF6PTA") + R("R_PGAMT") - R("R_PGAMT_Rev")  = 0.

Mgam6pe: -R("R_EX_gam6p_e_") + R("R_EX_gam6p_e__Rev") -
 R("R_GAMAN6Ptex") + R("R_GAMAN6Ptex_Rev")  = 0.

Mgam6pp: -R("R_GAM6Pt6_2pp") + R("R_GAMAN6Ptex") -
 R("R_GAMAN6Ptex_Rev")  = 0.

Mgamb: R("Mgamb_Transp") + R("R_EX_gam_e_") - R("R_EX_gam_e__Rev") 
 = 0.

Mgame: -R("R_EX_gam_e_") + R("R_EX_gam_e__Rev") - R("R_GAMtex") +
 R("R_GAMtex_Rev")  = 0.

Mgamp: -R("R_GAMptspp") + R("R_GAMtex") - R("R_GAMtex_Rev")  = 0.

Mgaragundc: -R("R_O16GALFT") + R("R_O16GLCT1")  = 0.

Mgarc: -R("R_GARFT") + R("R_GARFT_Rev") - R("R_GART") +
 R("R_PRAGSr") - R("R_PRAGSr_Rev")  = 0.

Mgbbtnb: R("Mgbbtnb_Transp") + R("R_EX_gbbtn_e_") -
 R("R_EX_gbbtn_e__Rev")  = 0.

Mgbbtnc: R("R_CRNBTCT") - R("R_CRNBTCT_Rev") - R("R_CRNt7pp")  = 0.

Mgbbtne: -R("R_EX_gbbtn_e_") + R("R_EX_gbbtn_e__Rev") -
 R("R_GBBTNtex") + R("R_GBBTNtex_Rev")  = 0.

Mgbbtnp: R("R_CRNt7pp") + R("R_GBBTNtex") - R("R_GBBTNtex_Rev")  =
 0.

Mgcaldc: R("R_DHNPA2") + R("R_FDMO") - R("R_GCALDD")  = 0.

Mgdpb: R("Mgdpb_Transp") + R("R_EX_gdp_e_") - R("R_EX_gdp_e__Rev") 
 = 0.

Mgdpc: R("R_ADK3") - R("R_ADK3_Rev") + R("R_ADSS") + R("R_CYTDK2") -
 R("R_GDPDPK") + R("R_GDPMNH") + R("R_GK1") - R("R_GK1_Rev") + 2
 R("R_GP4GH") - R("R_MAN1PT2") - R("R_NDPK1") + R("R_NDPK1_Rev") +
 R("R_NTP3") + R("R_PPGPPDP") - R("R_RNDR2") - R("R_RNDR2b") + R("R_SADT2") +
 R("R_URIK2")  = 0.

Mgdpddmanc: -R("R_GDMANE") + R("R_GMAND")  = 0.

Mgdpe: -R("R_EX_gdp_e_") + R("R_EX_gdp_e__Rev") - R("R_GDPtex") +
 R("R_GDPtex_Rev")  = 0.

Mgdpfucc: R("R_GOFUCR")  = 0.

Mgdpmannc: -R("R_GDPMNH") - R("R_GDPMNP") - R("R_GMAND") +
 R("R_MAN1PT2")  = 0.

Mgdpofucc: R("R_GDMANE") - R("R_GOFUCR")  = 0.

Mgdpp: R("R_GDPtex") - R("R_GDPtex_Rev") + R("R_NTP3pp")  = 0.

Mgdptpc: -R("R_GTPDPDP") + R("R_GTPDPK")  = 0.

Mgfgaragundc: R("R_O16GALFT") - R("R_O16GLCT2")  = 0.

Mgg4abutc: R("R_GGGABADr") - R("R_GGGABADr_Rev") - R("R_GGGABAH")  =
 0.

Mggagicolipac: R("R_GLCTR2") - R("R_GLCTR3")  = 0.

Mggbutalc: -R("R_GGGABADr") + R("R_GGGABADr_Rev") + R("R_GGPTRCO") 
 = 0.

Mgggagicolipac: R("R_GLCTR3") - R("R_HEPT4")  = 0.

Mggptrcc: -R("R_GGPTRCO") + R("R_GGPTRCS")  = 0.

Mgicolipac: -R("R_GALT1") + R("R_GLCTR1")  = 0.

MglcDb: R("MglcDb_Transp") + R("R_EX_glc_e_") - R("R_EX_glc_e__Rev")
  = 0.

MglcDc: R("R_AMALT1") + R("R_AMALT2") + R("R_AMALT3") +
 R("R_AMALT4") + R("R_G6PP") + R("R_GALS3") + R("R_GLCabcpp") - R("R_GLCATr") +
 R("R_GLCATr_Rev") + R("R_GLCt2pp") - R("R_HEX1") + R("R_LACZ") +
 R("R_MLTG1") + R("R_MLTG2") + R("R_MLTG3") + R("R_MLTG4") + R("R_MLTG5") +
 R("R_TRE6PH") + 2 R("R_TREH") - R("R_XYLI2") + R("R_XYLI2_Rev")  = 0.

MglcDe: -R("R_EX_glc_e_") + R("R_EX_glc_e__Rev") - R("R_GLCtex") -
 R("R_GLCtexi") + R("R_GLCtex_Rev")  = 0.

MglcDp: R("R_G1PPpp") - R("R_GLCabcpp") - R("R_GLCDpp") -
 R("R_GLCptspp") - R("R_GLCt2pp") + R("R_GLCtex") + R("R_GLCtexi") -
 R("R_GLCtex_Rev") + R("R_LACZpp") + 2 R("R_TREHpp")  = 0.

Mglcnb: R("Mglcnb_Transp") + R("R_EX_glcn_e_") -
 R("R_EX_glcn_e__Rev")  = 0.

Mglcnc: R("R_2DGLCNRx") + R("R_2DGLCNRy") + R("R_5DGLCNR") -
 R("R_5DGLCNR_Rev") + R("R_GLCNt2rpp") - R("R_GLCNt2rpp_Rev") - R("R_GNK")  =
 0.

Mglcne: -R("R_EX_glcn_e_") + R("R_EX_glcn_e__Rev") - R("R_GLCNtex")
 + R("R_GLCNtex_Rev")  = 0.

Mglcnp: R("R_GLCDpp") - R("R_GLCNt2rpp") + R("R_GLCNt2rpp_Rev") +
 R("R_GLCNtex") - R("R_GLCNtex_Rev")  = 0.

Mglcrb: R("Mglcrb_Transp") + R("R_EX_glcr_e_") -
 R("R_EX_glcr_e__Rev")  = 0.

Mglcrc: -R("R_GLCRD") + R("R_GLCRt2rpp") - R("R_GLCRt2rpp_Rev")  =
 0.

Mglcre: -R("R_EX_glcr_e_") + R("R_EX_glcr_e__Rev") - R("R_GLCRtex")
 + R("R_GLCRtex_Rev")  = 0.

Mglcrp: -R("R_GLCRt2rpp") + R("R_GLCRt2rpp_Rev") + R("R_GLCRtex") -
 R("R_GLCRtex_Rev")  = 0.

Mglcur1pb: R("Mglcur1pb_Transp") + R("R_EX_glcur1p_e_") -
 R("R_EX_glcur1p_e__Rev")  = 0.

Mglcur1pe: -R("R_EX_glcur1p_e_") + R("R_EX_glcur1p_e__Rev") -
 R("R_GLCUR1Ptex") + R("R_GLCUR1Ptex_Rev")  = 0.

Mglcur1pp: R("R_GLCUR1Ptex") - R("R_GLCUR1Ptex_Rev") -
 R("R_GUR1PPpp") + R("R_UGLCURPpp")  = 0.

Mglcurb: R("Mglcurb_Transp") + R("R_EX_glcur_e_") -
 R("R_EX_glcur_e__Rev")  = 0.

Mglcurc: R("R_GLCURt2rpp") - R("R_GLCURt2rpp_Rev") - R("R_GUI1") +
 R("R_GUI1_Rev")  = 0.

Mglcure: -R("R_EX_glcur_e_") + R("R_EX_glcur_e__Rev") -
 R("R_GLCURtex") + R("R_GLCURtex_Rev")  = 0.

Mglcurp: -R("R_GLCURt2rpp") + R("R_GLCURt2rpp_Rev") +
 R("R_GLCURtex") - R("R_GLCURtex_Rev") + R("R_GUR1PPpp")  = 0.

MglnLb: R("MglnLb_Transp") + R("R_EX_gln_L_e_") -
 R("R_EX_gln_L_e__Rev")  = 0.

MglnLc: -R("R_ADCS") - R("R_ANS") - R("R_ASNS1") - R("R_CBPS") -
 R("R_CTPS2") - 0.2632 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_GF6PTA") + R("R_GLNabcpp") + R("R_GLNS") - R("R_GLNTRS") - R("R_GLUN") -
 R("R_GLUPRT") - R("R_GLUSy") - R("R_GMPS2") - R("R_IG3PS") - R("R_PRFGS")
  = 0.

MglnLe: -R("R_EX_gln_L_e_") + R("R_EX_gln_L_e__Rev") - R("R_GLNtex")
 + R("R_GLNtex_Rev")  = 0.

MglnLp: -R("R_GLNabcpp") + R("R_GLNtex") - R("R_GLNtex_Rev") -
 R("R_GLUNpp")  = 0.

Mglntrnac: R("R_GLNTRS")  = 0.

Mglu1sac: -R("R_G1SAT") + R("R_G1SAT_Rev") + R("R_GLUTRR")  = 0.

Mglu5pc: -R("R_G5SD") + R("R_GLU5K")  = 0.

Mglu5sac: -R("R_G5SADs") + R("R_G5SD") + R("R_NACODA")  = 0.

Mglucysc: R("R_GLUCYS") - R("R_GTHS")  = 0.

MgluDc: -R("R_GLUR") + R("R_GLUR_Rev") - R("R_UAMAGS")  = 0.

MgluLb: R("MgluLb_Transp") + R("R_EX_glu_L_e_") -
 R("R_EX_glu_L_e__Rev")  = 0.

MgluLc: R("R_ABTA") - R("R_ACGS") + R("R_ACOTA") - R("R_ACOTA_Rev")
 + R("R_ADCS") + R("R_ALATA_L") - R("R_ALATA_L_Rev") + R("R_ANS") +
 R("R_ASNS1") + R("R_ASPTA") - R("R_ASPTA_Rev") + R("R_CBPS") +
 R("R_CTPS2") - R("R_DHFS") - 0.2632 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_GF6PTA") + R("R_GGGABAH") - R("R_GGPTRCS") - R("R_GLNS") -
 R("R_GLU5K") + R("R_GLUabcpp") + R("R_GLUABUTt7pp") - R("R_GLUABUTt7pp_Rev") -
 R("R_GLUCYS") - R("R_GLUDC") - R("R_GLUDy") + R("R_GLUDy_Rev") +
 R("R_GLUN") + R("R_GLUPRT") + R("R_GLUR") - R("R_GLUR_Rev") + 2 R("R_GLUSy")
 + R("R_GLUt2rpp") - R("R_GLUt2rpp_Rev") + R("R_GLUt4pp") -
 R("R_GLUTRS") + R("R_GMPS2") - R("R_HSTPT") + R("R_IG3PS") + R("R_ILETA") -
 R("R_ILETA_Rev") + R("R_LALGP") - R("R_LEUTAi") - R("R_OHPBAT") +
 R("R_OHPBAT_Rev") + R("R_P5CD") + R("R_PHETA1") - R("R_PHETA1_Rev") +
 R("R_PRFGS") - R("R_PSERT") + R("R_PTRCTA") + R("R_SDPTA") - R("R_SDPTA_Rev")
 + R("R_SGDS") + R("R_SOTA") - R("R_TDPAGTA") + R("R_TYRTA") -
 R("R_TYRTA_Rev") - R("R_UDPKAAT") + R("R_UDPKAAT_Rev") - R("R_UNK3") +
 R("R_VALTA") - R("R_VALTA_Rev")  = 0.

MgluLe: -R("R_EX_glu_L_e_") + R("R_EX_glu_L_e__Rev") - R("R_GLUtex")
 + R("R_GLUtex_Rev")  = 0.

MgluLp: -R("R_GLUabcpp") - R("R_GLUABUTt7pp") +
 R("R_GLUABUTt7pp_Rev") + R("R_GLUNpp") - R("R_GLUt2rpp") + R("R_GLUt2rpp_Rev") -
 R("R_GLUt4pp") + R("R_GLUtex") - R("R_GLUtex_Rev") + R("R_GTHRDHpp")  = 0.

Mglutrnac: -R("R_GLUTRR") + R("R_GLUTRS")  = 0.

Mglxc: R("R_FDMO6") - 2 R("R_GLXCL") - R("R_GLYCLTDx") -
 R("R_GLYCLTDy") + R("R_GLYCTO2") + R("R_GLYCTO3") + R("R_GLYCTO4") + R("R_ICL") -
 R("R_MALS") + R("R_UGLYCH")  = 0.

Mglyaldb: R("Mglyaldb_Transp") + R("R_EX_glyald_e_") -
 R("R_EX_glyald_e__Rev")  = 0.

Mglyaldc: -R("R_ALCD19") + R("R_ALCD19_Rev") + R("R_GLYALDtpp") -
 R("R_GLYALDtpp_Rev")  = 0.

Mglyalde: -R("R_EX_glyald_e_") + R("R_EX_glyald_e__Rev") -
 R("R_GLYALDtex") + R("R_GLYALDtex_Rev")  = 0.

Mglyaldp: R("R_GLYALDtex") - R("R_GLYALDtex_Rev") - R("R_GLYALDtpp")
 + R("R_GLYALDtpp_Rev")  = 0.

Mglyb: R("Mglyb_Transp") + R("R_EX_gly_e_") - R("R_EX_gly_e__Rev") 
 = 0.

Mglybb: R("Mglybb_Transp") + R("R_EX_glyb_e_") -
 R("R_EX_glyb_e__Rev")  = 0.

Mglybc: R("R_BETALDHx") + R("R_BETALDHy") + R("R_GLYBabcpp") +
 R("R_GLYBt2pp")  = 0.

Mglybe: -R("R_EX_glyb_e_") + R("R_EX_glyb_e__Rev") - R("R_GLYBtex")
 + R("R_GLYBtex_Rev")  = 0.

Mglybp: -R("R_GLYBabcpp") - R("R_GLYBt2pp") + R("R_GLYBtex") -
 R("R_GLYBtex_Rev")  = 0.

Mglyc: R("R_AMPTASECG") + R("R_AMPTASEPG") - 0.6126
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_GHMT2r") - R("R_GHMT2r_Rev") -
 R("R_GLYAT") + R("R_GLYAT_Rev") - R("R_GLYCL") + R("R_GLYt2pp") + R("R_GLYt4pp")
 - R("R_GLYTRS") - R("R_GTHS") - R("R_PRAGSr") + R("R_PRAGSr_Rev") +
 R("R_SARCOX") + R("R_THRA2i") + R("R_THRAi")  = 0.

Mglyc2pb: R("Mglyc2pb_Transp") + R("R_EX_glyc2p_e_") -
 R("R_EX_glyc2p_e__Rev")  = 0.

Mglyc2pe: -R("R_EX_glyc2p_e_") + R("R_EX_glyc2p_e__Rev") -
 R("R_GLYC2Ptex") + R("R_GLYC2Ptex_Rev")  = 0.

Mglyc2pp: -R("R_G2PPpp") + R("R_GLYC2Ptex") - R("R_GLYC2Ptex_Rev") 
 = 0.

Mglyc3pb: R("Mglyc3pb_Transp") + R("R_EX_glyc3p_e_") -
 R("R_EX_glyc3p_e__Rev")  = 0.

Mglyc3pc: -R("R_G3PAT120") - R("R_G3PAT140") - R("R_G3PAT141") -
 R("R_G3PAT160") - R("R_G3PAT161") - R("R_G3PAT180") - R("R_G3PAT181") -
 R("R_G3PD2") + R("R_G3PD2_Rev") - R("R_G3PD5") - R("R_G3PD6") -
 R("R_G3PD7") - R("R_G3PT") + R("R_GLYC3Pabcpp") + R("R_GLYC3Pt6pp") +
 R("R_GLYK") + R("R_GPDDA1") + R("R_GPDDA2") + R("R_GPDDA3") + R("R_GPDDA4") +
 R("R_GPDDA5") + R("R_LPLIPAL2A120") + R("R_LPLIPAL2A140") +
 R("R_LPLIPAL2A141") + R("R_LPLIPAL2A160") + R("R_LPLIPAL2A161") +
 R("R_LPLIPAL2A180") + R("R_LPLIPAL2A181") - R("R_PGSA120") - R("R_PGSA140") -
 R("R_PGSA141") - R("R_PGSA160") - R("R_PGSA161") - R("R_PGSA180") -
 R("R_PGSA181")  = 0.

Mglyc3pe: -R("R_EX_glyc3p_e_") + R("R_EX_glyc3p_e__Rev") -
 R("R_GLYC3Ptex") + R("R_GLYC3Ptex_Rev")  = 0.

Mglyc3pp: -R("R_GLYC3Pabcpp") - R("R_GLYC3Pt6pp") + R("R_GLYC3Ptex")
 - R("R_GLYC3Ptex_Rev") + R("R_GPDDA1pp") + R("R_GPDDA2pp") +
 R("R_GPDDA3pp") + R("R_GPDDA4pp") + R("R_GPDDA5pp") + R("R_LPLIPAL1A120pp") +
 R("R_LPLIPAL1A140pp") + R("R_LPLIPAL1A141pp") + R("R_LPLIPAL1A160pp")
 + R("R_LPLIPAL1A161pp") + R("R_LPLIPAL1A180pp") +
 R("R_LPLIPAL1A181pp")  = 0.

Mglycb: R("Mglycb_Transp") + R("R_EX_glyc_e_") -
 R("R_EX_glyc_e__Rev")  = 0.

Mglycc: R("R_ALCD19") - R("R_ALCD19_Rev") + R("R_G3PT") -
 R("R_GLYCDx") - R("R_GLYCtpp") + R("R_GLYCtpp_Rev") - R("R_GLYK") +
 R("R_GPDDA4")  = 0.

Mglyce: -R("R_EX_glyc_e_") + R("R_EX_glyc_e__Rev") - R("R_GLYCtex")
 + R("R_GLYCtex_Rev")  = 0.

Mglycltb: R("Mglycltb_Transp") + R("R_EX_glyclt_e_") -
 R("R_EX_glyclt_e__Rev")  = 0.

Mglycltc: R("R_GCALDD") + R("R_GLYCLTDx") + R("R_GLYCLTDy") +
 R("R_GLYCLTt2rpp") - R("R_GLYCLTt2rpp_Rev") + R("R_GLYCLTt4pp") -
 R("R_GLYCTO2") - R("R_GLYCTO3") - R("R_GLYCTO4") + R("R_PGLYCP")  = 0.

Mglyclte: -R("R_EX_glyclt_e_") + R("R_EX_glyclt_e__Rev") -
 R("R_GLYCLTtex") + R("R_GLYCLTtex_Rev")  = 0.

Mglycltp: -R("R_GLYCLTt2rpp") + R("R_GLYCLTt2rpp_Rev") -
 R("R_GLYCLTt4pp") + R("R_GLYCLTtex") - R("R_GLYCLTtex_Rev")  = 0.

Mglycogenc: -R("R_GLBRAN2") - R("R_GLCP") + R("R_GLCS1") +
 R("R_GLDBRAN2")  = 0.

Mglycp: R("R_CLPNS120pp") - R("R_CLPNS120pp_Rev") +
 R("R_CLPNS140pp") - R("R_CLPNS140pp_Rev") + R("R_CLPNS141pp") - R("R_CLPNS141pp_Rev")
 + R("R_CLPNS160pp") - R("R_CLPNS160pp_Rev") + R("R_CLPNS161pp") -
 R("R_CLPNS161pp_Rev") + R("R_CLPNS180pp") - R("R_CLPNS180pp_Rev") +
 R("R_CLPNS181pp") - R("R_CLPNS181pp_Rev") + R("R_G2PPpp") + R("R_GLYCtex") -
 R("R_GLYCtex_Rev") + R("R_GLYCtpp") - R("R_GLYCtpp_Rev") +
 R("R_GPDDA4pp")  = 0.

MglycRb: R("MglycRb_Transp") + R("R_EX_glyc_R_e_") -
 R("R_EX_glyc_R_e__Rev")  = 0.

MglycRc: R("R_GLYCAt2rpp") - R("R_GLYCAt2rpp_Rev") - R("R_GLYCK") -
 R("R_GLYCK2") + R("R_HPYRRx") + R("R_HPYRRy") + R("R_MANPGH") +
 R("R_TRSARr") - R("R_TRSARr_Rev")  = 0.

MglycRe: -R("R_EX_glyc_R_e_") + R("R_EX_glyc_R_e__Rev") -
 R("R_GLYCAtex") + R("R_GLYCAtex_Rev")  = 0.

MglycRp: -R("R_GLYCAt2rpp") + R("R_GLYCAt2rpp_Rev") +
 R("R_GLYCAtex") - R("R_GLYCAtex_Rev")  = 0.

Mglye: -R("R_EX_gly_e_") + R("R_EX_gly_e__Rev") - R("R_GLYtex") +
 R("R_GLYtex_Rev")  = 0.

Mglyp: -R("R_GLYt2pp") - R("R_GLYt4pp") + R("R_GLYtex") -
 R("R_GLYtex_Rev")  = 0.

Mglytrnac: R("R_GLYTRS")  = 0.

Mgmhep17bpc: R("R_GMHEPK") - R("R_GMHEPPA")  = 0.

Mgmhep1pc: -R("R_GMHEPAT") + R("R_GMHEPPA")  = 0.

Mgmhep7pc: -R("R_GMHEPK") + R("R_S7PI")  = 0.

Mgmpb: R("Mgmpb_Transp") + R("R_EX_gmp_e_") - R("R_EX_gmp_e__Rev") 
 = 0.

Mgmpc: R("R_ADOCBLS") + R("R_GDPMNP") - R("R_GK1") + R("R_GK1_Rev")
 - R("R_GMPR") + R("R_GMPS2") + R("R_GSNK") + R("R_GUAPRT") -
 R("R_NTD9") + R("R_NTPP2")  = 0.

Mgmpe: -R("R_EX_gmp_e_") + R("R_EX_gmp_e__Rev") - R("R_GMPtex") +
 R("R_GMPtex_Rev")  = 0.

Mgmpp: R("R_GMPtex") - R("R_GMPtex_Rev") - R("R_NTD9pp")  = 0.

Mgp4gc: -R("R_GP4GH")  = 0.

Mgrdpc: R("R_DMATT") - R("R_GRTT")  = 0.

Mgrxoxc: -R("R_GRXR") + R("R_PAPSR2") + R("R_RNDR1b") +
 R("R_RNDR2b") + R("R_RNDR3b") + R("R_RNDR4b")  = 0.

Mgrxrdc: R("R_GRXR") - R("R_PAPSR2") - R("R_RNDR1b") - R("R_RNDR2b")
 - R("R_RNDR3b") - R("R_RNDR4b")  = 0.

Mgsnb: R("Mgsnb_Transp") + R("R_EX_gsn_e_") - R("R_EX_gsn_e__Rev") 
 = 0.

Mgsnc: -R("R_GSNK") + R("R_GSNt2pp") + R("R_NTD9") + R("R_NTPTP2") -
 R("R_PUNP3") + R("R_PUNP3_Rev")  = 0.

Mgsne: -R("R_EX_gsn_e_") + R("R_EX_gsn_e__Rev") - R("R_GSNtex") +
 R("R_GSNtex_Rev")  = 0.

Mgsnp: R("R_3NTD9pp") - R("R_GSNt2pp") + R("R_GSNtex") -
 R("R_GSNtex_Rev") + R("R_NTD9pp")  = 0.

Mgthoxb: R("Mgthoxb_Transp") + R("R_EX_gthox_e_") -
 R("R_EX_gthox_e__Rev")  = 0.

Mgthoxc: R("R_ASR") + R("R_GRXR") - R("R_GTHOr") + R("R_GTHOr_Rev") 
 = 0.

Mgthoxe: -R("R_EX_gthox_e_") + R("R_EX_gthox_e__Rev") -
 R("R_GTHOXtex") + R("R_GTHOXtex_Rev")  = 0.

Mgthoxp: R("R_DSBCGT") + R("R_DSBGGT") + R("R_GTHOXtex") -
 R("R_GTHOXtex_Rev")  = 0.

Mgthrdb: R("Mgthrdb_Transp") + R("R_EX_gthrd_e_") -
 R("R_EX_gthrd_e__Rev")  = 0.

Mgthrdc: -2 R("R_ASR") - R("R_FALGTHLs") + R("R_FALGTHLs_Rev") +
 R("R_GLYOX") - 2 R("R_GRXR") + R("R_GSPMDA") - R("R_GSPMDS") + 2
 R("R_GTHOr") - 2 R("R_GTHOr_Rev") - R("R_GTHRDabc2pp") + R("R_GTHRDabcpp") +
 R("R_GTHS") - R("R_LGTHL") + R("R_SFGTHi")  = 0.

Mgthrde: -R("R_EX_gthrd_e_") + R("R_EX_gthrd_e__Rev") -
 R("R_GTHRDtex") + R("R_GTHRDtex_Rev")  = 0.

Mgthrdp: -2 R("R_DSBCGT") - 2 R("R_DSBGGT") + R("R_GTHRDabc2pp") -
 R("R_GTHRDabcpp") - R("R_GTHRDHpp") + R("R_GTHRDtex") -
 R("R_GTHRDtex_Rev")  = 0.

Mgtpb: R("Mgtpb_Transp") + R("R_EX_gtp_e_") - R("R_EX_gtp_e__Rev") 
 = 0.

Mgtpc: -R("R_ACBIPGT") - R("R_ADK3") + R("R_ADK3_Rev") - R("R_ADSS")
 - R("R_CYTDK2") - 0.2151 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_GTPCI") - R("R_GTPCII2") - R("R_GTPDPK") - R("R_GTPHs") +
 R("R_NDPK1") - R("R_NDPK1_Rev") - R("R_NTP3") - R("R_NTPP2") - R("R_NTPTP2") -
 R("R_RNTR2c") - R("R_SADT2") - R("R_URIK2")  = 0.

Mgtpe: -R("R_EX_gtp_e_") + R("R_EX_gtp_e__Rev") - R("R_GTPtex") +
 R("R_GTPtex_Rev")  = 0.

Mgtpp: R("R_GTPtex") - R("R_GTPtex_Rev") - R("R_NTP3pp")  = 0.

Mgtspmdc: -R("R_GSPMDA") + R("R_GSPMDS")  = 0.

Mguab: R("Mguab_Transp") + R("R_EX_gua_e_") - R("R_EX_gua_e__Rev") 
 = 0.

Mguac: -R("R_GUAD") - R("R_GUAPRT") + R("R_GUAt2pp") + R("R_GUAtpp")
 - R("R_GUAtpp_Rev") + R("R_PUNP3") - R("R_PUNP3_Rev") + R("R_PUNP4")
 - R("R_PUNP4_Rev")  = 0.

Mguae: -R("R_EX_gua_e_") + R("R_EX_gua_e__Rev") - R("R_GUAtex") +
 R("R_GUAtex_Rev")  = 0.

Mguap: -R("R_GUAt2pp") + R("R_GUAtex") - R("R_GUAtex_Rev") -
 R("R_GUAtpp") + R("R_GUAtpp_Rev")  = 0.

Mh2b: R("Mh2b_Transp") + R("R_EX_h2_e_") - R("R_EX_h2_e__Rev")  = 0.

Mh2c: R("R_FHL") + R("R_H2tpp") - R("R_H2tpp_Rev") - R("R_HYD1pp") -
 R("R_HYD2pp") - R("R_HYD3pp")  = 0.

Mh2e: -R("R_EX_h2_e_") + R("R_EX_h2_e__Rev") - R("R_H2tex") +
 R("R_H2tex_Rev")  = 0.

Mh2mb4pc: -R("R_DMPPS") - R("R_IPDPS") + R("R_MECDPDH2")  = 0.

Mh2o2b: R("Mh2o2b_Transp") + R("R_EX_h2o2_e_") -
 R("R_EX_h2o2_e__Rev")  = 0.

Mh2o2c: R("R_AACTOOR") + R("R_ASPO6") - 2 R("R_CAT") +
 R("R_GGPTRCO") - R("R_METOX1s") - R("R_METOX2s") + R("R_MTRPOX") + R("R_PDX5POi") +
 R("R_PYAM5PO") + R("R_SARCOX") + R("R_SPODM") - R("R_THIORDXi") +
 R("R_URIC")  = 0.

Mh2o2e: -R("R_EX_h2o2_e_") + R("R_EX_h2o2_e__Rev") - R("R_H2O2tex")
 + R("R_H2O2tex_Rev")  = 0.

Mh2o2p: R("R_42A12BOOXpp") + R("R_H2O2tex") - R("R_H2O2tex_Rev") +
 R("R_PEAMNOpp") + R("R_SPODMpp") + R("R_TYROXDApp")  = 0.

Mh2ob: R("Mh2ob_Transp") + R("R_EX_h2o_e_") - R("R_EX_h2o_e__Rev") 
 = 0.

Mh2oc: -R("R_14GLUCANabcpp") - R("R_2MAHMP") + R("R_3HAD100") +
 R("R_3HAD120") + R("R_3HAD121") + R("R_3HAD140") + R("R_3HAD141") +
 R("R_3HAD160") + R("R_3HAD161") + R("R_3HAD180") + R("R_3HAD181") +
 R("R_3HAD40") + R("R_3HAD60") + R("R_3HAD80") + R("R_3HCINNMH") +
 R("R_3HPPPNH") - R("R_3PEPTabcpp") - R("R_4HTHRS") - R("R_4PCP") -
 R("R_4PEPTabcpp") - R("R_5DOAN") - R("R_AACTOOR") - R("R_AB6PGH") - R("R_ABUTD") -
 R("R_ACM6PH") - R("R_ACODA") + R("R_ACONTa") - R("R_ACONTa_Rev") -
 R("R_ACONTb") + R("R_ACONTb_Rev") - R("R_ADA") - R("R_ADD") - R("R_ADNUC")
 - R("R_ADOCBLabcpp") - R("R_ADPRDP") - R("R_AGDC") - R("R_AGM3PA") -
 R("R_AGM3PH") - R("R_AGM4PA") - R("R_AGM4PCP") - R("R_AGM4PH") -
 R("R_AGMH") - R("R_AGMT") - R("R_AHCYSNS") - R("R_ALAabcpp") -
 R("R_ALAALAabcpp") - R("R_ALAALAD") - R("R_ALDD19x") - R("R_ALDD2x") -
 R("R_ALDD2y") - R("R_ALDD3y") - R("R_ALDD4") - R("R_ALLabcpp") - 2 R("R_ALLTAMH")
 - R("R_ALLTN") + R("R_ALTRH") - R("R_AM3PA") - R("R_AM4PA") -
 R("R_AM4PCP") - R("R_AMPMS2") - R("R_AMPN") - R("R_AMPTASECG") -
 R("R_AMPTASEPG") - R("R_ANHMK") - R("R_AP4AH") - R("R_AP5AH") - R("R_ARBabcpp") -
 R("R_ARBTNabcpp") - R("R_ARGabcpp") - R("R_ASCBPL") - R("R_ASNabcpp") -
 R("R_ASNN") - R("R_ASNS1") - R("R_ASO3t8pp") - R("R_ASPabcpp") +
 R("R_ASR") - R("R_ATPHs") - R("R_ATPM") + R("R_ATPS4rpp") -
 R("R_ATPS4rpp_Rev") - R("R_BETALDHx") - R("R_BETALDHy") - R("R_BPNT") + R("R_BSORx")
 + R("R_BSORy") - R("R_BUTSO3abcpp") + 2 R("R_CAT") - R("R_CBIuabcpp")
 - R("R_CBL1abcpp") - R("R_CBPS") - R("R_CD2abcpp") - R("R_CDAPPA120")
 - R("R_CDAPPA140") - R("R_CDAPPA141") - R("R_CDAPPA160") -
 R("R_CDAPPA161") - R("R_CDAPPA180") - R("R_CDAPPA181") - R("R_CGLYabcpp") -
 R("R_CHLabcpp") - R("R_CMPN") - R("R_COBALT2abcpp") - R("R_COLIPAabcpp") -
 R("R_CPGNabcpp") + 2 R("R_CPPPGO") - R("R_CRNabcpp") + R("R_CRNCDH")
 - R("R_CRNCDH_Rev") - R("R_CRNDabcpp") - R("R_CS") - R("R_CSND") -
 R("R_CTBTabcpp") - R("R_CTPS2") - R("R_CU1abcpp") - R("R_CU2abcpp") -
 R("R_CYSabc2pp") - R("R_CYSabcpp") - R("R_CYSDabcpp") - R("R_CYSDDS") -
 R("R_CYSDS") - R("R_CYSTL") + R("R_CYTBD2pp") + R("R_CYTBDpp") +
 R("R_CYTBO3_4pp") - R("R_CYTD") - R("R_CYTDH") - R("R_DAAD") - R("R_DADA") -
 R("R_DAPabcpp") - R("R_DAPAL") - R("R_DATPHs") - R("R_DCTPD") -
 R("R_DCYTD") - R("R_DDPA") + R("R_DHAD1") + R("R_DHAD2") - R("R_DHBSH") + 2
 R("R_DHDPS") - R("R_DHORTS") + R("R_DHORTS_Rev") - R("R_DHPPDA2") +
 R("R_DHPTDCs") + R("R_DHQTi") - R("R_DKMPPD") - 3 R("R_DKMPPD3") +
 R("R_DMPPS") + R("R_DMSOR1") + R("R_DMSOR2") - R("R_DNMPPA") -
 R("R_DNTPPA") - R("R_DUTPDP") - R("R_E4PD") + R("R_E4PD_Rev") - 54.462
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_ECOAH1") - R("R_ECOAH1_Rev") +
 R("R_ECOAH2") - R("R_ECOAH2_Rev") + R("R_ECOAH3") - R("R_ECOAH3_Rev") +
 R("R_ECOAH4") - R("R_ECOAH4_Rev") + R("R_ECOAH5") - R("R_ECOAH5_Rev") +
 R("R_ECOAH6") - R("R_ECOAH6_Rev") + R("R_ECOAH7") - R("R_ECOAH7_Rev")
 + R("R_ECOAH8") - R("R_ECOAH8_Rev") + R("R_EDD") + R("R_ENO") -
 R("R_ENO_Rev") - 3 R("R_ENTERES") - 3 R("R_ENTERES2") - R("R_ETHSO3abcpp") -
 R("R_F6PP") - R("R_FA100ACPHi") - R("R_FA120ACPHi") -
 R("R_FA140ACPHi") - R("R_FA141ACPHi") - R("R_FA160ACPHi") - R("R_FA161ACPHi") -
 R("R_FA80ACPHi") - R("R_FACOAE100") - R("R_FACOAE120") - R("R_FACOAE140") -
 R("R_FACOAE141") - R("R_FACOAE160") - R("R_FACOAE161") -
 R("R_FACOAE180") - R("R_FACOAE181") - R("R_FACOAE60") - R("R_FACOAE80") -
 R("R_FBP") + R("R_FDMO") + R("R_FDMO2") + R("R_FDMO3") + R("R_FDMO4") +
 R("R_FDMO6") - R("R_FE2abcpp") - R("R_FE3abcpp") - R("R_FE3DCITabcpp") -
 R("R_FE3DHBZSabcpp") - R("R_FE3HOXabcpp") - R("R_FECRMabcpp") -
 R("R_FEENTERabcpp") - R("R_FEOXAMabcpp") - R("R_FFSD") - R("R_FRULYSDG") +
 R("R_FRULYSDG_Rev") - R("R_FTHFD") - R("R_FUM") + R("R_FUM_Rev") -
 R("R_G3PCabcpp") - R("R_G3PEabcpp") - R("R_G3PGabcpp") - R("R_G3PIabcpp") -
 R("R_G3PSabcpp") - R("R_G3PT") + R("R_G5SADs") - R("R_G6PDA") -
 R("R_G6PP") - R("R_GALabcpp") + R("R_GALCTD") + R("R_GALCTND") - R("R_GALS3")
 - R("R_GCALDD") - R("R_GDPMNH") - R("R_GDPMNP") - R("R_GGGABADr") +
 R("R_GGGABADr_Rev") - R("R_GGGABAH") - R("R_GGPTRCO") + R("R_GHMT2r") -
 R("R_GHMT2r_Rev") - R("R_GLCabcpp") + R("R_GLCRD") - R("R_GLNabcpp")
 - R("R_GLUabcpp") - R("R_GLUDy") + R("R_GLUDy_Rev") - R("R_GLUN") -
 R("R_GLUPRT") - R("R_GLYBabcpp") - R("R_GLYC3Pabcpp") - R("R_GLYOX") -
 R("R_GLYOX3") + R("R_GMAND") - R("R_GMHEPPA") - R("R_GMPS2") -
 R("R_GP4GH") - R("R_GPDDA1") - R("R_GPDDA2") - R("R_GPDDA3") - R("R_GPDDA4") -
 R("R_GPDDA5") - R("R_GSPMDA") - R("R_GTHRDabc2pp") -
 R("R_GTHRDabcpp") - R("R_GTPCI") - 3 R("R_GTPCII2") - R("R_GTPDPDP") - R("R_GTPHs") -
 R("R_GUAD") + R("R_H2Otpp") - R("R_H2Otpp_Rev") - R("R_HCO3E") +
 R("R_HCO3E_Rev") - R("R_HEMEOS") - R("R_HG2abcpp") - R("R_HISabcpp") -
 R("R_HISTD") - R("R_HISTP") - R("R_HKNDDH") - R("R_HKNTDH") - R("R_HMBS")
 - R("R_HXAND") - R("R_HYPOE") - R("R_ICHORT") + R("R_IGPDH") +
 R("R_IGPS") - R("R_ILEabcpp") - R("R_IMPC") + R("R_IMPC_Rev") - R("R_IMPD") -
 R("R_INSH") + R("R_IPDPS") + R("R_IPPMIa") - R("R_IPPMIa_Rev") -
 R("R_IPPMIb") + R("R_IPPMIb_Rev") - R("R_IPPS") - R("R_ISETACabcpp") -
 R("R_K2L4Aabcpp") - R("R_Kabcpp") - R("R_KDOPP") - R("R_KDOPS") -
 R("R_LACZ") - R("R_LADGMDH") - R("R_LALGP") - R("R_LCADi") - R("R_LEUabcpp")
 - R("R_LIPAabcpp") - R("R_LIPACabcpp") - R("R_LPLIPAL2A120") -
 R("R_LPLIPAL2A140") - R("R_LPLIPAL2A141") - R("R_LPLIPAL2A160") -
 R("R_LPLIPAL2A161") - R("R_LPLIPAL2A180") - R("R_LPLIPAL2A181") -
 R("R_LPLIPAL2E120") - R("R_LPLIPAL2E140") - R("R_LPLIPAL2E141") - R("R_LPLIPAL2E160")
 - R("R_LPLIPAL2E161") - R("R_LPLIPAL2E180") - R("R_LPLIPAL2E181") -
 R("R_LPLIPAL2G120") - R("R_LPLIPAL2G140") - R("R_LPLIPAL2G141") -
 R("R_LPLIPAL2G160") - R("R_LPLIPAL2G161") - R("R_LPLIPAL2G180") -
 R("R_LPLIPAL2G181") - R("R_LYSabcpp") - R("R_MALS") - R("R_MALTabcpp") -
 R("R_MALTHXabcpp") - R("R_MALTPTabcpp") - R("R_MALTTRabcpp") -
 R("R_MALTTTRabcpp") - R("R_MANPGH") + R("R_MCITD") - R("R_MCITS") + R("R_MDRPD") +
 R("R_MECDPDH2") - R("R_METabcpp") - R("R_METAT") - R("R_METDabcpp") +
 R("R_METOX1s") + R("R_METOX2s") - R("R_METSOX1abcpp") -
 R("R_METSOX2abcpp") + R("R_METSOXR1") + R("R_METSOXR2") - R("R_MG2uabcpp") -
 R("R_MI1PP") - R("R_MICITD") - R("R_MLTG1") - R("R_MLTG2") - R("R_MLTG3") -
 R("R_MLTG4") - R("R_MLTG5") - R("R_MN6PP") + R("R_MNNH") -
 R("R_MOBDabcpp") - R("R_MOHMT") - R("R_MSO3abcpp") - R("R_MTAN") - R("R_MTHFC") +
 R("R_MTHFC_Rev") - R("R_MTRPOX") - R("R_NACODA") - R("R_NADDP") -
 R("R_NADN") - R("R_NADPPPS") - R("R_NAMNPP") + R("R_NHFRBO") -
 R("R_NI2abcpp") - R("R_NI2uabcpp") - R("R_NMNDA") - R("R_NMNN") - R("R_NMNt7pp")
 - R("R_NNAM") + R("R_NO3R1pp") + R("R_NO3R2pp") - R("R_NTD1") -
 R("R_NTD10") - R("R_NTD11") - R("R_NTD12") - R("R_NTD2") - R("R_NTD3") -
 R("R_NTD4") - R("R_NTD5") - R("R_NTD6") - R("R_NTD7") - R("R_NTD8") -
 R("R_NTD9") - R("R_NTP1") - R("R_NTP10") - R("R_NTP11") - R("R_NTP12") -
 R("R_NTP3") - R("R_NTP5") - R("R_NTPP1") - R("R_NTPP10") -
 R("R_NTPP11") - R("R_NTPP2") - R("R_NTPP3") - R("R_NTPP4") - R("R_NTPP5") -
 R("R_NTPP6") - R("R_NTPP7") - R("R_NTPP8") - R("R_NTPP9") - R("R_NTPTP1") -
 R("R_NTPTP2") + 2 R("R_NTRIR2x") - 3 R("R_OMMBLHX3") - 3
 R("R_OMPHHX3") - R("R_OP4ENH") - 3 R("R_OPHHX3") - R("R_ORNabcpp") - 2 R("R_P5CD")
 - R("R_PA120abcpp") - R("R_PA140abcpp") - R("R_PA141abcpp") -
 R("R_PA160abcpp") - R("R_PA161abcpp") - R("R_PA180abcpp") - R("R_PA181abcpp")
 - R("R_PAPA120") - R("R_PAPA140") - R("R_PAPA141") - R("R_PAPA160") -
 R("R_PAPA161") - R("R_PAPA180") - R("R_PAPA181") - R("R_PDE1") + 2
 R("R_PDX5PS") - R("R_PDXPP") - R("R_PE120abcpp") - R("R_PE140abcpp") -
 R("R_PE141abcpp") - R("R_PE160abcpp") - R("R_PE161abcpp") -
 R("R_PE180abcpp") - R("R_PE181abcpp") - R("R_PG120abcpp") - R("R_PG140abcpp") -
 R("R_PG141abcpp") - R("R_PG160abcpp") - R("R_PG161abcpp") -
 R("R_PG180abcpp") - R("R_PG181abcpp") - R("R_PGL") - R("R_PGLYCP") -
 R("R_PGP120abcpp") - R("R_PGP140abcpp") - R("R_PGP141abcpp") - R("R_PGP160abcpp")
 - R("R_PGP161abcpp") - R("R_PGP180abcpp") - R("R_PGP181abcpp") -
 R("R_PGPP120") - R("R_PGPP140") - R("R_PGPP141") - R("R_PGPP160") -
 R("R_PGPP161") - R("R_PGPP180") - R("R_PGPP181") - R("R_PHEMEabcpp") -
 R("R_PIuabcpp") - R("R_PMDPHT") - R("R_POX") - R("R_PPA") - R("R_PPA2") + 2
 R("R_PPBNGS") - R("R_PPC") - R("R_PPGPPDP") + R("R_PPNDH") + 3
 R("R_PPPGO") - R("R_PPS") - R("R_PRAMPC") - R("R_PRATPP") - R("R_PRFGS") -
 R("R_PROabcpp") - R("R_PROGLYabcpp") - R("R_PSP_L") - R("R_PTRCabcpp") -
 R("R_PYAM5PO") - R("R_PYDXPP") + 2 R("R_QULNS") - R("R_R5PP") + 2
 R("R_RBFSa") - R("R_RIBabcpp") + R("R_RNDR1") + R("R_RNDR1b") +
 R("R_RNDR2") + R("R_RNDR2b") + R("R_RNDR3") + R("R_RNDR3b") + R("R_RNDR4") +
 R("R_RNDR4b") + R("R_RNTR1c") + R("R_RNTR2c") + R("R_RNTR3c") +
 R("R_RNTR4c") - R("R_RZ5PP") - 2 R("R_SADH") - R("R_SADT2") - R("R_SARCOX") -
 R("R_SDPDS") - R("R_SELNPS") - R("R_SFGTHi") - R("R_SGDS") -
 R("R_SGSAD") - R("R_SPMDabcpp") - R("R_SSALx") - R("R_SSALy") + R("R_SUCBZS") -
 R("R_SULabcpp") - R("R_SULFACabcpp") + 3 R("R_SULRi") + R("R_TARTD")
 - R("R_TAURabcpp") - R("R_TDP") + R("R_TDPGDH") - R("R_THDPS") + 2
 R("R_THIORDXi") - R("R_THMabcpp") - R("R_THRabcpp") - R("R_THRS") +
 R("R_THZPSN") + R("R_TMAOR1") + R("R_TMAOR2") - R("R_TRE6PH") -
 R("R_TRE6PP") - R("R_TREH") - R("R_TRPAS2") + R("R_TRPAS2_Rev") + R("R_TRPS1") +
 R("R_TRPS2") - R("R_TSULabcpp") - R("R_TUNGSabcpp") - R("R_UACMAMO") -
 R("R_UDCPDP") - R("R_UDPGD") - R("R_UGLYCH") - R("R_UHGADA") -
 R("R_UM4PCP") - R("R_UPLA4FNF") + R("R_UPP3S") - 2 R("R_URIC") - R("R_URIH")
 - R("R_USHD") - R("R_VALabcpp") - R("R_XAND") - R("R_XTSNH") -
 R("R_XYLabcpp") - R("R_ZN2abcpp") - R("R_ZNabcpp")  = 0.

Mh2oe: -R("R_EX_h2o_e_") + R("R_EX_h2o_e__Rev") - R("R_H2Otex") +
 R("R_H2Otex_Rev") + R("R_LIPAHT2ex") + R("R_LIPAHTex")  = 0.

Mh2op: -R("R_23PDE2pp") - R("R_23PDE4pp") - R("R_23PDE7pp") -
 R("R_23PDE9pp") - R("R_3NTD2pp") - R("R_3NTD4pp") - R("R_3NTD7pp") -
 R("R_3NTD9pp") - R("R_42A12BOOXpp") - R("R_4PCPpp") - R("R_ACGAL1PPpp") -
 R("R_ACGAM1PPpp") - R("R_AGM3PApp") - R("R_AGM4PApp") - R("R_AGM4PCPpp")
 - R("R_ASNNpp") - R("R_CLPNH120pp") - R("R_CLPNH140pp") -
 R("R_CLPNH141pp") - R("R_CLPNH160pp") - R("R_CLPNH161pp") - R("R_CLPNH180pp") -
 R("R_CLPNH181pp") + 2 R("R_CU1Opp") + R("R_DMSOR1pp") + R("R_DMSOR2pp")
 + 2 R("R_FEROpp") - R("R_G1PPpp") - R("R_G2PPpp") - R("R_GAL1PPpp") -
 R("R_GLCDpp") - R("R_GLUNpp") - R("R_GPDDA1pp") - R("R_GPDDA2pp") -
 R("R_GPDDA3pp") - R("R_GPDDA4pp") - R("R_GPDDA5pp") - R("R_GTHRDHpp") -
 R("R_GUR1PPpp") + R("R_H2Otex") - R("R_H2Otex_Rev") - R("R_H2Otpp") +
 R("R_H2Otpp_Rev") - R("R_LACZpp") - R("R_LPLIPAL1A120pp") -
 R("R_LPLIPAL1A140pp") - R("R_LPLIPAL1A141pp") - R("R_LPLIPAL1A160pp") -
 R("R_LPLIPAL1A161pp") - R("R_LPLIPAL1A180pp") - R("R_LPLIPAL1A181pp") -
 R("R_LPLIPAL1E120pp") - R("R_LPLIPAL1E140pp") - R("R_LPLIPAL1E141pp") -
 R("R_LPLIPAL1E160pp") - R("R_LPLIPAL1E161pp") - R("R_LPLIPAL1E180pp") -
 R("R_LPLIPAL1E181pp") - R("R_LPLIPAL1G120pp") - R("R_LPLIPAL1G140pp") -
 R("R_LPLIPAL1G141pp") - R("R_LPLIPAL1G160pp") - R("R_LPLIPAL1G161pp")
 - R("R_LPLIPAL1G180pp") - R("R_LPLIPAL1G181pp") - R("R_MDDCP1pp") -
 R("R_MDDCP2pp") - R("R_MDDCP3pp") - R("R_MDDCP4pp") - R("R_MDDCP5pp") -
 R("R_MDDEP1pp") - R("R_MDDEP2pp") - R("R_MDDEP3pp") - R("R_MDDEP4pp")
 - R("R_MLDCP1App") - R("R_MLDCP1Bpp") - R("R_MLDCP2App") -
 R("R_MLDCP2Bpp") - R("R_MLDCP3App") - R("R_MLDEP1pp") - R("R_MLDEP2pp") +
 R("R_NO3R1bpp") + R("R_NO3R2bpp") - R("R_NTD10pp") - R("R_NTD11pp") -
 R("R_NTD12pp") - R("R_NTD1pp") - R("R_NTD2pp") - R("R_NTD3pp") -
 R("R_NTD4pp") - R("R_NTD5pp") - R("R_NTD6pp") - R("R_NTD7pp") - R("R_NTD8pp") -
 R("R_NTD9pp") - R("R_NTP3pp") + 2 R("R_NTRIR3pp") + 2 R("R_NTRIR4pp") -
 R("R_PAPA120pp") - R("R_PAPA140pp") - R("R_PAPA141pp") -
 R("R_PAPA160pp") - R("R_PAPA161pp") - R("R_PAPA180pp") - R("R_PAPA181pp") -
 R("R_PEAMNOpp") - R("R_PGPP120pp") - R("R_PGPP140pp") - R("R_PGPP141pp") -
 R("R_PGPP160pp") - R("R_PGPP161pp") - R("R_PGPP180pp") -
 R("R_PGPP181pp") - 6 R("R_PHYTSpp") - R("R_PLIPA1A120pp") - R("R_PLIPA1A140pp") -
 R("R_PLIPA1A141pp") - R("R_PLIPA1A160pp") - R("R_PLIPA1A161pp") -
 R("R_PLIPA1A180pp") - R("R_PLIPA1A181pp") - R("R_PLIPA1E120pp") -
 R("R_PLIPA1E140pp") - R("R_PLIPA1E141pp") - R("R_PLIPA1E160pp") -
 R("R_PLIPA1E161pp") - R("R_PLIPA1E180pp") - R("R_PLIPA1E181pp") - R("R_PLIPA1G120pp")
 - R("R_PLIPA1G140pp") - R("R_PLIPA1G141pp") - R("R_PLIPA1G160pp") -
 R("R_PLIPA1G161pp") - R("R_PLIPA1G180pp") - R("R_PLIPA1G181pp") -
 R("R_PLIPA2A120pp") - R("R_PLIPA2A140pp") - R("R_PLIPA2A141pp") -
 R("R_PLIPA2A160pp") - R("R_PLIPA2A161pp") - R("R_PLIPA2A180pp") -
 R("R_PLIPA2A181pp") - R("R_PLIPA2E120pp") - R("R_PLIPA2E140pp") -
 R("R_PLIPA2E141pp") - R("R_PLIPA2E160pp") - R("R_PLIPA2E161pp") - R("R_PLIPA2E180pp") -
 R("R_PLIPA2E181pp") - R("R_PLIPA2G120pp") - R("R_PLIPA2G140pp") -
 R("R_PLIPA2G141pp") - R("R_PLIPA2G160pp") - R("R_PLIPA2G161pp") -
 R("R_PLIPA2G180pp") - R("R_PLIPA2G181pp") - R("R_PPTHpp") - R("R_PSP_Lpp") -
 R("R_PTHRpp") - R("R_R5PPpp") + R("R_TMAOR1pp") + R("R_TMAOR2pp") -
 R("R_TREHpp") - R("R_TYROXDApp") - R("R_TYRPpp") - R("R_UACGALPpp") -
 R("R_UACGAMPpp") - R("R_UDCPDPpp") - R("R_UDPGALPpp") - R("R_UDPGPpp") -
 R("R_UGLCURPpp")  = 0.

Mh2p: R("R_H2tex") - R("R_H2tex_Rev") - R("R_H2tpp") +
 R("R_H2tpp_Rev") + R("R_PPTHpp")  = 0.

Mh2sb: R("Mh2sb_Transp") + R("R_EX_h2s_e_") - R("R_EX_h2s_e__Rev") 
 = 0.

Mh2sc: R("R_CYSDDS") + R("R_CYSDS") - R("R_CYSS") - R("R_H2SO") -
 R("R_H2St1pp") + R("R_SULRi")  = 0.

Mh2se: -R("R_EX_h2s_e_") + R("R_EX_h2s_e__Rev") - R("R_H2Stex") +
 R("R_H2Stex_Rev")  = 0.

Mh2sp: R("R_H2St1pp") + R("R_H2Stex") - R("R_H2Stex_Rev")  = 0.

Mhacolipab: R("Mhacolipab_Transp") + R("R_EX_hacolipa_e_") -
 R("R_EX_hacolipa_e__Rev")  = 0.

Mhacolipae: -R("R_EX_hacolipa_e_") + R("R_EX_hacolipa_e__Rev") +
 R("R_LIPAHT2ex")  = 0.

Mhalipab: R("Mhalipab_Transp") + R("R_EX_halipa_e_") -
 R("R_EX_halipa_e__Rev")  = 0.

Mhalipae: -R("R_EX_halipa_e_") + R("R_EX_halipa_e__Rev") +
 R("R_LIPAHTex")  = 0.

Mhb: R("Mhb_Transp") + R("R_EX_h_e_") - R("R_EX_h_e__Rev")  = 0.

Mhc: R("R_14GLUCANabcpp") + R("R_23DAPPAt2pp") - R("R_2DGLCNRx") -
 R("R_2DGLCNRy") - R("R_2DGULRx") - R("R_2DGULRy") + R("R_2MAHMP") -
 R("R_3HCINNMH") - R("R_3HPPPNH") + R("R_3KGK") - R("R_3OAR100") +
 R("R_3OAR100_Rev") - R("R_3OAR120") + R("R_3OAR120_Rev") - R("R_3OAR121") -
 R("R_3OAR140") + R("R_3OAR140_Rev") - R("R_3OAR141") - R("R_3OAR160") +
 R("R_3OAR160_Rev") - R("R_3OAR161") - R("R_3OAR180") +
 R("R_3OAR180_Rev") - R("R_3OAR181") - R("R_3OAR40") + R("R_3OAR40_Rev") -
 R("R_3OAR60") + R("R_3OAR60_Rev") - R("R_3OAR80") + R("R_3OAR80_Rev") -
 R("R_3OAS100") - R("R_3OAS120") - R("R_3OAS121") - R("R_3OAS140") -
 R("R_3OAS141") - R("R_3OAS160") - R("R_3OAS161") - R("R_3OAS180") -
 R("R_3OAS181") - R("R_3OAS60") - R("R_3OAS80") + R("R_3PEPTabcpp") +
 R("R_4PEPTabcpp") - R("R_5DGLCNR") + R("R_5DGLCNR_Rev") + R("R_5DGLCNt2rpp") -
 R("R_5DGLCNt2rpp_Rev") + R("R_AADDGT") + 2 R("R_ABUTD") + R("R_ABUTt2pp")
 + R("R_ACACt2pp") - R("R_ACACt2pp_Rev") + R("R_ACALD") -
 R("R_ACALD_Rev") - R("R_ACBIPGT") + R("R_ACCOAC") + R("R_ACGAMK") + R("R_ACGS") -
 R("R_ACHBS") - R("R_ACLS") + R("R_ACMAMUT") + R("R_ACNAMt2pp") +
 R("R_ACPS1") + R("R_ACt2rpp") - R("R_ACt2rpp_Rev") - R("R_ADA") +
 R("R_ADCL") - R("R_ADD") + R("R_ADEt2rpp") - R("R_ADEt2rpp_Rev") - R("R_ADMDC")
 + R("R_ADNK1") + R("R_ADNt2pp") + R("R_ADNt2rpp") -
 R("R_ADNt2rpp_Rev") + R("R_ADOCBIK") + R("R_ADOCBLabcpp") + R("R_ADOCBLS") +
 R("R_ADOCBLtonex") + 2 R("R_ADPRDP") + R("R_ADSK") + 2 R("R_ADSS") +
 R("R_AGM3Pt2pp") + R("R_AGM4Pt2pp") + R("R_AGMt2pp") + R("R_AGPR") -
 R("R_AGPR_Rev") + R("R_AGt3") + R("R_AIRC2") + R("R_AKGt2rpp") -
 R("R_AKGt2rpp_Rev") + R("R_ALAabcpp") + R("R_ALAALAabcpp") + R("R_ALAALAr") -
 R("R_ALAALAr_Rev") + R("R_ALAt2pp") - R("R_ALCD19") + R("R_ALCD19_Rev") +
 R("R_ALCD2x") - R("R_ALCD2x_Rev") + 2 R("R_ALDD19x") + 2 R("R_ALDD2x") + 2
 R("R_ALDD2y") + 2 R("R_ALDD3y") + 2 R("R_ALDD4") + R("R_ALLabcpp") +
 R("R_ALLK") - 2 R("R_ALLTAMH") + R("R_ALLTN") + R("R_ALLTNt2rpp") -
 R("R_ALLTNt2rpp_Rev") - R("R_ALR2") - R("R_ALR4x") + R("R_AMANK") +
 R("R_AMMQLT8") + 3 R("R_AMPMS2") + R("R_ANHMK") + R("R_ANS") -
 R("R_AOBUTDs") - R("R_AOXSr") + R("R_AOXSr_Rev") + 2 R("R_AP4AH") + 2 R("R_AP5AH")
 - R("R_APRAUR") + R("R_ARBabcpp") + R("R_ARBt2rpp") -
 R("R_ARBt2rpp_Rev") + R("R_ARBt3ipp") + R("R_ARBTNabcpp") + 2 R("R_ARBTNR1") + 2
 R("R_ARBTNR2") + 2 R("R_ARBTNR3") + R("R_ARBTNtex") + R("R_ARBTNtonex") +
 R("R_ARBTNtpp") + R("R_ARGabcpp") - R("R_ARGDC") + R("R_ARGSS") +
 R("R_ARGt3pp") + R("R_ASAD") - R("R_ASAD_Rev") + R("R_ASCBPL") +
 R("R_ASNabcpp") + R("R_ASNS1") + R("R_ASNS2") + R("R_ASNt2rpp") -
 R("R_ASNt2rpp_Rev") + R("R_ASO3t8pp") - R("R_ASP1DC") + R("R_ASPabcpp") +
 R("R_ASPCT") + R("R_ASPO3") + R("R_ASPO4") + R("R_ASPO5") + R("R_ASPO6") + 2
 R("R_ASPt2_2pp") + 3 R("R_ASPt2_3pp") + R("R_ASPt2pp") + R("R_AST") +
 R("R_ATHRDHr") - R("R_ATHRDHr_Rev") - R("R_ATPHs") + R("R_ATPM") + 3
 R("R_ATPS4rpp") - 3 R("R_ATPS4rpp_Rev") + R("R_BALAt2pp") + 2
 R("R_BETALDHx") + 2 R("R_BETALDHy") - R("R_BSORx") - R("R_BSORy") + R("R_BTS4") +
 R("R_BUTSO3abcpp") + R("R_BUTt2rpp") - R("R_BUTt2rpp_Rev") +
 R("R_CA2t3pp") + R("R_CADVtpp") - R("R_CBIAT") + R("R_CBIAT_Rev") +
 R("R_CBItonex") + R("R_CBIuabcpp") + R("R_CBL1abcpp") + R("R_CBL1tonex") -
 R("R_CBLAT") + R("R_CBLAT_Rev") + 2 R("R_CBMKr") - 2 R("R_CBMKr_Rev") + 2
 R("R_CBPS") + R("R_CD2abcpp") + R("R_CD2t3pp") + 2 R("R_CDAPPA120") + 2
 R("R_CDAPPA140") + 2 R("R_CDAPPA141") + 2 R("R_CDAPPA160") + 2
 R("R_CDAPPA161") + 2 R("R_CDAPPA180") + 2 R("R_CDAPPA181") + R("R_CDPMEK") +
 2 R("R_CFAS160E") + 2 R("R_CFAS160G") + 2 R("R_CFAS180E") + 2
 R("R_CFAS180G") + R("R_CGLYabcpp") + R("R_CHLabcpp") + R("R_CHLt2pp") -
 R("R_CINNDO") - R("R_CLt3_2pp") + R("R_COBALT2abcpp") + R("R_COBALT2t3pp") +
 R("R_COLIPAabcpp") + R("R_CPGNabcpp") + 2 R("R_CPGNR1") + 2
 R("R_CPGNR2") + 2 R("R_CPGNR3") + R("R_CPGNtonex") + R("R_CPGNUtex") +
 R("R_CPGNUtpp") - 2 R("R_CPPPGO") + R("R_CRNabcpp") + R("R_CRNDabcpp") +
 R("R_CRNDt2rpp") - R("R_CRNDt2rpp_Rev") + R("R_CRNt2rpp") -
 R("R_CRNt2rpp_Rev") + R("R_CS") - R("R_CSND") + R("R_CSNt2pp") + R("R_CTBTabcpp") +
 R("R_CTBTt2rpp") - R("R_CTBTt2rpp_Rev") + 2 R("R_CTPS2") +
 R("R_CU1abcpp") + R("R_CU2abcpp") + R("R_CUt3") + R("R_CYANST") - 3 R("R_CYNTAH")
 + R("R_CYNTt2pp") + R("R_CYSabc2pp") + R("R_CYSabcpp") +
 R("R_CYSDabcpp") + R("R_CYSS") - 2 R("R_CYSSADS") - 2 R("R_CYTBD2pp") - 2
 R("R_CYTBDpp") - 4 R("R_CYTBO3_4pp") - R("R_CYTD") + R("R_CYTDK2") +
 R("R_CYTDt2pp") + R("R_CYTDt2rpp") - R("R_CYTDt2rpp_Rev") - R("R_DADA") +
 R("R_DADNt2pp") + R("R_DAGK120") + R("R_DAGK140") + R("R_DAGK141") +
 R("R_DAGK160") + R("R_DAGK161") + R("R_DAGK180") + R("R_DAGK181") +
 R("R_DALAt2pp") + R("R_DAPabcpp") - R("R_DAPDC") - R("R_DASYN120") -
 R("R_DASYN140") - R("R_DASYN141") - R("R_DASYN160") - R("R_DASYN161") -
 R("R_DASYN180") - R("R_DASYN181") - R("R_DATPHs") + R("R_DB4PS") + 3
 R("R_DBTS") - R("R_DCTPD") - R("R_DCYTD") + R("R_DCYTt2pp") + R("R_DDGALK") +
 R("R_DDGLCNt2rpp") - R("R_DDGLCNt2rpp_Rev") + R("R_DDGLK") +
 R("R_DGSNt2pp") + R("R_DHBD") - R("R_DHBD_Rev") - R("R_DHBS") + R("R_DHCIND") +
 R("R_DHCINDO") - R("R_DHDPRy") + R("R_DHDPS") - R("R_DHFR") +
 R("R_DHFR_Rev") + R("R_DHFS") - R("R_DHNAOT4") + R("R_DHORTS") -
 R("R_DHORTS_Rev") + R("R_DHPPD") - R("R_DHPPDA2") + R("R_DINSt2pp") -
 R("R_DKGLCNR1") - R("R_DKGLCNR2x") - R("R_DKGLCNR2y") + 2 R("R_DKMPPD") + 4
 R("R_DKMPPD3") + R("R_D_LACt2pp") - R("R_D_LACt2pp_Rev") - R("R_DMPPS") +
 R("R_DMQMT") + R("R_DNTPPA") - R("R_DOGULNR") + R("R_DPCOAK") -
 R("R_DPR") + R("R_DSERDHr") - R("R_DSERDHr_Rev") + R("R_DSERt2pp") +
 R("R_DURIK1") + R("R_DURIt2pp") + R("R_DUTPDP") - R("R_DXPRIi") - R("R_DXPS") +
 R("R_DXYLK") + 2 R("R_E4PD") - 2 R("R_E4PD_Rev") - R("R_EAR100x") -
 R("R_EAR100y") - R("R_EAR120x") - R("R_EAR120y") - R("R_EAR121x") -
 R("R_EAR121y") - R("R_EAR140x") - R("R_EAR140y") - R("R_EAR141x") -
 R("R_EAR141y") - R("R_EAR160x") - R("R_EAR160y") - R("R_EAR161x") -
 R("R_EAR161y") - R("R_EAR180x") - R("R_EAR180y") - R("R_EAR181x") -
 R("R_EAR181y") - R("R_EAR40x") - R("R_EAR40y") - R("R_EAR60x") - R("R_EAR60y") -
 R("R_EAR80x") - R("R_EAR80y") + 59.81
 R("R_Ec_biomass_iAF1260_core_59p81M") + 9 R("R_ENTCS") + 3 R("R_ENTERES") + 3 R("R_ENTERES2") +
 R("R_ETHAt2pp") + R("R_ETHSO3abcpp") + R("R_ETOHt2rpp") -
 R("R_ETOHt2rpp_Rev") + R("R_FA100ACPHi") + R("R_FA120ACPHi") + R("R_FA140ACPHi") +
 R("R_FA141ACPHi") + R("R_FA160ACPHi") + R("R_FA161ACPHi") +
 R("R_FA80ACPHi") + R("R_FACOAE100") + R("R_FACOAE120") + R("R_FACOAE140") +
 R("R_FACOAE141") + R("R_FACOAE160") + R("R_FACOAE161") + R("R_FACOAE180") +
 R("R_FACOAE181") + R("R_FACOAE60") + R("R_FACOAE80") +
 R("R_FACOAL100t2pp") + R("R_FACOAL120t2pp") + R("R_FACOAL140t2pp") +
 R("R_FACOAL141t2pp") + R("R_FACOAL160t2pp") + R("R_FACOAL161t2pp") +
 R("R_FACOAL180t2pp") + R("R_FACOAL181t2pp") + R("R_FACOAL60t2pp") + R("R_FACOAL80t2pp") -
 R("R_FADRx") - R("R_FADRx2") + R("R_FALDH2") - R("R_FALDH2_Rev") +
 R("R_FCLK") + 2 R("R_FCLT") - 2 R("R_FDH4pp") - 2 R("R_FDH5pp") +
 R("R_FDMO") + R("R_FDMO2") + R("R_FDMO3") + R("R_FDMO4") + R("R_FDMO6") +
 R("R_FE2abcpp") + R("R_FE2t2pp") + R("R_FE2t3pp") + R("R_FE3abcpp") +
 R("R_FE3DCITabcpp") + R("R_FE3DCITtonex") + R("R_FE3DHBZSabcpp") +
 R("R_FE3DHBZStonex") + R("R_FE3HOXabcpp") + 2 R("R_FE3HOXR1") + 2
 R("R_FE3HOXR2") + 2 R("R_FE3HOXR3") + R("R_FE3HOXtonex") + R("R_FE3HOXUtex") +
 R("R_FE3HOXUtpp") + 2 R("R_FE3Ri") + R("R_FECRMabcpp") + 2
 R("R_FECRMR1") + 2 R("R_FECRMR2") + 2 R("R_FECRMR3") + R("R_FECRMtonex") +
 R("R_FECRMUtex") + R("R_FECRMUtpp") + R("R_FEENTERabcpp") + 2
 R("R_FEENTERR1") + 2 R("R_FEENTERR2") + 2 R("R_FEENTERR3") + R("R_FEENTERtex") +
 R("R_FEENTERtonex") + R("R_FEENTERtpp") + R("R_FEOXAMabcpp") + 2
 R("R_FEOXAMR1") + 2 R("R_FEOXAMR2") + 2 R("R_FEOXAMR3") + R("R_FEOXAMtonex") +
 R("R_FEOXAMUtex") + R("R_FEOXAMUtpp") - R("R_FHL") - R("R_FLDR") -
 R("R_FLVR") - R("R_FLVRx") + R("R_FMETTRS") - R("R_FMNAT") -
 R("R_FMNRx") - R("R_FMNRx2") + R("R_FORt2pp") + R("R_FRUK") + R("R_FRULYSK") +
 R("R_FRULYSt2pp") + R("R_FRUURt2rpp") - R("R_FRUURt2rpp_Rev") +
 R("R_FTHFD") + R("R_FUCtpp") - R("R_FUCtpp_Rev") + 2 R("R_FUMt2_2pp") + 3
 R("R_FUMt2_3pp") + R("R_G1PACT") - R("R_G1PTT") + R("R_G3PCabcpp") +
 R("R_G3PD2") - R("R_G3PD2_Rev") + R("R_G3PEabcpp") + R("R_G3PGabcpp") +
 R("R_G3PIabcpp") + R("R_G3PSabcpp") + R("R_G5SADs") - R("R_G5SD") +
 R("R_G6PDH2r") - R("R_G6PDH2r_Rev") + R("R_GALabcpp") + R("R_GALCTLO") +
 R("R_GALCTNLt2pp") + R("R_GALCTNt2pp") + R("R_GALCTt2rpp") -
 R("R_GALCTt2rpp_Rev") + R("R_GALKr") - R("R_GALKr_Rev") + R("R_GALT1") +
 R("R_GALt2pp") - R("R_GALUi") + R("R_GALURt2rpp") - R("R_GALURt2rpp_Rev") +
 R("R_GAPD") - R("R_GAPD_Rev") + R("R_GARFT") - R("R_GARFT_Rev") +
 R("R_GART") + 2 R("R_GCALDD") + R("R_GDPDPK") + R("R_GDPMNH") + 2
 R("R_GDPMNP") + 2 R("R_GGGABADr") - 2 R("R_GGGABADr_Rev") + R("R_GGPTRCS") +
 R("R_GLCabcpp") + R("R_GLCNt2rpp") - R("R_GLCNt2rpp_Rev") +
 R("R_GLCRt2rpp") - R("R_GLCRt2rpp_Rev") + R("R_GLCS1") + R("R_GLCt2pp") +
 R("R_GLCTR1") + R("R_GLCTR2") + R("R_GLCTR3") + R("R_GLCURt2rpp") -
 R("R_GLCURt2rpp_Rev") - R("R_GLGC") + R("R_GLNabcpp") + R("R_GLNS") +
 R("R_GLTPD") - R("R_GLTPD_Rev") + R("R_GLUabcpp") + R("R_GLUCYS") - R("R_GLUDC")
 + R("R_GLUDy") - R("R_GLUDy_Rev") - R("R_GLUSy") + R("R_GLUt2rpp") -
 R("R_GLUt2rpp_Rev") - R("R_GLUTRR") - R("R_GLXCL") + R("R_GLYBabcpp") +
 R("R_GLYBt2pp") + R("R_GLYC3Pabcpp") + R("R_GLYCAt2rpp") -
 R("R_GLYCAt2rpp_Rev") + R("R_GLYCDx") + R("R_GLYCK") + R("R_GLYCK2") -
 R("R_GLYCLTDx") - R("R_GLYCLTDy") + R("R_GLYCLTt2rpp") - R("R_GLYCLTt2rpp_Rev")
 + R("R_GLYK") + R("R_GLYOX") + R("R_GLYOX3") + R("R_GLYt2pp") -
 R("R_GMHEPAT") + R("R_GMHEPK") - 2 R("R_GMPR") + 2 R("R_GMPS2") +
 R("R_GNK") - R("R_GOFUCR") + 2 R("R_GP4GH") + R("R_GPDDA1") + R("R_GPDDA2") +
 R("R_GPDDA3") + R("R_GPDDA4") + R("R_GPDDA5") + R("R_GSNK") +
 R("R_GSNt2pp") + R("R_GSPMDS") - R("R_GTHOr") + R("R_GTHOr_Rev") +
 R("R_GTHRDabc2pp") + R("R_GTHRDabcpp") + R("R_GTHS") + R("R_GTPCI") + 2
 R("R_GTPCII2") + R("R_GTPDPDP") + R("R_GTPDPK") - R("R_GTPHs") - R("R_GUAD") +
 R("R_GUAt2pp") + 2 R("R_H2SO") + R("R_HACD1i") + R("R_HACD2i") +
 R("R_HACD3i") + R("R_HACD4i") + R("R_HACD5i") + R("R_HACD6i") +
 R("R_HACD7i") + R("R_HACD8i") + R("R_HCINNMt2rpp") - R("R_HCINNMt2rpp_Rev") +
 R("R_HCO3E") - R("R_HCO3E_Rev") + R("R_HCYSMT") + R("R_HCYSMT2") +
 R("R_HEPK1") + R("R_HEPK2") + R("R_HEPT1") + R("R_HEPT2") + R("R_HEPT3") +
 R("R_HEPT4") + R("R_HETZK") + R("R_HEX1") + R("R_HEX4") + R("R_HEX7") +
 R("R_HEXt2rpp") - R("R_HEXt2rpp_Rev") + R("R_HG2abcpp") +
 R("R_HG2t3pp") + R("R_HISabcpp") + R("R_HISt2rpp") - R("R_HISt2rpp_Rev") + 3
 R("R_HISTD") + R("R_HKNDDH") + R("R_HKNTDH") + R("R_HMPK1") +
 R("R_HOMt2pp") + R("R_HPPK2") + R("R_HPPPNDO") + R("R_HPPPNt2rpp") -
 R("R_HPPPNt2rpp_Rev") - R("R_HPYRRx") - R("R_HPYRRy") + R("R_HSDy") -
 R("R_HSDy_Rev") + R("R_HSK") + R("R_HXAND") - 2 R("R_HYD1pp") - 2 R("R_HYD2pp") - 2
 R("R_HYD3pp") - R("R_IDOND") - R("R_IDOND2") + R("R_IDOND_Rev") +
 R("R_IDONt2rpp") - R("R_IDONt2rpp_Rev") + R("R_IG3PS") - R("R_IGPS") +
 R("R_ILEabcpp") + R("R_ILEt2rpp") - R("R_ILEt2rpp_Rev") + R("R_IMPD") -
 R("R_INDOLEt2pp") + R("R_INDOLEt2rpp") - R("R_INDOLEt2rpp_Rev") +
 R("R_INSK") + R("R_INSt2pp") + R("R_INSt2rpp") - R("R_INSt2rpp_Rev") -
 R("R_IPDPS") + R("R_IPMD") + R("R_IPPS") + R("R_ISETACabcpp") +
 R("R_K2L4Aabcpp") + R("R_Kabcpp") + R("R_KARA1") - R("R_KARA1_Rev") -
 R("R_KARA2") + R("R_KARA2_Rev") - R("R_KAS14") - R("R_KAS15") - R("R_KG6PDC") +
 R("R_Kt2pp") + R("R_Kt3pp") - R("R_LALDO2x") + 2 R("R_LCADi") -
 R("R_LCARR") + R("R_LCARR_Rev") - R("R_LCARS") + R("R_LCARS_Rev") +
 R("R_LCTSt3ipp") + R("R_LCTStpp") - R("R_LCTStpp_Rev") + R("R_LDH_D") -
 R("R_LDH_D_Rev") + R("R_LEUabcpp") + R("R_LEUt2rpp") - R("R_LEUt2rpp_Rev") +
 R("R_LIPAabcpp") + R("R_LIPACabcpp") + R("R_L_LACt2rpp") -
 R("R_L_LACt2rpp_Rev") + R("R_LPADSS") + 2 R("R_LPLIPAL2A120") + 2
 R("R_LPLIPAL2A140") + 2 R("R_LPLIPAL2A141") + 2 R("R_LPLIPAL2A160") + 2
 R("R_LPLIPAL2A161") + 2 R("R_LPLIPAL2A180") + 2 R("R_LPLIPAL2A181") +
 R("R_LPLIPAL2E120") + R("R_LPLIPAL2E140") + R("R_LPLIPAL2E141") +
 R("R_LPLIPAL2E160") + R("R_LPLIPAL2E161") + R("R_LPLIPAL2E180") + R("R_LPLIPAL2E181")
 + R("R_LPLIPAL2G120") + R("R_LPLIPAL2G140") + R("R_LPLIPAL2G141") +
 R("R_LPLIPAL2G160") + R("R_LPLIPAL2G161") + R("R_LPLIPAL2G180") +
 R("R_LPLIPAL2G181") + R("R_LSERDHr") - R("R_LSERDHr_Rev") + R("R_LYSabcpp")
 - R("R_LYSDC") + R("R_LYSt2pp") + R("R_LYSt3pp") + R("R_LYXt2pp") +
 R("R_M1PD") - R("R_M1PD_Rev") - R("R_MACPD") + 2 R("R_MALDt2_2pp") +
 R("R_MALS") + 2 R("R_MALt2_2pp") + 3 R("R_MALt2_3pp") + R("R_MALTabcpp")
 + R("R_MALTHXabcpp") + R("R_MALTPTabcpp") + R("R_MALTTRabcpp") +
 R("R_MALTTTRabcpp") - R("R_MAN1PT2") + R("R_MANAO") - R("R_MANAO_Rev") +
 R("R_MCITS") + R("R_MCPST") + R("R_MDH") - R("R_MDH_Rev") +
 R("R_MELIBt2pp") + R("R_MELIBt3ipp") - R("R_MEPCT") + R("R_METabcpp") +
 R("R_METDabcpp") + R("R_METS") + R("R_METSOX1abcpp") + R("R_METSOX2abcpp") - 2
 R("R_MG2t3_2pp") + 2 R("R_MG2t3_2pp_Rev") + R("R_MG2uabcpp") -
 R("R_MMCD") + R("R_MMETt2pp") + R("R_MN2t3pp") + R("R_MNt2pp") + R("R_MOAT") +
 R("R_MOAT2") + R("R_MOAT3C") + R("R_MOBDabcpp") + 2 R("R_MPTG") +
 R("R_MPTG2") + R("R_MSO3abcpp") + R("R_MTHFC") - R("R_MTHFC_Rev") - 2
 R("R_MTHFR2") + R("R_MTRK") + 2 R("R_NADDP") - R("R_NADH10") - 4
 R("R_NADH16pp") - 4 R("R_NADH17pp") - 4 R("R_NADH18pp") - R("R_NADH5") -
 R("R_NADH9") + R("R_NADK") + R("R_NADN") - R("R_NADPHQR2") -
 R("R_NADPHQR3") - R("R_NADPHQR4") + R("R_NADS1") + 3 R("R_NAt3_1p5pp") + 2
 R("R_NAt3_2pp") + R("R_NAt3pp") - R("R_NHFRBO") + R("R_NI2abcpp") +
 R("R_NI2t3pp") + R("R_NI2uabcpp") - R("R_NMNAT") + R("R_NMNN") + R("R_NMNt7pp") -
 R("R_NNATr") + R("R_NNATr_Rev") + R("R_NNDMBRT") - 2 R("R_NNDPR") +
 R("R_NO2t2rpp") - R("R_NO2t2rpp_Rev") - 2 R("R_NO3R1pp") - 2
 R("R_NO3R2pp") + R("R_NODOx") + R("R_NODOy") + R("R_NTP1") + R("R_NTP10") +
 R("R_NTP11") + R("R_NTP12") + R("R_NTP3") + R("R_NTP5") + R("R_NTPP1") +
 R("R_NTPP10") + R("R_NTPP11") + R("R_NTPP2") + R("R_NTPP3") +
 R("R_NTPP4") + R("R_NTPP5") + R("R_NTPP6") + R("R_NTPP7") + R("R_NTPP8") +
 R("R_NTPP9") - 5 R("R_NTRIR2x") + R("R_O16GALFT") + R("R_O16GLCT1") +
 R("R_O16GLCT2") - R("R_OAADC") + R("R_OCBT") - R("R_OCBT_Rev") +
 R("R_OHPHM") + R("R_OMBZLM") - R("R_OMCDC") + 3 R("R_OMMBLHX3") - R("R_OMPDC")
 + 3 R("R_OMPHHX3") - R("R_OPHBDC") + 3 R("R_OPHHX3") + R("R_ORNabcpp")
 - R("R_ORNDC") + 2 R("R_OROTt2_2pp") - R("R_OXGDC2") + R("R_P5CD") -
 2 R("R_P5CR") + R("R_PA120abcpp") + R("R_PA140abcpp") +
 R("R_PA141abcpp") + R("R_PA160abcpp") + R("R_PA161abcpp") + R("R_PA180abcpp") +
 R("R_PA181abcpp") + R("R_PACALDt2rpp") - R("R_PACALDt2rpp_Rev") +
 R("R_PANTS") + 2 R("R_PAPSR") + 2 R("R_PAPSR2") + R("R_PDE1") + R("R_PDX5PO2")
 + R("R_PDX5PS") + R("R_PE120abcpp") + R("R_PE140abcpp") +
 R("R_PE141abcpp") + R("R_PE160abcpp") + R("R_PE161abcpp") + R("R_PE180abcpp") +
 R("R_PE181abcpp") + R("R_PERD") - R("R_PERD_Rev") + R("R_PFK") +
 R("R_PFK_2") + R("R_PG120abcpp") + R("R_PG140abcpp") + R("R_PG141abcpp") +
 R("R_PG160abcpp") + R("R_PG161abcpp") + R("R_PG180abcpp") +
 R("R_PG181abcpp") + R("R_PGCD") + R("R_PGL") + R("R_PGP120abcpp") +
 R("R_PGP140abcpp") + R("R_PGP141abcpp") + R("R_PGP160abcpp") + R("R_PGP161abcpp") +
 R("R_PGP180abcpp") + R("R_PGP181abcpp") + R("R_PGSA120") +
 R("R_PGSA140") + R("R_PGSA141") + R("R_PGSA160") + R("R_PGSA161") +
 R("R_PGSA180") + R("R_PGSA181") + R("R_PHEMEabcpp") + R("R_PHEt2rpp") -
 R("R_PHEt2rpp_Rev") + R("R_PIt2rpp") - R("R_PIt2rpp_Rev") + R("R_PIuabcpp") +
 R("R_PNTK") + R("R_PPA") + R("R_PPA2") + R("R_PPBNGS") + R("R_PPC") -
 R("R_PPCDC") + R("R_PPNCL2") - R("R_PPNDH") - R("R_PPPNDO") +
 R("R_PPPNt2rpp") - R("R_PPPNt2rpp_Rev") + 2 R("R_PPS") + R("R_PRAGSr") -
 R("R_PRAGSr_Rev") + 2 R("R_PRAIS") + R("R_PRASCSi") + R("R_PRATPP") +
 R("R_PRFGS") + R("R_PROabcpp") + R("R_PROD2") + R("R_PROGLYabcpp") +
 R("R_PROt2rpp") - R("R_PROt2rpp_Rev") + R("R_PRPPS") - R("R_PRPPS_Rev") +
 R("R_PSCLYSt2pp") - R("R_PSD120") - R("R_PSD140") - R("R_PSD141") -
 R("R_PSD160") - R("R_PSD161") - R("R_PSD180") - R("R_PSD181") +
 R("R_PSSA120") + R("R_PSSA140") + R("R_PSSA141") + R("R_PSSA160") + R("R_PSSA161")
 + R("R_PSSA180") + R("R_PSSA181") - R("R_PTPATi") + R("R_PTRCabcpp")
 + R("R_PTRCt2pp") + R("R_PYDAMK") + R("R_PYDXK") + R("R_PYDXNK") -
 R("R_PYK") + R("R_PYRt2rpp") - R("R_PYRt2rpp_Rev") + 2 R("R_QMO2") + 2
 R("R_QMO3") + R("R_R1PK") + R("R_RBFK") + R("R_RBK") + R("R_RBK_L1") +
 R("R_RHAT1") + R("R_RIBabcpp") + R("R_RMK") + R("R_RMNtpp") - 2
 R("R_SADH") + R("R_SBTPD") - R("R_SBTPD_Rev") + R("R_SELCYSS") -
 R("R_SERASr") + R("R_SERASr_Rev") + R("R_SERt2rpp") - R("R_SERt2rpp_Rev") +
 R("R_SFGTHi") + 2 R("R_SGSAD") + R("R_SHCHD2") + 3 R("R_SHCHF") -
 R("R_SHK3Dr") + R("R_SHK3Dr_Rev") + R("R_SHKK") + R("R_SHSL1") + R("R_SKMt2pp")
 + R("R_SPMDabcpp") + R("R_SPMDAT1") + R("R_SPMDAT2") + R("R_SPMS") - 2
 R("R_SPODM") + 2 R("R_SSALx") + 2 R("R_SSALy") + 2 R("R_SUCCt2_2pp")
 + 3 R("R_SUCCt2_3pp") + R("R_SUCCt3pp") + R("R_SULabcpp") +
 R("R_SULFACabcpp") - 5 R("R_SULRi") + R("R_TAGURr") - R("R_TAGURr_Rev") +
 R("R_TAUDO") + R("R_TAURabcpp") + R("R_TDP") + R("R_TDPADGAT") -
 R("R_TDPDRR") + R("R_TDSK") + 2 R("R_THD2pp") + R("R_THMabcpp") + R("R_THMDt2pp")
 + R("R_THMDt2rpp") - R("R_THMDt2rpp_Rev") + R("R_THRabcpp") +
 R("R_THRD") + R("R_THRt2pp") + R("R_THRt2rpp") - R("R_THRt2rpp_Rev") +
 R("R_THYMt3pp") + R("R_THZPSN") - R("R_TMAOR1") - R("R_TMAOR2") +
 R("R_TMDK1") + R("R_TMK") - R("R_TMPPP") - R("R_TRDR") + R("R_TRE6PS") +
 R("R_TRPt2rpp") - R("R_TRPt2rpp_Rev") - R("R_TRSARr") + R("R_TRSARr_Rev") +
 R("R_TSULabcpp") + R("R_TUNGSabcpp") + R("R_TYRt2rpp") -
 R("R_TYRt2rpp_Rev") + R("R_U23GAAT") + R("R_UAAGDS") + 3 R("R_UACMAMO") -
 R("R_UAGDP") + R("R_UAGPT3") + R("R_UAMAGS") + R("R_UAMAS") - R("R_UAPGR") +
 R("R_UDCPDP") + 3 R("R_UDPGD") - 2 R("R_UGLYCH") + R("R_UGMDDS") +
 R("R_ULA4NFT") + R("R_UM3PL") + R("R_UM4PL") + R("R_UPP3MT") - 4
 R("R_UPPDC1") + R("R_URAt2pp") + R("R_URAt2rpp") - R("R_URAt2rpp_Rev") +
 R("R_URDGLYCD") + R("R_URIK2") + R("R_URIt2pp") + R("R_URIt2rpp") -
 R("R_URIt2rpp_Rev") + 2 R("R_USHD") + R("R_VALabcpp") + R("R_VALt2rpp") -
 R("R_VALt2rpp_Rev") + R("R_XAND") + R("R_XANt2pp") + R("R_XTSNt2rpp") -
 R("R_XTSNt2rpp_Rev") + R("R_XYLabcpp") + R("R_XYLK") + R("R_XYLK2") +
 R("R_XYLt2pp") + R("R_XYLUt2pp") + R("R_ZN2abcpp") + R("R_ZN2t3pp") +
 R("R_ZNabcpp")  = 0.

Mhco3c: -R("R_ACCOAC") - R("R_AIRC2") - R("R_CBPS") - R("R_CYNTAH")
 + R("R_HCO3E") - R("R_HCO3E_Rev")  = 0.

MhcysLc: R("R_CYSTL") - R("R_HCYSMT") - R("R_HCYSMT2") - R("R_METS")
 + R("R_RHCCE")  = 0.

Mhdcab: R("Mhdcab_Transp") + R("R_EX_hdca_e_") -
 R("R_EX_hdca_e__Rev")  = 0.

Mhdcac: -R("R_2AGPEAT160") - R("R_2AGPGAT160") - R("R_AACPS3") +
 R("R_FA160ACPHi") + R("R_FACOAE160") + R("R_LPLIPAL2A160") +
 R("R_LPLIPAL2E160") + R("R_LPLIPAL2G160")  = 0.

Mhdcae: -R("R_EX_hdca_e_") + R("R_EX_hdca_e__Rev") - R("R_HDCAtexi")
 - R("R_LIPAHT2ex") - R("R_LIPAHTex")  = 0.

Mhdcap: -R("R_FACOAL160t2pp") + R("R_HDCAtexi") +
 R("R_LPLIPAL1A160pp") + R("R_LPLIPAL1E160pp") + R("R_LPLIPAL1G160pp") +
 R("R_PLIPA1A160pp") + R("R_PLIPA1E160pp") + R("R_PLIPA1G160pp") + R("R_PLIPA2A160pp")
 + R("R_PLIPA2E160pp") + R("R_PLIPA2G160pp")  = 0.

Mhdceab: R("Mhdceab_Transp") + R("R_EX_hdcea_e_") -
 R("R_EX_hdcea_e__Rev")  = 0.

Mhdceac: -R("R_2AGPEAT161") - R("R_2AGPGAT161") - R("R_AACPS4") +
 R("R_FA161ACPHi") + R("R_FACOAE161") + R("R_LPLIPAL2A161") +
 R("R_LPLIPAL2E161") + R("R_LPLIPAL2G161")  = 0.

Mhdceae: -R("R_EX_hdcea_e_") + R("R_EX_hdcea_e__Rev") -
 R("R_HDCEAtexi")  = 0.

Mhdceap: -R("R_FACOAL161t2pp") + R("R_HDCEAtexi") +
 R("R_LPLIPAL1A161pp") + R("R_LPLIPAL1E161pp") + R("R_LPLIPAL1G161pp") +
 R("R_PLIPA1A161pp") + R("R_PLIPA1E161pp") + R("R_PLIPA1G161pp") +
 R("R_PLIPA2A161pp") + R("R_PLIPA2E161pp") + R("R_PLIPA2G161pp")  = 0.

Mhdcoac: -R("R_FACOAE161") + R("R_FACOAL161t2pp") - R("R_HDCOAI")  =
 0.

Mhdd2coac: R("R_ACOAD7f") - R("R_ACOAD7f_Rev") + R("R_ECOAH7") -
 R("R_ECOAH7_Rev") + R("R_HDCOAI")  = 0.

MhdeACPc: -R("R_3OAS181") + R("R_AACPS4") - R("R_AGPAT161") +
 R("R_EAR161x") + R("R_EAR161y") - R("R_EDTXS3") - R("R_FA161ACPHi") -
 R("R_G3PAT161")  = 0.

Mhe: -R("R_AGt3") - R("R_CUt3") - R("R_EX_h_e_") +
 R("R_EX_h_e__Rev") - R("R_Htex") + R("R_Htex_Rev") - R("R_LIPAHT2ex") - R("R_LIPAHTex")
  = 0.

MhemeOc: R("R_HEMEOS")  = 0.

MhexACPc: -R("R_3OAS80") + R("R_EAR60x") + R("R_EAR60y")  = 0.

Mhg2b: R("Mhg2b_Transp") + R("R_EX_hg2_e_") - R("R_EX_hg2_e__Rev") 
 = 0.

Mhg2c: -R("R_HG2abcpp") - R("R_HG2t3pp")  = 0.

Mhg2e: -R("R_EX_hg2_e_") + R("R_EX_hg2_e__Rev") - R("R_HG2tex") +
 R("R_HG2tex_Rev")  = 0.

Mhg2p: R("R_HG2abcpp") + R("R_HG2t3pp") + R("R_HG2tex") -
 R("R_HG2tex_Rev")  = 0.

Mhhlipac: -R("R_HEPK1") + R("R_HEPT2")  = 0.

MhisLb: R("MhisLb_Transp") + R("R_EX_his_L_e_") -
 R("R_EX_his_L_e__Rev")  = 0.

MhisLc: -0.09474 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_HISabcpp") + R("R_HISt2rpp") - R("R_HISt2rpp_Rev") + R("R_HISTD") -
 R("R_HISTRS")  = 0.

MhisLe: -R("R_EX_his_L_e_") + R("R_EX_his_L_e__Rev") - R("R_HIStex")
 + R("R_HIStex_Rev")  = 0.

MhisLp: -R("R_HISabcpp") - R("R_HISt2rpp") + R("R_HISt2rpp_Rev") +
 R("R_HIStex") - R("R_HIStex_Rev")  = 0.

Mhispc: -R("R_HISTP") + R("R_HSTPT")  = 0.

Mhistdc: -R("R_HISTD") + R("R_HISTP")  = 0.

Mhistrnac: R("R_HISTRS")  = 0.

Mhknddc: -R("R_HKNDDH") + R("R_HPPPNDO")  = 0.

Mhkntdc: R("R_DHCINDO") - R("R_HKNTDH")  = 0.

Mhlipac: R("R_HEPT1") - R("R_HEPT2")  = 0.

Mhmbilc: R("R_HMBS") - R("R_UPP3S")  = 0.

Mhmfurnb: R("Mhmfurnb_Transp") + R("R_DM_HMFURN")  = 0.

Mhmfurnc: R("R_DHPTDCs") - R("R_DM_HMFURN")  = 0.

Mhmgthc: -R("R_FALDH2") + R("R_FALDH2_Rev") + R("R_FALGTHLs") -
 R("R_FALGTHLs_Rev")  = 0.

MhomLb: R("MhomLb_Transp") + R("R_EX_hom_L_e_") -
 R("R_EX_hom_L_e__Rev")  = 0.

MhomLc: -R("R_HOMt2pp") - R("R_HSDy") + R("R_HSDy_Rev") - R("R_HSK")
 - R("R_HSST")  = 0.

MhomLe: -R("R_EX_hom_L_e_") + R("R_EX_hom_L_e__Rev") - R("R_HOMtex")
 + R("R_HOMtex_Rev")  = 0.

MhomLp: R("R_HOMt2pp") + R("R_HOMtex") - R("R_HOMtex_Rev")  = 0.

Mhp: -R("R_23DAPPAt2pp") + R("R_23PDE2pp") + R("R_23PDE4pp") +
 R("R_23PDE7pp") + R("R_23PDE9pp") - R("R_5DGLCNt2rpp") +
 R("R_5DGLCNt2rpp_Rev") - R("R_ABUTt2pp") - R("R_ACACt2pp") + R("R_ACACt2pp_Rev") -
 R("R_ACNAMt2pp") - R("R_ACt2rpp") + R("R_ACt2rpp_Rev") - R("R_ADEt2rpp") +
 R("R_ADEt2rpp_Rev") - R("R_ADNt2pp") - R("R_ADNt2rpp") +
 R("R_ADNt2rpp_Rev") - R("R_ADOCBLtonex") - R("R_AGM3Pt2pp") - R("R_AGM4Pt2pp") -
 R("R_AGMt2pp") - R("R_AKGt2rpp") + R("R_AKGt2rpp_Rev") - R("R_ALAt2pp")
 - R("R_ALLTNt2rpp") + R("R_ALLTNt2rpp_Rev") - R("R_ARBt2rpp") +
 R("R_ARBt2rpp_Rev") - R("R_ARBt3ipp") - R("R_ARBTNtex") - R("R_ARBTNtonex")
 - R("R_ARBTNtpp") - R("R_ARGDCpp") - R("R_ARGt3pp") - R("R_ASNt2rpp")
 + R("R_ASNt2rpp_Rev") - 2 R("R_ASPt2_2pp") - 3 R("R_ASPt2_3pp") -
 R("R_ASPt2pp") - 4 R("R_ATPS4rpp") + 4 R("R_ATPS4rpp_Rev") -
 R("R_BALAt2pp") - R("R_BUTt2rpp") + R("R_BUTt2rpp_Rev") - R("R_CA2t3pp") -
 R("R_CADVtpp") - R("R_CBItonex") - R("R_CBL1tonex") - R("R_CD2t3pp") -
 R("R_CHLt2pp") + R("R_CLPNH120pp") + R("R_CLPNH140pp") + R("R_CLPNH141pp") +
 R("R_CLPNH160pp") + R("R_CLPNH161pp") + R("R_CLPNH180pp") +
 R("R_CLPNH181pp") + R("R_CLt3_2pp") - R("R_COBALT2t3pp") - R("R_CPGNtonex") -
 R("R_CPGNUtex") - R("R_CPGNUtpp") - R("R_CRNDt2rpp") +
 R("R_CRNDt2rpp_Rev") - R("R_CRNt2rpp") + R("R_CRNt2rpp_Rev") - R("R_CSNt2pp") -
 R("R_CTBTt2rpp") + R("R_CTBTt2rpp_Rev") - 4 R("R_CU1Opp") + R("R_CYANSTpp") -
 R("R_CYNTt2pp") + 2 R("R_CYTBD2pp") + 2 R("R_CYTBDpp") + 4
 R("R_CYTBO3_4pp") - R("R_CYTDt2pp") - R("R_CYTDt2rpp") + R("R_CYTDt2rpp_Rev") -
 R("R_DADNt2pp") - R("R_DALAt2pp") - R("R_DCYTt2pp") -
 R("R_DDGLCNt2rpp") + R("R_DDGLCNt2rpp_Rev") - R("R_DGSNt2pp") - R("R_DINSt2pp") -
 R("R_D_LACt2pp") + R("R_D_LACt2pp_Rev") - R("R_DSERt2pp") -
 R("R_DURIt2pp") + R("R_ECA4OALpp") + R("R_ECAP1pp") + R("R_ECAP2pp") +
 R("R_ECAP3pp") - R("R_ETHAt2pp") - R("R_ETOHt2rpp") + R("R_ETOHt2rpp_Rev") -
 R("R_FACOAL100t2pp") - R("R_FACOAL120t2pp") - R("R_FACOAL140t2pp") -
 R("R_FACOAL141t2pp") - R("R_FACOAL160t2pp") - R("R_FACOAL161t2pp") -
 R("R_FACOAL180t2pp") - R("R_FACOAL181t2pp") - R("R_FACOAL60t2pp") -
 R("R_FACOAL80t2pp") + R("R_FDH4pp") + R("R_FDH5pp") - R("R_FE2t2pp") -
 R("R_FE2t3pp") - R("R_FE3DCITtonex") - R("R_FE3DHBZStonex") -
 R("R_FE3HOXtonex") - R("R_FE3HOXUtex") - R("R_FE3HOXUtpp") - R("R_FECRMtonex") -
 R("R_FECRMUtex") - R("R_FECRMUtpp") - R("R_FEENTERtex") -
 R("R_FEENTERtonex") - R("R_FEENTERtpp") - R("R_FEOXAMtonex") - R("R_FEOXAMUtex") -
 R("R_FEOXAMUtpp") - 4 R("R_FEROpp") - R("R_FORt2pp") - R("R_FRULYSt2pp") -
 R("R_FRUURt2rpp") + R("R_FRUURt2rpp_Rev") - R("R_FUCtpp") +
 R("R_FUCtpp_Rev") - 2 R("R_FUMt2_2pp") - 3 R("R_FUMt2_3pp") - R("R_GALCTNLt2pp")
 - R("R_GALCTNt2pp") - R("R_GALCTt2rpp") + R("R_GALCTt2rpp_Rev") -
 R("R_GALt2pp") - R("R_GALURt2rpp") + R("R_GALURt2rpp_Rev") +
 R("R_GLCDpp") - R("R_GLCNt2rpp") + R("R_GLCNt2rpp_Rev") - R("R_GLCRt2rpp") +
 R("R_GLCRt2rpp_Rev") - R("R_GLCt2pp") - R("R_GLCURt2rpp") +
 R("R_GLCURt2rpp_Rev") - R("R_GLUt2rpp") + R("R_GLUt2rpp_Rev") - R("R_GLYBt2pp") -
 R("R_GLYCAt2rpp") + R("R_GLYCAt2rpp_Rev") - R("R_GLYCLTt2rpp") +
 R("R_GLYCLTt2rpp_Rev") - R("R_GLYt2pp") + R("R_GPDDA1pp") + R("R_GPDDA2pp") +
 R("R_GPDDA3pp") + R("R_GPDDA4pp") + R("R_GPDDA5pp") - R("R_GSNt2pp") -
 R("R_GUAt2pp") - R("R_HCINNMt2rpp") + R("R_HCINNMt2rpp_Rev") -
 R("R_HEXt2rpp") + R("R_HEXt2rpp_Rev") - R("R_HG2t3pp") - R("R_HISt2rpp") +
 R("R_HISt2rpp_Rev") - R("R_HOMt2pp") - R("R_HPPPNt2rpp") +
 R("R_HPPPNt2rpp_Rev") + R("R_Htex") - R("R_Htex_Rev") + 2 R("R_HYD1pp") + 2
 R("R_HYD2pp") + 2 R("R_HYD3pp") - R("R_IDONt2rpp") + R("R_IDONt2rpp_Rev") -
 R("R_ILEt2rpp") + R("R_ILEt2rpp_Rev") + R("R_INDOLEt2pp") -
 R("R_INDOLEt2rpp") + R("R_INDOLEt2rpp_Rev") - R("R_INSt2pp") - R("R_INSt2rpp") +
 R("R_INSt2rpp_Rev") - R("R_Kt2pp") - R("R_Kt3pp") - R("R_LCTSt3ipp")
 - R("R_LCTStpp") + R("R_LCTStpp_Rev") - R("R_LEUt2rpp") +
 R("R_LEUt2rpp_Rev") - R("R_L_LACt2rpp") + R("R_L_LACt2rpp_Rev") +
 R("R_LPLIPAL1A120pp") + R("R_LPLIPAL1A140pp") + R("R_LPLIPAL1A141pp") +
 R("R_LPLIPAL1A160pp") + R("R_LPLIPAL1A161pp") + R("R_LPLIPAL1A180pp") +
 R("R_LPLIPAL1A181pp") + R("R_LPLIPAL1E120pp") + R("R_LPLIPAL1E140pp") +
 R("R_LPLIPAL1E141pp") + R("R_LPLIPAL1E160pp") + R("R_LPLIPAL1E161pp") +
 R("R_LPLIPAL1E180pp") + R("R_LPLIPAL1E181pp") + R("R_LPLIPAL1G120pp") +
 R("R_LPLIPAL1G140pp") + R("R_LPLIPAL1G141pp") + R("R_LPLIPAL1G160pp") +
 R("R_LPLIPAL1G161pp") + R("R_LPLIPAL1G180pp") + R("R_LPLIPAL1G181pp") -
 R("R_LYSt2pp") - R("R_LYSt3pp") - R("R_LYXt2pp") - 2 R("R_MALDt2_2pp") -
 2 R("R_MALt2_2pp") - 3 R("R_MALt2_3pp") - R("R_MELIBt2pp") -
 R("R_MELIBt3ipp") + 2 R("R_MG2t3_2pp") - 2 R("R_MG2t3_2pp_Rev") -
 R("R_MMETt2pp") - R("R_MN2t3pp") - R("R_MNt2pp") + 3 R("R_NADH16pp") + 3
 R("R_NADH17pp") + 3 R("R_NADH18pp") - 3 R("R_NAt3_1p5pp") - 2 R("R_NAt3_2pp") -
 R("R_NAt3pp") - R("R_NI2t3pp") - R("R_NO2t2rpp") +
 R("R_NO2t2rpp_Rev") + 2 R("R_NO3R1pp") + 2 R("R_NO3R2pp") + R("R_NTP3pp") - 2
 R("R_NTRIR3pp") - 2 R("R_NTRIR4pp") + R("R_O16A4Lpp") + R("R_O16AP1pp") +
 R("R_O16AP2pp") + R("R_O16AP3pp") - 2 R("R_OROTt2_2pp") - R("R_PACALDt2rpp")
 + R("R_PACALDt2rpp_Rev") - R("R_PHEt2rpp") + R("R_PHEt2rpp_Rev") -
 R("R_PIt2rpp") + R("R_PIt2rpp_Rev") + R("R_PLIPA1E120pp") +
 R("R_PLIPA1E140pp") + R("R_PLIPA1E141pp") + R("R_PLIPA1E160pp") +
 R("R_PLIPA1E161pp") + R("R_PLIPA1E180pp") + R("R_PLIPA1E181pp") + R("R_PLIPA1G120pp")
 + R("R_PLIPA1G140pp") + R("R_PLIPA1G141pp") + R("R_PLIPA1G160pp") +
 R("R_PLIPA1G161pp") + R("R_PLIPA1G180pp") + R("R_PLIPA1G181pp") +
 R("R_PLIPA2A120pp") + R("R_PLIPA2A140pp") + R("R_PLIPA2A141pp") +
 R("R_PLIPA2A160pp") + R("R_PLIPA2A161pp") + R("R_PLIPA2A180pp") +
 R("R_PLIPA2A181pp") + R("R_PLIPA2E120pp") + R("R_PLIPA2E140pp") +
 R("R_PLIPA2E141pp") + R("R_PLIPA2E160pp") + R("R_PLIPA2E161pp") + R("R_PLIPA2E180pp") +
 R("R_PLIPA2E181pp") + R("R_PLIPA2G120pp") + R("R_PLIPA2G140pp") +
 R("R_PLIPA2G141pp") + R("R_PLIPA2G160pp") + R("R_PLIPA2G161pp") +
 R("R_PLIPA2G180pp") + R("R_PLIPA2G181pp") - R("R_PPPNt2rpp") +
 R("R_PPPNt2rpp_Rev") - R("R_PROt2rpp") + R("R_PROt2rpp_Rev") - R("R_PSCLYSt2pp") -
 R("R_PTRCt2pp") - R("R_PYRt2rpp") + R("R_PYRt2rpp_Rev") - R("R_RMNtpp")
 - R("R_SERt2rpp") + R("R_SERt2rpp_Rev") - R("R_SKMt2pp") - 2
 R("R_SPODMpp") - 2 R("R_SUCCt2_2pp") - 3 R("R_SUCCt2_3pp") - R("R_SUCCt3pp") -
 2 R("R_THD2pp") - R("R_THMDt2pp") - R("R_THMDt2rpp") +
 R("R_THMDt2rpp_Rev") - R("R_THRt2pp") - R("R_THRt2rpp") + R("R_THRt2rpp_Rev") -
 R("R_THYMt3pp") - R("R_TMAOR1pp") - R("R_TMAOR2pp") - R("R_TRPt2rpp") +
 R("R_TRPt2rpp_Rev") - R("R_TYRt2rpp") + R("R_TYRt2rpp_Rev") + 2
 R("R_UACGALPpp") + 2 R("R_UACGAMPpp") + R("R_UDCPDPpp") + 2 R("R_UDPGALPpp") +
 2 R("R_UDPGPpp") + 2 R("R_UGLCURPpp") - R("R_URAt2pp") -
 R("R_URAt2rpp") + R("R_URAt2rpp_Rev") - R("R_URIt2pp") - R("R_URIt2rpp") +
 R("R_URIt2rpp_Rev") - R("R_VALt2rpp") + R("R_VALt2rpp_Rev") - R("R_XANt2pp") -
 R("R_XTSNt2rpp") + R("R_XTSNt2rpp_Rev") - R("R_XYLt2pp") -
 R("R_XYLUt2pp") - R("R_ZN2t3pp")  = 0.

Mhphhlipac: -R("R_HEPK2") + R("R_HEPT3")  = 0.

Mhpyrc: -R("R_HPYRI") + R("R_HPYRI_Rev") - R("R_HPYRRx") -
 R("R_HPYRRy")  = 0.

Mhqnc: R("R_AB6PGH")  = 0.

Mhx2coac: R("R_ACOAD2f") - R("R_ACOAD2f_Rev") + R("R_ECOAH2") -
 R("R_ECOAH2_Rev")  = 0.

Mhxab: R("Mhxab_Transp") + R("R_EX_hxa_e_") - R("R_EX_hxa_e__Rev") 
 = 0.

Mhxac: R("R_FACOAE60") + R("R_HEXt2rpp") - R("R_HEXt2rpp_Rev") -
 R("R_HXCT")  = 0.

Mhxae: -R("R_EX_hxa_e_") + R("R_EX_hxa_e__Rev") - R("R_HXAtex") +
 R("R_HXAtex_Rev")  = 0.

Mhxanb: R("Mhxanb_Transp") + R("R_EX_hxan_e_") -
 R("R_EX_hxan_e__Rev")  = 0.

Mhxanc: R("R_ADD") - R("R_HXAND") - R("R_HXPRT") + R("R_HYXNtpp") -
 R("R_HYXNtpp_Rev") + R("R_INSH") + R("R_PUNP5") - R("R_PUNP5_Rev") +
 R("R_PUNP6") - R("R_PUNP6_Rev")  = 0.

Mhxane: -R("R_EX_hxan_e_") + R("R_EX_hxan_e__Rev") - R("R_HYXNtex")
 + R("R_HYXNtex_Rev")  = 0.

Mhxanp: R("R_HYXNtex") - R("R_HYXNtex_Rev") - R("R_HYXNtpp") +
 R("R_HYXNtpp_Rev")  = 0.

Mhxap: -R("R_FACOAL60t2pp") - R("R_HEXt2rpp") + R("R_HEXt2rpp_Rev")
 + R("R_HXAtex") - R("R_HXAtex_Rev")  = 0.

Mhxcoac: -R("R_ACOAD2f") + R("R_ACOAD2f_Rev") - R("R_FACOAE60") +
 R("R_FACOAL60t2pp") + R("R_HXCT") + R("R_KAT3")  = 0.

Miaspc: R("R_ASPO3") + R("R_ASPO4") + R("R_ASPO5") + R("R_ASPO6") -
 R("R_QULNS")  = 0.

Michorc: R("R_ICHORS") + R("R_ICHORSi") - R("R_ICHORS_Rev") -
 R("R_ICHORT") - R("R_SHCHCS2")  = 0.

Micitc: R("R_ACONTb") - R("R_ACONTb_Rev") - R("R_ICDHyr") +
 R("R_ICDHyr_Rev") - R("R_ICL")  = 0.

Micolipac: -R("R_GLCTR1") + R("R_RHAT1")  = 0.

MidonLb: R("MidonLb_Transp") + R("R_EX_idon_L_e_") -
 R("R_EX_idon_L_e__Rev")  = 0.

MidonLc: R("R_2DGULRx") + R("R_2DGULRy") + R("R_IDOND") +
 R("R_IDOND2") - R("R_IDOND_Rev") + R("R_IDONt2rpp") - R("R_IDONt2rpp_Rev")  =
 0.

MidonLe: -R("R_EX_idon_L_e_") + R("R_EX_idon_L_e__Rev") -
 R("R_IDONtex") + R("R_IDONtex_Rev")  = 0.

MidonLp: -R("R_IDONt2rpp") + R("R_IDONt2rpp_Rev") + R("R_IDONtex") -
 R("R_IDONtex_Rev")  = 0.

Midpc: R("R_ADK4") - R("R_ADK4_Rev") + R("R_NTP10")  = 0.

MileLb: R("MileLb_Transp") + R("R_EX_ile_L_e_") -
 R("R_EX_ile_L_e__Rev")  = 0.

MileLc: -0.2905 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_ILEabcpp") + R("R_ILEt2rpp") - R("R_ILEt2rpp_Rev") - R("R_ILETA") +
 R("R_ILETA_Rev") - R("R_ILETRS")  = 0.

MileLe: -R("R_EX_ile_L_e_") + R("R_EX_ile_L_e__Rev") - R("R_ILEtex")
 + R("R_ILEtex_Rev")  = 0.

MileLp: -R("R_ILEabcpp") - R("R_ILEt2rpp") + R("R_ILEt2rpp_Rev") +
 R("R_ILEtex") - R("R_ILEtex_Rev")  = 0.

Miletrnac: R("R_ILETRS")  = 0.

Mimacpc: -R("R_HSTPT") + R("R_IGPDH")  = 0.

Mimpb: R("Mimpb_Transp") + R("R_EX_imp_e_") - R("R_EX_imp_e__Rev") 
 = 0.

Mimpc: -R("R_ADSS") + R("R_GMPR") + R("R_HXPRT") - R("R_IMPC") +
 R("R_IMPC_Rev") - R("R_IMPD") + R("R_INSK") - R("R_NTD11") + R("R_NTPP9")
  = 0.

Mimpe: -R("R_EX_imp_e_") + R("R_EX_imp_e__Rev") - R("R_IMPtex") +
 R("R_IMPtex_Rev")  = 0.

Mimpp: R("R_IMPtex") - R("R_IMPtex_Rev") - R("R_NTD11pp")  = 0.

Mindoleb: R("Mindoleb_Transp") + R("R_EX_indole_e_") -
 R("R_EX_indole_e__Rev")  = 0.

Mindolec: -R("R_INDOLEt2pp") + R("R_INDOLEt2rpp") -
 R("R_INDOLEt2rpp_Rev") + R("R_TRPAS2") - R("R_TRPAS2_Rev") - R("R_TRPS2") +
 R("R_TRPS3")  = 0.

Mindolee: -R("R_EX_indole_e_") + R("R_EX_indole_e__Rev") -
 R("R_INDOLEtex") + R("R_INDOLEtex_Rev")  = 0.

Mindolep: R("R_INDOLEt2pp") - R("R_INDOLEt2rpp") +
 R("R_INDOLEt2rpp_Rev") + R("R_INDOLEtex") - R("R_INDOLEtex_Rev")  = 0.

Minostb: R("Minostb_Transp") + R("R_EX_inost_e_") -
 R("R_EX_inost_e__Rev")  = 0.

Minostc: R("R_GPDDA5") + R("R_INOSTt4pp") + R("R_MI1PP")  = 0.

Minoste: -R("R_EX_inost_e_") + R("R_EX_inost_e__Rev") -
 R("R_INSTtex") + R("R_INSTtex_Rev")  = 0.

Minostp: R("R_GPDDA5pp") - R("R_INOSTt4pp") + R("R_INSTtex") -
 R("R_INSTtex_Rev") + R("R_PHYTSpp")  = 0.

Minsb: R("Minsb_Transp") + R("R_EX_ins_e_") - R("R_EX_ins_e__Rev") 
 = 0.

Minsc: R("R_ADA") - R("R_INSH") - R("R_INSK") + R("R_INSt2pp") +
 R("R_INSt2rpp") - R("R_INSt2rpp_Rev") + R("R_NTD11") - R("R_PUNP5") +
 R("R_PUNP5_Rev")  = 0.

Minse: -R("R_EX_ins_e_") + R("R_EX_ins_e__Rev") - R("R_INStex") +
 R("R_INStex_Rev")  = 0.

Minsp: -R("R_INSt2pp") - R("R_INSt2rpp") + R("R_INSt2rpp_Rev") +
 R("R_INStex") - R("R_INStex_Rev") + R("R_NTD11pp")  = 0.

Mipdpc: -R("R_DMATT") - R("R_GRTT") - R("R_IPDDI") +
 R("R_IPDDI_Rev") + R("R_IPDPS") - 5 R("R_OCTDPS") - 8 R("R_UDCPDPS")  = 0.

Misetacb: R("Misetacb_Transp") + R("R_EX_isetac_e_") -
 R("R_EX_isetac_e__Rev")  = 0.

Misetacc: -R("R_FDMO") + R("R_ISETACabcpp")  = 0.

Misetace: -R("R_EX_isetac_e_") + R("R_EX_isetac_e__Rev") -
 R("R_ISETACtex") + R("R_ISETACtex_Rev")  = 0.

Misetacp: -R("R_ISETACabcpp") + R("R_ISETACtex") -
 R("R_ISETACtex_Rev")  = 0.

Mitpc: -R("R_ADK4") + R("R_ADK4_Rev") + R("R_ATPHs") - R("R_NTP10")
 - R("R_NTPP9")  = 0.

Mkb: R("Mkb_Transp") + R("R_EX_k_e_") - R("R_EX_k_e__Rev")  = 0.

Mkc: -0.1776 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_Kabcpp") +
 R("R_Kt2pp") - R("R_Kt3pp")  = 0.

Mkdo2lipid4b: R("Mkdo2lipid4b_Transp") + R("R_EX_kdo2lipid4_e_") -
 R("R_EX_kdo2lipid4_e__Rev")  = 0.

Mkdo2lipid4c: -R("R_EDTXS1") - R("R_EDTXS3") - R("R_K2L4Aabcpp") +
 R("R_MOAT2")  = 0.

Mkdo2lipid4e: -0.01945 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_EX_kdo2lipid4_e_") + R("R_EX_kdo2lipid4_e__Rev") + R("R_K2L4Atex")  =
 0.

Mkdo2lipid4Lc: R("R_EDTXS1") - R("R_EDTXS2")  = 0.

Mkdo2lipid4p: R("R_K2L4Aabcpp") - R("R_K2L4Atex")  = 0.

Mkdo2lipid4pc: R("R_EDTXS3") - R("R_EDTXS4")  = 0.

Mkdo8pc: -R("R_KDOPP") + R("R_KDOPS")  = 0.

Mkdoc: -R("R_KDOCT2") + R("R_KDOPP")  = 0.

Mkdolipid4c: R("R_MOAT") - R("R_MOAT2")  = 0.

Mke: -R("R_EX_k_e_") + R("R_EX_k_e__Rev") - R("R_Ktex") +
 R("R_Ktex_Rev")  = 0.

Mkp: -R("R_Kabcpp") - R("R_Kt2pp") + R("R_Kt3pp") + R("R_Ktex") -
 R("R_Ktex_Rev")  = 0.

Mkphphhlipac: R("R_MOAT3C") - R("R_RHAT1")  = 0.

MlacDb: R("MlacDb_Transp") + R("R_EX_lac_D_e_") -
 R("R_EX_lac_D_e__Rev")  = 0.

MlacDc: R("R_ACM6PH") + R("R_D_LACt2pp") - R("R_D_LACt2pp_Rev") +
 R("R_GLYOX") + R("R_GLYOX3") - R("R_LDH_D") - R("R_LDH_D2") +
 R("R_LDH_D_Rev")  = 0.

MlacDe: -R("R_D_LACtex") + R("R_D_LACtex_Rev") - R("R_EX_lac_D_e_")
 + R("R_EX_lac_D_e__Rev")  = 0.

MlacDp: -R("R_D_LACt2pp") + R("R_D_LACt2pp_Rev") + R("R_D_LACtex") -
 R("R_D_LACtex_Rev")  = 0.

MlacLb: R("MlacLb_Transp") + R("R_EX_lac_L_e_") -
 R("R_EX_lac_L_e__Rev")  = 0.

MlacLc: R("R_LCADi") - R("R_L_LACD2") - R("R_L_LACD3") +
 R("R_L_LACt2rpp") - R("R_L_LACt2rpp_Rev")  = 0.

MlacLe: -R("R_EX_lac_L_e_") + R("R_EX_lac_L_e__Rev") -
 R("R_L_LACtex") + R("R_L_LACtex_Rev")  = 0.

MlacLp: -R("R_L_LACt2rpp") + R("R_L_LACt2rpp_Rev") + R("R_L_LACtex")
 - R("R_L_LACtex_Rev")  = 0.

MLalaDgluc: -R("R_ALAGLUE") + R("R_ALAGLUE_Rev") + R("R_LADGMDH")  =
 0.

MLalaDgluMdapb: R("MLalaDgluMdapb_Transp") +
 R("R_EX_LalaDgluMdap_e_") - R("R_EX_LalaDgluMdap_e__Rev")  = 0.

MLalaDgluMdapc: R("R_3PEPTabcpp") + R("R_4PCP") + R("R_AGM3PA") +
 R("R_AM3PA") - R("R_LADGMDH") - R("R_UM3PL")  = 0.

MLalaDgluMdapDalab: R("MLalaDgluMdapDalab_Transp") +
 R("R_EX_LalaDgluMdapDala_e_") - R("R_EX_LalaDgluMdapDala_e__Rev")  = 0.

MLalaDgluMdapDalac: -R("R_4PCP") + R("R_4PEPTabcpp") + R("R_AGM4PA")
 + R("R_AM4PA") - R("R_UM4PL")  = 0.

MLalaDgluMdapDalae: -R("R_4PEPTtex") + R("R_4PEPTtex_Rev") -
 R("R_EX_LalaDgluMdapDala_e_") + R("R_EX_LalaDgluMdapDala_e__Rev")  = 0.

MLalaDgluMdapDalap: -R("R_4PCPpp") - R("R_4PEPTabcpp") +
 R("R_4PEPTtex") - R("R_4PEPTtex_Rev") + R("R_AGM4PApp")  = 0.

MLalaDgluMdape: -R("R_3PEPTtex") + R("R_3PEPTtex_Rev") -
 R("R_EX_LalaDgluMdap_e_") + R("R_EX_LalaDgluMdap_e__Rev")  = 0.

MLalaDgluMdapp: -R("R_3PEPTabcpp") + R("R_3PEPTtex") -
 R("R_3PEPTtex_Rev") + R("R_4PCPpp") + R("R_AGM3PApp")  = 0.

MLalaLgluc: R("R_ALAGLUE") - R("R_ALAGLUE_Rev") - R("R_LALGP")  = 0.

MlaldDc: R("R_LALDO2x") - R("R_LCARR") + R("R_LCARR_Rev")  = 0.

MlaldLc: R("R_FCLPA") - R("R_FCLPA_Rev") - R("R_LCADi") -
 R("R_LCARS") + R("R_LCARS_Rev") + R("R_RMPA") - R("R_RMPA_Rev")  = 0.

Mlctsb: R("Mlctsb_Transp") + R("R_EX_lcts_e_") -
 R("R_EX_lcts_e__Rev")  = 0.

Mlctsc: -R("R_LACZ") - R("R_LCTSt3ipp") + R("R_LCTStpp") -
 R("R_LCTStpp_Rev")  = 0.

Mlctse: -R("R_EX_lcts_e_") + R("R_EX_lcts_e__Rev") - R("R_LCTStex")
 + R("R_LCTStex_Rev")  = 0.

Mlctsp: -R("R_LACZpp") + R("R_LCTSt3ipp") + R("R_LCTStex") -
 R("R_LCTStex_Rev") - R("R_LCTStpp") + R("R_LCTStpp_Rev")  = 0.

MleuLb: R("MleuLb_Transp") + R("R_EX_leu_L_e_") -
 R("R_EX_leu_L_e__Rev")  = 0.

MleuLc: -0.4505 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_LEUabcpp") + R("R_LEUt2rpp") - R("R_LEUt2rpp_Rev") + R("R_LEUTAi") -
 R("R_LEUTRS")  = 0.

MleuLe: -R("R_EX_leu_L_e_") + R("R_EX_leu_L_e__Rev") - R("R_LEUtex")
 + R("R_LEUtex_Rev")  = 0.

MleuLp: -R("R_LEUabcpp") - R("R_LEUt2rpp") + R("R_LEUt2rpp_Rev") +
 R("R_LEUtex") - R("R_LEUtex_Rev")  = 0.

Mleutrnac: R("R_LEUTRS")  = 0.

MlgtSc: -R("R_GLYOX") + R("R_LGTHL")  = 0.

Mlipab: R("Mlipab_Transp") + R("R_EX_lipa_e_") -
 R("R_EX_lipa_e__Rev")  = 0.

Mlipac: R("R_EDTXS2") - R("R_HEPT1") - R("R_LIPAabcpp")  = 0.

Mlipacoldb: R("Mlipacoldb_Transp") + R("R_EX_lipa_cold_e_") -
 R("R_EX_lipa_cold_e__Rev")  = 0.

Mlipacoldc: R("R_EDTXS4") - R("R_LIPACabcpp")  = 0.

Mlipacolde: R("R_CLIPAtex") - R("R_EX_lipa_cold_e_") +
 R("R_EX_lipa_cold_e__Rev")  = 0.

Mlipacoldp: -R("R_CLIPAtex") + R("R_LIPACabcpp")  = 0.

Mlipae: -R("R_EX_lipa_e_") + R("R_EX_lipa_e__Rev") - R("R_LIPAHTex")
 + R("R_LIPAtex")  = 0.

Mlipap: R("R_LIPAabcpp") - R("R_LIPAtex") - R("R_PETNT161pp") -
 R("R_PETNT181pp")  = 0.

MlipidAc: -R("R_MOAT") + R("R_TDSK")  = 0.

MlipidAdsc: R("R_LPADSS") - R("R_TDSK")  = 0.

MlipidXc: -R("R_LPADSS") + R("R_USHD")  = 0.

Mlppp: R("R_ALPATE160pp") + R("R_ALPATG160pp")  = 0.

MlysLb: R("MlysLb_Transp") + R("R_EX_lys_L_e_") -
 R("R_EX_lys_L_e__Rev")  = 0.

MlysLc: R("R_CADVtpp") + R("R_DAPDC") - 0.3432
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_FRULYSDG") - R("R_FRULYSDG_Rev") +
 R("R_LYSabcpp") - R("R_LYSDC") + R("R_LYSt2pp") - R("R_LYSt3pp") - R("R_LYSTRS") 
 = 0.

MlysLe: -R("R_EX_lys_L_e_") + R("R_EX_lys_L_e__Rev") - R("R_LYStex")
 + R("R_LYStex_Rev")  = 0.

MlysLp: -R("R_CADVtpp") - R("R_LYSabcpp") - R("R_LYSt2pp") +
 R("R_LYSt3pp") + R("R_LYStex") - R("R_LYStex_Rev")  = 0.

Mlystrnac: R("R_LYSTRS")  = 0.

MlyxLb: R("MlyxLb_Transp") + R("R_EX_lyx_L_e_") -
 R("R_EX_lyx_L_e__Rev")  = 0.

MlyxLc: -R("R_LYXI") + R("R_LYXt2pp")  = 0.

MlyxLe: -R("R_EX_lyx_L_e_") + R("R_EX_lyx_L_e__Rev") - R("R_LYXtex")
 + R("R_LYXtex_Rev")  = 0.

MlyxLp: -R("R_LYXt2pp") + R("R_LYXtex") - R("R_LYXtex_Rev")  = 0.

MmalACPc: -R("R_3OAS100") - R("R_3OAS120") - R("R_3OAS121") -
 R("R_3OAS140") - R("R_3OAS141") - R("R_3OAS160") - R("R_3OAS161") -
 R("R_3OAS180") - R("R_3OAS181") - R("R_3OAS60") - R("R_3OAS80") - R("R_KAS14")
 - R("R_KAS15") - R("R_MACPD") + R("R_MCOATA") - R("R_MCOATA_Rev")  =
 0.

Mmalcoac: R("R_ACCOAC") - R("R_MCOATA") + R("R_MCOATA_Rev")  = 0.

MmalDb: R("MmalDb_Transp") + R("R_EX_mal_D_e_") -
 R("R_EX_mal_D_e__Rev")  = 0.

MmalDc: -R("R_MALDDH") + R("R_MALDt2_2pp")  = 0.

MmalDe: -R("R_EX_mal_D_e_") + R("R_EX_mal_D_e__Rev") -
 R("R_MALDtex") + R("R_MALDtex_Rev")  = 0.

MmalDp: -R("R_MALDt2_2pp") + R("R_MALDtex") - R("R_MALDtex_Rev")  =
 0.

MmalLb: R("MmalLb_Transp") + R("R_EX_mal_L_e_") -
 R("R_EX_mal_L_e__Rev")  = 0.

MmalLc: R("R_FUM") - R("R_FUM_Rev") + R("R_MALS") + R("R_MALt2_2pp")
 + R("R_MALt2_3pp") - R("R_MDH") - R("R_MDH2") - R("R_MDH3") +
 R("R_MDH_Rev") - R("R_ME1") - R("R_ME2")  = 0.

MmalLe: -R("R_EX_mal_L_e_") + R("R_EX_mal_L_e__Rev") - R("R_MALtex")
 + R("R_MALtex_Rev")  = 0.

MmalLp: -R("R_MALt2_2pp") - R("R_MALt2_3pp") + R("R_MALtex") -
 R("R_MALtex_Rev")  = 0.

Mmalt6pc: R("R_MALTptspp")  = 0.

Mmaltb: R("Mmaltb_Transp") + R("R_EX_malt_e_") -
 R("R_EX_malt_e__Rev")  = 0.

Mmaltc: -R("R_AMALT1") - R("R_AMALT2") - R("R_AMALT3") -
 R("R_AMALT4") + R("R_MALTabcpp") - R("R_MALTATr") + R("R_MALTATr_Rev") +
 R("R_MLTG1")  = 0.

Mmalte: -R("R_EX_malt_e_") + R("R_EX_malt_e__Rev") - R("R_MALTtexi")
  = 0.

Mmalthpc: R("R_AMALT4") - R("R_MLTG5") - R("R_MLTP3") +
 R("R_MLTP3_Rev")  = 0.

Mmalthxb: R("Mmalthxb_Transp") + R("R_EX_malthx_e_") -
 R("R_EX_malthx_e__Rev")  = 0.

Mmalthxc: R("R_AAMYL") + R("R_AMALT3") - R("R_AMALT4") +
 R("R_MALTHXabcpp") - R("R_MLTG4") + R("R_MLTG5") - R("R_MLTP2") +
 R("R_MLTP2_Rev") + R("R_MLTP3") - R("R_MLTP3_Rev")  = 0.

Mmalthxe: -R("R_EX_malthx_e_") + R("R_EX_malthx_e__Rev") -
 R("R_MALTHXtexi")  = 0.

Mmalthxp: R("R_AAMYLpp") - R("R_MALTHXabcpp") + R("R_MALTHXtexi")  =
 0.

Mmaltp: -R("R_MALTabcpp") - R("R_MALTptspp") + R("R_MALTtexi")  = 0.

Mmaltptb: R("Mmaltptb_Transp") + R("R_EX_maltpt_e_") -
 R("R_EX_maltpt_e__Rev")  = 0.

Mmaltptc: R("R_AMALT2") - R("R_AMALT3") + R("R_MALTPTabcpp") -
 R("R_MLTG3") + R("R_MLTG4") - R("R_MLTP1") + R("R_MLTP1_Rev") +
 R("R_MLTP2") - R("R_MLTP2_Rev")  = 0.

Mmaltpte: -R("R_EX_maltpt_e_") + R("R_EX_maltpt_e__Rev") -
 R("R_MALTPTtexi")  = 0.

Mmaltptp: -R("R_MALTPTabcpp") + R("R_MALTPTtexi")  = 0.

Mmalttrb: R("Mmalttrb_Transp") + R("R_EX_malttr_e_") -
 R("R_EX_malttr_e__Rev")  = 0.

Mmalttrc: -R("R_AMALT1") + R("R_MALTTRabcpp") - R("R_MLTG1") +
 R("R_MLTG2")  = 0.

Mmalttre: -R("R_EX_malttr_e_") + R("R_EX_malttr_e__Rev") -
 R("R_MALTTRtexi")  = 0.

Mmalttrp: -R("R_MALTTRabcpp") + R("R_MALTTRtexi")  = 0.

Mmaltttrb: R("Mmaltttrb_Transp") + R("R_EX_maltttr_e_") -
 R("R_EX_maltttr_e__Rev")  = 0.

Mmaltttrc: R("R_AMALT1") - R("R_AMALT2") + R("R_MALTTTRabcpp") -
 R("R_MLTG2") + R("R_MLTG3") + R("R_MLTP1") - R("R_MLTP1_Rev")  = 0.

Mmaltttre: -R("R_EX_maltttr_e_") + R("R_EX_maltttr_e__Rev") -
 R("R_MALTTTRtexi")  = 0.

Mmaltttrp: -R("R_MALTTTRabcpp") + R("R_MALTTTRtexi")  = 0.

Mman1pc: R("R_GDPMNP") - R("R_MAN1PT2") - R("R_PMANM") +
 R("R_PMANM_Rev")  = 0.

Mman6pb: R("Mman6pb_Transp") + R("R_EX_man6p_e_") -
 R("R_EX_man6p_e__Rev")  = 0.

Mman6pc: R("R_HEX4") - R("R_MAN6PI") + R("R_MAN6PI_Rev") +
 R("R_MAN6Pt6_2pp") + R("R_MANPGH") + R("R_MANptspp") - R("R_MN6PP") +
 R("R_PMANM") - R("R_PMANM_Rev")  = 0.

Mman6pe: -R("R_EX_man6p_e_") + R("R_EX_man6p_e__Rev") -
 R("R_MAN6Ptex") + R("R_MAN6Ptex_Rev")  = 0.

Mman6pglycc: R("R_MANGLYCptspp") - R("R_MANPGH")  = 0.

Mman6pp: -R("R_MAN6Pt6_2pp") + R("R_MAN6Ptex") - R("R_MAN6Ptex_Rev")
  = 0.

Mmanac: -R("R_MANAO") + R("R_MANAO_Rev") - R("R_MNNH")  = 0.

Mmanb: R("Mmanb_Transp") + R("R_EX_man_e_") - R("R_EX_man_e__Rev") 
 = 0.

Mmanc: R("R_GDPMNH") - R("R_HEX4") + R("R_MN6PP")  = 0.

Mmane: -R("R_EX_man_e_") + R("R_EX_man_e__Rev") - R("R_MANtex") +
 R("R_MANtex_Rev")  = 0.

Mmanglycb: R("Mmanglycb_Transp") + R("R_EX_manglyc_e_") -
 R("R_EX_manglyc_e__Rev")  = 0.

Mmanglyce: -R("R_EX_manglyc_e_") + R("R_EX_manglyc_e__Rev") -
 R("R_MANGLYCtex") + R("R_MANGLYCtex_Rev")  = 0.

Mmanglycp: -R("R_MANGLYCptspp") + R("R_MANGLYCtex") -
 R("R_MANGLYCtex_Rev")  = 0.

Mmanp: -R("R_MANptspp") + R("R_MANtex") - R("R_MANtex_Rev")  = 0.

Mmelibb: R("Mmelibb_Transp") + R("R_EX_melib_e_") -
 R("R_EX_melib_e__Rev")  = 0.

Mmelibc: -R("R_GALS3") + R("R_MELIBt2pp") - R("R_MELIBt3ipp")  = 0.

Mmelibe: -R("R_EX_melib_e_") + R("R_EX_melib_e__Rev") -
 R("R_MELIBtex") + R("R_MELIBtex_Rev")  = 0.

Mmelibp: -R("R_MELIBt2pp") + R("R_MELIBt3ipp") + R("R_MELIBtex") -
 R("R_MELIBtex_Rev")  = 0.

Mmercppyrc: -R("R_MCPST")  = 0.

MmetDb: R("MmetDb_Transp") + R("R_EX_met_D_e_") -
 R("R_EX_met_D_e__Rev")  = 0.

MmetDc: R("R_METDabcpp")  = 0.

MmetDe: -R("R_EX_met_D_e_") + R("R_EX_met_D_e__Rev") -
 R("R_METDtex") + R("R_METDtex_Rev")  = 0.

MmetDp: -R("R_METDabcpp") + R("R_METDtex") - R("R_METDtex_Rev")  =
 0.

Mmethfc: -R("R_MTHFC") + R("R_MTHFC_Rev") + R("R_MTHFD") -
 R("R_MTHFD_Rev")  = 0.

MmetLb: R("MmetLb_Transp") + R("R_EX_met_L_e_") -
 R("R_EX_met_L_e__Rev")  = 0.

MmetLc: R("R_BTS4") + 2 R("R_CPPPGO2") - 0.1537
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_HCYSMT") + 2 R("R_HCYSMT2") + R("R_METabcpp")
 - R("R_METAT") - R("R_METOX1s") - R("R_METOX2s") + R("R_METS") +
 R("R_METSOXR1") + R("R_METSOXR2") - R("R_METTRS") + R("R_UNK3")  = 0.

MmetLe: -R("R_EX_met_L_e_") + R("R_EX_met_L_e__Rev") - R("R_METtex")
 + R("R_METtex_Rev")  = 0.

MmetLp: -R("R_METabcpp") + R("R_METtex") - R("R_METtex_Rev")  = 0.

MmetsoxRLb: R("MmetsoxRLb_Transp") + R("R_EX_metsox_R_L_e_") -
 R("R_EX_metsox_R_L_e__Rev")  = 0.

MmetsoxRLc: R("R_METOX2s") + R("R_METSOX2abcpp") - R("R_METSOXR2") 
 = 0.

MmetsoxRLe: -R("R_EX_metsox_R_L_e_") + R("R_EX_metsox_R_L_e__Rev") -
 R("R_METSOX2tex") + R("R_METSOX2tex_Rev")  = 0.

MmetsoxRLp: -R("R_METSOX2abcpp") + R("R_METSOX2tex") -
 R("R_METSOX2tex_Rev")  = 0.

MmetsoxSLb: R("MmetsoxSLb_Transp") + R("R_EX_metsox_S_L_e_") -
 R("R_EX_metsox_S_L_e__Rev")  = 0.

MmetsoxSLc: R("R_METOX1s") + R("R_METSOX1abcpp") - R("R_METSOXR1") 
 = 0.

MmetsoxSLe: -R("R_EX_metsox_S_L_e_") + R("R_EX_metsox_S_L_e__Rev") -
 R("R_METSOX1tex") + R("R_METSOX1tex_Rev")  = 0.

MmetsoxSLp: -R("R_METSOX1abcpp") + R("R_METSOX1tex") -
 R("R_METSOX1tex_Rev")  = 0.

Mmettrnac: -R("R_FMETTRS") + R("R_METTRS")  = 0.

Mmg2b: R("Mmg2b_Transp") + R("R_EX_mg2_e_") - R("R_EX_mg2_e__Rev") 
 = 0.

Mmg2c: -0.007895 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_MG2t3_2pp") - R("R_MG2t3_2pp_Rev") + R("R_MG2tpp") + R("R_MG2uabcpp")  = 0.

Mmg2e: -R("R_EX_mg2_e_") + R("R_EX_mg2_e__Rev") - R("R_MG2tex") +
 R("R_MG2tex_Rev")  = 0.

Mmg2p: -R("R_MG2t3_2pp") + R("R_MG2t3_2pp_Rev") + R("R_MG2tex") -
 R("R_MG2tex_Rev") - R("R_MG2tpp") - R("R_MG2uabcpp")  = 0.

Mmi1pDc: -R("R_MI1PP")  = 0.

Mmicitc: -R("R_MCITL2") + R("R_MCITL2_Rev") + R("R_MICITD")  = 0.

Mminohpb: R("Mminohpb_Transp") + R("R_EX_minohp_e_") -
 R("R_EX_minohp_e__Rev")  = 0.

Mminohpe: -R("R_EX_minohp_e_") + R("R_EX_minohp_e__Rev") -
 R("R_MINOHPtexi")  = 0.

Mminohpp: R("R_MINOHPtexi") - R("R_PHYTSpp")  = 0.

Mmlthfc: -0.000223 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_GHMT2r") - R("R_GHMT2r_Rev") + R("R_GLYCL") - R("R_MOHMT") - R("R_MTHFD")
 + R("R_MTHFD_Rev") - R("R_MTHFR2") - R("R_TMDS")  = 0.

MmmcoaRc: -R("R_MME") + R("R_MME_Rev") + R("R_MMM2")  = 0.

MmmcoaSc: -R("R_MMCD") + R("R_MME") - R("R_MME_Rev")  = 0.

Mmmetb: R("Mmmetb_Transp") + R("R_EX_mmet_e_") -
 R("R_EX_mmet_e__Rev")  = 0.

Mmmetc: -R("R_HCYSMT2") + R("R_MMETt2pp")  = 0.

Mmmete: -R("R_EX_mmet_e_") + R("R_EX_mmet_e__Rev") - R("R_MMETtex")
 + R("R_MMETtex_Rev")  = 0.

Mmmetp: -R("R_MMETt2pp") + R("R_MMETtex") - R("R_MMETtex_Rev")  = 0.

Mmn2b: R("Mmn2b_Transp") + R("R_EX_mn2_e_") - R("R_EX_mn2_e__Rev") 
 = 0.

Mmn2c: -0.003158 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_MN2t3pp") + R("R_MN2tpp") + R("R_MNt2pp")  = 0.

Mmn2e: -R("R_EX_mn2_e_") + R("R_EX_mn2_e__Rev") - R("R_MNtex") +
 R("R_MNtex_Rev")  = 0.

Mmn2p: R("R_MN2t3pp") - R("R_MN2tpp") - R("R_MNt2pp") + R("R_MNtex")
 - R("R_MNtex_Rev")  = 0.

Mmnl1pc: -R("R_M1PD") + R("R_M1PD_Rev") + R("R_MNLptspp")  = 0.

Mmnlb: R("Mmnlb_Transp") + R("R_EX_mnl_e_") - R("R_EX_mnl_e__Rev") 
 = 0.

Mmnle: -R("R_EX_mnl_e_") + R("R_EX_mnl_e__Rev") - R("R_MNLtex") +
 R("R_MNLtex_Rev")  = 0.

Mmnlp: -R("R_MNLptspp") + R("R_MNLtex") - R("R_MNLtex_Rev")  = 0.

Mmobdb: R("Mmobdb_Transp") + R("R_EX_mobd_e_") -
 R("R_EX_mobd_e__Rev")  = 0.

Mmobdc: -0.003158 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_MOBDabcpp")  = 0.

Mmobde: -R("R_EX_mobd_e_") + R("R_EX_mobd_e__Rev") - R("R_MOBDtex")
 + R("R_MOBDtex_Rev")  = 0.

Mmobdp: -R("R_MOBDabcpp") + R("R_MOBDtex") - R("R_MOBDtex_Rev")  =
 0.

Mmql8c: R("R_AMMQLT8") + R("R_ASPO4") - R("R_CYTBD2pp") +
 R("R_DHORD5") - R("R_DMSOR1") - R("R_DMSOR1pp") + R("R_DSBAO2") + R("R_FDH5pp")
 - R("R_FRD2") + R("R_G3PD6") + R("R_GLYCTO3") + R("R_HYD2pp") +
 R("R_L_LACD3") + R("R_MDH3") + R("R_NADH10") + R("R_NADH17pp") +
 R("R_NADPHQR3") - R("R_NO3R2bpp") - R("R_NO3R2pp") - 3 R("R_NTRIR4pp") -
 R("R_QMO3") - R("R_TMAOR1") - R("R_TMAOR1pp")  = 0.

Mmqn8c: -R("R_ASPO4") + R("R_CYTBD2pp") - R("R_DHORD5") +
 R("R_DMSOR1") + R("R_DMSOR1pp") - R("R_DSBAO2") - R("R_FDH5pp") + R("R_FRD2") -
 R("R_G3PD6") - R("R_GLYCTO3") - R("R_HYD2pp") - R("R_L_LACD3") -
 R("R_MDH3") - R("R_NADH10") - R("R_NADH17pp") - R("R_NADPHQR3") +
 R("R_NO3R2bpp") + R("R_NO3R2pp") + 3 R("R_NTRIR4pp") + R("R_QMO3") +
 R("R_TMAOR1") + R("R_TMAOR1pp")  = 0.

Mmso3b: R("Mmso3b_Transp") + R("R_EX_mso3_e_") -
 R("R_EX_mso3_e__Rev")  = 0.

Mmso3c: -R("R_FDMO2") + R("R_MSO3abcpp")  = 0.

Mmso3e: -R("R_EX_mso3_e_") + R("R_EX_mso3_e__Rev") - R("R_MSO3tex")
 + R("R_MSO3tex_Rev")  = 0.

Mmso3p: -R("R_MSO3abcpp") + R("R_MSO3tex") - R("R_MSO3tex_Rev")  =
 0.

Mmthgxlc: R("R_AACTOOR") - R("R_ALR2") - R("R_GLYOX3") -
 R("R_LALDO2x") - R("R_LGTHL") + R("R_MGSA")  = 0.

Mmurein3p3pp: R("R_MLDCP2Bpp") + R("R_MLDEP1pp") - R("R_MLTGY3pp") 
 = 0.

Mmurein3px3pp: R("R_MLDCP3App") - R("R_MLDEP1pp")  = 0.

Mmurein3px4pp: -R("R_MDDEP2pp") + R("R_MLDCP1App")  = 0.

Mmurein4p3pp: R("R_MDDCP5pp") + R("R_MDDEP2pp") + R("R_MLDCP1Bpp") -
 R("R_MLDCP2Bpp") - R("R_MLTGY2pp")  = 0.

Mmurein4p4pp: R("R_MDDCP4pp") + R("R_MDDEP1pp") - R("R_MLDCP1Bpp") -
 R("R_MLTGY1pp")  = 0.

Mmurein4px4p4pp: R("R_MDDEP4pp") - R("R_MLTGY4pp")  = 0.

Mmurein4px4pp: R("R_MDDCP1pp") - R("R_MDDEP1pp") + R("R_MLTGY4pp") 
 = 0.

Mmurein4px4px4pp: R("R_MDDCP2pp") - R("R_MDDEP4pp")  = 0.

Mmurein5p3pp: -R("R_MDDCP5pp") + R("R_MLDCP2App") + R("R_MLDEP2pp") 
 = 0.

Mmurein5p4pp: R("R_MDDCP3pp") - R("R_MDDCP4pp") + R("R_MDDEP3pp")  =
 0.

Mmurein5p5p5pp: -R("R_MCTP2App") + R("R_MPTG2")  = 0.

Mmurein5p5pp: -R("R_MCTP1App") - R("R_MCTP1Bpp") - R("R_MDDCP3pp") -
 R("R_MLDCP2App") + R("R_MPTG") - R("R_MPTG2")  = 0.

Mmurein5px3pp: R("R_MCTP1Bpp") - R("R_MLDCP3App") - R("R_MLDEP2pp") 
 = 0.

Mmurein5px4pp: -0.01389 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_MCTP1App") - R("R_MDDCP1pp") - R("R_MDDEP3pp") - R("R_MLDCP1App")  =
 0.

Mmurein5px4px4pp: R("R_MCTP2App") - R("R_MDDCP2pp")  = 0.

MmyrsACPc: -R("R_3OAS160") + R("R_AACPS1") - R("R_AGPAT140") +
 R("R_EAR140x") + R("R_EAR140y") - R("R_EDTXS2") - R("R_EDTXS4") -
 R("R_FA140ACPHi") - R("R_G3PAT140")  = 0.

MN1aspmdc: R("R_SPMDAT1")  = 0.

Mn2ob: R("Mn2ob_Transp") + R("R_EX_n2o_e_") - R("R_EX_n2o_e__Rev") 
 = 0.

Mn2oc: R("R_N2Otpp") - R("R_N2Otpp_Rev") + R("R_NHFRBO")  = 0.

Mn2oe: -R("R_EX_n2o_e_") + R("R_EX_n2o_e__Rev") - R("R_N2Otex") +
 R("R_N2Otex_Rev")  = 0.

Mn2op: R("R_N2Otex") - R("R_N2Otex_Rev") - R("R_N2Otpp") +
 R("R_N2Otpp_Rev")  = 0.

Mn8aspmdc: R("R_SPMDAT2")  = 0.

Mna1b: R("Mna1b_Transp") + R("R_EX_na1_e_") - R("R_EX_na1_e__Rev") 
 = 0.

Mna1c: R("R_ACt4pp") + R("R_ALAt4pp") + R("R_CAt6pp") -
 R("R_CAt6pp_Rev") + R("R_GLUt4pp") + R("R_GLYCLTt4pp") + R("R_GLYt4pp") +
 R("R_INOSTt4pp") - 2 R("R_NAt3_1p5pp") - R("R_NAt3_2pp") - R("R_NAt3pp") +
 R("R_PNTOt4pp") + R("R_PPAt4pp") + R("R_PROt4pp") + R("R_SERt4pp") +
 R("R_THRt4pp")  = 0.

Mna1e: -R("R_EX_na1_e_") + R("R_EX_na1_e__Rev") - R("R_NAtex") +
 R("R_NAtex_Rev")  = 0.

Mna1p: -R("R_ACt4pp") - R("R_ALAt4pp") - R("R_CAt6pp") +
 R("R_CAt6pp_Rev") - R("R_GLUt4pp") - R("R_GLYCLTt4pp") - R("R_GLYt4pp") -
 R("R_INOSTt4pp") + 2 R("R_NAt3_1p5pp") + R("R_NAt3_2pp") + R("R_NAt3pp") +
 R("R_NAtex") - R("R_NAtex_Rev") - R("R_PNTOt4pp") - R("R_PPAt4pp") -
 R("R_PROt4pp") - R("R_SERt4pp") - R("R_THRt4pp")  = 0.

Mnacb: R("Mnacb_Transp") + R("R_EX_nac_e_") - R("R_EX_nac_e__Rev") 
 = 0.

Mnacc: R("R_NACtpp") - R("R_NAMNPP") + R("R_NNAM") + R("R_NNDMBRT") 
 = 0.

Mnace: -R("R_EX_nac_e_") + R("R_EX_nac_e__Rev") - R("R_NACtex") +
 R("R_NACtex_Rev")  = 0.

Mnacp: R("R_NACtex") - R("R_NACtex_Rev") - R("R_NACtpp")  = 0.

Mnadc: R("R_2DGLCNRx") + R("R_2DGULRx") + R("R_3HCINNMH") +
 R("R_3HPPPNH") - R("R_ABUTD") - R("R_ACALD") + R("R_ACALD_Rev") - R("R_AKGDH")
 + R("R_ALCD19") - R("R_ALCD19_Rev") - R("R_ALCD2x") +
 R("R_ALCD2x_Rev") - R("R_ALDD19x") - R("R_ALDD2x") - R("R_ALDD4") + R("R_ALR4x") -
 R("R_AMPMS2") - R("R_BETALDHx") + R("R_BSORx") + R("R_CINNDO") -
 R("R_DHBD") + R("R_DHBD_Rev") - R("R_DHCIND") - R("R_DHPPD") +
 R("R_DKGLCNR2x") - 2 R("R_DKMPPD3") + R("R_DMPPS") + R("R_DOGULNR") - R("R_E4PD") +
 R("R_E4PD_Rev") + R("R_EAR100x") + R("R_EAR120x") + R("R_EAR121x") +
 R("R_EAR140x") + R("R_EAR141x") + R("R_EAR160x") + R("R_EAR161x") +
 R("R_EAR180x") + R("R_EAR181x") + R("R_EAR40x") + R("R_EAR60x") +
 R("R_EAR80x") - 0.001831 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_FADRx")
 - R("R_FALDH2") + R("R_FALDH2_Rev") + R("R_FLVRx") + R("R_FMNRx") -
 R("R_GALCTLO") - R("R_GAPD") + R("R_GAPD_Rev") - R("R_GCALDD") -
 R("R_GLTPD") + R("R_GLTPD_Rev") - R("R_GLYCDx") - R("R_GLYCL") +
 R("R_GLYCLTDx") - R("R_HACD1i") - R("R_HACD2i") - R("R_HACD3i") - R("R_HACD4i") -
 R("R_HACD5i") - R("R_HACD6i") - R("R_HACD7i") - R("R_HACD8i") - 2
 R("R_HISTD") + R("R_HPYRRx") - R("R_HXAND") + R("R_IDOND") -
 R("R_IDOND_Rev") - R("R_IMPD") + R("R_IPDPS") - R("R_IPMD") + R("R_LALDO2x") -
 R("R_LCADi") + R("R_LCARR") - R("R_LCARR_Rev") + R("R_LCARS") -
 R("R_LCARS_Rev") - R("R_LDH_D") + R("R_LDH_D_Rev") - R("R_M1PD") +
 R("R_M1PD_Rev") - R("R_MALDDH") - R("R_MANAO") + R("R_MANAO_Rev") - R("R_MDH") +
 R("R_MDH_Rev") - R("R_ME1") + R("R_MECDPDH2") + R("R_MTHFR2") -
 R("R_NADDP") + R("R_NADH10") + R("R_NADH16pp") + R("R_NADH17pp") +
 R("R_NADH18pp") + R("R_NADH5") + R("R_NADH9") - R("R_NADK") - R("R_NADN") +
 R("R_NADPPPS") + R("R_NADS1") - R("R_NADTRHD") + R("R_NHFRBO") +
 R("R_NMNAT") + R("R_NODOx") + 3 R("R_NTRIR2x") - R("R_OMMBLHX3") -
 R("R_OMPHHX3") - R("R_OPHHX3") - R("R_P5CD") - R("R_PDH") - R("R_PDX5PO2") -
 R("R_PDX5PS") - R("R_PERD") + R("R_PERD_Rev") - R("R_PGCD") - R("R_PPND") +
 R("R_PPPNDO") - R("R_SBTPD") + R("R_SBTPD_Rev") - R("R_SGSAD") -
 R("R_SHCHD2") - R("R_SSALx") - R("R_TAGURr") + R("R_TAGURr_Rev") +
 R("R_THD2pp") - R("R_THRD") + R("R_TRSARr") - R("R_TRSARr_Rev") - 2
 R("R_UACMAMO") - 2 R("R_UDPGD") - R("R_UDPGDC") - R("R_URDGLYCD") - R("R_XAND")
  = 0.

Mnadhc: -R("R_2DGLCNRx") - R("R_2DGULRx") - R("R_3HCINNMH") -
 R("R_3HPPPNH") + R("R_ABUTD") + R("R_ACALD") - R("R_ACALD_Rev") +
 R("R_AKGDH") - R("R_ALCD19") + R("R_ALCD19_Rev") + R("R_ALCD2x") -
 R("R_ALCD2x_Rev") + R("R_ALDD19x") + R("R_ALDD2x") + R("R_ALDD4") - R("R_ALR4x") +
 R("R_AMPMS2") + R("R_BETALDHx") - R("R_BSORx") - R("R_CINNDO") +
 R("R_DHBD") - R("R_DHBD_Rev") + R("R_DHCIND") + R("R_DHPPD") -
 R("R_DKGLCNR2x") + 2 R("R_DKMPPD3") - R("R_DMPPS") - R("R_DOGULNR") + R("R_E4PD")
 - R("R_E4PD_Rev") - R("R_EAR100x") - R("R_EAR120x") - R("R_EAR121x") -
 R("R_EAR140x") - R("R_EAR141x") - R("R_EAR160x") - R("R_EAR161x") -
 R("R_EAR180x") - R("R_EAR181x") - R("R_EAR40x") - R("R_EAR60x") -
 R("R_EAR80x") - R("R_FADRx") + R("R_FALDH2") - R("R_FALDH2_Rev") -
 R("R_FLVRx") - R("R_FMNRx") + R("R_GALCTLO") + R("R_GAPD") - R("R_GAPD_Rev") +
 R("R_GCALDD") + R("R_GLTPD") - R("R_GLTPD_Rev") + R("R_GLYCDx") +
 R("R_GLYCL") - R("R_GLYCLTDx") + R("R_HACD1i") + R("R_HACD2i") +
 R("R_HACD3i") + R("R_HACD4i") + R("R_HACD5i") + R("R_HACD6i") + R("R_HACD7i")
 + R("R_HACD8i") + 2 R("R_HISTD") - R("R_HPYRRx") + R("R_HXAND") -
 R("R_IDOND") + R("R_IDOND_Rev") + R("R_IMPD") - R("R_IPDPS") + R("R_IPMD")
 - R("R_LALDO2x") + R("R_LCADi") - R("R_LCARR") + R("R_LCARR_Rev") -
 R("R_LCARS") + R("R_LCARS_Rev") + R("R_LDH_D") - R("R_LDH_D_Rev") +
 R("R_M1PD") - R("R_M1PD_Rev") + R("R_MALDDH") + R("R_MANAO") -
 R("R_MANAO_Rev") + R("R_MDH") - R("R_MDH_Rev") + R("R_ME1") - R("R_MECDPDH2") -
 R("R_MTHFR2") - R("R_NADH10") - R("R_NADH16pp") - R("R_NADH17pp") -
 R("R_NADH18pp") - R("R_NADH5") - R("R_NADH9") + R("R_NADTRHD") -
 R("R_NHFRBO") - R("R_NODOx") - 3 R("R_NTRIR2x") + R("R_OMMBLHX3") +
 R("R_OMPHHX3") + R("R_OPHHX3") + R("R_P5CD") + R("R_PDH") + R("R_PDX5PO2") +
 R("R_PDX5PS") + R("R_PERD") - R("R_PERD_Rev") + R("R_PGCD") +
 R("R_PPND") - R("R_PPPNDO") + R("R_SBTPD") - R("R_SBTPD_Rev") + R("R_SGSAD") +
 R("R_SHCHD2") + R("R_SSALx") + R("R_TAGURr") - R("R_TAGURr_Rev") -
 R("R_THD2pp") + R("R_THRD") - R("R_TRSARr") + R("R_TRSARr_Rev") + 2
 R("R_UACMAMO") + 2 R("R_UDPGD") + R("R_UDPGDC") + R("R_URDGLYCD") +
 R("R_XAND")  = 0.

Mnadpc: R("R_2DGLCNRy") + R("R_2DGULRy") + R("R_3OAR100") -
 R("R_3OAR100_Rev") + R("R_3OAR120") - R("R_3OAR120_Rev") + R("R_3OAR121") +
 R("R_3OAR140") - R("R_3OAR140_Rev") + R("R_3OAR141") + R("R_3OAR160") -
 R("R_3OAR160_Rev") + R("R_3OAR161") + R("R_3OAR180") -
 R("R_3OAR180_Rev") + R("R_3OAR181") + R("R_3OAR40") - R("R_3OAR40_Rev") +
 R("R_3OAR60") - R("R_3OAR60_Rev") + R("R_3OAR80") - R("R_3OAR80_Rev") +
 R("R_5DGLCNR") - R("R_5DGLCNR_Rev") - R("R_AGPR") + R("R_AGPR_Rev") -
 R("R_ALDD2y") - R("R_ALDD3y") + R("R_ALR2") + R("R_APRAUR") - R("R_ASAD") +
 R("R_ASAD_Rev") - R("R_ATHRDHr") + R("R_ATHRDHr_Rev") - R("R_BETALDHy") +
 R("R_BSORy") + R("R_DHDPRy") + R("R_DHFR") - R("R_DHFR_Rev") +
 R("R_DKGLCNR1") + R("R_DKGLCNR2y") + R("R_DPR") - R("R_DSERDHr") +
 R("R_DSERDHr_Rev") + R("R_DXPRIi") + R("R_EAR100y") + R("R_EAR120y") +
 R("R_EAR121y") + R("R_EAR140y") + R("R_EAR141y") + R("R_EAR160y") +
 R("R_EAR161y") + R("R_EAR180y") + R("R_EAR181y") + R("R_EAR40y") + R("R_EAR60y")
 + R("R_EAR80y") - 0.000447 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_FADRx2") + R("R_FLDR") + R("R_FLVR") + R("R_FMNRx2") - R("R_G3PD2")
 + R("R_G3PD2_Rev") + R("R_G5SD") - R("R_G6PDH2r") +
 R("R_G6PDH2r_Rev") - R("R_GGGABADr") + R("R_GGGABADr_Rev") - R("R_GLUDy") +
 R("R_GLUDy_Rev") + R("R_GLUSy") + R("R_GLUTRR") + R("R_GLYCLTDy") + R("R_GMPR") -
 R("R_GND") + R("R_GOFUCR") + R("R_GTHOr") - R("R_GTHOr_Rev") +
 R("R_HPYRRy") - R("R_HSDy") + R("R_HSDy_Rev") - R("R_ICDHyr") +
 R("R_ICDHyr_Rev") + R("R_IDOND2") - R("R_KARA1") + R("R_KARA1_Rev") + R("R_KARA2")
 - R("R_KARA2_Rev") - R("R_LSERDHr") + R("R_LSERDHr_Rev") - R("R_ME2")
 - R("R_MTHFD") + R("R_MTHFD_Rev") + R("R_NADK") + R("R_NADPHQR2") +
 R("R_NADPHQR3") + R("R_NADPHQR4") - R("R_NADPPPS") + R("R_NADTRHD") +
 R("R_NODOy") + R("R_P5CR") + R("R_SHK3Dr") - R("R_SHK3Dr_Rev") -
 R("R_SSALy") + 3 R("R_SULRi") + R("R_TDPDRR") - R("R_THD2pp") + R("R_TRDR")
 + R("R_UAPGR")  = 0.

Mnadphc: -R("R_2DGLCNRy") - R("R_2DGULRy") - R("R_3OAR100") +
 R("R_3OAR100_Rev") - R("R_3OAR120") + R("R_3OAR120_Rev") - R("R_3OAR121") -
 R("R_3OAR140") + R("R_3OAR140_Rev") - R("R_3OAR141") - R("R_3OAR160")
 + R("R_3OAR160_Rev") - R("R_3OAR161") - R("R_3OAR180") +
 R("R_3OAR180_Rev") - R("R_3OAR181") - R("R_3OAR40") + R("R_3OAR40_Rev") -
 R("R_3OAR60") + R("R_3OAR60_Rev") - R("R_3OAR80") + R("R_3OAR80_Rev") -
 R("R_5DGLCNR") + R("R_5DGLCNR_Rev") + R("R_AGPR") - R("R_AGPR_Rev") +
 R("R_ALDD2y") + R("R_ALDD3y") - R("R_ALR2") - R("R_APRAUR") + R("R_ASAD") -
 R("R_ASAD_Rev") + R("R_ATHRDHr") - R("R_ATHRDHr_Rev") + R("R_BETALDHy")
 - R("R_BSORy") - R("R_DHDPRy") - R("R_DHFR") + R("R_DHFR_Rev") -
 R("R_DKGLCNR1") - R("R_DKGLCNR2y") - R("R_DPR") + R("R_DSERDHr") -
 R("R_DSERDHr_Rev") - R("R_DXPRIi") - R("R_EAR100y") - R("R_EAR120y") -
 R("R_EAR121y") - R("R_EAR140y") - R("R_EAR141y") - R("R_EAR160y") -
 R("R_EAR161y") - R("R_EAR180y") - R("R_EAR181y") - R("R_EAR40y") -
 R("R_EAR60y") - R("R_EAR80y") - R("R_FADRx2") - R("R_FLDR") - R("R_FLVR") -
 R("R_FMNRx2") + R("R_G3PD2") - R("R_G3PD2_Rev") - R("R_G5SD") +
 R("R_G6PDH2r") - R("R_G6PDH2r_Rev") + R("R_GGGABADr") - R("R_GGGABADr_Rev") +
 R("R_GLUDy") - R("R_GLUDy_Rev") - R("R_GLUSy") - R("R_GLUTRR") -
 R("R_GLYCLTDy") - R("R_GMPR") + R("R_GND") - R("R_GOFUCR") - R("R_GTHOr") +
 R("R_GTHOr_Rev") - R("R_HPYRRy") + R("R_HSDy") - R("R_HSDy_Rev") +
 R("R_ICDHyr") - R("R_ICDHyr_Rev") - R("R_IDOND2") + R("R_KARA1") -
 R("R_KARA1_Rev") - R("R_KARA2") + R("R_KARA2_Rev") + R("R_LSERDHr") -
 R("R_LSERDHr_Rev") + R("R_ME2") + R("R_MTHFD") - R("R_MTHFD_Rev") -
 R("R_NADPHQR2") - R("R_NADPHQR3") - R("R_NADPHQR4") - R("R_NADTRHD") -
 R("R_NODOy") - R("R_P5CR") - R("R_SHK3Dr") + R("R_SHK3Dr_Rev") + R("R_SSALy") -
 3 R("R_SULRi") - R("R_TDPDRR") + R("R_THD2pp") - R("R_TRDR") -
 R("R_UAPGR")  = 0.

Mncamc: R("R_NADN") + R("R_NMNN") + R("R_NMNt7pp") - R("R_NNAM")  =
 0.

Mnh4b: R("Mnh4b_Transp") + R("R_EX_nh4_e_") - R("R_EX_nh4_e__Rev") 
 = 0.

Mnh4c: R("R_AACTOOR") + R("R_ADA") + R("R_ADD") + 2 R("R_ALLTAMH") +
 R("R_ASNN") - R("R_ASNS2") + R("R_ASPT") + R("R_ATPHs") -
 R("R_CBMKr") + R("R_CBMKr_Rev") + R("R_CSND") + R("R_CYNTAH") + R("R_CYSDDS") +
 R("R_CYSDS") + R("R_CYSTL") + R("R_CYTD") + R("R_DAAD") + R("R_DADA") +
 2 R("R_DAPAL") + R("R_DATPHs") + R("R_DCTPD") + R("R_DCYTD") +
 R("R_DHPPDA2") - 0.011843 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_ETHAAL") + R("R_G6PDA") + R("R_GGPTRCO") - R("R_GLNS") + R("R_GLUDy") -
 R("R_GLUDy_Rev") + R("R_GLUN") + R("R_GLYCL") + R("R_GMPR") +
 R("R_GTPHs") + R("R_GUAD") + 4 R("R_HMBS") - R("R_NADS1") + R("R_NH4tpp") -
 R("R_NH4tpp_Rev") + R("R_NMNDA") + R("R_NNAM") + R("R_NTRIR2x") +
 R("R_PYAM5PO") + 2 R("R_SADH") + R("R_SERD_D") + R("R_SERD_L") + R("R_THRD_L")
 + R("R_TRPAS2") - R("R_TRPAS2_Rev") + 2 R("R_UGLYCH")  = 0.

Mnh4e: -R("R_EX_nh4_e_") + R("R_EX_nh4_e__Rev") - R("R_NH4tex") +
 R("R_NH4tex_Rev")  = 0.

Mnh4p: R("R_42A12BOOXpp") + R("R_ASNNpp") + R("R_GLUNpp") +
 R("R_NH4tex") - R("R_NH4tex_Rev") - R("R_NH4tpp") + R("R_NH4tpp_Rev") +
 R("R_NTRIR3pp") + R("R_NTRIR4pp") + R("R_PEAMNOpp") + R("R_TYROXDApp")  =
 0.

Mni2b: R("Mni2b_Transp") + R("R_EX_ni2_e_") - R("R_EX_ni2_e__Rev") 
 = 0.

Mni2c: -R("R_NI2abcpp") - R("R_NI2t3pp") + R("R_NI2tpp") +
 R("R_NI2uabcpp")  = 0.

Mni2e: -R("R_EX_ni2_e_") + R("R_EX_ni2_e__Rev") - R("R_NI2tex") +
 R("R_NI2tex_Rev")  = 0.

Mni2p: R("R_NI2abcpp") + R("R_NI2t3pp") + R("R_NI2tex") -
 R("R_NI2tex_Rev") - R("R_NI2tpp") - R("R_NI2uabcpp")  = 0.

Mnicrntc: R("R_NAMNPP") + R("R_NMNDA") - R("R_NNATr") +
 R("R_NNATr_Rev") - R("R_NNDMBRT") + R("R_NNDPR")  = 0.

Mnmnb: R("Mnmnb_Transp") + R("R_EX_nmn_e_") - R("R_EX_nmn_e__Rev") 
 = 0.

Mnmnc: R("R_NADDP") - R("R_NMNAT") - R("R_NMNDA") - R("R_NMNN") +
 R("R_NMNPtpp")  = 0.

Mnmne: -R("R_EX_nmn_e_") + R("R_EX_nmn_e__Rev") - R("R_NMNtex") +
 R("R_NMNtex_Rev")  = 0.

Mnmnp: -R("R_NMNPtpp") - R("R_NMNt7pp") + R("R_NMNtex") -
 R("R_NMNtex_Rev")  = 0.

MNmtrpc: -R("R_MTRPOX")  = 0.

Mno2b: R("Mno2b_Transp") + R("R_EX_no2_e_") - R("R_EX_no2_e__Rev") 
 = 0.

Mno2c: R("R_NO2t2rpp") - R("R_NO2t2rpp_Rev") + R("R_NO3R1pp") +
 R("R_NO3R2pp") - R("R_NO3t7pp") - R("R_NTRIR2x")  = 0.

Mno2e: -R("R_EX_no2_e_") + R("R_EX_no2_e__Rev") - R("R_NO2tex") +
 R("R_NO2tex_Rev")  = 0.

Mno2p: -R("R_NO2t2rpp") + R("R_NO2t2rpp_Rev") + R("R_NO2tex") -
 R("R_NO2tex_Rev") + R("R_NO3R1bpp") + R("R_NO3R2bpp") + R("R_NO3t7pp") -
 R("R_NTRIR3pp") - R("R_NTRIR4pp")  = 0.

Mno3b: R("Mno3b_Transp") + R("R_EX_no3_e_") - R("R_EX_no3_e__Rev") 
 = 0.

Mno3c: -R("R_NO3R1pp") - R("R_NO3R2pp") + R("R_NO3t7pp") + 2
 R("R_NODOx") + 2 R("R_NODOy")  = 0.

Mno3e: -R("R_EX_no3_e_") + R("R_EX_no3_e__Rev") - R("R_NO3tex") +
 R("R_NO3tex_Rev")  = 0.

Mno3p: -R("R_NO3R1bpp") - R("R_NO3R2bpp") - R("R_NO3t7pp") +
 R("R_NO3tex") - R("R_NO3tex_Rev")  = 0.

Mnob: R("Mnob_Transp") + R("R_EX_no_e_") - R("R_EX_no_e__Rev")  = 0.

Mnoc: -2 R("R_NHFRBO") - 2 R("R_NODOx") - 2 R("R_NODOy") +
 R("R_NOtpp") - R("R_NOtpp_Rev")  = 0.

Mnoe: -R("R_EX_no_e_") + R("R_EX_no_e__Rev") - R("R_NOtex") +
 R("R_NOtex_Rev")  = 0.

Mnop: R("R_NOtex") - R("R_NOtex_Rev") - R("R_NOtpp") +
 R("R_NOtpp_Rev")  = 0.

Mo16a2undp: R("R_O16AP1pp") - R("R_O16AP2pp")  = 0.

Mo16a3undp: R("R_O16AP2pp") - R("R_O16AP3pp")  = 0.

Mo16a4colipab: R("Mo16a4colipab_Transp") + R("R_EX_o16a4colipa_e_")
 - R("R_EX_o16a4colipa_e__Rev")  = 0.

Mo16a4colipae: -R("R_EX_o16a4colipa_e_") +
 R("R_EX_o16a4colipa_e__Rev") + R("R_O16A4COLIPAtex")  = 0.

Mo16a4colipap: -R("R_O16A4COLIPAtex") + R("R_O16A4Lpp")  = 0.

Mo16a4undp: -R("R_O16A4Lpp") + R("R_O16AP3pp")  = 0.

Mo16aundc: -R("R_O16AUNDtpp") + R("R_O16GLCT2")  = 0.

Mo16aundp: -2 R("R_O16AP1pp") - R("R_O16AP2pp") - R("R_O16AP3pp") +
 R("R_O16AUNDtpp")  = 0.

Mo2b: R("Mo2b_Transp") + R("R_EX_o2_e_") - R("R_EX_o2_e__Rev")  = 0.

Mo2c: -R("R_3HCINNMH") - R("R_3HPPPNH") - R("R_AACTOOR") -
 R("R_ASPO6") + R("R_CAT") - R("R_CINNDO") - R("R_CPPPGO") - 0.5 R("R_CYTBD2pp")
 - 0.5 R("R_CYTBDpp") - 0.5 R("R_CYTBO3_4pp") - R("R_DHCINDO") -
 R("R_DKMPPD") - R("R_FDMO") - R("R_FDMO2") - R("R_FDMO3") - R("R_FDMO4") -
 R("R_FDMO6") - R("R_GGPTRCO") - 2 R("R_H2SO") - R("R_HPPPNDO") -
 R("R_MTRPOX") - 2 R("R_NODOx") - 2 R("R_NODOy") + R("R_O2tpp") -
 R("R_O2tpp_Rev") - 0.5 R("R_OMMBLHX") - 0.5 R("R_OMPHHX") - 0.5 R("R_OPHHX") -
 R("R_PDX5POi") - 1.5 R("R_PPPGO") - R("R_PPPNDO") - R("R_PYAM5PO") - 2
 R("R_QMO2") - 2 R("R_QMO3") - R("R_SARCOX") + R("R_SPODM") -
 R("R_TAUDO") - R("R_URIC")  = 0.

Mo2e: -R("R_EX_o2_e_") + R("R_EX_o2_e__Rev") - R("R_O2tex") +
 R("R_O2tex_Rev")  = 0.

Mo2p: -R("R_42A12BOOXpp") - R("R_CU1Opp") - R("R_FEROpp") +
 R("R_O2tex") - R("R_O2tex_Rev") - R("R_O2tpp") + R("R_O2tpp_Rev") -
 R("R_PEAMNOpp") + R("R_SPODMpp") - R("R_TYROXDApp")  = 0.

Mo2sb: R("Mo2sb_Transp") + R("R_EX_o2s_e_") - R("R_EX_o2s_e__Rev") 
 = 0.

Mo2sc: 2 R("R_QMO2") + 2 R("R_QMO3") - 2 R("R_SPODM")  = 0.

Mo2se: -R("R_EX_o2s_e_") + R("R_EX_o2s_e__Rev") - R("R_O2Stex") +
 R("R_O2Stex_Rev")  = 0.

Mo2sp: R("R_O2Stex") - R("R_O2Stex_Rev") - 2 R("R_SPODMpp")  = 0.

Moaac: R("R_ASPTA") - R("R_ASPTA_Rev") + R("R_CITL") - R("R_CS") -
 R("R_MCITS") + R("R_MDH") + R("R_MDH2") + R("R_MDH3") - R("R_MDH_Rev")
 - R("R_OAADC") + R("R_PPC") - R("R_PPCK") + R("R_TARTD")  = 0.

Moc2coac: R("R_ACOAD3f") - R("R_ACOAD3f_Rev") + R("R_ECOAH3") -
 R("R_ECOAH3_Rev")  = 0.

MocACPc: -R("R_3OAS100") + R("R_AACPS9") + R("R_EAR80x") +
 R("R_EAR80y") - R("R_FA80ACPHi")  = 0.

Moccoac: -R("R_ACOAD3f") + R("R_ACOAD3f_Rev") - R("R_FACOAE80") +
 R("R_FACOAL80t2pp") + R("R_KAT4")  = 0.

MocdcaACPc: R("R_AACPS6") - R("R_AGPAT180") + R("R_EAR180x") +
 R("R_EAR180y") - R("R_G3PAT180")  = 0.

Mocdcab: R("Mocdcab_Transp") + R("R_EX_ocdca_e_") -
 R("R_EX_ocdca_e__Rev")  = 0.

Mocdcac: -R("R_2AGPEAT180") - R("R_2AGPGAT180") - R("R_AACPS6") +
 R("R_FACOAE180") + R("R_LPLIPAL2A180") + R("R_LPLIPAL2E180") +
 R("R_LPLIPAL2G180")  = 0.

Mocdcae: -R("R_EX_ocdca_e_") + R("R_EX_ocdca_e__Rev") -
 R("R_OCDCAtexi")  = 0.

Mocdcap: -R("R_FACOAL180t2pp") + R("R_LPLIPAL1A180pp") +
 R("R_LPLIPAL1E180pp") + R("R_LPLIPAL1G180pp") + R("R_OCDCAtexi") +
 R("R_PLIPA1A180pp") + R("R_PLIPA1E180pp") + R("R_PLIPA1G180pp") +
 R("R_PLIPA2A180pp") + R("R_PLIPA2E180pp") + R("R_PLIPA2G180pp")  = 0.

Mocdceab: R("Mocdceab_Transp") + R("R_EX_ocdcea_e_") -
 R("R_EX_ocdcea_e__Rev")  = 0.

Mocdceac: -R("R_2AGPEAT181") - R("R_2AGPGAT181") - R("R_AACPS5") +
 R("R_FACOAE181") + R("R_LPLIPAL2A181") + R("R_LPLIPAL2E181") +
 R("R_LPLIPAL2G181")  = 0.

Mocdceae: -R("R_EX_ocdcea_e_") + R("R_EX_ocdcea_e__Rev") -
 R("R_OCDCEAtexi")  = 0.

Mocdceap: -R("R_FACOAL181t2pp") + R("R_LPLIPAL1A181pp") +
 R("R_LPLIPAL1E181pp") + R("R_LPLIPAL1G181pp") + R("R_OCDCEAtexi") +
 R("R_PLIPA1A181pp") + R("R_PLIPA1E181pp") + R("R_PLIPA1G181pp") +
 R("R_PLIPA2A181pp") + R("R_PLIPA2E181pp") + R("R_PLIPA2G181pp")  = 0.

Moctab: R("Moctab_Transp") + R("R_EX_octa_e_") -
 R("R_EX_octa_e__Rev")  = 0.

Moctac: -R("R_AACPS9") + R("R_FA80ACPHi") + R("R_FACOAE80")  = 0.

Moctae: -R("R_EX_octa_e_") + R("R_EX_octa_e__Rev") - R("R_OCTAtex")
 + R("R_OCTAtex_Rev")  = 0.

Moctap: -R("R_FACOAL80t2pp") + R("R_OCTAtex") - R("R_OCTAtex_Rev") 
 = 0.

Moctdpc: -R("R_DHNAOT4") - R("R_HBZOPT") + R("R_OCTDPS")  = 0.

MocteACPc: R("R_AACPS5") - R("R_AGPAT181") + R("R_EAR181x") +
 R("R_EAR181y") - R("R_G3PAT181")  = 0.

Mod2coac: R("R_ACOAD8f") - R("R_ACOAD8f_Rev") + R("R_ECOAH8") -
 R("R_ECOAH8_Rev") + R("R_ODECOAI")  = 0.

Modecoac: -R("R_FACOAE181") + R("R_FACOAL181t2pp") - R("R_ODECOAI") 
 = 0.

Mohpbc: -R("R_OHPBAT") + R("R_OHPBAT_Rev") + R("R_PERD") -
 R("R_PERD_Rev")  = 0.

Mop4enc: R("R_HKNDDH") + R("R_HKNTDH") - R("R_OP4ENH")  = 0.

Mornb: R("Mornb_Transp") + R("R_EX_orn_e_") - R("R_EX_orn_e__Rev") 
 = 0.

Mornc: R("R_ACODA") - R("R_ARGORNt7pp") + R("R_ARGORNt7pp_Rev") -
 R("R_OCBT") + R("R_OCBT_Rev") + R("R_ORNabcpp") - R("R_ORNDC") -
 R("R_PTRCORNt7pp") + R("R_PTRCORNt7pp_Rev")  = 0.

Morne: -R("R_EX_orn_e_") + R("R_EX_orn_e__Rev") - R("R_ORNtex") +
 R("R_ORNtex_Rev")  = 0.

Mornp: R("R_ARGORNt7pp") - R("R_ARGORNt7pp_Rev") - R("R_ORNabcpp") +
 R("R_ORNtex") - R("R_ORNtex_Rev") + R("R_PTRCORNt7pp") -
 R("R_PTRCORNt7pp_Rev")  = 0.

Morot5pc: -R("R_OMPDC") - R("R_ORPT") + R("R_ORPT_Rev")  = 0.

Morotb: R("Morotb_Transp") + R("R_EX_orot_e_") -
 R("R_EX_orot_e__Rev")  = 0.

Morotc: R("R_DHORD2") + R("R_DHORD5") + R("R_OROTt2_2pp") +
 R("R_ORPT") - R("R_ORPT_Rev")  = 0.

Morote: -R("R_EX_orot_e_") + R("R_EX_orot_e__Rev") - R("R_OROTtex")
 + R("R_OROTtex_Rev")  = 0.

Morotp: -R("R_OROTt2_2pp") + R("R_OROTtex") - R("R_OROTtex_Rev")  =
 0.

Moxamb: R("Moxamb_Transp") + R("R_DM_OXAM")  = 0.

Moxamc: -R("R_DM_OXAM") + R("R_OXAMTC")  = 0.

Moxurc: -R("R_OXAMTC") + R("R_URDGLYCD")  = 0.

Mpa120c: R("R_AGPAT120") + R("R_CDAPPA120") + R("R_DAGK120") -
 R("R_DASYN120") - R("R_PA120abcpp") - R("R_PAPA120")  = 0.

Mpa120p: R("R_CLPNH120pp") + R("R_PA120abcpp") - R("R_PAPA120pp") -
 R("R_PLIPA1A120pp") - R("R_PLIPA2A120pp")  = 0.

Mpa140c: R("R_AGPAT140") + R("R_CDAPPA140") + R("R_DAGK140") -
 R("R_DASYN140") - R("R_PA140abcpp") - R("R_PAPA140")  = 0.

Mpa140p: R("R_CLPNH140pp") + R("R_PA140abcpp") - R("R_PAPA140pp") -
 R("R_PLIPA1A140pp") - R("R_PLIPA2A140pp")  = 0.

Mpa141c: R("R_AGPAT141") + R("R_CDAPPA141") + R("R_DAGK141") -
 R("R_DASYN141") - R("R_PA141abcpp") - R("R_PAPA141")  = 0.

Mpa141p: R("R_CLPNH141pp") + R("R_PA141abcpp") - R("R_PAPA141pp") -
 R("R_PLIPA1A141pp") - R("R_PLIPA2A141pp")  = 0.

Mpa160c: R("R_AGPAT160") + R("R_CDAPPA160") + R("R_DAGK160") -
 R("R_DASYN160") - R("R_PA160abcpp") - R("R_PAPA160")  = 0.

Mpa160p: R("R_CLPNH160pp") + R("R_PA160abcpp") - R("R_PAPA160pp") -
 R("R_PLIPA1A160pp") - R("R_PLIPA2A160pp")  = 0.

Mpa161c: R("R_AGPAT161") + R("R_CDAPPA161") + R("R_DAGK161") -
 R("R_DASYN161") - R("R_PA161abcpp") - R("R_PAPA161")  = 0.

Mpa161p: R("R_CLPNH161pp") + R("R_PA161abcpp") - R("R_PAPA161pp") -
 R("R_PLIPA1A161pp") - R("R_PLIPA2A161pp")  = 0.

Mpa180c: R("R_AGPAT180") + R("R_CDAPPA180") + R("R_DAGK180") -
 R("R_DASYN180") - R("R_PA180abcpp") - R("R_PAPA180")  = 0.

Mpa180p: R("R_CLPNH180pp") + R("R_PA180abcpp") - R("R_PAPA180pp") -
 R("R_PLIPA1A180pp") - R("R_PLIPA2A180pp")  = 0.

Mpa181c: R("R_AGPAT181") + R("R_CDAPPA181") + R("R_DAGK181") -
 R("R_DASYN181") - R("R_PA181abcpp") - R("R_PAPA181")  = 0.

Mpa181p: R("R_CLPNH181pp") + R("R_PA181abcpp") - R("R_PAPA181pp") -
 R("R_PLIPA1A181pp") - R("R_PLIPA2A181pp")  = 0.

Mpacaldb: R("Mpacaldb_Transp") + R("R_EX_pacald_e_") -
 R("R_EX_pacald_e__Rev")  = 0.

Mpacaldc: -R("R_ALDD19x") + R("R_PACALDt2rpp") -
 R("R_PACALDt2rpp_Rev")  = 0.

Mpacalde: -R("R_EX_pacald_e_") + R("R_EX_pacald_e__Rev") -
 R("R_PACALDtex") + R("R_PACALDtex_Rev")  = 0.

Mpacaldp: -R("R_PACALDt2rpp") + R("R_PACALDt2rpp_Rev") +
 R("R_PACALDtex") - R("R_PACALDtex_Rev") + R("R_PEAMNOpp")  = 0.

Mpacc: R("R_ALDD19x") - R("R_PACCOAL")  = 0.

MpalmACPc: -R("R_3OAS180") + R("R_AACPS3") - R("R_AGPAT160") +
 R("R_EAR160x") + R("R_EAR160y") - R("R_FA160ACPHi") - R("R_G3PAT160")  =
 0.

Mpan4pc: R("R_PPCDC") - R("R_PTPATi")  = 0.

MpantRc: R("R_DPR") - R("R_PANTS")  = 0.

Mpapc: R("R_ACPS1") - R("R_BPNT") + R("R_PAPSR") + R("R_PAPSR2")  =
 0.

Mpapsc: R("R_ADSK") - R("R_PAPSR") - R("R_PAPSR2")  = 0.

Mpdx5pc: -R("R_PDX5PO2") - R("R_PDX5POi") + R("R_PDX5PS") -
 R("R_PDXPP") + R("R_PYDXNK")  = 0.

Mpe120c: R("R_2AGPEAT120") - R("R_PE120abcpp") + R("R_PSD120")  = 0.

Mpe120p: R("R_PE120abcpp") - R("R_PLIPA1E120pp") -
 R("R_PLIPA2E120pp")  = 0.

Mpe140c: R("R_2AGPEAT140") - R("R_PE140abcpp") + R("R_PSD140")  = 0.

Mpe140p: R("R_PE140abcpp") - R("R_PLIPA1E140pp") -
 R("R_PLIPA2E140pp")  = 0.

Mpe141c: R("R_2AGPEAT141") - R("R_PE141abcpp") + R("R_PSD141")  = 0.

Mpe141p: R("R_PE141abcpp") - R("R_PLIPA1E141pp") -
 R("R_PLIPA2E141pp")  = 0.

Mpe160c: R("R_2AGPEAT160") - 0.02233
 R("R_Ec_biomass_iAF1260_core_59p81M") - R("R_PE160abcpp") + R("R_PSD160")  = 0.

Mpe160p: -R("R_ALPATE160pp") - 0.04148
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_PE160abcpp") - R("R_PLIPA1E160pp") -
 R("R_PLIPA2E160pp")  = 0.

Mpe161c: R("R_2AGPEAT161") - R("R_CFAS160E") - 0.02632
 R("R_Ec_biomass_iAF1260_core_59p81M") - R("R_PE161abcpp") + R("R_PSD161")  = 0.

Mpe161p: -0.04889 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_PE161abcpp") - R("R_PETNT161pp") - R("R_PLIPA1E161pp") -
 R("R_PLIPA2E161pp")  = 0.

Mpe180c: R("R_2AGPEAT180") - R("R_PE180abcpp") + R("R_PSD180")  = 0.

Mpe180p: R("R_PE180abcpp") - R("R_PLIPA1E180pp") -
 R("R_PLIPA2E180pp")  = 0.

Mpe181c: R("R_2AGPEAT181") - R("R_CFAS180E") - R("R_PE181abcpp") +
 R("R_PSD181")  = 0.

Mpe181p: R("R_PE181abcpp") - R("R_PETNT181pp") - R("R_PLIPA1E181pp")
 - R("R_PLIPA2E181pp")  = 0.

Mpeamnb: R("Mpeamnb_Transp") + R("R_EX_peamn_e_") -
 R("R_EX_peamn_e__Rev")  = 0.

Mpeamne: -R("R_EX_peamn_e_") + R("R_EX_peamn_e__Rev") -
 R("R_PEAMNtex") + R("R_PEAMNtex_Rev")  = 0.

Mpeamnp: -R("R_PEAMNOpp") + R("R_PEAMNtex") - R("R_PEAMNtex_Rev")  =
 0.

Mpepc: -R("R_ACGAptspp") - R("R_ACMANAptspp") - R("R_ACMUMptspp") -
 R("R_ASCBptspp") - R("R_DDPA") - R("R_DHAPT") + R("R_ENO") -
 R("R_ENO_Rev") - R("R_FRUpts2pp") - R("R_FRUptspp") - R("R_GALTptspp") -
 R("R_GAMptspp") - R("R_GLCptspp") - R("R_KDOPS") - R("R_MALTptspp") -
 R("R_MANGLYCptspp") - R("R_MANptspp") - R("R_MNLptspp") - R("R_PPC") +
 R("R_PPCK") + R("R_PPS") - R("R_PSCVT") + R("R_PSCVT_Rev") - R("R_PYK") -
 R("R_SBTptspp") - R("R_SUCptspp") - R("R_TREptspp") - R("R_UAGCVT")  =
 0.

Mpg120c: R("R_2AGPGAT120") - R("R_LPLIPAL2ATE120") -
 R("R_LPLIPAL2ATG120") - R("R_PG120abcpp") + R("R_PGPP120")  = 0.

Mpg120p: R("R_CLPNH120pp") - 2 R("R_CLPNS120pp") + 2
 R("R_CLPNS120pp_Rev") + R("R_PG120abcpp") + R("R_PGPP120pp") - R("R_PLIPA1G120pp") -
 R("R_PLIPA2G120pp")  = 0.

Mpg140c: R("R_2AGPGAT140") - R("R_LPLIPAL2ATE140") -
 R("R_LPLIPAL2ATG140") - R("R_PG140abcpp") + R("R_PGPP140")  = 0.

Mpg140p: R("R_CLPNH140pp") - 2 R("R_CLPNS140pp") + 2
 R("R_CLPNS140pp_Rev") + R("R_PG140abcpp") + R("R_PGPP140pp") - R("R_PLIPA1G140pp") -
 R("R_PLIPA2G140pp")  = 0.

Mpg141c: R("R_2AGPGAT141") - R("R_LPLIPAL2ATE141") -
 R("R_LPLIPAL2ATG141") - R("R_PG141abcpp") + R("R_PGPP141")  = 0.

Mpg141p: R("R_CLPNH141pp") - 2 R("R_CLPNS141pp") + 2
 R("R_CLPNS141pp_Rev") + R("R_PG141abcpp") + R("R_PGPP141pp") - R("R_PLIPA1G141pp") -
 R("R_PLIPA2G141pp")  = 0.

Mpg160c: R("R_2AGPGAT160") - R("R_LPLIPAL2ATE160") -
 R("R_LPLIPAL2ATG160") - R("R_PG160abcpp") + R("R_PGPP160")  = 0.

Mpg160p: -R("R_ALPATG160pp") + R("R_CLPNH160pp") - 2
 R("R_CLPNS160pp") + 2 R("R_CLPNS160pp_Rev") + R("R_PG160abcpp") + R("R_PGPP160pp") -
 R("R_PLIPA1G160pp") - R("R_PLIPA2G160pp")  = 0.

Mpg161c: R("R_2AGPGAT161") - R("R_CFAS160G") - R("R_LPLIPAL2ATE161")
 - R("R_LPLIPAL2ATG161") - R("R_PG161abcpp") + R("R_PGPP161")  = 0.

Mpg161p: R("R_CLPNH161pp") - 2 R("R_CLPNS161pp") + 2
 R("R_CLPNS161pp_Rev") + R("R_PG161abcpp") + R("R_PGPP161pp") - R("R_PLIPA1G161pp") -
 R("R_PLIPA2G161pp")  = 0.

Mpg180c: R("R_2AGPGAT180") - R("R_LPLIPAL2ATE180") -
 R("R_LPLIPAL2ATG180") - R("R_PG180abcpp") + R("R_PGPP180")  = 0.

Mpg180p: R("R_CLPNH180pp") - 2 R("R_CLPNS180pp") + 2
 R("R_CLPNS180pp_Rev") + R("R_PG180abcpp") + R("R_PGPP180pp") - R("R_PLIPA1G180pp") -
 R("R_PLIPA2G180pp")  = 0.

Mpg181c: R("R_2AGPGAT181") - R("R_CFAS180G") - R("R_LPLIPAL2ATE181")
 - R("R_LPLIPAL2ATG181") - R("R_PG181abcpp") + R("R_PGPP181")  = 0.

Mpg181p: R("R_CLPNH181pp") - 2 R("R_CLPNS181pp") + 2
 R("R_CLPNS181pp_Rev") + R("R_PG181abcpp") + R("R_PGPP181pp") - R("R_PLIPA1G181pp") -
 R("R_PLIPA2G181pp")  = 0.

Mpgp120c: -R("R_PGP120abcpp") - R("R_PGPP120") + R("R_PGSA120")  =
 0.

Mpgp120p: R("R_PGP120abcpp") - R("R_PGPP120pp")  = 0.

Mpgp140c: -R("R_PGP140abcpp") - R("R_PGPP140") + R("R_PGSA140")  =
 0.

Mpgp140p: R("R_PGP140abcpp") - R("R_PGPP140pp")  = 0.

Mpgp141c: -R("R_PGP141abcpp") - R("R_PGPP141") + R("R_PGSA141")  =
 0.

Mpgp141p: R("R_PGP141abcpp") - R("R_PGPP141pp")  = 0.

Mpgp160c: -R("R_PGP160abcpp") - R("R_PGPP160") + R("R_PGSA160")  =
 0.

Mpgp160p: R("R_PGP160abcpp") - R("R_PGPP160pp")  = 0.

Mpgp161c: -R("R_PGP161abcpp") - R("R_PGPP161") + R("R_PGSA161")  =
 0.

Mpgp161p: R("R_PGP161abcpp") - R("R_PGPP161pp")  = 0.

Mpgp180c: -R("R_PGP180abcpp") - R("R_PGPP180") + R("R_PGSA180")  =
 0.

Mpgp180p: R("R_PGP180abcpp") - R("R_PGPP180pp")  = 0.

Mpgp181c: -R("R_PGP181abcpp") - R("R_PGPP181") + R("R_PGSA181")  =
 0.

Mpgp181p: R("R_PGP181abcpp") - R("R_PGPP181pp")  = 0.

Mphaccoac: R("R_PACCOAL")  = 0.

MpheLb: R("MpheLb_Transp") + R("R_EX_phe_L_e_") -
 R("R_EX_phe_L_e__Rev")  = 0.

MpheLc: -0.1759 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_PHEt2rpp") - R("R_PHEt2rpp_Rev") - R("R_PHETA1") + R("R_PHETA1_Rev") -
 R("R_PHETRS")  = 0.

MpheLe: -R("R_EX_phe_L_e_") + R("R_EX_phe_L_e__Rev") - R("R_PHEtex")
 + R("R_PHEtex_Rev")  = 0.

MpheLp: -R("R_PHEt2rpp") + R("R_PHEt2rpp_Rev") + R("R_PHEtex") -
 R("R_PHEtex_Rev")  = 0.

Mphemeb: R("Mphemeb_Transp") + R("R_EX_pheme_e_") -
 R("R_EX_pheme_e__Rev")  = 0.

Mphemec: -0.000223 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_FCLT") - R("R_HEMEOS") - R("R_PHEMEabcpp")  = 0.

Mphemee: -R("R_EX_pheme_e_") + R("R_EX_pheme_e__Rev") +
 R("R_PHEMEtiex")  = 0.

Mphemep: R("R_PHEMEabcpp") - R("R_PHEMEtiex")  = 0.

Mphetrnac: R("R_PHETRS")  = 0.

Mphhlipac: R("R_HEPK1") - R("R_HEPT3")  = 0.

Mphomc: R("R_HSK") - R("R_THRS")  = 0.

Mphphhlipac: R("R_HEPK2") - R("R_MOAT3C")  = 0.

Mphpyrc: R("R_PHETA1") - R("R_PHETA1_Rev") + R("R_PPNDH")  = 0.

Mphthrc: -R("R_4HTHRS") + R("R_OHPBAT") - R("R_OHPBAT_Rev") -
 R("R_PDX5PS")  = 0.

Mpib: R("Mpib_Transp") + R("R_EX_pi_e_") - R("R_EX_pi_e__Rev")  = 0.

Mpic: R("R_14GLUCANabcpp") + R("R_2MAHMP") + R("R_3PEPTabcpp") +
 R("R_4HTHRS") + R("R_4PEPTabcpp") + R("R_ACCOAC") + R("R_ACCOAL") +
 R("R_ADOCBLabcpp") + R("R_ADSS") - R("R_AGPR") + R("R_AGPR_Rev") +
 R("R_AIRC2") + R("R_ALAabcpp") + R("R_ALAALAabcpp") + R("R_ALAALAr") -
 R("R_ALAALAr_Rev") + R("R_ALLabcpp") + R("R_ARBabcpp") + R("R_ARBTNabcpp") +
 R("R_ARGabcpp") - R("R_ASAD") + R("R_ASAD_Rev") + R("R_ASNabcpp") +
 R("R_ASO3t8pp") + R("R_ASPabcpp") + R("R_ASPCT") + R("R_ATPM") -
 R("R_ATPS4rpp") + R("R_ATPS4rpp_Rev") + R("R_BPNT") + R("R_BUTSO3abcpp") +
 R("R_CBIuabcpp") + R("R_CBL1abcpp") + R("R_CBPS") + R("R_CD2abcpp") +
 R("R_CGLYabcpp") + R("R_CHLabcpp") + R("R_CHORS") + R("R_COBALT2abcpp")
 + R("R_COLIPAabcpp") + R("R_CPGNabcpp") + R("R_CRNabcpp") +
 R("R_CRNCAL2") + R("R_CRNDabcpp") + R("R_CRNDCAL2") + R("R_CTBTabcpp") +
 R("R_CTBTCAL2") + R("R_CTPS2") + R("R_CU1abcpp") + R("R_CU2abcpp") +
 R("R_CYSabc2pp") + R("R_CYSabcpp") + R("R_CYSDabcpp") + R("R_DAPabcpp") +
 R("R_DBTS") + R("R_DDPA") + R("R_DHFS") + R("R_DHQS") + R("R_DKMPPD") +
 R("R_DKMPPD3") + R("R_DNMPPA") - R("R_DURIPP") + R("R_DURIPP_Rev") +
 59.806 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_ETHSO3abcpp") +
 R("R_F6PP") - 2 R("R_F6Pt6_2pp") + R("R_FBP") + R("R_FE2abcpp") +
 R("R_FE3abcpp") + R("R_FE3DCITabcpp") + R("R_FE3DHBZSabcpp") +
 R("R_FE3HOXabcpp") + R("R_FECRMabcpp") + R("R_FEENTERabcpp") + R("R_FEOXAMabcpp") +
 R("R_G3PCabcpp") + R("R_G3PEabcpp") + R("R_G3PGabcpp") +
 R("R_G3PIabcpp") + R("R_G3PSabcpp") + R("R_G3PT") + R("R_G5SD") + R("R_G6PP") - 2
 R("R_G6Pt6_2pp") + R("R_GALabcpp") - 2 R("R_GAM6Pt6_2pp") - R("R_GAPD") +
 R("R_GAPD_Rev") + R("R_GART") + R("R_GGPTRCS") + R("R_GLCabcpp") -
 R("R_GLCP") - R("R_GLCP2") + R("R_GLNabcpp") + R("R_GLNS") +
 R("R_GLUabcpp") + R("R_GLUCYS") + R("R_GLYBabcpp") + R("R_GLYC3Pabcpp") -
 R("R_GLYC3Pt6pp") + R("R_GMHEPPA") + R("R_GSPMDS") + R("R_GTHRDabc2pp") +
 R("R_GTHRDabcpp") + R("R_GTHS") + R("R_GTPDPDP") + R("R_HG2abcpp") +
 R("R_HISabcpp") + R("R_HISTP") + R("R_HYPOE") + R("R_ILEabcpp") +
 R("R_ISETACabcpp") + R("R_K2L4Aabcpp") + R("R_Kabcpp") + R("R_KDOPP") +
 R("R_KDOPS") + R("R_LEUabcpp") + R("R_LIPAabcpp") + R("R_LIPACabcpp") +
 R("R_LYSabcpp") + R("R_MALTabcpp") + R("R_MALTHXabcpp") +
 R("R_MALTPTabcpp") + R("R_MALTTRabcpp") + R("R_MALTTTRabcpp") + R("R_MAN1PT2") - 2
 R("R_MAN6Pt6_2pp") + R("R_METabcpp") + R("R_METAT") + R("R_METDabcpp") +
 R("R_METSOX1abcpp") + R("R_METSOX2abcpp") + R("R_MG2uabcpp") +
 R("R_MGSA") + R("R_MI1PP") - R("R_MLTP1") + R("R_MLTP1_Rev") - R("R_MLTP2") +
 R("R_MLTP2_Rev") - R("R_MLTP3") + R("R_MLTP3_Rev") + R("R_MN6PP") +
 R("R_MOBDabcpp") + R("R_MSO3abcpp") + R("R_NADPPPS") + R("R_NAMNPP") +
 R("R_NI2abcpp") + R("R_NI2uabcpp") + R("R_NTD1") + R("R_NTD10") +
 R("R_NTD11") + R("R_NTD12") + R("R_NTD2") + R("R_NTD3") + R("R_NTD4") +
 R("R_NTD5") + R("R_NTD6") + R("R_NTD7") + R("R_NTD8") + R("R_NTD9") +
 R("R_NTP1") + R("R_NTP10") + R("R_NTP11") + R("R_NTP12") + R("R_NTP3") +
 R("R_NTP5") + R("R_OCBT") - R("R_OCBT_Rev") + 2 R("R_OMMBLHX3") + 2
 R("R_OMPHHX3") + 2 R("R_OPHHX3") + R("R_ORNabcpp") - R("R_OXAMTC") +
 R("R_PA120abcpp") + R("R_PA140abcpp") + R("R_PA141abcpp") +
 R("R_PA160abcpp") + R("R_PA161abcpp") + R("R_PA180abcpp") + R("R_PA181abcpp") +
 R("R_PAPA120") + R("R_PAPA140") + R("R_PAPA141") + R("R_PAPA160") +
 R("R_PAPA161") + R("R_PAPA180") + R("R_PAPA181") + R("R_PDX5PS") +
 R("R_PDXPP") + R("R_PE120abcpp") + R("R_PE140abcpp") + R("R_PE141abcpp") +
 R("R_PE160abcpp") + R("R_PE161abcpp") + R("R_PE180abcpp") +
 R("R_PE181abcpp") + R("R_PG120abcpp") + R("R_PG140abcpp") + R("R_PG141abcpp") +
 R("R_PG160abcpp") + R("R_PG161abcpp") + R("R_PG180abcpp") +
 R("R_PG181abcpp") + R("R_PGLYCP") + R("R_PGP120abcpp") + R("R_PGP140abcpp") +
 R("R_PGP141abcpp") + R("R_PGP160abcpp") + R("R_PGP161abcpp") +
 R("R_PGP180abcpp") + R("R_PGP181abcpp") + R("R_PGPP120") + R("R_PGPP140") +
 R("R_PGPP141") + R("R_PGPP160") + R("R_PGPP161") + R("R_PGPP180") +
 R("R_PGPP181") + R("R_PHEMEabcpp") + R("R_PIt2rpp") - R("R_PIt2rpp_Rev") + 2
 R("R_PIuabcpp") + R("R_PMDPHT") + 2 R("R_PPA") + R("R_PPA2") +
 R("R_PPC") - R("R_PPKr") + R("R_PPKr_Rev") + R("R_PPS") + R("R_PRAGSr") -
 R("R_PRAGSr_Rev") + R("R_PRAIS") + R("R_PRASCSi") + R("R_PRFGS") +
 R("R_PROabcpp") + R("R_PROGLYabcpp") + R("R_PSCVT") - R("R_PSCVT_Rev") +
 R("R_PSP_L") - R("R_PTA2") - R("R_PTAr") + R("R_PTAr_Rev") +
 R("R_PTRCabcpp") - R("R_PUNP1") + R("R_PUNP1_Rev") - R("R_PUNP2") +
 R("R_PUNP2_Rev") - R("R_PUNP3") + R("R_PUNP3_Rev") - R("R_PUNP4") + R("R_PUNP4_Rev")
 - R("R_PUNP5") + R("R_PUNP5_Rev") - R("R_PUNP6") + R("R_PUNP6_Rev") -
 R("R_PUNP7") + R("R_PUNP7_Rev") + R("R_PYDXPP") - R("R_PYNP2r") +
 R("R_PYNP2r_Rev") + R("R_QULNS") + R("R_R5PP") + R("R_RBFSa") +
 R("R_RIBabcpp") + R("R_RZ5PP") + R("R_SADT2") + R("R_SELCYSS") + R("R_SELNPS")
 + R("R_SPMDabcpp") + R("R_SUCOAS") - R("R_SUCOAS_Rev") +
 R("R_SULabcpp") + R("R_SULFACabcpp") + R("R_TAURabcpp") + R("R_TDP") +
 R("R_THMabcpp") + R("R_THRabcpp") + R("R_THRS") - R("R_TMDPP") + R("R_TMDPP_Rev")
 + R("R_TRE6PP") + R("R_TSULabcpp") + R("R_TUNGSabcpp") + R("R_UAAGDS")
 + R("R_UAGCVT") + R("R_UAMAGS") + R("R_UAMAS") + R("R_UDCPDP") +
 R("R_UGMDDS") + R("R_UM3PL") + R("R_UM4PL") + R("R_VALabcpp") +
 R("R_XYLabcpp") + R("R_ZN2abcpp") + R("R_ZNabcpp")  = 0.

Mpie: -R("R_EX_pi_e_") + R("R_EX_pi_e__Rev") - R("R_PItex") +
 R("R_PItex_Rev")  = 0.

Mpip: R("R_3NTD2pp") + R("R_3NTD4pp") + R("R_3NTD7pp") +
 R("R_3NTD9pp") + R("R_ACGAL1PPpp") + R("R_ACGAM1PPpp") + 2 R("R_F6Pt6_2pp") +
 R("R_G1PPpp") + R("R_G2PPpp") + 2 R("R_G6Pt6_2pp") + R("R_GAL1PPpp") + 2
 R("R_GAM6Pt6_2pp") + R("R_GLYC3Pt6pp") + R("R_GUR1PPpp") + 2
 R("R_MAN6Pt6_2pp") + R("R_NTD10pp") + R("R_NTD11pp") + R("R_NTD12pp") +
 R("R_NTD1pp") + R("R_NTD2pp") + R("R_NTD3pp") + R("R_NTD4pp") + R("R_NTD5pp")
 + R("R_NTD6pp") + R("R_NTD7pp") + R("R_NTD8pp") + R("R_NTD9pp") +
 R("R_NTP3pp") + R("R_PAPA120pp") + R("R_PAPA140pp") + R("R_PAPA141pp") +
 R("R_PAPA160pp") + R("R_PAPA161pp") + R("R_PAPA180pp") +
 R("R_PAPA181pp") + R("R_PGPP120pp") + R("R_PGPP140pp") + R("R_PGPP141pp") +
 R("R_PGPP160pp") + R("R_PGPP161pp") + R("R_PGPP180pp") + R("R_PGPP181pp") + 6
 R("R_PHYTSpp") - R("R_PIt2rpp") + R("R_PIt2rpp_Rev") + R("R_PItex") -
 R("R_PItex_Rev") - R("R_PIuabcpp") + R("R_PPTHpp") + R("R_PSP_Lpp") +
 R("R_PTHRpp") + R("R_R5PPpp") + R("R_TYRPpp") + R("R_UDCPDPpp")  =
 0.

Mpmcoac: -R("R_AOXSr") + R("R_AOXSr_Rev")  = 0.

Mpmtcoac: -R("R_ACOAD7f") + R("R_ACOAD7f_Rev") - R("R_FACOAE160") +
 R("R_FACOAL160t2pp") + R("R_KAT8")  = 0.

MpntoRb: R("MpntoRb_Transp") + R("R_EX_pnto_R_e_") -
 R("R_EX_pnto_R_e__Rev")  = 0.

MpntoRc: R("R_PANTS") - R("R_PNTK") + R("R_PNTOt4pp")  = 0.

MpntoRe: -R("R_EX_pnto_R_e_") + R("R_EX_pnto_R_e__Rev") -
 R("R_PNTOtex") + R("R_PNTOtex_Rev")  = 0.

MpntoRp: -R("R_PNTOt4pp") + R("R_PNTOtex") - R("R_PNTOtex_Rev")  =
 0.

Mppab: R("Mppab_Transp") + R("R_EX_ppa_e_") - R("R_EX_ppa_e__Rev") 
 = 0.

Mppac: -R("R_ACCOAL") + R("R_ALDD3y") + R("R_PPAKr") -
 R("R_PPAKr_Rev") + R("R_PPAt4pp") + R("R_PPCSCT")  = 0.

Mppae: -R("R_EX_ppa_e_") + R("R_EX_ppa_e__Rev") - R("R_PPAtex") +
 R("R_PPAtex_Rev")  = 0.

Mppalb: R("Mppalb_Transp") + R("R_EX_ppal_e_") -
 R("R_EX_ppal_e__Rev")  = 0.

Mppalc: -R("R_ALDD3y") - R("R_PPALtpp") + R("R_PPALtpp_Rev")  = 0.

Mppale: -R("R_EX_ppal_e_") + R("R_EX_ppal_e__Rev") - R("R_PPALtex")
 + R("R_PPALtex_Rev")  = 0.

Mppalp: R("R_PPALtex") - R("R_PPALtex_Rev") + R("R_PPALtpp") -
 R("R_PPALtpp_Rev")  = 0.

Mppap: -R("R_PPAt4pp") + R("R_PPAtex") - R("R_PPAtex_Rev")  = 0.

Mppapc: -R("R_PPAKr") + R("R_PPAKr_Rev") + R("R_PTA2")  = 0.

Mppbngc: -4 R("R_HMBS") + R("R_PPBNGS")  = 0.

Mppcoac: R("R_ACCOAL") - R("R_MCITS") + R("R_MMCD") + R("R_OBTFL") -
 R("R_PPCSCT") - R("R_PTA2")  = 0.

Mppgppc: R("R_GDPDPK") + R("R_GTPDPDP") - R("R_PPGPPDP")  = 0.

Mpphnc: R("R_CHORM") - R("R_PPND") - R("R_PPNDH")  = 0.

Mppic: R("R_2AGPEAT120") + R("R_2AGPEAT140") + R("R_2AGPEAT141") +
 R("R_2AGPEAT160") + R("R_2AGPEAT161") + R("R_2AGPEAT180") +
 R("R_2AGPEAT181") + R("R_2AGPGAT120") + R("R_2AGPGAT140") + R("R_2AGPGAT141") +
 R("R_2AGPGAT160") + R("R_2AGPGAT161") + R("R_2AGPGAT180") +
 R("R_2AGPGAT181") + R("R_AACPS1") + R("R_AACPS2") + R("R_AACPS3") + R("R_AACPS4")
 + R("R_AACPS5") + R("R_AACPS6") + R("R_AACPS7") + R("R_AACPS8") +
 R("R_AACPS9") + R("R_ACBIPGT") + R("R_ACS") + R("R_ADNCYC") + R("R_ADPT")
 + R("R_ALATRS") + R("R_ANPRT") + R("R_ARGSS") + R("R_ARGTRS") +
 R("R_ASNS1") + R("R_ASNS2") + R("R_ASNTRS") + R("R_ASPTRS") + R("R_ATPPRT")
 + R("R_CYSTRS") + R("R_DASYN120") + R("R_DASYN140") + R("R_DASYN141")
 + R("R_DASYN160") + R("R_DASYN161") + R("R_DASYN180") +
 R("R_DASYN181") + R("R_DHBS") + R("R_DHNAOT4") + R("R_DHPS2") + R("R_DMATT") +
 R("R_DNTPPA") + R("R_DUTPDP") + 0.7739
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_FACOAL100t2pp") + R("R_FACOAL120t2pp") +
 R("R_FACOAL140t2pp") + R("R_FACOAL141t2pp") + R("R_FACOAL160t2pp") + R("R_FACOAL161t2pp")
 + R("R_FACOAL180t2pp") + R("R_FACOAL181t2pp") + R("R_FACOAL60t2pp") +
 R("R_FACOAL80t2pp") + R("R_FMNAT") + R("R_G1PTT") + R("R_GALUi") +
 R("R_GLGC") + R("R_GLNTRS") + R("R_GLUPRT") + R("R_GLUTRS") +
 R("R_GLYTRS") + R("R_GMHEPAT") + R("R_GMPS2") + R("R_GRTT") + R("R_GTPCII2") +
 R("R_GUAPRT") + R("R_HBZOPT") + R("R_HEMEOS") + R("R_HISTRS") +
 R("R_HXPRT") + R("R_ILETRS") + R("R_KDOCT2") + R("R_LEUTRS") + R("R_LYSTRS")
 + R("R_MEPCT") + R("R_METAT") + R("R_METTRS") + R("R_NADS1") +
 R("R_NAMNPP") + R("R_NMNAT") + R("R_NNATr") - R("R_NNATr_Rev") + R("R_NNDPR")
 + R("R_NTPP1") + R("R_NTPP10") + R("R_NTPP11") + R("R_NTPP2") +
 R("R_NTPP3") + R("R_NTPP4") + R("R_NTPP5") + R("R_NTPP6") + R("R_NTPP7") +
 R("R_NTPP8") + R("R_NTPP9") + 5 R("R_OCTDPS") - R("R_ORPT") +
 R("R_ORPT_Rev") + R("R_PACCOAL") + R("R_PANTS") + R("R_PHETRS") - R("R_PPA") +
 R("R_PPA2") + R("R_PPGPPDP") - R("R_PPK2r") + R("R_PPK2r_Rev") +
 R("R_PPKr") - R("R_PPKr_Rev") + R("R_PPNCL2") + R("R_PRATPP") +
 R("R_PROTRS") + R("R_PTPATi") + R("R_SADT2") + R("R_SERASr") - R("R_SERASr_Rev")
 + R("R_SERTRS") + R("R_SERTRS2") + R("R_SUCBZL") + R("R_THRTRS") +
 R("R_THZPSN") + R("R_TMPPP") + R("R_TRPTRS") + R("R_TYRTRS") +
 R("R_UAGDP") + 8 R("R_UDCPDPS") + R("R_UPPRT") + R("R_VALTRS") + R("R_XPPT")  =
 0.

Mppp9c: -R("R_FCLT") + R("R_PPPGO") + R("R_PPPGO3")  = 0.

Mpppg9c: R("R_CPPPGO") + R("R_CPPPGO2") - R("R_PPPGO") -
 R("R_PPPGO3")  = 0.

Mpppic: R("R_CBIAT") - R("R_CBIAT_Rev") + R("R_CBLAT") -
 R("R_CBLAT_Rev") + R("R_NTPTP1") + R("R_NTPTP2") - R("R_PPA2") + R("R_PPK2r") -
 R("R_PPK2r_Rev")  = 0.

Mpppnb: R("Mpppnb_Transp") + R("R_EX_pppn_e_") -
 R("R_EX_pppn_e__Rev")  = 0.

Mpppnc: -R("R_PPPNDO") + R("R_PPPNt2rpp") - R("R_PPPNt2rpp_Rev")  =
 0.

Mpppne: -R("R_EX_pppn_e_") + R("R_EX_pppn_e__Rev") - R("R_PPPNtex")
 + R("R_PPPNtex_Rev")  = 0.

Mpppnp: -R("R_PPPNt2rpp") + R("R_PPPNt2rpp_Rev") + R("R_PPPNtex") -
 R("R_PPPNtex_Rev")  = 0.

Mpptb: R("Mpptb_Transp") + R("R_EX_ppt_e_") - R("R_EX_ppt_e__Rev") 
 = 0.

Mppte: -R("R_EX_ppt_e_") + R("R_EX_ppt_e__Rev") - R("R_PPTtex") +
 R("R_PPTtex_Rev")  = 0.

Mpptp: -R("R_PPTHpp") + R("R_PPTtex") - R("R_PPTtex_Rev")  = 0.

Mpramc: R("R_GLUPRT") - R("R_PRAGSr") + R("R_PRAGSr_Rev")  = 0.

Mpranc: R("R_ANPRT") - R("R_PRAIi")  = 0.

Mprbampc: -R("R_PRAMPC") + R("R_PRATPP")  = 0.

Mprbatpc: R("R_ATPPRT") - R("R_PRATPP")  = 0.

Mprfpc: R("R_PRAMPC") - R("R_PRMICI") + R("R_PRMICI_Rev")  = 0.

Mprlpc: -R("R_IG3PS") + R("R_PRMICI") - R("R_PRMICI_Rev")  = 0.

Mproglyb: R("Mproglyb_Transp") + R("R_EX_progly_e_") -
 R("R_EX_progly_e__Rev")  = 0.

Mproglyc: -R("R_AMPTASEPG") + R("R_PROGLYabcpp")  = 0.

Mproglye: -R("R_EX_progly_e_") + R("R_EX_progly_e__Rev") -
 R("R_PROGLYtex") + R("R_PROGLYtex_Rev")  = 0.

Mproglyp: -R("R_PROGLYabcpp") + R("R_PROGLYtex") -
 R("R_PROGLYtex_Rev")  = 0.

MproLb: R("MproLb_Transp") + R("R_EX_pro_L_e_") -
 R("R_EX_pro_L_e__Rev")  = 0.

MproLc: R("R_AMPTASEPG") - 0.2211
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_P5CR") + R("R_PROabcpp") - R("R_PROD2") + R("R_PROt2rpp") -
 R("R_PROt2rpp_Rev") + R("R_PROt4pp") - R("R_PROTRS")  = 0.

MproLe: -R("R_EX_pro_L_e_") + R("R_EX_pro_L_e__Rev") - R("R_PROtex")
 + R("R_PROtex_Rev")  = 0.

MproLp: -R("R_PROabcpp") - R("R_PROt2rpp") + R("R_PROt2rpp_Rev") -
 R("R_PROt4pp") + R("R_PROtex") - R("R_PROtex_Rev")  = 0.

Mprotrnac: R("R_PROTRS")  = 0.

Mprppc: -R("R_ADPT") - R("R_ANPRT") - R("R_ATPPRT") - R("R_GLUPRT")
 - R("R_GUAPRT") - R("R_HXPRT") - R("R_NAMNPP") - R("R_NNDPR") +
 R("R_ORPT") - R("R_ORPT_Rev") + R("R_PRPPS") - R("R_PRPPS_Rev") +
 R("R_R15BPK") - R("R_UPPRT") - R("R_XPPT")  = 0.

Mps120c: -R("R_PSD120") + R("R_PSSA120")  = 0.

Mps140c: -R("R_PSD140") + R("R_PSSA140")  = 0.

Mps141c: -R("R_PSD141") + R("R_PSSA141")  = 0.

Mps160c: -R("R_PSD160") + R("R_PSSA160")  = 0.

Mps161c: -R("R_PSD161") + R("R_PSSA161")  = 0.

Mps180c: -R("R_PSD180") + R("R_PSSA180")  = 0.

Mps181c: -R("R_PSD181") + R("R_PSSA181")  = 0.

Mpsclysb: R("Mpsclysb_Transp") + R("R_EX_psclys_e_") -
 R("R_EX_psclys_e__Rev")  = 0.

Mpsclysc: -R("R_FRULYSE") + R("R_FRULYSE_Rev") + R("R_PSCLYSt2pp") 
 = 0.

Mpsclyse: -R("R_EX_psclys_e_") + R("R_EX_psclys_e__Rev") -
 R("R_PSCLYStex") + R("R_PSCLYStex_Rev")  = 0.

Mpsclysp: -R("R_PSCLYSt2pp") + R("R_PSCLYStex") -
 R("R_PSCLYStex_Rev")  = 0.

MpserLb: R("MpserLb_Transp") + R("R_EX_pser_L_e_") -
 R("R_EX_pser_L_e__Rev")  = 0.

MpserLc: R("R_PSERT") - R("R_PSP_L")  = 0.

MpserLe: -R("R_EX_pser_L_e_") + R("R_EX_pser_L_e__Rev") -
 R("R_PSERtex") + R("R_PSERtex_Rev")  = 0.

MpserLp: R("R_PSERtex") - R("R_PSERtex_Rev") - R("R_PSP_Lpp")  = 0.

Mptrcb: R("Mptrcb_Transp") + R("R_EX_ptrc_e_") -
 R("R_EX_ptrc_e__Rev")  = 0.

Mptrcc: R("R_AGMT") - R("R_GGPTRCS") + R("R_ORNDC") +
 R("R_PTRCabcpp") + R("R_PTRCORNt7pp") - R("R_PTRCORNt7pp_Rev") + R("R_PTRCt2pp") -
 R("R_PTRCTA") - R("R_SPMS")  = 0.

Mptrce: -R("R_EX_ptrc_e_") + R("R_EX_ptrc_e__Rev") - R("R_PTRCtex")
 + R("R_PTRCtex_Rev")  = 0.

Mptrcp: -R("R_PTRCabcpp") - R("R_PTRCORNt7pp") +
 R("R_PTRCORNt7pp_Rev") - R("R_PTRCt2pp") + R("R_PTRCtex") - R("R_PTRCtex_Rev")  = 0.

Mpyam5pc: R("R_ALATA_D2") + R("R_ALATA_L2") - R("R_HYPOE") -
 R("R_PYAM5PO") + R("R_PYDAMK")  = 0.

Mpydamc: R("R_HYPOE") - R("R_PYDAMK")  = 0.

Mpydx5pc: -R("R_ALATA_D2") - R("R_ALATA_L2") - 0.000223
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_PDX5PO2") + R("R_PDX5POi") +
 R("R_PYAM5PO") + R("R_PYDXK") - R("R_PYDXPP")  = 0.

Mpydxc: -R("R_PYDXK") + R("R_PYDXPP")  = 0.

Mpydxnc: R("R_PDXPP") - R("R_PYDXNK")  = 0.

Mpyrb: R("Mpyrb_Transp") + R("R_EX_pyr_e_") - R("R_EX_pyr_e__Rev") 
 = 0.

Mpyrc: R("R_ACGAptspp") - R("R_ACHBS") - 2 R("R_ACLS") +
 R("R_ACMANAptspp") + R("R_ACMUMptspp") + R("R_ACNML") + R("R_ADCL") +
 R("R_ALATA_D2") + R("R_ALATA_L") + R("R_ALATA_L2") - R("R_ALATA_L_Rev") +
 R("R_ANS") + R("R_ASCBptspp") + R("R_CHRPL") + R("R_CYSDDS") + R("R_CYSDS") +
 R("R_CYSTL") + R("R_DAAD") + R("R_DAPAL") + R("R_DDPGALA") -
 R("R_DDPGALA_Rev") + R("R_DHAPT") - R("R_DHDPS") - R("R_DXPS") + R("R_EDA") +
 R("R_FRUpts2pp") + R("R_FRUptspp") + R("R_GALTptspp") +
 R("R_GAMptspp") + R("R_GLCptspp") + R("R_GLCRAL") + R("R_HOPNTAL") + R("R_ICHORT") +
 R("R_LDH_D") + R("R_LDH_D2") - R("R_LDH_D_Rev") + R("R_L_LACD2") +
 R("R_L_LACD3") + R("R_MALDDH") + R("R_MALTptspp") + R("R_MANGLYCptspp")
 + R("R_MANptspp") + R("R_MCITL2") - R("R_MCITL2_Rev") + R("R_MCPST") +
 R("R_ME1") + R("R_ME2") + R("R_MNLptspp") + R("R_OAADC") - R("R_PDH")
 - R("R_PFL") - R("R_POX") - R("R_PPS") + R("R_PYK") + R("R_PYRt2rpp")
 - R("R_PYRt2rpp_Rev") + R("R_SBTptspp") + R("R_SERD_D") +
 R("R_SERD_L") + R("R_SHCHCS2") + R("R_SUCptspp") + R("R_TREptspp") +
 R("R_TRPAS2") - R("R_TRPAS2_Rev") + R("R_VPAMT")  = 0.

Mpyre: -R("R_EX_pyr_e_") + R("R_EX_pyr_e__Rev") - R("R_PYRtex") +
 R("R_PYRtex_Rev")  = 0.

Mpyrp: -R("R_PYRt2rpp") + R("R_PYRt2rpp_Rev") + R("R_PYRtex") -
 R("R_PYRtex_Rev")  = 0.

Mq8c: -R("R_ASPO3") + R("R_CYTBDpp") + R("R_CYTBO3_4pp") -
 R("R_DHORD2") - R("R_DSBAO1") - R("R_FDH4pp") - R("R_G3PD5") - R("R_GLCDpp") -
 R("R_GLYCTO2") - R("R_HYD1pp") - R("R_LDH_D2") - R("R_L_LACD2") -
 R("R_MDH2") - R("R_NADH16pp") - R("R_NADH5") - R("R_NADPHQR2") +
 R("R_NO3R1bpp") + R("R_NO3R1pp") + 3 R("R_NTRIR3pp") - R("R_POX") + R("R_QMO2")
 - R("R_SUCDi")  = 0.

Mq8h2c: R("R_ASPO3") - R("R_CYTBDpp") - R("R_CYTBO3_4pp") +
 R("R_DHORD2") + R("R_DMQMT") + R("R_DSBAO1") + R("R_FDH4pp") + R("R_G3PD5") +
 R("R_GLCDpp") + R("R_GLYCTO2") + R("R_HYD1pp") + R("R_LDH_D2") +
 R("R_L_LACD2") + R("R_MDH2") + R("R_NADH16pp") + R("R_NADH5") +
 R("R_NADPHQR2") - R("R_NO3R1bpp") - R("R_NO3R1pp") - 3 R("R_NTRIR3pp") +
 R("R_POX") - R("R_QMO2") + R("R_SUCDi")  = 0.

Mqulnc: -R("R_NNDPR") + R("R_QULNS")  = 0.

Mr15bpc: -R("R_R15BPK") + R("R_R1PK")  = 0.

Mr1pc: -R("R_PPM") + R("R_PPM_Rev") + R("R_PUNP1") -
 R("R_PUNP1_Rev") + R("R_PUNP3") - R("R_PUNP3_Rev") + R("R_PUNP5") - R("R_PUNP5_Rev")
 + R("R_PUNP7") - R("R_PUNP7_Rev") + R("R_PYNP2r") - R("R_PYNP2r_Rev")
 - R("R_R1PK")  = 0.

Mr5pb: R("Mr5pb_Transp") + R("R_EX_r5p_e_") - R("R_EX_r5p_e__Rev") 
 = 0.

Mr5pc: R("R_ADPRDP") + R("R_AMPN") + R("R_CMPN") + R("R_NMNN") +
 R("R_NMNt7pp") + R("R_PPM") - R("R_PPM_Rev") - R("R_PRPPS") +
 R("R_PRPPS_Rev") - R("R_R5PP") + R("R_RBK") - R("R_RPI") + R("R_RPI_Rev") -
 R("R_TKT1") + R("R_TKT1_Rev")  = 0.

Mr5pe: -R("R_EX_r5p_e_") + R("R_EX_r5p_e__Rev") - R("R_R5Ptex") +
 R("R_R5Ptex_Rev")  = 0.

Mr5pp: -R("R_R5PPpp") + R("R_R5Ptex") - R("R_R5Ptex_Rev")  = 0.

Mragundc: -R("R_O16AT")  = 0.

Mrbflvrdc: -R("R_ARBTNR3") - R("R_CPGNR3") - R("R_FE3HOXR3") -
 R("R_FECRMR3") - R("R_FEENTERR3") - R("R_FEOXAMR3") + R("R_FLVR") +
 R("R_FLVRx")  = 0.

MrblLc: R("R_ARAI") - R("R_ARAI_Rev") - R("R_RBK_L1")  = 0.

Mrdmbzic: -R("R_ADOCBLS") + R("R_RZ5PP")  = 0.

Mrhcysc: R("R_AHCYSNS") - R("R_RHCCE")  = 0.

MribDb: R("MribDb_Transp") + R("R_EX_rib_D_e_") -
 R("R_EX_rib_D_e__Rev")  = 0.

MribDc: R("R_ADNUC") + R("R_CYTDH") + R("R_INSH") + R("R_R5PP") -
 R("R_RBK") + R("R_RIBabcpp") + R("R_URIH") + R("R_XTSNH")  = 0.

MribDe: -R("R_EX_rib_D_e_") + R("R_EX_rib_D_e__Rev") - R("R_RIBtex")
 + R("R_RIBtex_Rev")  = 0.

MribDp: R("R_R5PPpp") - R("R_RIBabcpp") + R("R_RIBtex") -
 R("R_RIBtex_Rev")  = 0.

Mribflvc: R("R_ARBTNR3") + R("R_CPGNR3") - 0.000223
 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_FE3HOXR3") + R("R_FECRMR3") +
 R("R_FEENTERR3") + R("R_FEOXAMR3") - R("R_FLVR") - R("R_FLVRx") - R("R_RBFK") +
 R("R_RBFSb")  = 0.

Mrml1pc: R("R_RMK") - R("R_RMPA") + R("R_RMPA_Rev")  = 0.

Mrmlc: R("R_RMI") - R("R_RMI_Rev") - R("R_RMK")  = 0.

Mrmnb: R("Mrmnb_Transp") + R("R_EX_rmn_e_") - R("R_EX_rmn_e__Rev") 
 = 0.

Mrmnc: -R("R_RMI") + R("R_RMI_Rev") + R("R_RMNtpp")  = 0.

Mrmne: -R("R_EX_rmn_e_") + R("R_EX_rmn_e__Rev") - R("R_RMNtex") +
 R("R_RMNtex_Rev")  = 0.

Mrmnp: R("R_RMNtex") - R("R_RMNtex_Rev") - R("R_RMNtpp")  = 0.

Mru5pDc: -R("R_A5PISO") + R("R_A5PISO_Rev") - R("R_DB4PS") +
 R("R_GND") - R("R_RPE") + R("R_RPE_Rev") + R("R_RPI") - R("R_RPI_Rev")  = 0.

Mru5pLc: R("R_RBK_L1") - R("R_RBP4E") + R("R_RBP4E_Rev") +
 R("R_X5PL3E")  = 0.

Ms7pc: -R("R_S7PI") - R("R_TALA") + R("R_TALA_Rev") + R("R_TKT1") -
 R("R_TKT1_Rev")  = 0.

Msarcsc: -R("R_SARCOX")  = 0.

Msbt6pc: -R("R_SBTPD") + R("R_SBTPD_Rev") + R("R_SBTptspp")  = 0.

MsbtDb: R("MsbtDb_Transp") + R("R_EX_sbt_D_e_") -
 R("R_EX_sbt_D_e__Rev")  = 0.

MsbtDe: -R("R_EX_sbt_D_e_") + R("R_EX_sbt_D_e__Rev") - R("R_SBTtex")
 + R("R_SBTtex_Rev")  = 0.

MsbtDp: -R("R_SBTptspp") + R("R_SBTtex") - R("R_SBTtex_Rev")  = 0.

Msbzcoac: -R("R_NPHS") + R("R_SUCBZL")  = 0.

Msc: -R("R_BTS4")  = 0.

Msclc: R("R_SHCHD2") - R("R_SHCHF")  = 0.

Msectrnac: R("R_SELCYSS")  = 0.

Mselnc: -R("R_SELNPS")  = 0.

Mselnpc: -R("R_SELCYSS") + R("R_SELNPS")  = 0.

Mserampc: -3 R("R_ENTCS") + R("R_SERASr") - R("R_SERASr_Rev")  = 0.

MserDb: R("MserDb_Transp") + R("R_EX_ser_D_e_") -
 R("R_EX_ser_D_e__Rev")  = 0.

MserDc: -R("R_DSERDHr") + R("R_DSERDHr_Rev") + R("R_DSERt2pp") -
 R("R_SERD_D")  = 0.

MserDe: -R("R_DSERtex") + R("R_DSERtex_Rev") - R("R_EX_ser_D_e_") +
 R("R_EX_ser_D_e__Rev")  = 0.

MserDp: -R("R_DSERt2pp") + R("R_DSERtex") - R("R_DSERtex_Rev")  = 0.

MserLb: R("MserLb_Transp") + R("R_EX_ser_L_e_") -
 R("R_EX_ser_L_e__Rev")  = 0.

MserLc: R("R_DHBSH") - 0.2158 R("R_Ec_biomass_iAF1260_core_59p81M")
 - R("R_GHMT2r") + R("R_GHMT2r_Rev") + R("R_GPDDA3") - R("R_LSERDHr") +
 R("R_LSERDHr_Rev") + R("R_PSP_L") - R("R_PSSA120") - R("R_PSSA140") -
 R("R_PSSA141") - R("R_PSSA160") - R("R_PSSA161") - R("R_PSSA180") -
 R("R_PSSA181") - R("R_SERASr") + R("R_SERASr_Rev") - R("R_SERAT") +
 R("R_SERAT_Rev") - R("R_SERD_L") + R("R_SERt2rpp") - R("R_SERt2rpp_Rev")
 + R("R_SERt4pp") - R("R_SERTRS") - R("R_SERTRS2") - R("R_TRPS1") -
 R("R_TRPS2")  = 0.

MserLe: -R("R_EX_ser_L_e_") + R("R_EX_ser_L_e__Rev") - R("R_SERtex")
 + R("R_SERtex_Rev")  = 0.

MserLp: R("R_GPDDA3pp") + R("R_PSP_Lpp") - R("R_SERt2rpp") +
 R("R_SERt2rpp_Rev") - R("R_SERt4pp") + R("R_SERtex") - R("R_SERtex_Rev")  =
 0.

Msertrnac: R("R_SERTRS")  = 0.

Msertrnasecoooc: -R("R_SELCYSS") + R("R_SERTRS2")  = 0.

MSfglutthc: R("R_FALDH2") - R("R_FALDH2_Rev") - R("R_SFGTHi")  = 0.

Mshemec: -0.000223 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_SHCHF")  = 0.

Mskm5pc: -R("R_PSCVT") + R("R_PSCVT_Rev") + R("R_SHKK")  = 0.

Mskmb: R("Mskmb_Transp") + R("R_EX_skm_e_") - R("R_EX_skm_e__Rev") 
 = 0.

Mskmc: R("R_SHK3Dr") - R("R_SHK3Dr_Rev") - R("R_SHKK") +
 R("R_SKMt2pp")  = 0.

Mskme: -R("R_EX_skm_e_") + R("R_EX_skm_e__Rev") - R("R_SKMtex") +
 R("R_SKMtex_Rev")  = 0.

Mskmp: -R("R_SKMt2pp") + R("R_SKMtex") - R("R_SKMtex_Rev")  = 0.

Msl26dac: -R("R_SDPDS") - R("R_SDPTA") + R("R_SDPTA_Rev")  = 0.

Msl2a6oc: R("R_SDPTA") - R("R_SDPTA_Rev") + R("R_THDPS")  = 0.

Mso2b: R("Mso2b_Transp") + R("R_EX_so2_e_") - R("R_EX_so2_e__Rev") 
 = 0.

Mso2c: R("R_CYSSADS") + R("R_SO2tpp") - R("R_SO2tpp_Rev")  = 0.

Mso2e: -R("R_EX_so2_e_") + R("R_EX_so2_e__Rev") - R("R_SO2tex") +
 R("R_SO2tex_Rev")  = 0.

Mso2p: R("R_SO2tex") - R("R_SO2tex_Rev") - R("R_SO2tpp") +
 R("R_SO2tpp_Rev")  = 0.

Mso3b: R("Mso3b_Transp") + R("R_EX_so3_e_") - R("R_EX_so3_e__Rev") 
 = 0.

Mso3c: R("R_CYANST") + R("R_FDMO") + R("R_FDMO2") + R("R_FDMO3") +
 R("R_FDMO4") + R("R_FDMO6") + R("R_PAPSR") + R("R_PAPSR2") -
 R("R_SULRi") + R("R_TAUDO")  = 0.

Mso3e: -R("R_EX_so3_e_") + R("R_EX_so3_e__Rev") - R("R_SO3tex") +
 R("R_SO3tex_Rev")  = 0.

Mso3p: R("R_CYANSTpp") + R("R_SO3tex") - R("R_SO3tex_Rev")  = 0.

Mso4b: R("Mso4b_Transp") + R("R_EX_so4_e_") - R("R_EX_so4_e__Rev") 
 = 0.

Mso4c: -0.003948 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_H2SO")
 - R("R_SADT2") + R("R_SULabcpp")  = 0.

Mso4e: -R("R_EX_so4_e_") + R("R_EX_so4_e__Rev") - R("R_SO4tex") +
 R("R_SO4tex_Rev")  = 0.

Mso4p: R("R_SO4tex") - R("R_SO4tex_Rev") - R("R_SULabcpp")  = 0.

Mspmdb: R("Mspmdb_Transp") + R("R_EX_spmd_e_") -
 R("R_EX_spmd_e__Rev")  = 0.

Mspmdc: R("R_GSPMDA") - R("R_GSPMDS") + R("R_SPMDabcpp") -
 R("R_SPMDAT1") - R("R_SPMDAT2") + R("R_SPMS")  = 0.

Mspmde: -R("R_EX_spmd_e_") + R("R_EX_spmd_e__Rev") - R("R_SPMDtex")
 + R("R_SPMDtex_Rev")  = 0.

Mspmdp: -R("R_SPMDabcpp") + R("R_SPMDtex") - R("R_SPMDtex_Rev")  =
 0.

Mssaltppc: R("R_OXGDC2") - R("R_SHCHCS2")  = 0.

Mstcoac: -R("R_ACOAD8f") + R("R_ACOAD8f_Rev") - R("R_FACOAE180") +
 R("R_FACOAL180t2pp")  = 0.

Msuc6pc: -R("R_FFSD") + R("R_SUCptspp")  = 0.

Msucargc: R("R_AST") - R("R_SADH")  = 0.

Msucbzc: -R("R_SUCBZL") + R("R_SUCBZS")  = 0.

Msuccb: R("Msuccb_Transp") + R("R_EX_succ_e_") -
 R("R_EX_succ_e__Rev")  = 0.

Msuccc: R("R_ASPO5") - R("R_CITt7pp") + R("R_FRD2") + R("R_FRD3") +
 R("R_HKNDDH") + R("R_ICL") + R("R_MCITL2") - R("R_MCITL2_Rev") -
 R("R_PPCSCT") + 3 R("R_PPPGO3") + R("R_SDPDS") + R("R_SGDS") + R("R_SHSL1")
 + R("R_SSALx") + R("R_SSALy") + R("R_SUCCt2_2pp") + R("R_SUCCt2_3pp")
 - R("R_SUCCt3pp") - R("R_SUCDi") - R("R_SUCFUMtpp") +
 R("R_SUCFUMtpp_Rev") - R("R_SUCOAS") + R("R_SUCOAS_Rev") - R("R_TARTRt7pp") +
 R("R_TARTRt7pp_Rev") + R("R_TAUDO")  = 0.

Msucce: -R("R_EX_succ_e_") + R("R_EX_succ_e__Rev") - R("R_SUCCtex")
 + R("R_SUCCtex_Rev")  = 0.

Msuccoac: R("R_AKGDH") - R("R_AST") - R("R_HSST") - R("R_MMM2") +
 R("R_PPCSCT") + R("R_SUCOAS") - R("R_SUCOAS_Rev") - R("R_THDPS")  = 0.

Msuccp: R("R_CITt7pp") - R("R_SUCCt2_2pp") - R("R_SUCCt2_3pp") +
 R("R_SUCCt3pp") + R("R_SUCCtex") - R("R_SUCCtex_Rev") + R("R_SUCFUMtpp")
 - R("R_SUCFUMtpp_Rev") + R("R_TARTRt7pp") - R("R_TARTRt7pp_Rev")  =
 0.

Msucgluc: -R("R_SGDS") + R("R_SGSAD")  = 0.

Msucgsac: -R("R_SGSAD") + R("R_SOTA")  = 0.

Msuchmsc: R("R_HSST") - R("R_SHSL1")  = 0.

Msucornc: R("R_SADH") - R("R_SOTA")  = 0.

Msucrb: R("Msucrb_Transp") + R("R_EX_sucr_e_") -
 R("R_EX_sucr_e__Rev")  = 0.

Msucre: -R("R_EX_sucr_e_") + R("R_EX_sucr_e__Rev") - R("R_SUCRtex")
 + R("R_SUCRtex_Rev")  = 0.

Msucrp: -R("R_SUCptspp") + R("R_SUCRtex") - R("R_SUCRtex_Rev")  = 0.

Msucsalc: R("R_ABTA") - R("R_SSALx") - R("R_SSALy")  = 0.

Msulfacb: R("Msulfacb_Transp") + R("R_EX_sulfac_e_") -
 R("R_EX_sulfac_e__Rev")  = 0.

Msulfacc: -R("R_FDMO6") + R("R_SULFACabcpp")  = 0.

Msulface: -R("R_EX_sulfac_e_") + R("R_EX_sulfac_e__Rev") -
 R("R_SULFACtex") + R("R_SULFACtex_Rev")  = 0.

Msulfacp: -R("R_SULFACabcpp") + R("R_SULFACtex") -
 R("R_SULFACtex_Rev")  = 0.

Mt3c11vaceACPc: R("R_3HAD181") - R("R_EAR181x") - R("R_EAR181y")  =
 0.

Mt3c5ddeceACPc: R("R_3HAD121") - R("R_EAR121x") - R("R_EAR121y")  =
 0.

Mt3c7mrseACPc: R("R_3HAD141") - R("R_EAR141x") - R("R_EAR141y")  =
 0.

Mt3c9palmeACPc: R("R_3HAD161") - R("R_EAR161x") - R("R_EAR161y")  =
 0.

Mtag6pDc: R("R_GLTPD") - R("R_GLTPD_Rev") - R("R_PFK_2")  = 0.

MtagdpDc: R("R_PFK_2") - R("R_TGBPA") + R("R_TGBPA_Rev")  = 0.

Mtagurc: R("R_GALCTLO") + R("R_GUI2") - R("R_GUI2_Rev") +
 R("R_TAGURr") - R("R_TAGURr_Rev")  = 0.

MtartrLb: R("MtartrLb_Transp") + R("R_EX_tartr_L_e_") -
 R("R_EX_tartr_L_e__Rev")  = 0.

MtartrLc: -R("R_TARTD") + R("R_TARTRt7pp") - R("R_TARTRt7pp_Rev")  =
 0.

MtartrLe: -R("R_EX_tartr_L_e_") + R("R_EX_tartr_L_e__Rev") -
 R("R_TARTRtex") + R("R_TARTRtex_Rev")  = 0.

MtartrLp: -R("R_TARTRt7pp") + R("R_TARTRt7pp_Rev") + R("R_TARTRtex")
 - R("R_TARTRtex_Rev")  = 0.

Mtaurb: R("Mtaurb_Transp") + R("R_EX_taur_e_") -
 R("R_EX_taur_e__Rev")  = 0.

Mtaurc: -R("R_TAUDO") + R("R_TAURabcpp")  = 0.

Mtaure: -R("R_EX_taur_e_") + R("R_EX_taur_e__Rev") - R("R_TAURtex")
 + R("R_TAURtex_Rev")  = 0.

Mtaurp: -R("R_TAURabcpp") + R("R_TAURtex") - R("R_TAURtex_Rev")  =
 0.

Mtcyntb: R("Mtcyntb_Transp") + R("R_EX_tcynt_e_") -
 R("R_EX_tcynt_e__Rev")  = 0.

Mtcyntc: R("R_CYANST") + R("R_MCPST")  = 0.

Mtcynte: -R("R_EX_tcynt_e_") + R("R_EX_tcynt_e__Rev") -
 R("R_TCYNTtex") + R("R_TCYNTtex_Rev")  = 0.

Mtcyntp: R("R_CYANSTpp") + R("R_TCYNTtex") - R("R_TCYNTtex_Rev")  =
 0.

Mtd2coac: R("R_ACOAD6f") - R("R_ACOAD6f_Rev") + R("R_ECOAH6") -
 R("R_ECOAH6_Rev") + R("R_TDECOAI")  = 0.

Mtdcoac: -R("R_ACOAD6f") + R("R_ACOAD6f_Rev") - R("R_FACOAE140") +
 R("R_FACOAL140t2pp") + R("R_KAT7")  = 0.

Mtddec2eACPc: R("R_3HAD120") - R("R_EAR120x") - R("R_EAR120y")  = 0.

MtdeACPc: -R("R_3OAS161") + R("R_AACPS2") - R("R_AGPAT141") +
 R("R_EAR141x") + R("R_EAR141y") - R("R_FA141ACPHi") - R("R_G3PAT141")  = 0.

Mtdec2eACPc: R("R_3HAD100") - R("R_EAR100x") - R("R_EAR100y") -
 R("R_T2DECAI") + R("R_T2DECAI_Rev")  = 0.

Mtdecoac: -R("R_FACOAE141") + R("R_FACOAL141t2pp") - R("R_TDECOAI") 
 = 0.

Mthdpc: R("R_DHDPRy") - R("R_THDPS")  = 0.

Mthex2eACPc: R("R_3HAD60") - R("R_EAR60x") - R("R_EAR60y")  = 0.

Mthfc: R("R_AICART") - R("R_AICART_Rev") + R("R_DHFR") -
 R("R_DHFR_Rev") - 0.000223 R("R_Ec_biomass_iAF1260_core_59p81M") + R("R_FMETTRS")
 + R("R_FTHFD") + R("R_GARFT") - R("R_GARFT_Rev") - R("R_GHMT2r") +
 R("R_GHMT2r_Rev") - R("R_GLYCL") + R("R_METS") + R("R_MOHMT") +
 R("R_ULA4NFT")  = 0.

Mthmb: R("Mthmb_Transp") + R("R_EX_thm_e_") - R("R_EX_thm_e__Rev") 
 = 0.

Mthmc: R("R_THMabcpp") - R("R_TMK")  = 0.

Mthme: -R("R_EX_thm_e_") + R("R_EX_thm_e__Rev") - R("R_THMtex") +
 R("R_THMtex_Rev")  = 0.

Mthmmpc: R("R_TDP") + R("R_TMK") - R("R_TMPK") + R("R_TMPPP")  = 0.

Mthmp: -R("R_THMabcpp") + R("R_THMtex") - R("R_THMtex_Rev")  = 0.

Mthmppc: -0.000223 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_OXGDC2") + R("R_SHCHCS2") - R("R_TDP") + R("R_TMPK")  = 0.

MthrLb: R("MthrLb_Transp") + R("R_EX_thr_L_e_") -
 R("R_EX_thr_L_e__Rev")  = 0.

MthrLc: -0.2537 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_THRabcpp") - R("R_THRAi") - R("R_THRD") - R("R_THRD_L") + R("R_THRS") -
 R("R_THRt2pp") + R("R_THRt2rpp") - R("R_THRt2rpp_Rev") + R("R_THRt4pp") -
 R("R_THRTRS")  = 0.

MthrLe: -R("R_EX_thr_L_e_") + R("R_EX_thr_L_e__Rev") - R("R_THRtex")
 + R("R_THRtex_Rev")  = 0.

MthrLp: R("R_PTHRpp") - R("R_THRabcpp") + R("R_THRt2pp") -
 R("R_THRt2rpp") + R("R_THRt2rpp_Rev") - R("R_THRt4pp") + R("R_THRtex") -
 R("R_THRtex_Rev")  = 0.

Mthrpb: R("Mthrpb_Transp") + R("R_EX_thrp_e_") -
 R("R_EX_thrp_e__Rev")  = 0.

Mthrpe: -R("R_EX_thrp_e_") + R("R_EX_thrp_e__Rev") - R("R_THRPtex")
 + R("R_THRPtex_Rev")  = 0.

Mthrpp: -R("R_PTHRpp") + R("R_THRPtex") - R("R_THRPtex_Rev")  = 0.

Mthrtrnac: R("R_THRTRS")  = 0.

Mthymb: R("Mthymb_Transp") + R("R_EX_thym_e_") -
 R("R_EX_thym_e__Rev")  = 0.

Mthymc: -R("R_THYMt3pp") + R("R_TMDPP") - R("R_TMDPP_Rev")  = 0.

Mthymdb: R("Mthymdb_Transp") + R("R_EX_thymd_e_") -
 R("R_EX_thymd_e__Rev")  = 0.

Mthymdc: R("R_NTD5") + R("R_THMDt2pp") + R("R_THMDt2rpp") -
 R("R_THMDt2rpp_Rev") - R("R_TMDK1") - R("R_TMDPP") + R("R_TMDPP_Rev")  = 0.

Mthymde: -R("R_EX_thymd_e_") + R("R_EX_thymd_e__Rev") -
 R("R_THMDtex") + R("R_THMDtex_Rev")  = 0.

Mthymdp: R("R_NTD5pp") - R("R_THMDt2pp") - R("R_THMDt2rpp") +
 R("R_THMDt2rpp_Rev") + R("R_THMDtex") - R("R_THMDtex_Rev")  = 0.

Mthyme: -R("R_EX_thym_e_") + R("R_EX_thym_e__Rev") - R("R_THYMtex")
 + R("R_THYMtex_Rev")  = 0.

Mthymp: R("R_THYMt3pp") + R("R_THYMtex") - R("R_THYMtex_Rev")  = 0.

Mtmab: R("Mtmab_Transp") + R("R_EX_tma_e_") - R("R_EX_tma_e__Rev") 
 = 0.

Mtmac: R("R_TMAOR1") + R("R_TMAOR2")  = 0.

Mtmae: -R("R_EX_tma_e_") + R("R_EX_tma_e__Rev") - R("R_TMAtex") +
 R("R_TMAtex_Rev")  = 0.

Mtmaob: R("Mtmaob_Transp") + R("R_EX_tmao_e_") -
 R("R_EX_tmao_e__Rev")  = 0.

Mtmaoc: -R("R_TMAOR1") - R("R_TMAOR2")  = 0.

Mtmaoe: -R("R_EX_tmao_e_") + R("R_EX_tmao_e__Rev") - R("R_TMAOtex")
 + R("R_TMAOtex_Rev")  = 0.

Mtmaop: -R("R_TMAOR1pp") - R("R_TMAOR2pp") + R("R_TMAOtex") -
 R("R_TMAOtex_Rev")  = 0.

Mtmap: R("R_TMAOR1pp") + R("R_TMAOR2pp") + R("R_TMAtex") -
 R("R_TMAtex_Rev")  = 0.

Mtmrs2eACPc: R("R_3HAD140") - R("R_EAR140x") - R("R_EAR140y")  = 0.

Mtoct2eACPc: R("R_3HAD80") - R("R_EAR80x") - R("R_EAR80y")  = 0.

Mtoctd2eACPc: R("R_3HAD180") - R("R_EAR180x") - R("R_EAR180y")  = 0.

Mtpalm2eACPc: R("R_3HAD160") - R("R_EAR160x") - R("R_EAR160y")  = 0.

Mtrdoxc: R("R_DSBDR") + R("R_METSOXR1") + R("R_METSOXR2") +
 R("R_PAPSR") + R("R_RNDR1") + R("R_RNDR2") + R("R_RNDR3") + R("R_RNDR4") +
 R("R_THIORDXi") - R("R_TRDR")  = 0.

Mtrdrdc: -R("R_DSBDR") - R("R_METSOXR1") - R("R_METSOXR2") -
 R("R_PAPSR") - R("R_RNDR1") - R("R_RNDR2") - R("R_RNDR3") - R("R_RNDR4") -
 R("R_THIORDXi") + R("R_TRDR")  = 0.

Mtre6pc: -R("R_TRE6PH") - R("R_TRE6PP") + R("R_TRE6PS") +
 R("R_TREptspp")  = 0.

Mtreb: R("Mtreb_Transp") + R("R_EX_tre_e_") - R("R_EX_tre_e__Rev") 
 = 0.

Mtrec: R("R_TRE6PP") - R("R_TREH")  = 0.

Mtree: -R("R_EX_tre_e_") + R("R_EX_tre_e__Rev") - R("R_TREtex") +
 R("R_TREtex_Rev")  = 0.

Mtrep: -R("R_TREHpp") - R("R_TREptspp") + R("R_TREtex") -
 R("R_TREtex_Rev")  = 0.

Mtrnaalac: -R("R_ALATRS")  = 0.

Mtrnaargc: -R("R_ARGTRS")  = 0.

Mtrnaasnc: -R("R_ASNTRS")  = 0.

Mtrnaaspc: -R("R_ASPTRS")  = 0.

Mtrnacysc: -R("R_CYSTRS")  = 0.

Mtrnaglnc: -R("R_GLNTRS")  = 0.

Mtrnagluc: R("R_GLUTRR") - R("R_GLUTRS")  = 0.

Mtrnaglyc: -R("R_GLYTRS")  = 0.

Mtrnahisc: -R("R_HISTRS")  = 0.

Mtrnailec: -R("R_ILETRS")  = 0.

Mtrnaleuc: -R("R_LEUTRS")  = 0.

Mtrnalysc: -R("R_LYSTRS")  = 0.

Mtrnametc: -R("R_METTRS")  = 0.

Mtrnaphec: -R("R_PHETRS")  = 0.

Mtrnaproc: -R("R_PROTRS")  = 0.

Mtrnasecysc: -R("R_SERTRS2")  = 0.

Mtrnaserc: -R("R_SERTRS")  = 0.

Mtrnathrc: -R("R_THRTRS")  = 0.

Mtrnatrpc: -R("R_TRPTRS")  = 0.

Mtrnatyrc: -R("R_TYRTRS")  = 0.

Mtrnavalc: -R("R_VALTRS")  = 0.

MtrpLb: R("MtrpLb_Transp") + R("R_EX_trp_L_e_") -
 R("R_EX_trp_L_e__Rev")  = 0.

MtrpLc: -0.05684 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_MTRPOX") - R("R_TRPAS2") + R("R_TRPAS2_Rev") + R("R_TRPS1") + R("R_TRPS2") +
 R("R_TRPt2rpp") - R("R_TRPt2rpp_Rev") - R("R_TRPTRS")  = 0.

MtrpLe: -R("R_EX_trp_L_e_") + R("R_EX_trp_L_e__Rev") - R("R_TRPtex")
 + R("R_TRPtex_Rev")  = 0.

MtrpLp: -R("R_TRPt2rpp") + R("R_TRPt2rpp_Rev") + R("R_TRPtex") -
 R("R_TRPtex_Rev")  = 0.

Mtrptrnac: R("R_TRPTRS")  = 0.

Mtsulb: R("Mtsulb_Transp") + R("R_EX_tsul_e_") -
 R("R_EX_tsul_e__Rev")  = 0.

Mtsulc: -R("R_CYANST") + R("R_TSULabcpp")  = 0.

Mtsule: -R("R_EX_tsul_e_") + R("R_EX_tsul_e__Rev") - R("R_TSULtex")
 + R("R_TSULtex_Rev")  = 0.

Mtsulp: -R("R_CYANSTpp") - R("R_TSULabcpp") + R("R_TSULtex") -
 R("R_TSULtex_Rev")  = 0.

Mttdcab: R("Mttdcab_Transp") + R("R_EX_ttdca_e_") -
 R("R_EX_ttdca_e__Rev")  = 0.

Mttdcac: -R("R_2AGPEAT140") - R("R_2AGPGAT140") - R("R_AACPS1") +
 R("R_FA140ACPHi") + R("R_FACOAE140") + R("R_LPLIPAL2A140") +
 R("R_LPLIPAL2E140") + R("R_LPLIPAL2G140")  = 0.

Mttdcae: -R("R_EX_ttdca_e_") + R("R_EX_ttdca_e__Rev") -
 R("R_TTDCAtexi")  = 0.

Mttdcap: -R("R_FACOAL140t2pp") + R("R_LPLIPAL1A140pp") +
 R("R_LPLIPAL1E140pp") + R("R_LPLIPAL1G140pp") + R("R_PLIPA1A140pp") +
 R("R_PLIPA1E140pp") + R("R_PLIPA1G140pp") + R("R_PLIPA2A140pp") +
 R("R_PLIPA2E140pp") + R("R_PLIPA2G140pp") + R("R_TTDCAtexi")  = 0.

Mttdceab: R("Mttdceab_Transp") + R("R_EX_ttdcea_e_") -
 R("R_EX_ttdcea_e__Rev")  = 0.

Mttdceac: -R("R_2AGPEAT141") - R("R_2AGPGAT141") - R("R_AACPS2") +
 R("R_FA141ACPHi") + R("R_FACOAE141") + R("R_LPLIPAL2A141") +
 R("R_LPLIPAL2E141") + R("R_LPLIPAL2G141")  = 0.

Mttdceae: -R("R_EX_ttdcea_e_") + R("R_EX_ttdcea_e__Rev") -
 R("R_TTDCEAtexi")  = 0.

Mttdceap: -R("R_FACOAL141t2pp") + R("R_LPLIPAL1A141pp") +
 R("R_LPLIPAL1E141pp") + R("R_LPLIPAL1G141pp") + R("R_PLIPA1A141pp") +
 R("R_PLIPA1E141pp") + R("R_PLIPA1G141pp") + R("R_PLIPA2A141pp") +
 R("R_PLIPA2E141pp") + R("R_PLIPA2G141pp") + R("R_TTDCEAtexi")  = 0.

Mtungsb: R("Mtungsb_Transp") + R("R_EX_tungs_e_") -
 R("R_EX_tungs_e__Rev")  = 0.

Mtungsc: R("R_TUNGSabcpp")  = 0.

Mtungse: -R("R_EX_tungs_e_") + R("R_EX_tungs_e__Rev") -
 R("R_TUNGStex") + R("R_TUNGStex_Rev")  = 0.

Mtungsp: -R("R_TUNGSabcpp") + R("R_TUNGStex") - R("R_TUNGStex_Rev") 
 = 0.

Mtymb: R("Mtymb_Transp") + R("R_EX_tym_e_") - R("R_EX_tym_e__Rev") 
 = 0.

Mtyme: -R("R_EX_tym_e_") + R("R_EX_tym_e__Rev") - R("R_TYMtex") +
 R("R_TYMtex_Rev")  = 0.

Mtymp: R("R_TYMtex") - R("R_TYMtex_Rev") - R("R_TYROXDApp")  = 0.

MtyrLb: R("MtyrLb_Transp") + R("R_EX_tyr_L_e_") -
 R("R_EX_tyr_L_e__Rev")  = 0.

MtyrLc: -0.1379 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_THZPSN") + R("R_TYRt2rpp") - R("R_TYRt2rpp_Rev") - R("R_TYRTA") +
 R("R_TYRTA_Rev") - R("R_TYRTRS")  = 0.

MtyrLe: -R("R_EX_tyr_L_e_") + R("R_EX_tyr_L_e__Rev") - R("R_TYRtex")
 + R("R_TYRtex_Rev")  = 0.

MtyrLp: R("R_TYRPpp") - R("R_TYRt2rpp") + R("R_TYRt2rpp_Rev") +
 R("R_TYRtex") - R("R_TYRtex_Rev")  = 0.

Mtyrpb: R("Mtyrpb_Transp") + R("R_EX_tyrp_e_") -
 R("R_EX_tyrp_e__Rev")  = 0.

Mtyrpe: -R("R_EX_tyrp_e_") + R("R_EX_tyrp_e__Rev") - R("R_TYRPtex")
 + R("R_TYRPtex_Rev")  = 0.

Mtyrpp: -R("R_TYRPpp") + R("R_TYRPtex") - R("R_TYRPtex_Rev")  = 0.

Mtyrtrnac: R("R_TYRTRS")  = 0.

Mu23gac: -R("R_LPADSS") + R("R_U23GAAT") - R("R_USHD")  = 0.

Mu3agac: R("R_UAGAAT") - R("R_UAGAAT_Rev") - R("R_UHGADA")  = 0.

Mu3hgac: -R("R_U23GAAT") + R("R_UHGADA")  = 0.

Muaagmdac: -2 R("R_MPTG") - R("R_MPTG2") + R("R_UAGPT3")  = 0.

Muaccgc: R("R_UAGCVT") - R("R_UAPGR")  = 0.

Muacgamb: R("Muacgamb_Transp") + R("R_EX_uacgam_e_") -
 R("R_EX_uacgam_e__Rev")  = 0.

Muacgamc: -R("R_ACGAMT") - R("R_UAG2E") + R("R_UAG2E_Rev") -
 R("R_UAGAAT") + R("R_UAGAAT_Rev") - R("R_UAGCVT") + R("R_UAGDP") -
 R("R_UAGPT3")  = 0.

Muacgame: -R("R_EX_uacgam_e_") + R("R_EX_uacgam_e__Rev") -
 R("R_UACGAMtex") + R("R_UACGAMtex_Rev")  = 0.

Muacgamp: -R("R_UACGAMPpp") + R("R_UACGAMtex") -
 R("R_UACGAMtex_Rev")  = 0.

Muacmamc: -R("R_UACMAMO") + R("R_UAG2E") - R("R_UAG2E_Rev")  = 0.

Muacmamuc: -R("R_ACMAMUT") + R("R_UACMAMO")  = 0.

Muagmdac: R("R_PAPPT3") - R("R_UAGPT3")  = 0.

Muamac: -R("R_UAMAGS") + R("R_UAMAS")  = 0.

Muamagc: -R("R_UAAGDS") + R("R_UAMAGS")  = 0.

Muamrc: -R("R_UAMAS") + R("R_UAPGR") - R("R_UM3PL") - R("R_UM4PL") 
 = 0.

Mudcpdpc: -0.000055 R("R_Ec_biomass_iAF1260_core_59p81M") + 2
 R("R_MPTG") + R("R_MPTG2") - R("R_UDCPDP") + R("R_UDCPDPS")  = 0.

Mudcpdpp: R("R_ECA4OALpp") + R("R_ECAP1pp") + R("R_ECAP2pp") +
 R("R_ECAP3pp") + R("R_O16A4Lpp") + R("R_O16AP1pp") + R("R_O16AP2pp") +
 R("R_O16AP3pp") - R("R_UDCPDPpp")  = 0.

Mudcppc: -R("R_ACGAMT") - R("R_PAPPT3") + R("R_UDCPDP") +
 R("R_UDCPPtppi") - R("R_UPLA4FNT")  = 0.

Mudcppp: R("R_LA4NTpp") + R("R_UDCPDPpp") - R("R_UDCPPtppi")  = 0.

Mudpacgalb: R("Mudpacgalb_Transp") + R("R_EX_udpacgal_e_") -
 R("R_EX_udpacgal_e__Rev")  = 0.

Mudpacgale: -R("R_EX_udpacgal_e_") + R("R_EX_udpacgal_e__Rev") -
 R("R_UDPACGALtex") + R("R_UDPACGALtex_Rev")  = 0.

Mudpacgalp: -R("R_UACGALPpp") + R("R_UDPACGALtex") -
 R("R_UDPACGALtex_Rev")  = 0.

Mudpc: R("R_ACMAMUT") + R("R_GALT1") + R("R_GLCTR1") + R("R_GLCTR2")
 + R("R_GLCTR3") + R("R_LPADSS") - R("R_NDPK2") + R("R_NDPK2_Rev") +
 R("R_O16GALFT") + R("R_O16GLCT1") + R("R_O16GLCT2") - R("R_RNDR4") -
 R("R_RNDR4b") + R("R_TRE6PS") + R("R_UAGPT3") + R("R_UMPK") -
 R("R_UMPK_Rev") + R("R_UPLA4FNT")  = 0.

Mudpgalb: R("Mudpgalb_Transp") + R("R_EX_udpgal_e_") -
 R("R_EX_udpgal_e__Rev")  = 0.

Mudpgalc: R("R_UDPG4E") - R("R_UDPG4E_Rev") - R("R_UDPGALM") +
 R("R_UGLT") - R("R_UGLT_Rev")  = 0.

Mudpgale: -R("R_EX_udpgal_e_") + R("R_EX_udpgal_e__Rev") -
 R("R_UDPGALtex") + R("R_UDPGALtex_Rev")  = 0.

Mudpgalfurc: -R("R_O16GALFT") + R("R_UDPGALM")  = 0.

Mudpgalp: -R("R_UDPGALPpp") + R("R_UDPGALtex") -
 R("R_UDPGALtex_Rev")  = 0.

Mudpgb: R("Mudpgb_Transp") + R("R_EX_udpg_e_") -
 R("R_EX_udpg_e__Rev")  = 0.

Mudpgc: -R("R_GALT1") + R("R_GALUi") - R("R_GLCTR1") - R("R_GLCTR2")
 - R("R_GLCTR3") - R("R_O16GLCT1") - R("R_O16GLCT2") - R("R_TRE6PS") -
 R("R_UDPG4E") + R("R_UDPG4E_Rev") - R("R_UDPGD") - R("R_UGLT") +
 R("R_UGLT_Rev")  = 0.

Mudpge: -R("R_EX_udpg_e_") + R("R_EX_udpg_e__Rev") - R("R_UDPGtex")
 + R("R_UDPGtex_Rev")  = 0.

Mudpglcurb: R("Mudpglcurb_Transp") + R("R_EX_udpglcur_e_") -
 R("R_EX_udpglcur_e__Rev")  = 0.

Mudpglcurc: R("R_UDPGD") - R("R_UDPGDC")  = 0.

Mudpglcure: -R("R_EX_udpglcur_e_") + R("R_EX_udpglcur_e__Rev") -
 R("R_UDPGLCURtex") + R("R_UDPGLCURtex_Rev")  = 0.

Mudpglcurp: R("R_UDPGLCURtex") - R("R_UDPGLCURtex_Rev") -
 R("R_UGLCURPpp")  = 0.

Mudpgp: -R("R_UDPGPpp") + R("R_UDPGtex") - R("R_UDPGtex_Rev")  = 0.

MudpLa4fnc: R("R_ULA4NFT") - R("R_UPLA4FNT")  = 0.

MudpLa4nc: R("R_UDPKAAT") - R("R_UDPKAAT_Rev") - R("R_ULA4NFT")  =
 0.

MudpLa4oc: R("R_UDPGDC") - R("R_UDPKAAT") + R("R_UDPKAAT_Rev")  = 0.

Mugmdac: -R("R_PAPPT3") + R("R_UGMDDS")  = 0.

Mugmdc: R("R_UAAGDS") - R("R_UGMDDS") + R("R_UM3PL") + R("R_UM4PCP")
  = 0.

MuLa4fnc: -R("R_UPLA4FNF") + R("R_UPLA4FNT")  = 0.

MuLa4nc: -R("R_ULA4Ntppi") + R("R_UPLA4FNF")  = 0.

MuLa4np: -R("R_LA4NTpp") + R("R_ULA4Ntppi")  = 0.

Mum4pc: -R("R_UM4PCP") + R("R_UM4PL")  = 0.

Mumpb: R("Mumpb_Transp") + R("R_EX_ump_e_") - R("R_EX_ump_e__Rev") 
 = 0.

Mumpc: R("R_ACGAMT") - R("R_NTD2") + R("R_NTPP8") + R("R_OMPDC") +
 R("R_PAPPT3") - R("R_UMPK") + R("R_UMPK_Rev") + R("R_UPPRT") +
 R("R_URIK2") + R("R_USHD")  = 0.

Mumpe: -R("R_EX_ump_e_") + R("R_EX_ump_e__Rev") - R("R_UMPtex") +
 R("R_UMPtex_Rev")  = 0.

Mumpp: -R("R_NTD2pp") + R("R_UACGALPpp") + R("R_UACGAMPpp") +
 R("R_UDPGALPpp") + R("R_UDPGPpp") + R("R_UGLCURPpp") + R("R_UMPtex") -
 R("R_UMPtex_Rev")  = 0.

Munagac: R("R_ACGAMT") - R("R_ACMAMUT")  = 0.

Munagamuc: -R("R_AADDGT") + R("R_ACMAMUT")  = 0.

Munagamufc: R("R_AADDGT") - R("R_ECAtpp")  = 0.

Munagamufp: -2 R("R_ECAP1pp") - R("R_ECAP2pp") - R("R_ECAP3pp") +
 R("R_ECAtpp")  = 0.

Muppg3c: -R("R_UPP3MT") + R("R_UPP3S") - R("R_UPPDC1")  = 0.

Murab: R("Murab_Transp") + R("R_EX_ura_e_") - R("R_EX_ura_e__Rev") 
 = 0.

Murac: R("R_CSND") + R("R_DURIPP") - R("R_DURIPP_Rev") +
 R("R_PYNP2r") - R("R_PYNP2r_Rev") - R("R_UPPRT") + R("R_URAt2pp") +
 R("R_URAt2rpp") - R("R_URAt2rpp_Rev") + R("R_URIH")  = 0.

Murae: -R("R_EX_ura_e_") + R("R_EX_ura_e__Rev") - R("R_URAtex") +
 R("R_URAtex_Rev")  = 0.

Murap: -R("R_URAt2pp") - R("R_URAt2rpp") + R("R_URAt2rpp_Rev") +
 R("R_URAtex") - R("R_URAtex_Rev")  = 0.

Muratec: -R("R_URIC") + R("R_XAND")  = 0.

Murdglycc: R("R_ALLTAMH") - R("R_UGLYCH") - R("R_URDGLYCD")  = 0.

Mureab: R("Mureab_Transp") + R("R_EX_urea_e_") -
 R("R_EX_urea_e__Rev")  = 0.

Mureac: R("R_AGMT") + R("R_UREAtpp") - R("R_UREAtpp_Rev")  = 0.

Mureae: -R("R_EX_urea_e_") + R("R_EX_urea_e__Rev") - R("R_UREAtex")
 + R("R_UREAtex_Rev")  = 0.

Mureap: R("R_UREAtex") - R("R_UREAtex_Rev") - R("R_UREAtpp") +
 R("R_UREAtpp_Rev")  = 0.

Murib: R("Murib_Transp") + R("R_EX_uri_e_") - R("R_EX_uri_e__Rev") 
 = 0.

Muric: R("R_CYTD") + R("R_NTD2") - R("R_PYNP2r") + R("R_PYNP2r_Rev")
 - R("R_URIH") - R("R_URIK2") + R("R_URIt2pp") + R("R_URIt2rpp") -
 R("R_URIt2rpp_Rev")  = 0.

Murie: -R("R_EX_uri_e_") + R("R_EX_uri_e__Rev") - R("R_URItex") +
 R("R_URItex_Rev")  = 0.

Murip: R("R_3NTD2pp") + R("R_NTD2pp") - R("R_URIt2pp") -
 R("R_URIt2rpp") + R("R_URIt2rpp_Rev") + R("R_URItex") - R("R_URItex_Rev")  = 0.

Mutpc: -R("R_CTPS2") - 0.1441 R("R_Ec_biomass_iAF1260_core_59p81M")
 - R("R_GALUi") + R("R_NDPK2") - R("R_NDPK2_Rev") - R("R_NTPP8") -
 R("R_RNTR4c") - R("R_UAGDP")  = 0.

MvalLb: R("MvalLb_Transp") + R("R_EX_val_L_e_") -
 R("R_EX_val_L_e__Rev")  = 0.

MvalLc: -0.4232 R("R_Ec_biomass_iAF1260_core_59p81M") +
 R("R_VALabcpp") + R("R_VALt2rpp") - R("R_VALt2rpp_Rev") - R("R_VALTA") +
 R("R_VALTA_Rev") - R("R_VALTRS") + R("R_VPAMT")  = 0.

MvalLe: -R("R_EX_val_L_e_") + R("R_EX_val_L_e__Rev") - R("R_VALtex")
 + R("R_VALtex_Rev")  = 0.

MvalLp: -R("R_VALabcpp") - R("R_VALt2rpp") + R("R_VALt2rpp_Rev") +
 R("R_VALtex") - R("R_VALtex_Rev")  = 0.

Mvaltrnac: R("R_VALTRS")  = 0.

Mxanb: R("Mxanb_Transp") + R("R_EX_xan_e_") - R("R_EX_xan_e__Rev") 
 = 0.

Mxanc: R("R_GUAD") + R("R_HXAND") + R("R_PUNP7") - R("R_PUNP7_Rev")
 - R("R_XAND") + R("R_XANt2pp") + R("R_XANtpp") - R("R_XANtpp_Rev") -
 R("R_XPPT") + R("R_XTSNH")  = 0.

Mxane: -R("R_EX_xan_e_") + R("R_EX_xan_e__Rev") - R("R_XANtex") +
 R("R_XANtex_Rev")  = 0.

Mxanp: -R("R_XANt2pp") + R("R_XANtex") - R("R_XANtex_Rev") -
 R("R_XANtpp") + R("R_XANtpp_Rev")  = 0.

Mxdpc: R("R_NTP12")  = 0.

Mxmpb: R("Mxmpb_Transp") + R("R_EX_xmp_e_") - R("R_EX_xmp_e__Rev") 
 = 0.

Mxmpc: -R("R_GMPS2") + R("R_IMPD") - R("R_NTD10") + R("R_NTPP11") +
 R("R_XPPT")  = 0.

Mxmpe: -R("R_EX_xmp_e_") + R("R_EX_xmp_e__Rev") - R("R_XMPtex") +
 R("R_XMPtex_Rev")  = 0.

Mxmpp: -R("R_NTD10pp") + R("R_XMPtex") - R("R_XMPtex_Rev")  = 0.

Mxtpc: R("R_GTPHs") - R("R_NTP12") - R("R_NTPP11")  = 0.

Mxtsnb: R("Mxtsnb_Transp") + R("R_EX_xtsn_e_") -
 R("R_EX_xtsn_e__Rev")  = 0.

Mxtsnc: R("R_NTD10") - R("R_PUNP7") + R("R_PUNP7_Rev") -
 R("R_XTSNH") + R("R_XTSNt2rpp") - R("R_XTSNt2rpp_Rev")  = 0.

Mxtsne: -R("R_EX_xtsn_e_") + R("R_EX_xtsn_e__Rev") - R("R_XTSNtex")
 + R("R_XTSNtex_Rev")  = 0.

Mxtsnp: R("R_NTD10pp") - R("R_XTSNt2rpp") + R("R_XTSNt2rpp_Rev") +
 R("R_XTSNtex") - R("R_XTSNtex_Rev")  = 0.

Mxu5pDc: R("R_RBP4E") - R("R_RBP4E_Rev") + R("R_RPE") -
 R("R_RPE_Rev") - R("R_TKT1") + R("R_TKT1_Rev") - R("R_TKT2") + R("R_TKT2_Rev") +
 R("R_XYLK")  = 0.

Mxu5pLc: R("R_KG6PDC") - R("R_X5PL3E") + R("R_XYLK2")  = 0.

MxylDb: R("MxylDb_Transp") + R("R_EX_xyl_D_e_") -
 R("R_EX_xyl_D_e__Rev")  = 0.

MxylDc: R("R_XYLabcpp") - R("R_XYLI1") + R("R_XYLI1_Rev") +
 R("R_XYLt2pp")  = 0.

MxylDe: -R("R_EX_xyl_D_e_") + R("R_EX_xyl_D_e__Rev") - R("R_XYLtex")
 + R("R_XYLtex_Rev")  = 0.

MxylDp: -R("R_XYLabcpp") - R("R_XYLt2pp") + R("R_XYLtex") -
 R("R_XYLtex_Rev")  = 0.

MxyluDc: R("R_XYLI1") - R("R_XYLI1_Rev") - R("R_XYLK")  = 0.

MxyluLb: R("MxyluLb_Transp") + R("R_EX_xylu_L_e_") -
 R("R_EX_xylu_L_e__Rev")  = 0.

MxyluLc: R("R_LYXI") - R("R_XYLK2") + R("R_XYLUt2pp")  = 0.

MxyluLe: -R("R_EX_xylu_L_e_") + R("R_EX_xylu_L_e__Rev") -
 R("R_XYLUtex") + R("R_XYLUtex_Rev")  = 0.

MxyluLp: -R("R_XYLUt2pp") + R("R_XYLUtex") - R("R_XYLUtex_Rev")  =
 0.

Mzn2b: R("Mzn2b_Transp") + R("R_EX_zn2_e_") - R("R_EX_zn2_e__Rev") 
 = 0.

Mzn2c: -0.003158 R("R_Ec_biomass_iAF1260_core_59p81M") -
 R("R_ZN2abcpp") - R("R_ZN2t3pp") + R("R_ZN2tpp") + R("R_ZNabcpp")  = 0.

Mzn2e: -R("R_EX_zn2_e_") + R("R_EX_zn2_e__Rev") - R("R_Zn2tex") +
 R("R_Zn2tex_Rev")  = 0.

Mzn2p: R("R_ZN2abcpp") + R("R_ZN2t3pp") + R("R_Zn2tex") -
 R("R_Zn2tex_Rev") - R("R_ZN2tpp") - R("R_ZNabcpp")  = 0.

Bounds

0 <= R("R_EX_12ppd_R_e_") <= 999999
0 <= R("R_EX_12ppd_R_e__Rev") <= 999999
0 <= R("R_12PPDRtex") <= 999999
0 <= R("R_12PPDRtex_Rev") <= 999999
0 <= R("R_12PPDRtpp") <= 999999
0 <= R("R_12PPDRtpp_Rev") <= 999999
0 <= R("R_EX_12ppd_S_e_") <= 999999
0 <= R("R_EX_12ppd_S_e__Rev") <= 999999
0 <= R("R_12PPDStex") <= 999999
0 <= R("R_12PPDStex_Rev") <= 999999
0 <= R("R_12PPDStpp") <= 999999
0 <= R("R_12PPDStpp_Rev") <= 999999
0 <= R("R_EX_14glucan_e_") <= 999999
0 <= R("R_EX_14glucan_e__Rev") <= 999999
0 <= R("R_EX_15dap_e_") <= 999999
0 <= R("R_EX_15dap_e__Rev") <= 999999
0 <= R("R_DAPtex") <= 999999
0 <= R("R_DAPtex_Rev") <= 999999
0 <= R("R_EX_23camp_e_") <= 999999
0 <= R("R_EX_23camp_e__Rev") <= 999999
0 <= R("R_23CAMPtex") <= 999999
0 <= R("R_23CAMPtex_Rev") <= 999999
0 <= R("R_EX_23ccmp_e_") <= 999999
0 <= R("R_EX_23ccmp_e__Rev") <= 999999
0 <= R("R_23CCMPtex") <= 999999
0 <= R("R_23CCMPtex_Rev") <= 999999
0 <= R("R_EX_23cgmp_e_") <= 999999
0 <= R("R_EX_23cgmp_e__Rev") <= 999999
0 <= R("R_23CGMPtex") <= 999999
0 <= R("R_23CGMPtex_Rev") <= 999999
0 <= R("R_EX_23cump_e_") <= 999999
0 <= R("R_EX_23cump_e__Rev") <= 999999
0 <= R("R_23CUMPtex") <= 999999
0 <= R("R_23CUMPtex_Rev") <= 999999
0 <= R("R_EX_23dappa_e_") <= 999999
0 <= R("R_EX_23dappa_e__Rev") <= 999999
0 <= R("R_23DAPPAtex") <= 999999
0 <= R("R_23DAPPAtex_Rev") <= 999999
0 <= R("R_ADSL2r") <= 999999
0 <= R("R_ADSL2r_Rev") <= 999999
0 <= R("R_DAPE") <= 999999
0 <= R("R_DAPE_Rev") <= 999999
0 <= R("R_EX_26dap_M_e_") <= 999999
0 <= R("R_EX_26dap_M_e__Rev") <= 999999
0 <= R("R_26DAHtex") <= 999999
0 <= R("R_26DAHtex_Rev") <= 999999
0 <= R("R_EX_2ddglcn_e_") <= 999999
0 <= R("R_EX_2ddglcn_e__Rev") <= 999999
0 <= R("R_DDGLCNtex") <= 999999
0 <= R("R_DDGLCNtex_Rev") <= 999999
0 <= R("R_DDPGALA") <= 999999
0 <= R("R_DDPGALA_Rev") <= 999999
0 <= R("R_PPM2") <= 999999
0 <= R("R_PPM2_Rev") <= 999999
0 <= R("R_PGM") <= 999999
0 <= R("R_PGM_Rev") <= 999999
0 <= R("R_ENO") <= 999999
0 <= R("R_ENO_Rev") <= 999999
0 <= R("R_EX_34dhpac_e_") <= 999999
0 <= R("R_EX_34dhpac_e__Rev") <= 999999
0 <= R("R_34dhpactex") <= 999999
0 <= R("R_34dhpactex_Rev") <= 999999
0 <= R("R_EX_3amp_e_") <= 999999
0 <= R("R_EX_3amp_e__Rev") <= 999999
0 <= R("R_3AMPtex") <= 999999
0 <= R("R_3AMPtex_Rev") <= 999999
0 <= R("R_IPPMIa") <= 999999
0 <= R("R_IPPMIa_Rev") <= 999999
0 <= R("R_EX_3cmp_e_") <= 999999
0 <= R("R_EX_3cmp_e__Rev") <= 999999
0 <= R("R_3CMPtex") <= 999999
0 <= R("R_3CMPtex_Rev") <= 999999
0 <= R("R_EX_3gmp_e_") <= 999999
0 <= R("R_EX_3gmp_e__Rev") <= 999999
0 <= R("R_3GMPtex") <= 999999
0 <= R("R_3GMPtex_Rev") <= 999999
0 <= R("R_ECOAH1") <= 999999
0 <= R("R_ECOAH1_Rev") <= 999999
0 <= R("R_EX_3hcinnm_e_") <= 999999
0 <= R("R_EX_3hcinnm_e__Rev") <= 999999
0 <= R("R_HCINNMtex") <= 999999
0 <= R("R_HCINNMtex_Rev") <= 999999
0 <= R("R_ECOAH4") <= 999999
0 <= R("R_ECOAH4_Rev") <= 999999
0 <= R("R_ECOAH5") <= 999999
0 <= R("R_ECOAH5_Rev") <= 999999
0 <= R("R_ECOAH2") <= 999999
0 <= R("R_ECOAH2_Rev") <= 999999
0 <= R("R_ECOAH7") <= 999999
0 <= R("R_ECOAH7_Rev") <= 999999
0 <= R("R_ECOAH3") <= 999999
0 <= R("R_ECOAH3_Rev") <= 999999
0 <= R("R_ECOAH8") <= 999999
0 <= R("R_ECOAH8_Rev") <= 999999
0 <= R("R_EX_3hpppn_e_") <= 999999
0 <= R("R_EX_3hpppn_e__Rev") <= 999999
0 <= R("R_HPPPNtex") <= 999999
0 <= R("R_HPPPNtex_Rev") <= 999999
0 <= R("R_ECOAH6") <= 999999
0 <= R("R_ECOAH6_Rev") <= 999999
0 <= R("R_EX_3ump_e_") <= 999999
0 <= R("R_EX_3ump_e__Rev") <= 999999
0 <= R("R_3UMPtex") <= 999999
0 <= R("R_3UMPtex_Rev") <= 999999
0 <= R("R_EX_4abut_e_") <= 999999
0 <= R("R_EX_4abut_e__Rev") <= 999999
0 <= R("R_ABUTtex") <= 999999
0 <= R("R_ABUTtex_Rev") <= 999999
0 <= R("R_EX_4hoxpacd_e_") <= 999999
0 <= R("R_EX_4hoxpacd_e__Rev") <= 999999
0 <= R("R_4HOXPACDtex") <= 999999
0 <= R("R_4HOXPACDtex_Rev") <= 999999
0 <= R("R_AIRC3") <= 999999
0 <= R("R_AIRC3_Rev") <= 999999
0 <= R("R_EX_5dglcn_e_") <= 999999
0 <= R("R_EX_5dglcn_e__Rev") <= 999999
0 <= R("R_5DGLCNtex") <= 999999
0 <= R("R_5DGLCNtex_Rev") <= 999999
0 <= R("R_MTRI") <= 999999
0 <= R("R_MTRI_Rev") <= 999999
0 <= R("R_EX_acac_e_") <= 999999
0 <= R("R_EX_acac_e__Rev") <= 999999
0 <= R("R_ACACtex") <= 999999
0 <= R("R_ACACtex_Rev") <= 999999
0 <= R("R_EX_acald_e_") <= 999999
0 <= R("R_EX_acald_e__Rev") <= 999999
0 <= R("R_ACALDtex") <= 999999
0 <= R("R_ACALDtex_Rev") <= 999999
0 <= R("R_ACALDtpp") <= 999999
0 <= R("R_ACALDtpp_Rev") <= 999999
0 <= R("R_ACACT1r") <= 999999
0 <= R("R_ACACT1r_Rev") <= 999999
0 <= R("R_EX_ac_e_") <= 999999
0 <= R("R_EX_ac_e__Rev") <= 999999
0 <= R("R_ACtex") <= 999999
0 <= R("R_ACtex_Rev") <= 999999
0 <= R("R_EX_acgal1p_e_") <= 999999
0 <= R("R_EX_acgal1p_e__Rev") <= 999999
0 <= R("R_ACGAL1Ptex") <= 999999
0 <= R("R_ACGAL1Ptex_Rev") <= 999999
0 <= R("R_EX_acgal_e_") <= 999999
0 <= R("R_EX_acgal_e__Rev") <= 999999
0 <= R("R_ACGALtex") <= 999999
0 <= R("R_ACGALtex_Rev") <= 999999
0 <= R("R_EX_acgam1p_e_") <= 999999
0 <= R("R_EX_acgam1p_e__Rev") <= 999999
0 <= R("R_ACGAM1Ptex") <= 999999
0 <= R("R_ACGAM1Ptex_Rev") <= 999999
0 <= R("R_EX_acgam_e_") <= 999999
0 <= R("R_EX_acgam_e__Rev") <= 999999
0 <= R("R_ACGAtex") <= 999999
0 <= R("R_ACGAtex_Rev") <= 999999
0 <= R("R_EX_acmana_e_") <= 999999
0 <= R("R_EX_acmana_e__Rev") <= 999999
0 <= R("R_ACMANAtex") <= 999999
0 <= R("R_ACMANAtex_Rev") <= 999999
0 <= R("R_AMANAPEr") <= 999999
0 <= R("R_AMANAPEr_Rev") <= 999999
0 <= R("R_EX_acmum_e_") <= 999999
0 <= R("R_EX_acmum_e__Rev") <= 999999
0 <= R("R_EX_acnam_e_") <= 999999
0 <= R("R_EX_acnam_e__Rev") <= 999999
0 <= R("R_ACNAMtex") <= 999999
0 <= R("R_ACNAMtex_Rev") <= 999999
0 <= R("R_EX_acolipa_e_") <= 999999
0 <= R("R_EX_acolipa_e__Rev") <= 999999
0 <= R("R_ACONIs") <= 999999
0 <= R("R_ACONIs_Rev") <= 999999
0 <= R("R_ACOATA") <= 999999
0 <= R("R_ACOATA_Rev") <= 999999
0 <= R("R_EX_acser_e_") <= 999999
0 <= R("R_EX_acser_e__Rev") <= 999999
0 <= R("R_ACSERtex") <= 999999
0 <= R("R_ACSERtex_Rev") <= 999999
0 <= R("R_EX_ade_e_") <= 999999
0 <= R("R_EX_ade_e__Rev") <= 999999
0 <= R("R_ADEtex") <= 999999
0 <= R("R_ADEtex_Rev") <= 999999
0 <= R("R_EX_adn_e_") <= 999999
0 <= R("R_EX_adn_e__Rev") <= 999999
0 <= R("R_ADNtex") <= 999999
0 <= R("R_ADNtex_Rev") <= 999999
0 <= R("R_EX_adocbl_e_") <= 999999
0 <= R("R_EX_adocbl_e__Rev") <= 999999
0 <= R("R_EX_ag_e_") <= 999999
0 <= R("R_EX_ag_e__Rev") <= 999999
0 <= R("R_EX_agm_e_") <= 999999
0 <= R("R_EX_agm_e__Rev") <= 999999
0 <= R("R_AGMtex") <= 999999
0 <= R("R_AGMtex_Rev") <= 999999
0 <= R("R_DHPTPE") <= 999999
0 <= R("R_DHPTPE_Rev") <= 999999
0 <= R("R_AICART") <= 999999
0 <= R("R_AICART_Rev") <= 999999
0 <= R("R_ACOTA") <= 999999
0 <= R("R_ACOTA_Rev") <= 999999
0 <= R("R_EX_akg_e_") <= 999999
0 <= R("R_EX_akg_e__Rev") <= 999999
0 <= R("R_AKGtex") <= 999999
0 <= R("R_AKGtex_Rev") <= 999999
0 <= R("R_EX_alaala_e_") <= 999999
0 <= R("R_EX_alaala_e__Rev") <= 999999
0 <= R("R_ALAALAtex") <= 999999
0 <= R("R_ALAALAtex_Rev") <= 999999
0 <= R("R_EX_ala_B_e_") <= 999999
0 <= R("R_EX_ala_B_e__Rev") <= 999999
0 <= R("R_BALAtex") <= 999999
0 <= R("R_BALAtex_Rev") <= 999999
0 <= R("R_EX_ala_D_e_") <= 999999
0 <= R("R_EX_ala_D_e__Rev") <= 999999
0 <= R("R_DALAtex") <= 999999
0 <= R("R_DALAtex_Rev") <= 999999
0 <= R("R_ALAR") <= 999999
0 <= R("R_ALAR_Rev") <= 999999
0 <= R("R_ALATA_L") <= 999999
0 <= R("R_ALATA_L_Rev") <= 999999
0 <= R("R_EX_ala_L_e_") <= 999999
0 <= R("R_EX_ala_L_e__Rev") <= 999999
0 <= R("R_ALAtex") <= 999999
0 <= R("R_ALAtex_Rev") <= 999999
0 <= R("R_ALLPI") <= 999999
0 <= R("R_ALLPI_Rev") <= 999999
0 <= R("R_EX_all_D_e_") <= 999999
0 <= R("R_EX_all_D_e__Rev") <= 999999
0 <= R("R_ALLtex") <= 999999
0 <= R("R_ALLtex_Rev") <= 999999
0 <= R("R_EX_alltn_e_") <= 999999
0 <= R("R_EX_alltn_e__Rev") <= 999999
0 <= R("R_ALLTNtex") <= 999999
0 <= R("R_ALLTNtex_Rev") <= 999999
0 <= R("R_ALLULPE") <= 999999
0 <= R("R_ALLULPE_Rev") <= 999999
0 <= R("R_AMAOTr") <= 999999
0 <= R("R_AMAOTr_Rev") <= 999999
0 <= R("R_EX_amp_e_") <= 999999
0 <= R("R_EX_amp_e__Rev") <= 999999
0 <= R("R_AMPtex") <= 999999
0 <= R("R_AMPtex_Rev") <= 999999
0 <= R("R_EX_anhgm_e_") <= 999999
0 <= R("R_EX_anhgm_e__Rev") <= 999999
0 <= R("R_ANHGMtex") <= 999999
0 <= R("R_ANHGMtex_Rev") <= 999999
0 <= R("R_ARAI") <= 999999
0 <= R("R_ARAI_Rev") <= 999999
0 <= R("R_EX_arab_L_e_") <= 999999
0 <= R("R_EX_arab_L_e__Rev") <= 999999
0 <= R("R_ARBtex") <= 999999
0 <= R("R_ARBtex_Rev") <= 999999
0 <= R("R_EX_arbtn_e_") <= 999999
0 <= R("R_EX_arbtn_e__Rev") <= 999999
0 <= R("R_EX_arbtn_fe3_e_") <= 999999
0 <= R("R_EX_arbtn_fe3_e__Rev") <= 999999
0 <= R("R_EX_arg_L_e_") <= 999999
0 <= R("R_EX_arg_L_e__Rev") <= 999999
0 <= R("R_ARGtex") <= 999999
0 <= R("R_ARGtex_Rev") <= 999999
0 <= R("R_ARGAGMt7pp") <= 999999
0 <= R("R_ARGAGMt7pp_Rev") <= 999999
0 <= R("R_ARGSL") <= 999999
0 <= R("R_ARGSL_Rev") <= 999999
0 <= R("R_EX_ascb_L_e_") <= 999999
0 <= R("R_EX_ascb_L_e__Rev") <= 999999
0 <= R("R_ASCBtex") <= 999999
0 <= R("R_ASCBtex_Rev") <= 999999
0 <= R("R_EX_asn_L_e_") <= 999999
0 <= R("R_EX_asn_L_e__Rev") <= 999999
0 <= R("R_ASNtex") <= 999999
0 <= R("R_ASNtex_Rev") <= 999999
0 <= R("R_EX_aso3_e_") <= 999999
0 <= R("R_EX_aso3_e__Rev") <= 999999
0 <= R("R_ASO3tex") <= 999999
0 <= R("R_ASO3tex_Rev") <= 999999
0 <= R("R_ASPTA") <= 999999
0 <= R("R_ASPTA_Rev") <= 999999
0 <= R("R_EX_asp_L_e_") <= 999999
0 <= R("R_EX_asp_L_e__Rev") <= 999999
0 <= R("R_ASPtex") <= 999999
0 <= R("R_ASPtex_Rev") <= 999999
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
0 <= R("R_EX_but_e_") <= 999999
0 <= R("R_EX_but_e__Rev") <= 999999
0 <= R("R_BUTtex") <= 999999
0 <= R("R_BUTtex_Rev") <= 999999
0 <= R("R_EX_butso3_e_") <= 999999
0 <= R("R_EX_butso3_e__Rev") <= 999999
0 <= R("R_BUTSO3tex") <= 999999
0 <= R("R_BUTSO3tex_Rev") <= 999999
0 <= R("R_EX_ca2_e_") <= 999999
0 <= R("R_EX_ca2_e__Rev") <= 999999
0 <= R("R_CA2tex") <= 999999
0 <= R("R_CA2tex_Rev") <= 999999
0 <= R("R_EX_cbi_e_") <= 999999
0 <= R("R_EX_cbi_e__Rev") <= 999999
0 <= R("R_EX_cbl1_e_") <= 999999
0 <= R("R_EX_cbl1_e__Rev") <= 999999
0 <= R("R_EX_cd2_e_") <= 999999
0 <= R("R_EX_cd2_e__Rev") <= 999999
0 <= R("R_CD2tex") <= 999999
0 <= R("R_CD2tex_Rev") <= 999999
0 <= R("R_NDPK3") <= 999999
0 <= R("R_NDPK3_Rev") <= 999999
0 <= R("R_EX_cgly_e_") <= 999999
0 <= R("R_EX_cgly_e__Rev") <= 999999
0 <= R("R_CGLYtex") <= 999999
0 <= R("R_CGLYtex_Rev") <= 999999
0 <= R("R_EX_chol_e_") <= 999999
0 <= R("R_EX_chol_e__Rev") <= 999999
0 <= R("R_CHLtex") <= 999999
0 <= R("R_CHLtex_Rev") <= 999999
0 <= R("R_ICHORS") <= 999999
0 <= R("R_ICHORS_Rev") <= 999999
0 <= R("R_ACONTa") <= 999999
0 <= R("R_ACONTa_Rev") <= 999999
0 <= R("R_EX_cit_e_") <= 999999
0 <= R("R_EX_cit_e__Rev") <= 999999
0 <= R("R_CITtex") <= 999999
0 <= R("R_CITtex_Rev") <= 999999
0 <= R("R_EX_cl_e_") <= 999999
0 <= R("R_EX_cl_e__Rev") <= 999999
0 <= R("R_CLtex") <= 999999
0 <= R("R_CLtex_Rev") <= 999999
0 <= R("R_CYTK1") <= 999999
0 <= R("R_CYTK1_Rev") <= 999999
0 <= R("R_EX_cmp_e_") <= 999999
0 <= R("R_EX_cmp_e__Rev") <= 999999
0 <= R("R_CMPtex") <= 999999
0 <= R("R_CMPtex_Rev") <= 999999
0 <= R("R_EX_co2_e_") <= 999999
0 <= R("R_EX_co2_e__Rev") <= 999999
0 <= R("R_CO2tex") <= 999999
0 <= R("R_CO2tex_Rev") <= 999999
0 <= R("R_CO2tpp") <= 999999
0 <= R("R_CO2tpp_Rev") <= 999999
0 <= R("R_EX_cobalt2_e_") <= 999999
0 <= R("R_EX_cobalt2_e__Rev") <= 999999
0 <= R("R_COBALT2tex") <= 999999
0 <= R("R_COBALT2tex_Rev") <= 999999
0 <= R("R_EX_colipa_e_") <= 999999
0 <= R("R_EX_colipa_e__Rev") <= 999999
0 <= R("R_EX_cpgn_e_") <= 999999
0 <= R("R_EX_cpgn_e__Rev") <= 999999
0 <= R("R_EX_cpgn_un_e_") <= 999999
0 <= R("R_EX_cpgn_un_e__Rev") <= 999999
0 <= R("R_CRNBTCT") <= 999999
0 <= R("R_CRNBTCT_Rev") <= 999999
0 <= R("R_CRNCAR") <= 999999
0 <= R("R_CRNCAR_Rev") <= 999999
0 <= R("R_CRNCDH") <= 999999
0 <= R("R_CRNCDH_Rev") <= 999999
0 <= R("R_EX_crn_e_") <= 999999
0 <= R("R_EX_crn_e__Rev") <= 999999
0 <= R("R_CRNtex") <= 999999
0 <= R("R_CRNtex_Rev") <= 999999
0 <= R("R_EX_csn_e_") <= 999999
0 <= R("R_EX_csn_e__Rev") <= 999999
0 <= R("R_CSNtex") <= 999999
0 <= R("R_CSNtex_Rev") <= 999999
0 <= R("R_CRNCBCT") <= 999999
0 <= R("R_CRNCBCT_Rev") <= 999999
0 <= R("R_EX_cu2_e_") <= 999999
0 <= R("R_EX_cu2_e__Rev") <= 999999
0 <= R("R_CU2tex") <= 999999
0 <= R("R_CU2tex_Rev") <= 999999
0 <= R("R_EX_cu_e_") <= 999999
0 <= R("R_EX_cu_e__Rev") <= 999999
0 <= R("R_CUtex") <= 999999
0 <= R("R_CUtex_Rev") <= 999999
0 <= R("R_EX_cyan_e_") <= 999999
0 <= R("R_EX_cyan_e__Rev") <= 999999
0 <= R("R_CYANtex") <= 999999
0 <= R("R_CYANtex_Rev") <= 999999
0 <= R("R_EX_cynt_e_") <= 999999
0 <= R("R_EX_cynt_e__Rev") <= 999999
0 <= R("R_CYNTtex") <= 999999
0 <= R("R_CYNTtex_Rev") <= 999999
0 <= R("R_EX_cys_D_e_") <= 999999
0 <= R("R_EX_cys_D_e__Rev") <= 999999
0 <= R("R_CYSDtex") <= 999999
0 <= R("R_CYSDtex_Rev") <= 999999
0 <= R("R_EX_cys_L_e_") <= 999999
0 <= R("R_EX_cys_L_e__Rev") <= 999999
0 <= R("R_CYStex") <= 999999
0 <= R("R_CYStex_Rev") <= 999999
0 <= R("R_EX_cytd_e_") <= 999999
0 <= R("R_EX_cytd_e__Rev") <= 999999
0 <= R("R_CYTDtex") <= 999999
0 <= R("R_CYTDtex_Rev") <= 999999
0 <= R("R_EX_dad_2_e_") <= 999999
0 <= R("R_EX_dad_2_e__Rev") <= 999999
0 <= R("R_DADNtex") <= 999999
0 <= R("R_DADNtex_Rev") <= 999999
0 <= R("R_NDPK8") <= 999999
0 <= R("R_NDPK8_Rev") <= 999999
0 <= R("R_DADK") <= 999999
0 <= R("R_DADK_Rev") <= 999999
0 <= R("R_EX_damp_e_") <= 999999
0 <= R("R_EX_damp_e__Rev") <= 999999
0 <= R("R_DAMPtex") <= 999999
0 <= R("R_DAMPtex_Rev") <= 999999
0 <= R("R_EX_dca_e_") <= 999999
0 <= R("R_EX_dca_e__Rev") <= 999999
0 <= R("R_DCAtex") <= 999999
0 <= R("R_DCAtex_Rev") <= 999999
0 <= R("R_ADSL1r") <= 999999
0 <= R("R_ADSL1r_Rev") <= 999999
0 <= R("R_NDPK7") <= 999999
0 <= R("R_NDPK7_Rev") <= 999999
0 <= R("R_CYTK2") <= 999999
0 <= R("R_CYTK2_Rev") <= 999999
0 <= R("R_EX_dcmp_e_") <= 999999
0 <= R("R_EX_dcmp_e__Rev") <= 999999
0 <= R("R_DCMPtex") <= 999999
0 <= R("R_DCMPtex_Rev") <= 999999
0 <= R("R_EX_dcyt_e_") <= 999999
0 <= R("R_EX_dcyt_e__Rev") <= 999999
0 <= R("R_DCYTtex") <= 999999
0 <= R("R_DCYTtex_Rev") <= 999999
0 <= R("R_EX_ddca_e_") <= 999999
0 <= R("R_EX_ddca_e__Rev") <= 999999
0 <= R("R_NDPK5") <= 999999
0 <= R("R_NDPK5_Rev") <= 999999
0 <= R("R_DGK1") <= 999999
0 <= R("R_DGK1_Rev") <= 999999
0 <= R("R_EX_dgmp_e_") <= 999999
0 <= R("R_EX_dgmp_e__Rev") <= 999999
0 <= R("R_DGMPtex") <= 999999
0 <= R("R_DGMPtex_Rev") <= 999999
0 <= R("R_EX_dgsn_e_") <= 999999
0 <= R("R_EX_dgsn_e__Rev") <= 999999
0 <= R("R_DGSNtex") <= 999999
0 <= R("R_DGSNtex_Rev") <= 999999
0 <= R("R_EX_dha_e_") <= 999999
0 <= R("R_EX_dha_e__Rev") <= 999999
0 <= R("R_DHAtex") <= 999999
0 <= R("R_DHAtex_Rev") <= 999999
0 <= R("R_DHAtpp") <= 999999
0 <= R("R_DHAtpp_Rev") <= 999999
0 <= R("R_TPI") <= 999999
0 <= R("R_TPI_Rev") <= 999999
0 <= R("R_EX_dimp_e_") <= 999999
0 <= R("R_EX_dimp_e__Rev") <= 999999
0 <= R("R_DIMPtex") <= 999999
0 <= R("R_DIMPtex_Rev") <= 999999
0 <= R("R_EX_din_e_") <= 999999
0 <= R("R_EX_din_e__Rev") <= 999999
0 <= R("R_DINStex") <= 999999
0 <= R("R_DINStex_Rev") <= 999999
0 <= R("R_EX_dms_e_") <= 999999
0 <= R("R_EX_dms_e__Rev") <= 999999
0 <= R("R_DMStex") <= 999999
0 <= R("R_DMStex_Rev") <= 999999
0 <= R("R_EX_dmso_e_") <= 999999
0 <= R("R_EX_dmso_e__Rev") <= 999999
0 <= R("R_DMSOtex") <= 999999
0 <= R("R_DMSOtex_Rev") <= 999999
0 <= R("R_EX_dopa_e_") <= 999999
0 <= R("R_EX_dopa_e__Rev") <= 999999
0 <= R("R_DOPAtex") <= 999999
0 <= R("R_DOPAtex_Rev") <= 999999
0 <= R("R_NDPK4") <= 999999
0 <= R("R_NDPK4_Rev") <= 999999
0 <= R("R_DTMPK") <= 999999
0 <= R("R_DTMPK_Rev") <= 999999
0 <= R("R_EX_dtmp_e_") <= 999999
0 <= R("R_EX_dtmp_e__Rev") <= 999999
0 <= R("R_DTMPtex") <= 999999
0 <= R("R_DTMPtex_Rev") <= 999999
0 <= R("R_NDPK6") <= 999999
0 <= R("R_NDPK6_Rev") <= 999999
0 <= R("R_URIDK2r") <= 999999
0 <= R("R_URIDK2r_Rev") <= 999999
0 <= R("R_EX_dump_e_") <= 999999
0 <= R("R_EX_dump_e__Rev") <= 999999
0 <= R("R_DUMPtex") <= 999999
0 <= R("R_DUMPtex_Rev") <= 999999
0 <= R("R_EX_duri_e_") <= 999999
0 <= R("R_EX_duri_e__Rev") <= 999999
0 <= R("R_DURItex") <= 999999
0 <= R("R_DURItex_Rev") <= 999999
0 <= R("R_EX_eca4colipa_e_") <= 999999
0 <= R("R_EX_eca4colipa_e__Rev") <= 999999
0 <= R("R_EX_enlipa_e_") <= 999999
0 <= R("R_EX_enlipa_e__Rev") <= 999999
0 <= R("R_EX_enter_e_") <= 999999
0 <= R("R_EX_enter_e__Rev") <= 999999
0 <= R("R_EX_etha_e_") <= 999999
0 <= R("R_EX_etha_e__Rev") <= 999999
0 <= R("R_ETHAtex") <= 999999
0 <= R("R_ETHAtex_Rev") <= 999999
0 <= R("R_EX_ethso3_e_") <= 999999
0 <= R("R_EX_ethso3_e__Rev") <= 999999
0 <= R("R_ETHSO3tex") <= 999999
0 <= R("R_ETHSO3tex_Rev") <= 999999
0 <= R("R_EX_etoh_e_") <= 999999
0 <= R("R_EX_etoh_e__Rev") <= 999999
0 <= R("R_ETOHtex") <= 999999
0 <= R("R_ETOHtex_Rev") <= 999999
0 <= R("R_F6PA") <= 999999
0 <= R("R_F6PA_Rev") <= 999999
0 <= R("R_EX_f6p_e_") <= 999999
0 <= R("R_EX_f6p_e__Rev") <= 999999
0 <= R("R_F6Ptex") <= 999999
0 <= R("R_F6Ptex_Rev") <= 999999
0 <= R("R_ACOAD1f") <= 999999
0 <= R("R_ACOAD1f_Rev") <= 999999
0 <= R("R_ACOAD4f") <= 999999
0 <= R("R_ACOAD4f_Rev") <= 999999
0 <= R("R_ACOAD5f") <= 999999
0 <= R("R_ACOAD5f_Rev") <= 999999
0 <= R("R_EX_fald_e_") <= 999999
0 <= R("R_EX_fald_e__Rev") <= 999999
0 <= R("R_FALDtex") <= 999999
0 <= R("R_FALDtex_Rev") <= 999999
0 <= R("R_FALDtpp") <= 999999
0 <= R("R_FALDtpp_Rev") <= 999999
0 <= R("R_FCLPA") <= 999999
0 <= R("R_FCLPA_Rev") <= 999999
0 <= R("R_FBA") <= 999999
0 <= R("R_FBA_Rev") <= 999999
0 <= R("R_EX_fe2_e_") <= 999999
0 <= R("R_EX_fe2_e__Rev") <= 999999
0 <= R("R_FE2tex") <= 999999
0 <= R("R_FE2tex_Rev") <= 999999
0 <= R("R_EX_fe3dcit_e_") <= 999999
0 <= R("R_EX_fe3dcit_e__Rev") <= 999999
0 <= R("R_EX_fe3dhbzs_e_") <= 999999
0 <= R("R_EX_fe3dhbzs_e__Rev") <= 999999
0 <= R("R_EX_fe3_e_") <= 999999
0 <= R("R_EX_fe3_e__Rev") <= 999999
0 <= R("R_FE3tex") <= 999999
0 <= R("R_FE3tex_Rev") <= 999999
0 <= R("R_EX_fe3hox_e_") <= 999999
0 <= R("R_EX_fe3hox_e__Rev") <= 999999
0 <= R("R_EX_fe3hox_un_e_") <= 999999
0 <= R("R_EX_fe3hox_un_e__Rev") <= 999999
0 <= R("R_EX_fecrm_e_") <= 999999
0 <= R("R_EX_fecrm_e__Rev") <= 999999
0 <= R("R_EX_fecrm_un_e_") <= 999999
0 <= R("R_EX_fecrm_un_e__Rev") <= 999999
0 <= R("R_EX_feenter_e_") <= 999999
0 <= R("R_EX_feenter_e__Rev") <= 999999
0 <= R("R_EX_feoxam_e_") <= 999999
0 <= R("R_EX_feoxam_e__Rev") <= 999999
0 <= R("R_EX_feoxam_un_e_") <= 999999
0 <= R("R_EX_feoxam_un_e__Rev") <= 999999
0 <= R("R_EX_for_e_") <= 999999
0 <= R("R_EX_for_e__Rev") <= 999999
0 <= R("R_FORtex") <= 999999
0 <= R("R_FORtex_Rev") <= 999999
0 <= R("R_EX_fru_e_") <= 999999
0 <= R("R_EX_fru_e__Rev") <= 999999
0 <= R("R_FRUtex") <= 999999
0 <= R("R_FRUtex_Rev") <= 999999
0 <= R("R_EX_frulys_e_") <= 999999
0 <= R("R_EX_frulys_e__Rev") <= 999999
0 <= R("R_FRULYStex") <= 999999
0 <= R("R_FRULYStex_Rev") <= 999999
0 <= R("R_EX_fruur_e_") <= 999999
0 <= R("R_EX_fruur_e__Rev") <= 999999
0 <= R("R_FRUURtex") <= 999999
0 <= R("R_FRUURtex_Rev") <= 999999
0 <= R("R_FCI") <= 999999
0 <= R("R_FCI_Rev") <= 999999
0 <= R("R_EX_fuc_L_e_") <= 999999
0 <= R("R_EX_fuc_L_e__Rev") <= 999999
0 <= R("R_FUCtex") <= 999999
0 <= R("R_FUCtex_Rev") <= 999999
0 <= R("R_EX_fum_e_") <= 999999
0 <= R("R_EX_fum_e__Rev") <= 999999
0 <= R("R_FUMtex") <= 999999
0 <= R("R_FUMtex_Rev") <= 999999
0 <= R("R_PGMT") <= 999999
0 <= R("R_PGMT_Rev") <= 999999
0 <= R("R_EX_g1p_e_") <= 999999
0 <= R("R_EX_g1p_e__Rev") <= 999999
0 <= R("R_G1Ptex") <= 999999
0 <= R("R_G1Ptex_Rev") <= 999999
0 <= R("R_EX_g3pc_e_") <= 999999
0 <= R("R_EX_g3pc_e__Rev") <= 999999
0 <= R("R_G3PCtex") <= 999999
0 <= R("R_G3PCtex_Rev") <= 999999
0 <= R("R_EX_g3pe_e_") <= 999999
0 <= R("R_EX_g3pe_e__Rev") <= 999999
0 <= R("R_G3PEtex") <= 999999
0 <= R("R_G3PEtex_Rev") <= 999999
0 <= R("R_EX_g3pg_e_") <= 999999
0 <= R("R_EX_g3pg_e__Rev") <= 999999
0 <= R("R_G3PGtex") <= 999999
0 <= R("R_G3PGtex_Rev") <= 999999
0 <= R("R_EX_g3pi_e_") <= 999999
0 <= R("R_EX_g3pi_e__Rev") <= 999999
0 <= R("R_G3PItex") <= 999999
0 <= R("R_G3PItex_Rev") <= 999999
0 <= R("R_EX_g3ps_e_") <= 999999
0 <= R("R_EX_g3ps_e__Rev") <= 999999
0 <= R("R_G3PStex") <= 999999
0 <= R("R_G3PStex_Rev") <= 999999
0 <= R("R_PGI") <= 999999
0 <= R("R_PGI_Rev") <= 999999
0 <= R("R_EX_g6p_e_") <= 999999
0 <= R("R_EX_g6p_e__Rev") <= 999999
0 <= R("R_G6Ptex") <= 999999
0 <= R("R_G6Ptex_Rev") <= 999999
0 <= R("R_EX_gal1p_e_") <= 999999
0 <= R("R_EX_gal1p_e__Rev") <= 999999
0 <= R("R_GAL1Ptex") <= 999999
0 <= R("R_GAL1Ptex_Rev") <= 999999
0 <= R("R_EX_gal_bD_e_") <= 999999
0 <= R("R_EX_gal_bD_e__Rev") <= 999999
0 <= R("R_GALBDtex") <= 999999
0 <= R("R_GALBDtex_Rev") <= 999999
0 <= R("R_GALKr") <= 999999
0 <= R("R_GALKr_Rev") <= 999999
0 <= R("R_EX_galct_D_e_") <= 999999
0 <= R("R_EX_galct_D_e__Rev") <= 999999
0 <= R("R_GALCTtex") <= 999999
0 <= R("R_GALCTtex_Rev") <= 999999
0 <= R("R_EX_galctn_D_e_") <= 999999
0 <= R("R_EX_galctn_D_e__Rev") <= 999999
0 <= R("R_GALCTNtex") <= 999999
0 <= R("R_GALCTNtex_Rev") <= 999999
0 <= R("R_EX_galctn_L_e_") <= 999999
0 <= R("R_EX_galctn_L_e__Rev") <= 999999
0 <= R("R_GALCTNLtex") <= 999999
0 <= R("R_GALCTNLtex_Rev") <= 999999
0 <= R("R_EX_gal_e_") <= 999999
0 <= R("R_EX_gal_e__Rev") <= 999999
0 <= R("R_GALtex") <= 999999
0 <= R("R_GALtex_Rev") <= 999999
0 <= R("R_EX_galt_e_") <= 999999
0 <= R("R_EX_galt_e__Rev") <= 999999
0 <= R("R_GALTtex") <= 999999
0 <= R("R_GALTtex_Rev") <= 999999
0 <= R("R_GUI2") <= 999999
0 <= R("R_GUI2_Rev") <= 999999
0 <= R("R_EX_galur_e_") <= 999999
0 <= R("R_EX_galur_e__Rev") <= 999999
0 <= R("R_GALURtex") <= 999999
0 <= R("R_GALURtex_Rev") <= 999999
0 <= R("R_PGAMT") <= 999999
0 <= R("R_PGAMT_Rev") <= 999999
0 <= R("R_EX_gam6p_e_") <= 999999
0 <= R("R_EX_gam6p_e__Rev") <= 999999
0 <= R("R_GAMAN6Ptex") <= 999999
0 <= R("R_GAMAN6Ptex_Rev") <= 999999
0 <= R("R_EX_gam_e_") <= 999999
0 <= R("R_EX_gam_e__Rev") <= 999999
0 <= R("R_GAMtex") <= 999999
0 <= R("R_GAMtex_Rev") <= 999999
0 <= R("R_GARFT") <= 999999
0 <= R("R_GARFT_Rev") <= 999999
0 <= R("R_EX_gbbtn_e_") <= 999999
0 <= R("R_EX_gbbtn_e__Rev") <= 999999
0 <= R("R_GBBTNtex") <= 999999
0 <= R("R_GBBTNtex_Rev") <= 999999
0 <= R("R_NDPK1") <= 999999
0 <= R("R_NDPK1_Rev") <= 999999
0 <= R("R_EX_gdp_e_") <= 999999
0 <= R("R_EX_gdp_e__Rev") <= 999999
0 <= R("R_GDPtex") <= 999999
0 <= R("R_GDPtex_Rev") <= 999999
0 <= R("R_XYLI2") <= 999999
0 <= R("R_XYLI2_Rev") <= 999999
0 <= R("R_GLCATr") <= 999999
0 <= R("R_GLCATr_Rev") <= 999999
0 <= R("R_EX_glc_e_") <= 999999
0 <= R("R_EX_glc_e__Rev") <= 999999
0 <= R("R_GLCtex") <= 999999
0 <= R("R_GLCtex_Rev") <= 999999
0 <= R("R_EX_glcn_e_") <= 999999
0 <= R("R_EX_glcn_e__Rev") <= 999999
0 <= R("R_GLCNtex") <= 999999
0 <= R("R_GLCNtex_Rev") <= 999999
0 <= R("R_EX_glcr_e_") <= 999999
0 <= R("R_EX_glcr_e__Rev") <= 999999
0 <= R("R_GLCRtex") <= 999999
0 <= R("R_GLCRtex_Rev") <= 999999
0 <= R("R_EX_glcur1p_e_") <= 999999
0 <= R("R_EX_glcur1p_e__Rev") <= 999999
0 <= R("R_GLCUR1Ptex") <= 999999
0 <= R("R_GLCUR1Ptex_Rev") <= 999999
0 <= R("R_GUI1") <= 999999
0 <= R("R_GUI1_Rev") <= 999999
0 <= R("R_EX_glcur_e_") <= 999999
0 <= R("R_EX_glcur_e__Rev") <= 999999
0 <= R("R_GLCURtex") <= 999999
0 <= R("R_GLCURtex_Rev") <= 999999
0 <= R("R_EX_gln_L_e_") <= 999999
0 <= R("R_EX_gln_L_e__Rev") <= 999999
0 <= R("R_GLNtex") <= 999999
0 <= R("R_GLNtex_Rev") <= 999999
0 <= R("R_G1SAT") <= 999999
0 <= R("R_G1SAT_Rev") <= 999999
0 <= R("R_GLUR") <= 999999
0 <= R("R_GLUR_Rev") <= 999999
0 <= R("R_EX_glu_L_e_") <= 999999
0 <= R("R_EX_glu_L_e__Rev") <= 999999
0 <= R("R_GLUtex") <= 999999
0 <= R("R_GLUtex_Rev") <= 999999
0 <= R("R_GLUABUTt7pp") <= 999999
0 <= R("R_GLUABUTt7pp_Rev") <= 999999
0 <= R("R_EX_glyald_e_") <= 999999
0 <= R("R_EX_glyald_e__Rev") <= 999999
0 <= R("R_GLYALDtex") <= 999999
0 <= R("R_GLYALDtex_Rev") <= 999999
0 <= R("R_GLYALDtpp") <= 999999
0 <= R("R_GLYALDtpp_Rev") <= 999999
0 <= R("R_EX_glyb_e_") <= 999999
0 <= R("R_EX_glyb_e__Rev") <= 999999
0 <= R("R_GLYBtex") <= 999999
0 <= R("R_GLYBtex_Rev") <= 999999
0 <= R("R_GLYAT") <= 999999
0 <= R("R_GLYAT_Rev") <= 999999
0 <= R("R_EX_glyc2p_e_") <= 999999
0 <= R("R_EX_glyc2p_e__Rev") <= 999999
0 <= R("R_GLYC2Ptex") <= 999999
0 <= R("R_GLYC2Ptex_Rev") <= 999999
0 <= R("R_EX_glyc3p_e_") <= 999999
0 <= R("R_EX_glyc3p_e__Rev") <= 999999
0 <= R("R_GLYC3Ptex") <= 999999
0 <= R("R_GLYC3Ptex_Rev") <= 999999
0 <= R("R_GLYCtpp") <= 999999
0 <= R("R_GLYCtpp_Rev") <= 999999
0 <= R("R_EX_glyc_e_") <= 999999
0 <= R("R_EX_glyc_e__Rev") <= 999999
0 <= R("R_GLYCtex") <= 999999
0 <= R("R_GLYCtex_Rev") <= 999999
0 <= R("R_EX_glyclt_e_") <= 999999
0 <= R("R_EX_glyclt_e__Rev") <= 999999
0 <= R("R_GLYCLTtex") <= 999999
0 <= R("R_GLYCLTtex_Rev") <= 999999
0 <= R("R_EX_glyc_R_e_") <= 999999
0 <= R("R_EX_glyc_R_e__Rev") <= 999999
0 <= R("R_GLYCAtex") <= 999999
0 <= R("R_GLYCAtex_Rev") <= 999999
0 <= R("R_EX_gly_e_") <= 999999
0 <= R("R_EX_gly_e__Rev") <= 999999
0 <= R("R_GLYtex") <= 999999
0 <= R("R_GLYtex_Rev") <= 999999
0 <= R("R_GK1") <= 999999
0 <= R("R_GK1_Rev") <= 999999
0 <= R("R_EX_gmp_e_") <= 999999
0 <= R("R_EX_gmp_e__Rev") <= 999999
0 <= R("R_GMPtex") <= 999999
0 <= R("R_GMPtex_Rev") <= 999999
0 <= R("R_EX_gsn_e_") <= 999999
0 <= R("R_EX_gsn_e__Rev") <= 999999
0 <= R("R_GSNtex") <= 999999
0 <= R("R_GSNtex_Rev") <= 999999
0 <= R("R_EX_gthox_e_") <= 999999
0 <= R("R_EX_gthox_e__Rev") <= 999999
0 <= R("R_GTHOXtex") <= 999999
0 <= R("R_GTHOXtex_Rev") <= 999999
0 <= R("R_FALGTHLs") <= 999999
0 <= R("R_FALGTHLs_Rev") <= 999999
0 <= R("R_EX_gthrd_e_") <= 999999
0 <= R("R_EX_gthrd_e__Rev") <= 999999
0 <= R("R_GTHRDtex") <= 999999
0 <= R("R_GTHRDtex_Rev") <= 999999
0 <= R("R_ADK3") <= 999999
0 <= R("R_ADK3_Rev") <= 999999
0 <= R("R_EX_gtp_e_") <= 999999
0 <= R("R_EX_gtp_e__Rev") <= 999999
0 <= R("R_GTPtex") <= 999999
0 <= R("R_GTPtex_Rev") <= 999999
0 <= R("R_EX_gua_e_") <= 999999
0 <= R("R_EX_gua_e__Rev") <= 999999
0 <= R("R_GUAtex") <= 999999
0 <= R("R_GUAtex_Rev") <= 999999
0 <= R("R_GUAtpp") <= 999999
0 <= R("R_GUAtpp_Rev") <= 999999
0 <= R("R_EX_h2_e_") <= 999999
0 <= R("R_EX_h2_e__Rev") <= 999999
0 <= R("R_H2tex") <= 999999
0 <= R("R_H2tex_Rev") <= 999999
0 <= R("R_EX_h2o2_e_") <= 999999
0 <= R("R_EX_h2o2_e__Rev") <= 999999
0 <= R("R_H2O2tex") <= 999999
0 <= R("R_H2O2tex_Rev") <= 999999
0 <= R("R_IPPMIb") <= 999999
0 <= R("R_IPPMIb_Rev") <= 999999
0 <= R("R_ACONTb") <= 999999
0 <= R("R_ACONTb_Rev") <= 999999
0 <= R("R_HCO3E") <= 999999
0 <= R("R_HCO3E_Rev") <= 999999
0 <= R("R_DHORTS") <= 999999
0 <= R("R_DHORTS_Rev") <= 999999
0 <= R("R_FRULYSDG") <= 999999
0 <= R("R_FRULYSDG_Rev") <= 999999
0 <= R("R_FUM") <= 999999
0 <= R("R_FUM_Rev") <= 999999
0 <= R("R_EX_h2o_e_") <= 999999
0 <= R("R_EX_h2o_e__Rev") <= 999999
0 <= R("R_H2Otex") <= 999999
0 <= R("R_H2Otex_Rev") <= 999999
0 <= R("R_H2Otpp") <= 999999
0 <= R("R_H2Otpp_Rev") <= 999999
0 <= R("R_H2tpp") <= 999999
0 <= R("R_H2tpp_Rev") <= 999999
0 <= R("R_EX_h2s_e_") <= 999999
0 <= R("R_EX_h2s_e__Rev") <= 999999
0 <= R("R_H2Stex") <= 999999
0 <= R("R_H2Stex_Rev") <= 999999
0 <= R("R_EX_hacolipa_e_") <= 999999
0 <= R("R_EX_hacolipa_e__Rev") <= 999999
0 <= R("R_EX_halipa_e_") <= 999999
0 <= R("R_EX_halipa_e__Rev") <= 999999
0 <= R("R_CBIAT") <= 999999
0 <= R("R_CBIAT_Rev") <= 999999
0 <= R("R_CBLAT") <= 999999
0 <= R("R_CBLAT_Rev") <= 999999
0 <= R("R_EX_hdca_e_") <= 999999
0 <= R("R_EX_hdca_e__Rev") <= 999999
0 <= R("R_EX_hdcea_e_") <= 999999
0 <= R("R_EX_hdcea_e__Rev") <= 999999
0 <= R("R_EX_h_e_") <= 999999
0 <= R("R_EX_h_e__Rev") <= 999999
0 <= R("R_Htex") <= 999999
0 <= R("R_Htex_Rev") <= 999999
0 <= R("R_EX_hg2_e_") <= 999999
0 <= R("R_EX_hg2_e__Rev") <= 999999
0 <= R("R_HG2tex") <= 999999
0 <= R("R_HG2tex_Rev") <= 999999
0 <= R("R_EX_his_L_e_") <= 999999
0 <= R("R_EX_his_L_e__Rev") <= 999999
0 <= R("R_HIStex") <= 999999
0 <= R("R_HIStex_Rev") <= 999999
0 <= R("R_EX_hom_L_e_") <= 999999
0 <= R("R_EX_hom_L_e__Rev") <= 999999
0 <= R("R_HOMtex") <= 999999
0 <= R("R_HOMtex_Rev") <= 999999
0 <= R("R_DDGLCNt2rpp") <= 999999
0 <= R("R_DDGLCNt2rpp_Rev") <= 999999
0 <= R("R_HCINNMt2rpp") <= 999999
0 <= R("R_HCINNMt2rpp_Rev") <= 999999
0 <= R("R_HPPPNt2rpp") <= 999999
0 <= R("R_HPPPNt2rpp_Rev") <= 999999
0 <= R("R_5DGLCNt2rpp") <= 999999
0 <= R("R_5DGLCNt2rpp_Rev") <= 999999
0 <= R("R_ACACt2pp") <= 999999
0 <= R("R_ACACt2pp_Rev") <= 999999
0 <= R("R_ACt2rpp") <= 999999
0 <= R("R_ACt2rpp_Rev") <= 999999
0 <= R("R_ADEt2rpp") <= 999999
0 <= R("R_ADEt2rpp_Rev") <= 999999
0 <= R("R_ADNt2rpp") <= 999999
0 <= R("R_ADNt2rpp_Rev") <= 999999
0 <= R("R_AKGt2rpp") <= 999999
0 <= R("R_AKGt2rpp_Rev") <= 999999
0 <= R("R_ALLTNt2rpp") <= 999999
0 <= R("R_ALLTNt2rpp_Rev") <= 999999
0 <= R("R_ARBt2rpp") <= 999999
0 <= R("R_ARBt2rpp_Rev") <= 999999
0 <= R("R_ASNt2rpp") <= 999999
0 <= R("R_ASNt2rpp_Rev") <= 999999
0 <= R("R_BUTt2rpp") <= 999999
0 <= R("R_BUTt2rpp_Rev") <= 999999
0 <= R("R_CRNDt2rpp") <= 999999
0 <= R("R_CRNDt2rpp_Rev") <= 999999
0 <= R("R_CRNt2rpp") <= 999999
0 <= R("R_CRNt2rpp_Rev") <= 999999
0 <= R("R_CTBTt2rpp") <= 999999
0 <= R("R_CTBTt2rpp_Rev") <= 999999
0 <= R("R_CYTDt2rpp") <= 999999
0 <= R("R_CYTDt2rpp_Rev") <= 999999
0 <= R("R_ETOHt2rpp") <= 999999
0 <= R("R_ETOHt2rpp_Rev") <= 999999
0 <= R("R_FRUURt2rpp") <= 999999
0 <= R("R_FRUURt2rpp_Rev") <= 999999
0 <= R("R_FUCtpp") <= 999999
0 <= R("R_FUCtpp_Rev") <= 999999
0 <= R("R_GALCTt2rpp") <= 999999
0 <= R("R_GALCTt2rpp_Rev") <= 999999
0 <= R("R_GALURt2rpp") <= 999999
0 <= R("R_GALURt2rpp_Rev") <= 999999
0 <= R("R_GLCNt2rpp") <= 999999
0 <= R("R_GLCNt2rpp_Rev") <= 999999
0 <= R("R_GLCRt2rpp") <= 999999
0 <= R("R_GLCRt2rpp_Rev") <= 999999
0 <= R("R_GLCURt2rpp") <= 999999
0 <= R("R_GLCURt2rpp_Rev") <= 999999
0 <= R("R_GLUt2rpp") <= 999999
0 <= R("R_GLUt2rpp_Rev") <= 999999
0 <= R("R_GLYCLTt2rpp") <= 999999
0 <= R("R_GLYCLTt2rpp_Rev") <= 999999
0 <= R("R_GLYCAt2rpp") <= 999999
0 <= R("R_GLYCAt2rpp_Rev") <= 999999
0 <= R("R_HISt2rpp") <= 999999
0 <= R("R_HISt2rpp_Rev") <= 999999
0 <= R("R_HPYRI") <= 999999
0 <= R("R_HPYRI_Rev") <= 999999
0 <= R("R_EX_hxa_e_") <= 999999
0 <= R("R_EX_hxa_e__Rev") <= 999999
0 <= R("R_HXAtex") <= 999999
0 <= R("R_HXAtex_Rev") <= 999999
0 <= R("R_EX_hxan_e_") <= 999999
0 <= R("R_EX_hxan_e__Rev") <= 999999
0 <= R("R_HYXNtex") <= 999999
0 <= R("R_HYXNtex_Rev") <= 999999
0 <= R("R_HYXNtpp") <= 999999
0 <= R("R_HYXNtpp_Rev") <= 999999
0 <= R("R_HEXt2rpp") <= 999999
0 <= R("R_HEXt2rpp_Rev") <= 999999
0 <= R("R_ACOAD2f") <= 999999
0 <= R("R_ACOAD2f_Rev") <= 999999
0 <= R("R_EX_idon_L_e_") <= 999999
0 <= R("R_EX_idon_L_e__Rev") <= 999999
0 <= R("R_IDONtex") <= 999999
0 <= R("R_IDONtex_Rev") <= 999999
0 <= R("R_IDONt2rpp") <= 999999
0 <= R("R_IDONt2rpp_Rev") <= 999999
0 <= R("R_ILETA") <= 999999
0 <= R("R_ILETA_Rev") <= 999999
0 <= R("R_EX_ile_L_e_") <= 999999
0 <= R("R_EX_ile_L_e__Rev") <= 999999
0 <= R("R_ILEtex") <= 999999
0 <= R("R_ILEtex_Rev") <= 999999
0 <= R("R_ILEt2rpp") <= 999999
0 <= R("R_ILEt2rpp_Rev") <= 999999
0 <= R("R_IMPC") <= 999999
0 <= R("R_IMPC_Rev") <= 999999
0 <= R("R_EX_imp_e_") <= 999999
0 <= R("R_EX_imp_e__Rev") <= 999999
0 <= R("R_IMPtex") <= 999999
0 <= R("R_IMPtex_Rev") <= 999999
0 <= R("R_EX_indole_e_") <= 999999
0 <= R("R_EX_indole_e__Rev") <= 999999
0 <= R("R_INDOLEtex") <= 999999
0 <= R("R_INDOLEtex_Rev") <= 999999
0 <= R("R_INDOLEt2rpp") <= 999999
0 <= R("R_INDOLEt2rpp_Rev") <= 999999
0 <= R("R_EX_inost_e_") <= 999999
0 <= R("R_EX_inost_e__Rev") <= 999999
0 <= R("R_INSTtex") <= 999999
0 <= R("R_INSTtex_Rev") <= 999999
0 <= R("R_EX_ins_e_") <= 999999
0 <= R("R_EX_ins_e__Rev") <= 999999
0 <= R("R_INStex") <= 999999
0 <= R("R_INStex_Rev") <= 999999
0 <= R("R_INSt2rpp") <= 999999
0 <= R("R_INSt2rpp_Rev") <= 999999
0 <= R("R_IPDDI") <= 999999
0 <= R("R_IPDDI_Rev") <= 999999
0 <= R("R_EX_isetac_e_") <= 999999
0 <= R("R_EX_isetac_e__Rev") <= 999999
0 <= R("R_ISETACtex") <= 999999
0 <= R("R_ISETACtex_Rev") <= 999999
0 <= R("R_ADK4") <= 999999
0 <= R("R_ADK4_Rev") <= 999999
0 <= R("R_EX_kdo2lipid4_e_") <= 999999
0 <= R("R_EX_kdo2lipid4_e__Rev") <= 999999
0 <= R("R_EX_k_e_") <= 999999
0 <= R("R_EX_k_e__Rev") <= 999999
0 <= R("R_Ktex") <= 999999
0 <= R("R_Ktex_Rev") <= 999999
0 <= R("R_EX_lac_D_e_") <= 999999
0 <= R("R_EX_lac_D_e__Rev") <= 999999
0 <= R("R_D_LACtex") <= 999999
0 <= R("R_D_LACtex_Rev") <= 999999
0 <= R("R_D_LACt2pp") <= 999999
0 <= R("R_D_LACt2pp_Rev") <= 999999
0 <= R("R_EX_lac_L_e_") <= 999999
0 <= R("R_EX_lac_L_e__Rev") <= 999999
0 <= R("R_L_LACtex") <= 999999
0 <= R("R_L_LACtex_Rev") <= 999999
0 <= R("R_L_LACt2rpp") <= 999999
0 <= R("R_L_LACt2rpp_Rev") <= 999999
0 <= R("R_ALAGLUE") <= 999999
0 <= R("R_ALAGLUE_Rev") <= 999999
0 <= R("R_EX_LalaDgluMdapDala_e_") <= 999999
0 <= R("R_EX_LalaDgluMdapDala_e__Rev") <= 999999
0 <= R("R_4PEPTtex") <= 999999
0 <= R("R_4PEPTtex_Rev") <= 999999
0 <= R("R_EX_LalaDgluMdap_e_") <= 999999
0 <= R("R_EX_LalaDgluMdap_e__Rev") <= 999999
0 <= R("R_3PEPTtex") <= 999999
0 <= R("R_3PEPTtex_Rev") <= 999999
0 <= R("R_EX_lcts_e_") <= 999999
0 <= R("R_EX_lcts_e__Rev") <= 999999
0 <= R("R_LCTStex") <= 999999
0 <= R("R_LCTStex_Rev") <= 999999
0 <= R("R_LCTStpp") <= 999999
0 <= R("R_LCTStpp_Rev") <= 999999
0 <= R("R_EX_leu_L_e_") <= 999999
0 <= R("R_EX_leu_L_e__Rev") <= 999999
0 <= R("R_LEUtex") <= 999999
0 <= R("R_LEUtex_Rev") <= 999999
0 <= R("R_LEUt2rpp") <= 999999
0 <= R("R_LEUt2rpp_Rev") <= 999999
0 <= R("R_EX_lipa_cold_e_") <= 999999
0 <= R("R_EX_lipa_cold_e__Rev") <= 999999
0 <= R("R_EX_lipa_e_") <= 999999
0 <= R("R_EX_lipa_e__Rev") <= 999999
0 <= R("R_EX_lys_L_e_") <= 999999
0 <= R("R_EX_lys_L_e__Rev") <= 999999
0 <= R("R_LYStex") <= 999999
0 <= R("R_LYStex_Rev") <= 999999
0 <= R("R_EX_lyx_L_e_") <= 999999
0 <= R("R_EX_lyx_L_e__Rev") <= 999999
0 <= R("R_LYXtex") <= 999999
0 <= R("R_LYXtex_Rev") <= 999999
0 <= R("R_MCOATA") <= 999999
0 <= R("R_MCOATA_Rev") <= 999999
0 <= R("R_EX_mal_D_e_") <= 999999
0 <= R("R_EX_mal_D_e__Rev") <= 999999
0 <= R("R_MALDtex") <= 999999
0 <= R("R_MALDtex_Rev") <= 999999
0 <= R("R_EX_mal_L_e_") <= 999999
0 <= R("R_EX_mal_L_e__Rev") <= 999999
0 <= R("R_MALtex") <= 999999
0 <= R("R_MALtex_Rev") <= 999999
0 <= R("R_MALTATr") <= 999999
0 <= R("R_MALTATr_Rev") <= 999999
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
0 <= R("R_MAN6Ptex") <= 999999
0 <= R("R_MAN6Ptex_Rev") <= 999999
0 <= R("R_EX_man_e_") <= 999999
0 <= R("R_EX_man_e__Rev") <= 999999
0 <= R("R_MANtex") <= 999999
0 <= R("R_MANtex_Rev") <= 999999
0 <= R("R_EX_manglyc_e_") <= 999999
0 <= R("R_EX_manglyc_e__Rev") <= 999999
0 <= R("R_MANGLYCtex") <= 999999
0 <= R("R_MANGLYCtex_Rev") <= 999999
0 <= R("R_EX_melib_e_") <= 999999
0 <= R("R_EX_melib_e__Rev") <= 999999
0 <= R("R_MELIBtex") <= 999999
0 <= R("R_MELIBtex_Rev") <= 999999
0 <= R("R_EX_met_D_e_") <= 999999
0 <= R("R_EX_met_D_e__Rev") <= 999999
0 <= R("R_METDtex") <= 999999
0 <= R("R_METDtex_Rev") <= 999999
0 <= R("R_MTHFC") <= 999999
0 <= R("R_MTHFC_Rev") <= 999999
0 <= R("R_EX_met_L_e_") <= 999999
0 <= R("R_EX_met_L_e__Rev") <= 999999
0 <= R("R_METtex") <= 999999
0 <= R("R_METtex_Rev") <= 999999
0 <= R("R_EX_metsox_R_L_e_") <= 999999
0 <= R("R_EX_metsox_R_L_e__Rev") <= 999999
0 <= R("R_METSOX2tex") <= 999999
0 <= R("R_METSOX2tex_Rev") <= 999999
0 <= R("R_EX_metsox_S_L_e_") <= 999999
0 <= R("R_EX_metsox_S_L_e__Rev") <= 999999
0 <= R("R_METSOX1tex") <= 999999
0 <= R("R_METSOX1tex_Rev") <= 999999
0 <= R("R_EX_mg2_e_") <= 999999
0 <= R("R_EX_mg2_e__Rev") <= 999999
0 <= R("R_MG2tex") <= 999999
0 <= R("R_MG2tex_Rev") <= 999999
0 <= R("R_MG2t3_2pp") <= 999999
0 <= R("R_MG2t3_2pp_Rev") <= 999999
0 <= R("R_MCITL2") <= 999999
0 <= R("R_MCITL2_Rev") <= 999999
0 <= R("R_EX_minohp_e_") <= 999999
0 <= R("R_EX_minohp_e__Rev") <= 999999
0 <= R("R_MME") <= 999999
0 <= R("R_MME_Rev") <= 999999
0 <= R("R_EX_mmet_e_") <= 999999
0 <= R("R_EX_mmet_e__Rev") <= 999999
0 <= R("R_MMETtex") <= 999999
0 <= R("R_MMETtex_Rev") <= 999999
0 <= R("R_EX_mn2_e_") <= 999999
0 <= R("R_EX_mn2_e__Rev") <= 999999
0 <= R("R_MNtex") <= 999999
0 <= R("R_MNtex_Rev") <= 999999
0 <= R("R_EX_mnl_e_") <= 999999
0 <= R("R_EX_mnl_e__Rev") <= 999999
0 <= R("R_MNLtex") <= 999999
0 <= R("R_MNLtex_Rev") <= 999999
0 <= R("R_EX_mobd_e_") <= 999999
0 <= R("R_EX_mobd_e__Rev") <= 999999
0 <= R("R_MOBDtex") <= 999999
0 <= R("R_MOBDtex_Rev") <= 999999
0 <= R("R_EX_mso3_e_") <= 999999
0 <= R("R_EX_mso3_e__Rev") <= 999999
0 <= R("R_MSO3tex") <= 999999
0 <= R("R_MSO3tex_Rev") <= 999999
0 <= R("R_EX_n2o_e_") <= 999999
0 <= R("R_EX_n2o_e__Rev") <= 999999
0 <= R("R_N2Otex") <= 999999
0 <= R("R_N2Otex_Rev") <= 999999
0 <= R("R_N2Otpp") <= 999999
0 <= R("R_N2Otpp_Rev") <= 999999
0 <= R("R_EX_na1_e_") <= 999999
0 <= R("R_EX_na1_e__Rev") <= 999999
0 <= R("R_NAtex") <= 999999
0 <= R("R_NAtex_Rev") <= 999999
0 <= R("R_CAt6pp") <= 999999
0 <= R("R_CAt6pp_Rev") <= 999999
0 <= R("R_EX_nac_e_") <= 999999
0 <= R("R_EX_nac_e__Rev") <= 999999
0 <= R("R_NACtex") <= 999999
0 <= R("R_NACtex_Rev") <= 999999
0 <= R("R_DHBD") <= 999999
0 <= R("R_DHBD_Rev") <= 999999
0 <= R("R_PERD") <= 999999
0 <= R("R_PERD_Rev") <= 999999
0 <= R("R_TAGURr") <= 999999
0 <= R("R_TAGURr_Rev") <= 999999
0 <= R("R_ACALD") <= 999999
0 <= R("R_ACALD_Rev") <= 999999
0 <= R("R_ALCD2x") <= 999999
0 <= R("R_ALCD2x_Rev") <= 999999
0 <= R("R_GLTPD") <= 999999
0 <= R("R_GLTPD_Rev") <= 999999
0 <= R("R_E4PD") <= 999999
0 <= R("R_E4PD_Rev") <= 999999
0 <= R("R_FALDH2") <= 999999
0 <= R("R_FALDH2_Rev") <= 999999
0 <= R("R_LDH_D") <= 999999
0 <= R("R_LDH_D_Rev") <= 999999
0 <= R("R_MDH") <= 999999
0 <= R("R_MDH_Rev") <= 999999
0 <= R("R_MANAO") <= 999999
0 <= R("R_MANAO_Rev") <= 999999
0 <= R("R_M1PD") <= 999999
0 <= R("R_M1PD_Rev") <= 999999
0 <= R("R_TRSARr") <= 999999
0 <= R("R_TRSARr_Rev") <= 999999
0 <= R("R_IDOND") <= 999999
0 <= R("R_IDOND_Rev") <= 999999
0 <= R("R_ALCD19") <= 999999
0 <= R("R_ALCD19_Rev") <= 999999
0 <= R("R_LCARR") <= 999999
0 <= R("R_LCARR_Rev") <= 999999
0 <= R("R_LCARS") <= 999999
0 <= R("R_LCARS_Rev") <= 999999
0 <= R("R_KARA1") <= 999999
0 <= R("R_KARA1_Rev") <= 999999
0 <= R("R_ATHRDHr") <= 999999
0 <= R("R_ATHRDHr_Rev") <= 999999
0 <= R("R_G6PDH2r") <= 999999
0 <= R("R_G6PDH2r_Rev") <= 999999
0 <= R("R_G3PD2") <= 999999
0 <= R("R_G3PD2_Rev") <= 999999
0 <= R("R_GGGABADr") <= 999999
0 <= R("R_GGGABADr_Rev") <= 999999
0 <= R("R_GLUDy") <= 999999
0 <= R("R_GLUDy_Rev") <= 999999
0 <= R("R_HSDy") <= 999999
0 <= R("R_HSDy_Rev") <= 999999
0 <= R("R_ICDHyr") <= 999999
0 <= R("R_ICDHyr_Rev") <= 999999
0 <= R("R_MTHFD") <= 999999
0 <= R("R_MTHFD_Rev") <= 999999
0 <= R("R_KARA2") <= 999999
0 <= R("R_KARA2_Rev") <= 999999
0 <= R("R_SHK3Dr") <= 999999
0 <= R("R_SHK3Dr_Rev") <= 999999
0 <= R("R_3OAR120") <= 999999
0 <= R("R_3OAR120_Rev") <= 999999
0 <= R("R_3OAR100") <= 999999
0 <= R("R_3OAR100_Rev") <= 999999
0 <= R("R_3OAR60") <= 999999
0 <= R("R_3OAR60_Rev") <= 999999
0 <= R("R_3OAR140") <= 999999
0 <= R("R_3OAR140_Rev") <= 999999
0 <= R("R_3OAR80") <= 999999
0 <= R("R_3OAR80_Rev") <= 999999
0 <= R("R_3OAR180") <= 999999
0 <= R("R_3OAR180_Rev") <= 999999
0 <= R("R_3OAR160") <= 999999
0 <= R("R_3OAR160_Rev") <= 999999
0 <= R("R_5DGLCNR") <= 999999
0 <= R("R_5DGLCNR_Rev") <= 999999
0 <= R("R_3OAR40") <= 999999
0 <= R("R_3OAR40_Rev") <= 999999
0 <= R("R_DHFR") <= 999999
0 <= R("R_DHFR_Rev") <= 999999
0 <= R("R_GTHOr") <= 999999
0 <= R("R_GTHOr_Rev") <= 999999
0 <= R("R_CBMKr") <= 999999
0 <= R("R_CBMKr_Rev") <= 999999
0 <= R("R_EX_nh4_e_") <= 999999
0 <= R("R_EX_nh4_e__Rev") <= 999999
0 <= R("R_NH4tex") <= 999999
0 <= R("R_NH4tex_Rev") <= 999999
0 <= R("R_NH4tpp") <= 999999
0 <= R("R_NH4tpp_Rev") <= 999999
0 <= R("R_EX_ni2_e_") <= 999999
0 <= R("R_EX_ni2_e__Rev") <= 999999
0 <= R("R_NI2tex") <= 999999
0 <= R("R_NI2tex_Rev") <= 999999
0 <= R("R_NNATr") <= 999999
0 <= R("R_NNATr_Rev") <= 999999
0 <= R("R_EX_nmn_e_") <= 999999
0 <= R("R_EX_nmn_e__Rev") <= 999999
0 <= R("R_NMNtex") <= 999999
0 <= R("R_NMNtex_Rev") <= 999999
0 <= R("R_EX_no2_e_") <= 999999
0 <= R("R_EX_no2_e__Rev") <= 999999
0 <= R("R_NO2tex") <= 999999
0 <= R("R_NO2tex_Rev") <= 999999
0 <= R("R_NO2t2rpp") <= 999999
0 <= R("R_NO2t2rpp_Rev") <= 999999
0 <= R("R_EX_no3_e_") <= 999999
0 <= R("R_EX_no3_e__Rev") <= 999999
0 <= R("R_NO3tex") <= 999999
0 <= R("R_NO3tex_Rev") <= 999999
0 <= R("R_EX_no_e_") <= 999999
0 <= R("R_EX_no_e__Rev") <= 999999
0 <= R("R_NOtex") <= 999999
0 <= R("R_NOtex_Rev") <= 999999
0 <= R("R_NOtpp") <= 999999
0 <= R("R_NOtpp_Rev") <= 999999
0 <= R("R_EX_o16a4colipa_e_") <= 999999
0 <= R("R_EX_o16a4colipa_e__Rev") <= 999999
0 <= R("R_EX_o2_e_") <= 999999
0 <= R("R_EX_o2_e__Rev") <= 999999
0 <= R("R_O2tex") <= 999999
0 <= R("R_O2tex_Rev") <= 999999
0 <= R("R_O2tpp") <= 999999
0 <= R("R_O2tpp_Rev") <= 999999
0 <= R("R_EX_o2s_e_") <= 999999
0 <= R("R_EX_o2s_e__Rev") <= 999999
0 <= R("R_O2Stex") <= 999999
0 <= R("R_O2Stex_Rev") <= 999999
0 <= R("R_ACOAD3f") <= 999999
0 <= R("R_ACOAD3f_Rev") <= 999999
0 <= R("R_EX_ocdca_e_") <= 999999
0 <= R("R_EX_ocdca_e__Rev") <= 999999
0 <= R("R_EX_ocdcea_e_") <= 999999
0 <= R("R_EX_ocdcea_e__Rev") <= 999999
0 <= R("R_EX_octa_e_") <= 999999
0 <= R("R_EX_octa_e__Rev") <= 999999
0 <= R("R_OCTAtex") <= 999999
0 <= R("R_OCTAtex_Rev") <= 999999
0 <= R("R_OHPBAT") <= 999999
0 <= R("R_OHPBAT_Rev") <= 999999
0 <= R("R_ARGORNt7pp") <= 999999
0 <= R("R_ARGORNt7pp_Rev") <= 999999
0 <= R("R_OCBT") <= 999999
0 <= R("R_OCBT_Rev") <= 999999
0 <= R("R_EX_orn_e_") <= 999999
0 <= R("R_EX_orn_e__Rev") <= 999999
0 <= R("R_ORNtex") <= 999999
0 <= R("R_ORNtex_Rev") <= 999999
0 <= R("R_EX_orot_e_") <= 999999
0 <= R("R_EX_orot_e__Rev") <= 999999
0 <= R("R_OROTtex") <= 999999
0 <= R("R_OROTtex_Rev") <= 999999
0 <= R("R_EX_pacald_e_") <= 999999
0 <= R("R_EX_pacald_e__Rev") <= 999999
0 <= R("R_PACALDtex") <= 999999
0 <= R("R_PACALDtex_Rev") <= 999999
0 <= R("R_PACALDt2rpp") <= 999999
0 <= R("R_PACALDt2rpp_Rev") <= 999999
0 <= R("R_EX_peamn_e_") <= 999999
0 <= R("R_EX_peamn_e__Rev") <= 999999
0 <= R("R_PEAMNtex") <= 999999
0 <= R("R_PEAMNtex_Rev") <= 999999
0 <= R("R_CLPNS120pp") <= 999999
0 <= R("R_CLPNS120pp_Rev") <= 999999
0 <= R("R_CLPNS140pp") <= 999999
0 <= R("R_CLPNS140pp_Rev") <= 999999
0 <= R("R_CLPNS141pp") <= 999999
0 <= R("R_CLPNS141pp_Rev") <= 999999
0 <= R("R_CLPNS160pp") <= 999999
0 <= R("R_CLPNS160pp_Rev") <= 999999
0 <= R("R_CLPNS161pp") <= 999999
0 <= R("R_CLPNS161pp_Rev") <= 999999
0 <= R("R_CLPNS180pp") <= 999999
0 <= R("R_CLPNS180pp_Rev") <= 999999
0 <= R("R_CLPNS181pp") <= 999999
0 <= R("R_CLPNS181pp_Rev") <= 999999
0 <= R("R_PHETA1") <= 999999
0 <= R("R_PHETA1_Rev") <= 999999
0 <= R("R_EX_phe_L_e_") <= 999999
0 <= R("R_EX_phe_L_e__Rev") <= 999999
0 <= R("R_PHEtex") <= 999999
0 <= R("R_PHEtex_Rev") <= 999999
0 <= R("R_PHEt2rpp") <= 999999
0 <= R("R_PHEt2rpp_Rev") <= 999999
0 <= R("R_EX_pheme_e_") <= 999999
0 <= R("R_EX_pheme_e__Rev") <= 999999
0 <= R("R_PTAr") <= 999999
0 <= R("R_PTAr_Rev") <= 999999
0 <= R("R_PUNP1") <= 999999
0 <= R("R_PUNP1_Rev") <= 999999
0 <= R("R_PPKr") <= 999999
0 <= R("R_PPKr_Rev") <= 999999
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
0 <= R("R_ATPS4rpp") <= 999999
0 <= R("R_ATPS4rpp_Rev") <= 999999
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
0 <= R("R_PItex") <= 999999
0 <= R("R_PItex_Rev") <= 999999
0 <= R("R_PIt2rpp") <= 999999
0 <= R("R_PIt2rpp_Rev") <= 999999
0 <= R("R_AOXSr") <= 999999
0 <= R("R_AOXSr_Rev") <= 999999
0 <= R("R_ACOAD7f") <= 999999
0 <= R("R_ACOAD7f_Rev") <= 999999
0 <= R("R_EX_pnto_R_e_") <= 999999
0 <= R("R_EX_pnto_R_e__Rev") <= 999999
0 <= R("R_PNTOtex") <= 999999
0 <= R("R_PNTOtex_Rev") <= 999999
0 <= R("R_EX_ppa_e_") <= 999999
0 <= R("R_EX_ppa_e__Rev") <= 999999
0 <= R("R_PPAtex") <= 999999
0 <= R("R_PPAtex_Rev") <= 999999
0 <= R("R_PPALtpp") <= 999999
0 <= R("R_PPALtpp_Rev") <= 999999
0 <= R("R_EX_ppal_e_") <= 999999
0 <= R("R_EX_ppal_e__Rev") <= 999999
0 <= R("R_PPALtex") <= 999999
0 <= R("R_PPALtex_Rev") <= 999999
0 <= R("R_PPAKr") <= 999999
0 <= R("R_PPAKr_Rev") <= 999999
0 <= R("R_PPK2r") <= 999999
0 <= R("R_PPK2r_Rev") <= 999999
0 <= R("R_ORPT") <= 999999
0 <= R("R_ORPT_Rev") <= 999999
0 <= R("R_EX_pppn_e_") <= 999999
0 <= R("R_EX_pppn_e__Rev") <= 999999
0 <= R("R_PPPNtex") <= 999999
0 <= R("R_PPPNtex_Rev") <= 999999
0 <= R("R_PPPNt2rpp") <= 999999
0 <= R("R_PPPNt2rpp_Rev") <= 999999
0 <= R("R_EX_ppt_e_") <= 999999
0 <= R("R_EX_ppt_e__Rev") <= 999999
0 <= R("R_PPTtex") <= 999999
0 <= R("R_PPTtex_Rev") <= 999999
0 <= R("R_PRAGSr") <= 999999
0 <= R("R_PRAGSr_Rev") <= 999999
0 <= R("R_PRMICI") <= 999999
0 <= R("R_PRMICI_Rev") <= 999999
0 <= R("R_EX_progly_e_") <= 999999
0 <= R("R_EX_progly_e__Rev") <= 999999
0 <= R("R_PROGLYtex") <= 999999
0 <= R("R_PROGLYtex_Rev") <= 999999
0 <= R("R_EX_pro_L_e_") <= 999999
0 <= R("R_EX_pro_L_e__Rev") <= 999999
0 <= R("R_PROtex") <= 999999
0 <= R("R_PROtex_Rev") <= 999999
0 <= R("R_PROt2rpp") <= 999999
0 <= R("R_PROt2rpp_Rev") <= 999999
0 <= R("R_FRULYSE") <= 999999
0 <= R("R_FRULYSE_Rev") <= 999999
0 <= R("R_EX_psclys_e_") <= 999999
0 <= R("R_EX_psclys_e__Rev") <= 999999
0 <= R("R_PSCLYStex") <= 999999
0 <= R("R_PSCLYStex_Rev") <= 999999
0 <= R("R_EX_pser_L_e_") <= 999999
0 <= R("R_EX_pser_L_e__Rev") <= 999999
0 <= R("R_PSERtex") <= 999999
0 <= R("R_PSERtex_Rev") <= 999999
0 <= R("R_EX_ptrc_e_") <= 999999
0 <= R("R_EX_ptrc_e__Rev") <= 999999
0 <= R("R_PTRCtex") <= 999999
0 <= R("R_PTRCtex_Rev") <= 999999
0 <= R("R_PTRCORNt7pp") <= 999999
0 <= R("R_PTRCORNt7pp_Rev") <= 999999
0 <= R("R_EX_pyr_e_") <= 999999
0 <= R("R_EX_pyr_e__Rev") <= 999999
0 <= R("R_PYRtex") <= 999999
0 <= R("R_PYRtex_Rev") <= 999999
0 <= R("R_PYRt2rpp") <= 999999
0 <= R("R_PYRt2rpp_Rev") <= 999999
0 <= R("R_PPM") <= 999999
0 <= R("R_PPM_Rev") <= 999999
0 <= R("R_RPI") <= 999999
0 <= R("R_RPI_Rev") <= 999999
0 <= R("R_PRPPS") <= 999999
0 <= R("R_PRPPS_Rev") <= 999999
0 <= R("R_EX_r5p_e_") <= 999999
0 <= R("R_EX_r5p_e__Rev") <= 999999
0 <= R("R_R5Ptex") <= 999999
0 <= R("R_R5Ptex_Rev") <= 999999
0 <= R("R_EX_rib_D_e_") <= 999999
0 <= R("R_EX_rib_D_e__Rev") <= 999999
0 <= R("R_RIBtex") <= 999999
0 <= R("R_RIBtex_Rev") <= 999999
0 <= R("R_RMPA") <= 999999
0 <= R("R_RMPA_Rev") <= 999999
0 <= R("R_RMI") <= 999999
0 <= R("R_RMI_Rev") <= 999999
0 <= R("R_EX_rmn_e_") <= 999999
0 <= R("R_EX_rmn_e__Rev") <= 999999
0 <= R("R_RMNtex") <= 999999
0 <= R("R_RMNtex_Rev") <= 999999
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
0 <= R("R_SBTtex") <= 999999
0 <= R("R_SBTtex_Rev") <= 999999
0 <= R("R_DSERDHr") <= 999999
0 <= R("R_DSERDHr_Rev") <= 999999
0 <= R("R_EX_ser_D_e_") <= 999999
0 <= R("R_EX_ser_D_e__Rev") <= 999999
0 <= R("R_DSERtex") <= 999999
0 <= R("R_DSERtex_Rev") <= 999999
0 <= R("R_SERAT") <= 999999
0 <= R("R_SERAT_Rev") <= 999999
0 <= R("R_SERASr") <= 999999
0 <= R("R_SERASr_Rev") <= 999999
0 <= R("R_LSERDHr") <= 999999
0 <= R("R_LSERDHr_Rev") <= 999999
0 <= R("R_EX_ser_L_e_") <= 999999
0 <= R("R_EX_ser_L_e__Rev") <= 999999
0 <= R("R_SERtex") <= 999999
0 <= R("R_SERtex_Rev") <= 999999
0 <= R("R_SERt2rpp") <= 999999
0 <= R("R_SERt2rpp_Rev") <= 999999
0 <= R("R_PSCVT") <= 999999
0 <= R("R_PSCVT_Rev") <= 999999
0 <= R("R_EX_skm_e_") <= 999999
0 <= R("R_EX_skm_e__Rev") <= 999999
0 <= R("R_SKMtex") <= 999999
0 <= R("R_SKMtex_Rev") <= 999999
0 <= R("R_SDPTA") <= 999999
0 <= R("R_SDPTA_Rev") <= 999999
0 <= R("R_EX_so2_e_") <= 999999
0 <= R("R_EX_so2_e__Rev") <= 999999
0 <= R("R_SO2tex") <= 999999
0 <= R("R_SO2tex_Rev") <= 999999
0 <= R("R_SO2tpp") <= 999999
0 <= R("R_SO2tpp_Rev") <= 999999
0 <= R("R_EX_so3_e_") <= 999999
0 <= R("R_EX_so3_e__Rev") <= 999999
0 <= R("R_SO3tex") <= 999999
0 <= R("R_SO3tex_Rev") <= 999999
0 <= R("R_EX_so4_e_") <= 999999
0 <= R("R_EX_so4_e__Rev") <= 999999
0 <= R("R_SO4tex") <= 999999
0 <= R("R_SO4tex_Rev") <= 999999
0 <= R("R_EX_spmd_e_") <= 999999
0 <= R("R_EX_spmd_e__Rev") <= 999999
0 <= R("R_SPMDtex") <= 999999
0 <= R("R_SPMDtex_Rev") <= 999999
0 <= R("R_ACOAD8f") <= 999999
0 <= R("R_ACOAD8f_Rev") <= 999999
0 <= R("R_SUCOAS") <= 999999
0 <= R("R_SUCOAS_Rev") <= 999999
0 <= R("R_SUCFUMtpp") <= 999999
0 <= R("R_SUCFUMtpp_Rev") <= 999999
0 <= R("R_EX_succ_e_") <= 999999
0 <= R("R_EX_succ_e__Rev") <= 999999
0 <= R("R_SUCCtex") <= 999999
0 <= R("R_SUCCtex_Rev") <= 999999
0 <= R("R_EX_sucr_e_") <= 999999
0 <= R("R_EX_sucr_e__Rev") <= 999999
0 <= R("R_SUCRtex") <= 999999
0 <= R("R_SUCRtex_Rev") <= 999999
0 <= R("R_EX_sulfac_e_") <= 999999
0 <= R("R_EX_sulfac_e__Rev") <= 999999
0 <= R("R_SULFACtex") <= 999999
0 <= R("R_SULFACtex_Rev") <= 999999
0 <= R("R_TGBPA") <= 999999
0 <= R("R_TGBPA_Rev") <= 999999
0 <= R("R_EX_tartr_L_e_") <= 999999
0 <= R("R_EX_tartr_L_e__Rev") <= 999999
0 <= R("R_TARTRtex") <= 999999
0 <= R("R_TARTRtex_Rev") <= 999999
0 <= R("R_TARTRt7pp") <= 999999
0 <= R("R_TARTRt7pp_Rev") <= 999999
0 <= R("R_EX_taur_e_") <= 999999
0 <= R("R_EX_taur_e__Rev") <= 999999
0 <= R("R_TAURtex") <= 999999
0 <= R("R_TAURtex_Rev") <= 999999
0 <= R("R_EX_tcynt_e_") <= 999999
0 <= R("R_EX_tcynt_e__Rev") <= 999999
0 <= R("R_TCYNTtex") <= 999999
0 <= R("R_TCYNTtex_Rev") <= 999999
0 <= R("R_ACOAD6f") <= 999999
0 <= R("R_ACOAD6f_Rev") <= 999999
0 <= R("R_T2DECAI") <= 999999
0 <= R("R_T2DECAI_Rev") <= 999999
0 <= R("R_GHMT2r") <= 999999
0 <= R("R_GHMT2r_Rev") <= 999999
0 <= R("R_EX_thm_e_") <= 999999
0 <= R("R_EX_thm_e__Rev") <= 999999
0 <= R("R_THMtex") <= 999999
0 <= R("R_THMtex_Rev") <= 999999
0 <= R("R_EX_thr_L_e_") <= 999999
0 <= R("R_EX_thr_L_e__Rev") <= 999999
0 <= R("R_THRtex") <= 999999
0 <= R("R_THRtex_Rev") <= 999999
0 <= R("R_THRt2rpp") <= 999999
0 <= R("R_THRt2rpp_Rev") <= 999999
0 <= R("R_EX_thrp_e_") <= 999999
0 <= R("R_EX_thrp_e__Rev") <= 999999
0 <= R("R_THRPtex") <= 999999
0 <= R("R_THRPtex_Rev") <= 999999
0 <= R("R_TMDPP") <= 999999
0 <= R("R_TMDPP_Rev") <= 999999
0 <= R("R_EX_thymd_e_") <= 999999
0 <= R("R_EX_thymd_e__Rev") <= 999999
0 <= R("R_THMDtex") <= 999999
0 <= R("R_THMDtex_Rev") <= 999999
0 <= R("R_THMDt2rpp") <= 999999
0 <= R("R_THMDt2rpp_Rev") <= 999999
0 <= R("R_EX_thym_e_") <= 999999
0 <= R("R_EX_thym_e__Rev") <= 999999
0 <= R("R_THYMtex") <= 999999
0 <= R("R_THYMtex_Rev") <= 999999
0 <= R("R_EX_tma_e_") <= 999999
0 <= R("R_EX_tma_e__Rev") <= 999999
0 <= R("R_TMAtex") <= 999999
0 <= R("R_TMAtex_Rev") <= 999999
0 <= R("R_EX_tmao_e_") <= 999999
0 <= R("R_EX_tmao_e__Rev") <= 999999
0 <= R("R_TMAOtex") <= 999999
0 <= R("R_TMAOtex_Rev") <= 999999
0 <= R("R_EX_tre_e_") <= 999999
0 <= R("R_EX_tre_e__Rev") <= 999999
0 <= R("R_TREtex") <= 999999
0 <= R("R_TREtex_Rev") <= 999999
0 <= R("R_TRPAS2") <= 999999
0 <= R("R_TRPAS2_Rev") <= 999999
0 <= R("R_EX_trp_L_e_") <= 999999
0 <= R("R_EX_trp_L_e__Rev") <= 999999
0 <= R("R_TRPtex") <= 999999
0 <= R("R_TRPtex_Rev") <= 999999
0 <= R("R_TRPt2rpp") <= 999999
0 <= R("R_TRPt2rpp_Rev") <= 999999
0 <= R("R_EX_tsul_e_") <= 999999
0 <= R("R_EX_tsul_e__Rev") <= 999999
0 <= R("R_TSULtex") <= 999999
0 <= R("R_TSULtex_Rev") <= 999999
0 <= R("R_EX_ttdca_e_") <= 999999
0 <= R("R_EX_ttdca_e__Rev") <= 999999
0 <= R("R_EX_ttdcea_e_") <= 999999
0 <= R("R_EX_ttdcea_e__Rev") <= 999999
0 <= R("R_EX_tungs_e_") <= 999999
0 <= R("R_EX_tungs_e__Rev") <= 999999
0 <= R("R_TUNGStex") <= 999999
0 <= R("R_TUNGStex_Rev") <= 999999
0 <= R("R_EX_tym_e_") <= 999999
0 <= R("R_EX_tym_e__Rev") <= 999999
0 <= R("R_TYMtex") <= 999999
0 <= R("R_TYMtex_Rev") <= 999999
0 <= R("R_TYRTA") <= 999999
0 <= R("R_TYRTA_Rev") <= 999999
0 <= R("R_EX_tyr_L_e_") <= 999999
0 <= R("R_EX_tyr_L_e__Rev") <= 999999
0 <= R("R_TYRtex") <= 999999
0 <= R("R_TYRtex_Rev") <= 999999
0 <= R("R_TYRt2rpp") <= 999999
0 <= R("R_TYRt2rpp_Rev") <= 999999
0 <= R("R_EX_tyrp_e_") <= 999999
0 <= R("R_EX_tyrp_e__Rev") <= 999999
0 <= R("R_TYRPtex") <= 999999
0 <= R("R_TYRPtex_Rev") <= 999999
0 <= R("R_UAG2E") <= 999999
0 <= R("R_UAG2E_Rev") <= 999999
0 <= R("R_UAGAAT") <= 999999
0 <= R("R_UAGAAT_Rev") <= 999999
0 <= R("R_EX_uacgam_e_") <= 999999
0 <= R("R_EX_uacgam_e__Rev") <= 999999
0 <= R("R_UACGAMtex") <= 999999
0 <= R("R_UACGAMtex_Rev") <= 999999
0 <= R("R_EX_udpacgal_e_") <= 999999
0 <= R("R_EX_udpacgal_e__Rev") <= 999999
0 <= R("R_UDPACGALtex") <= 999999
0 <= R("R_UDPACGALtex_Rev") <= 999999
0 <= R("R_NDPK2") <= 999999
0 <= R("R_NDPK2_Rev") <= 999999
0 <= R("R_EX_udpgal_e_") <= 999999
0 <= R("R_EX_udpgal_e__Rev") <= 999999
0 <= R("R_UDPGALtex") <= 999999
0 <= R("R_UDPGALtex_Rev") <= 999999
0 <= R("R_UDPG4E") <= 999999
0 <= R("R_UDPG4E_Rev") <= 999999
0 <= R("R_UGLT") <= 999999
0 <= R("R_UGLT_Rev") <= 999999
0 <= R("R_EX_udpg_e_") <= 999999
0 <= R("R_EX_udpg_e__Rev") <= 999999
0 <= R("R_UDPGtex") <= 999999
0 <= R("R_UDPGtex_Rev") <= 999999
0 <= R("R_EX_udpglcur_e_") <= 999999
0 <= R("R_EX_udpglcur_e__Rev") <= 999999
0 <= R("R_UDPGLCURtex") <= 999999
0 <= R("R_UDPGLCURtex_Rev") <= 999999
0 <= R("R_UDPKAAT") <= 999999
0 <= R("R_UDPKAAT_Rev") <= 999999
0 <= R("R_UMPK") <= 999999
0 <= R("R_UMPK_Rev") <= 999999
0 <= R("R_EX_ump_e_") <= 999999
0 <= R("R_EX_ump_e__Rev") <= 999999
0 <= R("R_UMPtex") <= 999999
0 <= R("R_UMPtex_Rev") <= 999999
0 <= R("R_EX_ura_e_") <= 999999
0 <= R("R_EX_ura_e__Rev") <= 999999
0 <= R("R_URAtex") <= 999999
0 <= R("R_URAtex_Rev") <= 999999
0 <= R("R_URAt2rpp") <= 999999
0 <= R("R_URAt2rpp_Rev") <= 999999
0 <= R("R_EX_urea_e_") <= 999999
0 <= R("R_EX_urea_e__Rev") <= 999999
0 <= R("R_UREAtex") <= 999999
0 <= R("R_UREAtex_Rev") <= 999999
0 <= R("R_UREAtpp") <= 999999
0 <= R("R_UREAtpp_Rev") <= 999999
0 <= R("R_PYNP2r") <= 999999
0 <= R("R_PYNP2r_Rev") <= 999999
0 <= R("R_EX_uri_e_") <= 999999
0 <= R("R_EX_uri_e__Rev") <= 999999
0 <= R("R_URItex") <= 999999
0 <= R("R_URItex_Rev") <= 999999
0 <= R("R_URIt2rpp") <= 999999
0 <= R("R_URIt2rpp_Rev") <= 999999
0 <= R("R_VALTA") <= 999999
0 <= R("R_VALTA_Rev") <= 999999
0 <= R("R_EX_val_L_e_") <= 999999
0 <= R("R_EX_val_L_e__Rev") <= 999999
0 <= R("R_VALtex") <= 999999
0 <= R("R_VALtex_Rev") <= 999999
0 <= R("R_VALt2rpp") <= 999999
0 <= R("R_VALt2rpp_Rev") <= 999999
0 <= R("R_EX_xan_e_") <= 999999
0 <= R("R_EX_xan_e__Rev") <= 999999
0 <= R("R_XANtex") <= 999999
0 <= R("R_XANtex_Rev") <= 999999
0 <= R("R_XANtpp") <= 999999
0 <= R("R_XANtpp_Rev") <= 999999
0 <= R("R_EX_xmp_e_") <= 999999
0 <= R("R_EX_xmp_e__Rev") <= 999999
0 <= R("R_XMPtex") <= 999999
0 <= R("R_XMPtex_Rev") <= 999999
0 <= R("R_PUNP7") <= 999999
0 <= R("R_PUNP7_Rev") <= 999999
0 <= R("R_EX_xtsn_e_") <= 999999
0 <= R("R_EX_xtsn_e__Rev") <= 999999
0 <= R("R_XTSNtex") <= 999999
0 <= R("R_XTSNtex_Rev") <= 999999
0 <= R("R_XTSNt2rpp") <= 999999
0 <= R("R_XTSNt2rpp_Rev") <= 999999
0 <= R("R_TKT2") <= 999999
0 <= R("R_TKT2_Rev") <= 999999
0 <= R("R_TKT1") <= 999999
0 <= R("R_TKT1_Rev") <= 999999
0 <= R("R_XYLI1") <= 999999
0 <= R("R_XYLI1_Rev") <= 999999
0 <= R("R_EX_xyl_D_e_") <= 999999
0 <= R("R_EX_xyl_D_e__Rev") <= 999999
0 <= R("R_XYLtex") <= 999999
0 <= R("R_XYLtex_Rev") <= 999999
0 <= R("R_EX_xylu_L_e_") <= 999999
0 <= R("R_EX_xylu_L_e__Rev") <= 999999
0 <= R("R_XYLUtex") <= 999999
0 <= R("R_XYLUtex_Rev") <= 999999
0 <= R("R_EX_zn2_e_") <= 999999
0 <= R("R_EX_zn2_e__Rev") <= 999999
0 <= R("R_Zn2tex") <= 999999
0 <= R("R_Zn2tex_Rev") <= 999999
0 <= R("R_12DGR120tipp") <= 999999
0 <= R("R_12DGR140tipp") <= 999999
0 <= R("R_12DGR141tipp") <= 999999
0 <= R("R_12DGR160tipp") <= 999999
0 <= R("R_12DGR161tipp") <= 999999
0 <= R("R_12DGR180tipp") <= 999999
0 <= R("R_12DGR181tipp") <= 999999
0 <= R("R_AAMYL") <= 999999
0 <= R("R_14GLUCANtexi") <= 999999
0 <= R("R_AAMYLpp") <= 999999
0 <= R("R_DHAD1") <= 999999
0 <= R("R_DHAD2") <= 999999
0 <= R("R_2AGPE120tipp") <= 999999
0 <= R("R_2AGPE140tipp") <= 999999
0 <= R("R_2AGPE141tipp") <= 999999
0 <= R("R_2AGPE160tipp") <= 999999
0 <= R("R_2AGPE161tipp") <= 999999
0 <= R("R_2AGPE180tipp") <= 999999
0 <= R("R_2AGPE181tipp") <= 999999
0 <= R("R_2AGPG120tipp") <= 999999
0 <= R("R_2AGPG140tipp") <= 999999
0 <= R("R_2AGPG141tipp") <= 999999
0 <= R("R_2AGPG160tipp") <= 999999
0 <= R("R_2AGPG161tipp") <= 999999
0 <= R("R_2AGPG180tipp") <= 999999
0 <= R("R_2AGPG181tipp") <= 999999
0 <= R("R_DHQS") <= 999999
0 <= R("R_2AGPA120tipp") <= 999999
0 <= R("R_EDA") <= 999999
0 <= R("R_DRPA") <= 999999
0 <= R("R_2AGPA161tipp") <= 999999
0 <= R("R_2AGPA160tipp") <= 999999
0 <= R("R_MCITD") <= 999999
0 <= R("R_2AGPA181tipp") <= 999999
0 <= R("R_2AGPA180tipp") <= 999999
0 <= R("R_MECDPS") <= 999999
0 <= R("R_SUCBZS") <= 999999
0 <= R("R_2AGPA141tipp") <= 999999
0 <= R("R_2AGPA140tipp") <= 999999
0 <= R("R_DHQTi") <= 999999
0 <= R("R_3HAD40") <= 999999
0 <= R("R_3HAD121") <= 999999
0 <= R("R_3HAD141") <= 999999
0 <= R("R_3HAD161") <= 999999
0 <= R("R_3HAD181") <= 999999
0 <= R("R_3HAD120") <= 999999
0 <= R("R_3HAD100") <= 999999
0 <= R("R_3HAD60") <= 999999
0 <= R("R_3HAD140") <= 999999
0 <= R("R_3HAD180") <= 999999
0 <= R("R_3HAD80") <= 999999
0 <= R("R_3HAD160") <= 999999
0 <= R("R_TRPS3") <= 999999
0 <= R("R_CHORS") <= 999999
0 <= R("R_ADCL") <= 999999
0 <= R("R_HOPNTAL") <= 999999
0 <= R("R_DM_4HBA") <= 999999
0 <= R("R_PPBNGS") <= 999999
0 <= R("R_GLCRAL") <= 999999
0 <= R("R_DM_5DRIB") <= 999999
0 <= R("R_MDRPD") <= 999999
0 <= R("R_DHPS2") <= 999999
0 <= R("R_EDD") <= 999999
0 <= R("R_DM_AACALD") <= 999999
0 <= R("R_ACACCT") <= 999999
0 <= R("R_ACMUMtex") <= 999999
0 <= R("R_ACNML") <= 999999
0 <= R("R_ACOLIPAtex") <= 999999
0 <= R("R_ACSERtpp") <= 999999
0 <= R("R_GLCS1") <= 999999
0 <= R("R_AGMHE") <= 999999
0 <= R("R_ABTA") <= 999999
0 <= R("R_VPAMT") <= 999999
0 <= R("R_ALTRH") <= 999999
0 <= R("R_AMMQLT8") <= 999999
0 <= R("R_OHPHM") <= 999999
0 <= R("R_OMBZLM") <= 999999
0 <= R("R_DMQMT") <= 999999
0 <= R("R_ACONMT") <= 999999
0 <= R("R_ACANTHAT") <= 999999
0 <= R("R_ASPT") <= 999999
0 <= R("R_THRA2i") <= 999999
0 <= R("R_ADNCYC") <= 999999
0 <= R("R_DAGK120") <= 999999
0 <= R("R_DAGK140") <= 999999
0 <= R("R_DAGK141") <= 999999
0 <= R("R_DAGK160") <= 999999
0 <= R("R_DAGK161") <= 999999
0 <= R("R_DAGK180") <= 999999
0 <= R("R_DAGK181") <= 999999
0 <= R("R_DDGLK") <= 999999
0 <= R("R_DDGALK") <= 999999
0 <= R("R_3KGK") <= 999999
0 <= R("R_HMPK1") <= 999999
0 <= R("R_PMPK") <= 999999
0 <= R("R_CDPMEK") <= 999999
0 <= R("R_HETZK") <= 999999
0 <= R("R_MTRK") <= 999999
0 <= R("R_HPPK2") <= 999999
0 <= R("R_ACGAMK") <= 999999
0 <= R("R_ACGK") <= 999999
0 <= R("R_AMANK") <= 999999
0 <= R("R_ADNK1") <= 999999
0 <= R("R_ADOCBIK") <= 999999
0 <= R("R_ALLK") <= 999999
0 <= R("R_ADSK") <= 999999
0 <= R("R_PRASCSi") <= 999999
0 <= R("R_GLDBRAN2") <= 999999
0 <= R("R_BUTCT") <= 999999
0 <= R("R_ASPCT") <= 999999
0 <= R("R_CD2tpp") <= 999999
0 <= R("R_ICHORSi") <= 999999
0 <= R("R_CHORM") <= 999999
0 <= R("R_CHRPL") <= 999999
0 <= R("R_CITL") <= 999999
0 <= R("R_ARGSS") <= 999999
0 <= R("R_OBTFL") <= 999999
0 <= R("R_KAT4") <= 999999
0 <= R("R_KAT5") <= 999999
0 <= R("R_KAT2") <= 999999
0 <= R("R_KAT7") <= 999999
0 <= R("R_KAT3") <= 999999
0 <= R("R_KAT8") <= 999999
0 <= R("R_KAT6") <= 999999
0 <= R("R_KAT1") <= 999999
0 <= R("R_ACPS1") <= 999999
0 <= R("R_ACS") <= 999999
0 <= R("R_COBALT2tpp") <= 999999
0 <= R("R_COLIPAtex") <= 999999
0 <= R("R_CPPPGO2") <= 999999
0 <= R("R_CRNCAL2") <= 999999
0 <= R("R_CRNDCAL2") <= 999999
0 <= R("R_CRNt8pp") <= 999999
0 <= R("R_CTBTCAL2") <= 999999
0 <= R("R_CU2tpp") <= 999999
0 <= R("R_CYStpp") <= 999999
0 <= R("R_PPNCL2") <= 999999
0 <= R("R_DBTS") <= 999999
0 <= R("R_RBFSa") <= 999999
0 <= R("R_AACPS8") <= 999999
0 <= R("R_AGPAT120") <= 999999
0 <= R("R_2AGPEAT120") <= 999999
0 <= R("R_2AGPGAT120") <= 999999
0 <= R("R_AACPS7") <= 999999
0 <= R("R_DDCAtexi") <= 999999
0 <= R("R_MGSA") <= 999999
0 <= R("R_DHNPA2") <= 999999
0 <= R("R_DHPTDCs") <= 999999
0 <= R("R_RBFSb") <= 999999
0 <= R("R_DMSOR2") <= 999999
0 <= R("R_DMSOR2pp") <= 999999
0 <= R("R_TPRDCOAS") <= 999999
0 <= R("R_DPCOAK") <= 999999
0 <= R("R_TDSR1") <= 999999
0 <= R("R_TDSR2") <= 999999
0 <= R("R_TDPADGAT") <= 999999
0 <= R("R_TDPDRE") <= 999999
0 <= R("R_TDPGDH") <= 999999
0 <= R("R_DURIK1") <= 999999
0 <= R("R_DXYLK") <= 999999
0 <= R("R_ECA4COLIPAtex") <= 999999
0 <= R("R_ECA4OALpp") <= 999999
0 <= R("R_IGPDH") <= 999999
0 <= R("R_ENLIPAtex") <= 999999
0 <= R("R_ETHAAL") <= 999999
0 <= R("R_FRUK") <= 999999
0 <= R("R_PFK") <= 999999
0 <= R("R_ARBTNR1") <= 999999
0 <= R("R_CPGNR1") <= 999999
0 <= R("R_FCLK") <= 999999
0 <= R("R_FE2tpp") <= 999999
0 <= R("R_FE3Ri") <= 999999
0 <= R("R_FE3DHBZR") <= 999999
0 <= R("R_ARBTNexs") <= 999999
0 <= R("R_CPGNexs") <= 999999
0 <= R("R_FEENTERexs") <= 999999
0 <= R("R_FE3HOXR1") <= 999999
0 <= R("R_FE3HOXexs") <= 999999
0 <= R("R_FECRMR1") <= 999999
0 <= R("R_FECRMexs") <= 999999
0 <= R("R_FEENTERR1") <= 999999
0 <= R("R_FEOXAMR1") <= 999999
0 <= R("R_FEOXAMexs") <= 999999
0 <= R("R_RNTR1c") <= 999999
0 <= R("R_RNTR3c") <= 999999
0 <= R("R_ARBTNR2") <= 999999
0 <= R("R_CPGNR2") <= 999999
0 <= R("R_FE3HOXR2") <= 999999
0 <= R("R_FECRMR2") <= 999999
0 <= R("R_FEENTERR2") <= 999999
0 <= R("R_FEOXAMR2") <= 999999
0 <= R("R_FORtppi") <= 999999
0 <= R("R_PRAIS") <= 999999
0 <= R("R_HEX7") <= 999999
0 <= R("R_FRULYSK") <= 999999
0 <= R("R_FRD3") <= 999999
0 <= R("R_ASPO5") <= 999999
0 <= R("R_GALM2pp") <= 999999
0 <= R("R_GALCTD") <= 999999
0 <= R("R_GALCTND") <= 999999
0 <= R("R_G1PACT") <= 999999
0 <= R("R_GART") <= 999999
0 <= R("R_CRNt7pp") <= 999999
0 <= R("R_GDPDPK") <= 999999
0 <= R("R_GDMANE") <= 999999
0 <= R("R_GMAND") <= 999999
0 <= R("R_HEPT4") <= 999999
0 <= R("R_HEX1") <= 999999
0 <= R("R_GLCtexi") <= 999999
0 <= R("R_GNK") <= 999999
0 <= R("R_GLCRD") <= 999999
0 <= R("R_ADCS") <= 999999
0 <= R("R_ANS") <= 999999
0 <= R("R_GF6PTA") <= 999999
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
0 <= R("R_PGSA120") <= 999999
0 <= R("R_PGSA161") <= 999999
0 <= R("R_PGSA160") <= 999999
0 <= R("R_PGSA181") <= 999999
0 <= R("R_PGSA180") <= 999999
0 <= R("R_PGSA141") <= 999999
0 <= R("R_PGSA140") <= 999999
0 <= R("R_G3PAT120") <= 999999
0 <= R("R_GLYK") <= 999999
0 <= R("R_GLYCTO4") <= 999999
0 <= R("R_GLBRAN2") <= 999999
0 <= R("R_GLYCK2") <= 999999
0 <= R("R_GLYCK") <= 999999
0 <= R("R_GMHEPK") <= 999999
0 <= R("R_RNDR1b") <= 999999
0 <= R("R_RNDR3b") <= 999999
0 <= R("R_RNDR2b") <= 999999
0 <= R("R_GSNK") <= 999999
0 <= R("R_ASR") <= 999999
0 <= R("R_GRXR") <= 999999
0 <= R("R_DSBCGT") <= 999999
0 <= R("R_DSBGGT") <= 999999
0 <= R("R_GTPDPK") <= 999999
0 <= R("R_CYTDK2") <= 999999
0 <= R("R_RNTR2c") <= 999999
0 <= R("R_CAT") <= 999999
0 <= R("R_FTHFD") <= 999999
0 <= R("R_DAPAL") <= 999999
0 <= R("R_DHBSH") <= 999999
0 <= R("R_LPLIPAL2E120") <= 999999
0 <= R("R_LPLIPAL2E140") <= 999999
0 <= R("R_LPLIPAL2E141") <= 999999
0 <= R("R_LPLIPAL2E160") <= 999999
0 <= R("R_LPLIPAL2E161") <= 999999
0 <= R("R_LPLIPAL2E180") <= 999999
0 <= R("R_LPLIPAL2E181") <= 999999
0 <= R("R_LPLIPAL2G120") <= 999999
0 <= R("R_LPLIPAL2G140") <= 999999
0 <= R("R_LPLIPAL2G141") <= 999999
0 <= R("R_LPLIPAL2G160") <= 999999
0 <= R("R_LPLIPAL2G161") <= 999999
0 <= R("R_LPLIPAL2G180") <= 999999
0 <= R("R_LPLIPAL2G181") <= 999999
0 <= R("R_LPLIPAL2A120") <= 999999
0 <= R("R_LPLIPAL2A161") <= 999999
0 <= R("R_LPLIPAL2A160") <= 999999
0 <= R("R_2MAHMP") <= 999999
0 <= R("R_MICITD") <= 999999
0 <= R("R_LPLIPAL2A181") <= 999999
0 <= R("R_LPLIPAL2A180") <= 999999
0 <= R("R_PGLYCP") <= 999999
0 <= R("R_LPLIPAL2A141") <= 999999
0 <= R("R_LPLIPAL2A140") <= 999999
0 <= R("R_PMDPHT") <= 999999
0 <= R("R_MTAN") <= 999999
0 <= R("R_RZ5PP") <= 999999
0 <= R("R_PGL") <= 999999
0 <= R("R_IPPS") <= 999999
0 <= R("R_NACODA") <= 999999
0 <= R("R_AGDC") <= 999999
0 <= R("R_ACM6PH") <= 999999
0 <= R("R_ACODA") <= 999999
0 <= R("R_ADNUC") <= 999999
0 <= R("R_ADPRDP") <= 999999
0 <= R("R_AGMT") <= 999999
0 <= R("R_AHCYSNS") <= 999999
0 <= R("R_DNTPPA") <= 999999
0 <= R("R_ALAALAD") <= 999999
0 <= R("R_ALLTN") <= 999999
0 <= R("R_NTD7") <= 999999
0 <= R("R_AMPN") <= 999999
0 <= R("R_AGM3PH") <= 999999
0 <= R("R_AGM3PA") <= 999999
0 <= R("R_AGM4PCP") <= 999999
0 <= R("R_AGM4PH") <= 999999
0 <= R("R_AGM4PA") <= 999999
0 <= R("R_AGMH") <= 999999
0 <= R("R_AM3PA") <= 999999
0 <= R("R_AM4PCP") <= 999999
0 <= R("R_AM4PA") <= 999999
0 <= R("R_AP4AH") <= 999999
0 <= R("R_AP5AH") <= 999999
0 <= R("R_AB6PGH") <= 999999
0 <= R("R_ASCBPL") <= 999999
0 <= R("R_ASNN") <= 999999
0 <= R("R_ATPM") <= 999999
0 <= R("R_NTP1") <= 999999
0 <= R("R_NTPP6") <= 999999
0 <= R("R_14GLUCANabcpp") <= 999999
0 <= R("R_DAPabcpp") <= 999999
0 <= R("R_ADOCBLabcpp") <= 999999
0 <= R("R_ALAALAabcpp") <= 999999
0 <= R("R_ALAabcpp") <= 999999
0 <= R("R_ALLabcpp") <= 999999
0 <= R("R_ANHMK") <= 999999
0 <= R("R_ARBabcpp") <= 999999
0 <= R("R_ARBTNabcpp") <= 999999
0 <= R("R_ARGabcpp") <= 999999
0 <= R("R_ASNabcpp") <= 999999
0 <= R("R_ASO3t8pp") <= 999999
0 <= R("R_ASPabcpp") <= 999999
0 <= R("R_BUTSO3abcpp") <= 999999
0 <= R("R_PDE1") <= 999999
0 <= R("R_CBIuabcpp") <= 999999
0 <= R("R_CBL1abcpp") <= 999999
0 <= R("R_CD2abcpp") <= 999999
0 <= R("R_CDAPPA120") <= 999999
0 <= R("R_CDAPPA161") <= 999999
0 <= R("R_CDAPPA160") <= 999999
0 <= R("R_CDAPPA181") <= 999999
0 <= R("R_CDAPPA180") <= 999999
0 <= R("R_CDAPPA141") <= 999999
0 <= R("R_CDAPPA140") <= 999999
0 <= R("R_AMPTASECG") <= 999999
0 <= R("R_CGLYabcpp") <= 999999
0 <= R("R_CHLabcpp") <= 999999
0 <= R("R_NTD4") <= 999999
0 <= R("R_CMPN") <= 999999
0 <= R("R_COBALT2abcpp") <= 999999
0 <= R("R_COLIPAabcpp") <= 999999
0 <= R("R_CPGNabcpp") <= 999999
0 <= R("R_CRNDabcpp") <= 999999
0 <= R("R_CRNabcpp") <= 999999
0 <= R("R_CTBTabcpp") <= 999999
0 <= R("R_NTP5") <= 999999
0 <= R("R_NTPP4") <= 999999
0 <= R("R_CU2abcpp") <= 999999
0 <= R("R_CU1abcpp") <= 999999
0 <= R("R_CYSDDS") <= 999999
0 <= R("R_CYSDabcpp") <= 999999
0 <= R("R_CYSDS") <= 999999
0 <= R("R_CYSabc2pp") <= 999999
0 <= R("R_CYSabcpp") <= 999999
0 <= R("R_CYSTL") <= 999999
0 <= R("R_CYTDH") <= 999999
0 <= R("R_5DOAN") <= 999999
0 <= R("R_NTD6") <= 999999
0 <= R("R_NTPP5") <= 999999
0 <= R("R_FA100ACPHi") <= 999999
0 <= R("R_FACOAE100") <= 999999
0 <= R("R_NTD3") <= 999999
0 <= R("R_NTPP3") <= 999999
0 <= R("R_FA120ACPHi") <= 999999
0 <= R("R_FACOAE120") <= 999999
0 <= R("R_NTD8") <= 999999
0 <= R("R_NTPP1") <= 999999
0 <= R("R_NTPTP1") <= 999999
0 <= R("R_DNMPPA") <= 999999
0 <= R("R_NTD12") <= 999999
0 <= R("R_NTP11") <= 999999
0 <= R("R_NTPP10") <= 999999
0 <= R("R_NTD5") <= 999999
0 <= R("R_NTPP7") <= 999999
0 <= R("R_NTD1") <= 999999
0 <= R("R_DUTPDP") <= 999999
0 <= R("R_ETHSO3abcpp") <= 999999
0 <= R("R_F6PP") <= 999999
0 <= R("R_DAAD") <= 999999
0 <= R("R_FBP") <= 999999
0 <= R("R_FE2abcpp") <= 999999
0 <= R("R_FE3DCITabcpp") <= 999999
0 <= R("R_FE3DHBZSabcpp") <= 999999
0 <= R("R_FE3HOXabcpp") <= 999999
0 <= R("R_FE3abcpp") <= 999999
0 <= R("R_FECRMabcpp") <= 999999
0 <= R("R_FEENTERabcpp") <= 999999
0 <= R("R_FEOXAMabcpp") <= 999999
0 <= R("R_GPDDA1") <= 999999
0 <= R("R_G3PCabcpp") <= 999999
0 <= R("R_GPDDA2") <= 999999
0 <= R("R_G3PEabcpp") <= 999999
0 <= R("R_GPDDA4") <= 999999
0 <= R("R_G3PGabcpp") <= 999999
0 <= R("R_GPDDA5") <= 999999
0 <= R("R_G3PIabcpp") <= 999999
0 <= R("R_GPDDA3") <= 999999
0 <= R("R_G3PSabcpp") <= 999999
0 <= R("R_G6PP") <= 999999
0 <= R("R_GALabcpp") <= 999999
0 <= R("R_G6PDA") <= 999999
0 <= R("R_GDPMNP") <= 999999
0 <= R("R_GDPMNH") <= 999999
0 <= R("R_GTPDPDP") <= 999999
0 <= R("R_GGGABAH") <= 999999
0 <= R("R_GLCabcpp") <= 999999
0 <= R("R_GLUN") <= 999999
0 <= R("R_ASNS1") <= 999999
0 <= R("R_PRFGS") <= 999999
0 <= R("R_GLNabcpp") <= 999999
0 <= R("R_GLUabcpp") <= 999999
0 <= R("R_MALS") <= 999999
0 <= R("R_GLYBabcpp") <= 999999
0 <= R("R_G3PT") <= 999999
0 <= R("R_GLYC3Pabcpp") <= 999999
0 <= R("R_GMHEPPA") <= 999999
0 <= R("R_NTD9") <= 999999
0 <= R("R_GP4GH") <= 999999
0 <= R("R_GTHRDabc2pp") <= 999999
0 <= R("R_GTHRDabcpp") <= 999999
0 <= R("R_GTPCI") <= 999999
0 <= R("R_NTP3") <= 999999
0 <= R("R_NTPP2") <= 999999
0 <= R("R_NTPTP2") <= 999999
0 <= R("R_GSPMDA") <= 999999
0 <= R("R_ENTERES") <= 999999
0 <= R("R_ENTERES2") <= 999999
0 <= R("R_GTPCII2") <= 999999
0 <= R("R_LPLIPAL1E120pp") <= 999999
0 <= R("R_LPLIPAL1E140pp") <= 999999
0 <= R("R_LPLIPAL1E141pp") <= 999999
0 <= R("R_LPLIPAL1E160pp") <= 999999
0 <= R("R_LPLIPAL1E161pp") <= 999999
0 <= R("R_LPLIPAL1E180pp") <= 999999
0 <= R("R_LPLIPAL1E181pp") <= 999999
0 <= R("R_LPLIPAL1G120pp") <= 999999
0 <= R("R_LPLIPAL1G140pp") <= 999999
0 <= R("R_LPLIPAL1G141pp") <= 999999
0 <= R("R_LPLIPAL1G160pp") <= 999999
0 <= R("R_LPLIPAL1G161pp") <= 999999
0 <= R("R_LPLIPAL1G180pp") <= 999999
0 <= R("R_LPLIPAL1G181pp") <= 999999
0 <= R("R_LPLIPAL1A120pp") <= 999999
0 <= R("R_LPLIPAL1A161pp") <= 999999
0 <= R("R_LPLIPAL1A160pp") <= 999999
0 <= R("R_LPLIPAL1A181pp") <= 999999
0 <= R("R_LPLIPAL1A180pp") <= 999999
0 <= R("R_LPLIPAL1A141pp") <= 999999
0 <= R("R_LPLIPAL1A140pp") <= 999999
0 <= R("R_23PDE7pp") <= 999999
0 <= R("R_23PDE4pp") <= 999999
0 <= R("R_23PDE9pp") <= 999999
0 <= R("R_23PDE2pp") <= 999999
0 <= R("R_3NTD7pp") <= 999999
0 <= R("R_3NTD4pp") <= 999999
0 <= R("R_3NTD9pp") <= 999999
0 <= R("R_3NTD2pp") <= 999999
0 <= R("R_ACGAL1PPpp") <= 999999
0 <= R("R_ACGAM1PPpp") <= 999999
0 <= R("R_NTD7pp") <= 999999
0 <= R("R_AGM3PApp") <= 999999
0 <= R("R_AGM4PCPpp") <= 999999
0 <= R("R_AGM4PApp") <= 999999
0 <= R("R_ASNNpp") <= 999999
0 <= R("R_CLPNH120pp") <= 999999
0 <= R("R_CLPNH140pp") <= 999999
0 <= R("R_CLPNH141pp") <= 999999
0 <= R("R_CLPNH160pp") <= 999999
0 <= R("R_CLPNH161pp") <= 999999
0 <= R("R_CLPNH180pp") <= 999999
0 <= R("R_CLPNH181pp") <= 999999
0 <= R("R_NTD4pp") <= 999999
0 <= R("R_NTD6pp") <= 999999
0 <= R("R_NTD3pp") <= 999999
0 <= R("R_NTD8pp") <= 999999
0 <= R("R_NTD12pp") <= 999999
0 <= R("R_NTD5pp") <= 999999
0 <= R("R_NTD1pp") <= 999999
0 <= R("R_G1PPpp") <= 999999
0 <= R("R_GPDDA1pp") <= 999999
0 <= R("R_GPDDA2pp") <= 999999
0 <= R("R_GPDDA4pp") <= 999999
0 <= R("R_GPDDA5pp") <= 999999
0 <= R("R_GPDDA3pp") <= 999999
0 <= R("R_GAL1PPpp") <= 999999
0 <= R("R_GUR1PPpp") <= 999999
0 <= R("R_GLUNpp") <= 999999
0 <= R("R_G2PPpp") <= 999999
0 <= R("R_NTD9pp") <= 999999
0 <= R("R_GTHRDHpp") <= 999999
0 <= R("R_NTP3pp") <= 999999
0 <= R("R_H2St1pp") <= 999999
0 <= R("R_CYSS") <= 999999
0 <= R("R_DAPDC") <= 999999
0 <= R("R_AOBUTDs") <= 999999
0 <= R("R_IGPS") <= 999999
0 <= R("R_OMCDC") <= 999999
0 <= R("R_KG6PDC") <= 999999
0 <= R("R_OPHBDC") <= 999999
0 <= R("R_TMPPP") <= 999999
0 <= R("R_PPCDC") <= 999999
0 <= R("R_ADMDC") <= 999999
0 <= R("R_ARGDC") <= 999999
0 <= R("R_ASP1DC") <= 999999
0 <= R("R_DHBS") <= 999999
0 <= R("R_CLt3_2pp") <= 999999
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
0 <= R("R_ATPHs") <= 999999
0 <= R("R_CSND") <= 999999
0 <= R("R_CYTD") <= 999999
0 <= R("R_DADA") <= 999999
0 <= R("R_DATPHs") <= 999999
0 <= R("R_DCTPD") <= 999999
0 <= R("R_DCYTD") <= 999999
0 <= R("R_GTPHs") <= 999999
0 <= R("R_GUAD") <= 999999
0 <= R("R_CYSSADS") <= 999999
0 <= R("R_HYD3pp") <= 999999
0 <= R("R_ALLTAMH") <= 999999
0 <= R("R_ACCOAC") <= 999999
0 <= R("R_AIRC2") <= 999999
0 <= R("R_CBPS") <= 999999
0 <= R("R_CYNTAH") <= 999999
0 <= R("R_METS") <= 999999
0 <= R("R_HCYSMT") <= 999999
0 <= R("R_2AGPEAT160") <= 999999
0 <= R("R_2AGPGAT160") <= 999999
0 <= R("R_AACPS3") <= 999999
0 <= R("R_HDCAtexi") <= 999999
0 <= R("R_2AGPEAT161") <= 999999
0 <= R("R_2AGPGAT161") <= 999999
0 <= R("R_AACPS4") <= 999999
0 <= R("R_HDCEAtexi") <= 999999
0 <= R("R_HDCOAI") <= 999999
0 <= R("R_FACOAE161") <= 999999
0 <= R("R_AGPAT161") <= 999999
0 <= R("R_G3PAT161") <= 999999
0 <= R("R_FA161ACPHi") <= 999999
0 <= R("R_AGt3") <= 999999
0 <= R("R_CUt3") <= 999999
0 <= R("R_LIPAHT2ex") <= 999999
0 <= R("R_HG2abcpp") <= 999999
0 <= R("R_HEPK1") <= 999999
0 <= R("R_HISabcpp") <= 999999
0 <= R("R_HISTP") <= 999999
0 <= R("R_HKNDDH") <= 999999
0 <= R("R_HKNTDH") <= 999999
0 <= R("R_HEPT2") <= 999999
0 <= R("R_UPP3S") <= 999999
0 <= R("R_DM_HMFURN") <= 999999
0 <= R("R_HSK") <= 999999
0 <= R("R_23DAPPAt2pp") <= 999999
0 <= R("R_ABUTt2pp") <= 999999
0 <= R("R_ACNAMt2pp") <= 999999
0 <= R("R_ADNt2pp") <= 999999
0 <= R("R_ADOCBLtonex") <= 999999
0 <= R("R_BALAt2pp") <= 999999
0 <= R("R_DALAt2pp") <= 999999
0 <= R("R_ALAt2pp") <= 999999
0 <= R("R_AGM3Pt2pp") <= 999999
0 <= R("R_AGM4Pt2pp") <= 999999
0 <= R("R_AGMt2pp") <= 999999
0 <= R("R_ARBt3ipp") <= 999999
0 <= R("R_ARBTNtpp") <= 999999
0 <= R("R_ARBTNtonex") <= 999999
0 <= R("R_ARBTNtex") <= 999999
0 <= R("R_ARGt3pp") <= 999999
0 <= R("R_ARGDCpp") <= 999999
0 <= R("R_ASPt2pp") <= 999999
0 <= R("R_CA2t3pp") <= 999999
0 <= R("R_CBItonex") <= 999999
0 <= R("R_CBL1tonex") <= 999999
0 <= R("R_CD2t3pp") <= 999999
0 <= R("R_CHLt2pp") <= 999999
0 <= R("R_COBALT2t3pp") <= 999999
0 <= R("R_CPGNtonex") <= 999999
0 <= R("R_CPGNUtpp") <= 999999
0 <= R("R_CPGNUtex") <= 999999
0 <= R("R_CSNt2pp") <= 999999
0 <= R("R_CYNTt2pp") <= 999999
0 <= R("R_CYTDt2pp") <= 999999
0 <= R("R_DADNt2pp") <= 999999
0 <= R("R_FACOAL100t2pp") <= 999999
0 <= R("R_DCYTt2pp") <= 999999
0 <= R("R_FACOAL120t2pp") <= 999999
0 <= R("R_DGSNt2pp") <= 999999
0 <= R("R_DINSt2pp") <= 999999
0 <= R("R_DURIt2pp") <= 999999
0 <= R("R_FEENTERtpp") <= 999999
0 <= R("R_FEENTERtex") <= 999999
0 <= R("R_ETHAt2pp") <= 999999
0 <= R("R_FE2t3pp") <= 999999
0 <= R("R_FE2t2pp") <= 999999
0 <= R("R_FE3DCITtonex") <= 999999
0 <= R("R_FE3DHBZStonex") <= 999999
0 <= R("R_FE3HOXtonex") <= 999999
0 <= R("R_FE3HOXUtex") <= 999999
0 <= R("R_FE3HOXUtpp") <= 999999
0 <= R("R_FECRMtonex") <= 999999
0 <= R("R_FECRMUtpp") <= 999999
0 <= R("R_FECRMUtex") <= 999999
0 <= R("R_FEENTERtonex") <= 999999
0 <= R("R_FEOXAMtonex") <= 999999
0 <= R("R_FEOXAMUtpp") <= 999999
0 <= R("R_FEOXAMUtex") <= 999999
0 <= R("R_FORt2pp") <= 999999
0 <= R("R_FRULYSt2pp") <= 999999
0 <= R("R_GALCTNt2pp") <= 999999
0 <= R("R_GALCTNLt2pp") <= 999999
0 <= R("R_GALt2pp") <= 999999
0 <= R("R_GLCt2pp") <= 999999
0 <= R("R_GLYBt2pp") <= 999999
0 <= R("R_GLYt2pp") <= 999999
0 <= R("R_GSNt2pp") <= 999999
0 <= R("R_GUAt2pp") <= 999999
0 <= R("R_FACOAL160t2pp") <= 999999
0 <= R("R_FACOAL161t2pp") <= 999999
0 <= R("R_HG2t3pp") <= 999999
0 <= R("R_HOMt2pp") <= 999999
0 <= R("R_ASPt2_2pp") <= 999999
0 <= R("R_FUMt2_2pp") <= 999999
0 <= R("R_ASPt2_3pp") <= 999999
0 <= R("R_FUMt2_3pp") <= 999999
0 <= R("R_HEPK2") <= 999999
0 <= R("R_HXCT") <= 999999
0 <= R("R_FACOAL60t2pp") <= 999999
0 <= R("R_FACOAE60") <= 999999
0 <= R("R_QULNS") <= 999999
0 <= R("R_ICHORT") <= 999999
0 <= R("R_ICL") <= 999999
0 <= R("R_ILEabcpp") <= 999999
0 <= R("R_HSTPT") <= 999999
0 <= R("R_ADSS") <= 999999
0 <= R("R_NTD11") <= 999999
0 <= R("R_NTD11pp") <= 999999
0 <= R("R_INDOLEt2pp") <= 999999
0 <= R("R_INSK") <= 999999
0 <= R("R_INSH") <= 999999
0 <= R("R_INSt2pp") <= 999999
0 <= R("R_DMATT") <= 999999
0 <= R("R_GRTT") <= 999999
0 <= R("R_OCTDPS") <= 999999
0 <= R("R_UDCPDPS") <= 999999
0 <= R("R_ISETACabcpp") <= 999999
0 <= R("R_NTP10") <= 999999
0 <= R("R_NTPP9") <= 999999
0 <= R("R_Kt3pp") <= 999999
0 <= R("R_EDTXS1") <= 999999
0 <= R("R_K2L4Aabcpp") <= 999999
0 <= R("R_EDTXS3") <= 999999
0 <= R("R_K2L4Atex") <= 999999
0 <= R("R_KDOPP") <= 999999
0 <= R("R_KDOCT2") <= 999999
0 <= R("R_MOAT2") <= 999999
0 <= R("R_Kabcpp") <= 999999
0 <= R("R_Kt2pp") <= 999999
0 <= R("R_RHAT1") <= 999999
0 <= R("R_LADGMDH") <= 999999
0 <= R("R_4PCP") <= 999999
0 <= R("R_4PEPTabcpp") <= 999999
0 <= R("R_4PCPpp") <= 999999
0 <= R("R_3PEPTabcpp") <= 999999
0 <= R("R_LALGP") <= 999999
0 <= R("R_LACZ") <= 999999
0 <= R("R_LCTSt3ipp") <= 999999
0 <= R("R_LACZpp") <= 999999
0 <= R("R_LEUabcpp") <= 999999
0 <= R("R_GLYOX") <= 999999
0 <= R("R_HEPT1") <= 999999
0 <= R("R_LIPAabcpp") <= 999999
0 <= R("R_LIPACabcpp") <= 999999
0 <= R("R_CLIPAtex") <= 999999
0 <= R("R_LIPAHTex") <= 999999
0 <= R("R_LIPAtex") <= 999999
0 <= R("R_MOAT") <= 999999
0 <= R("R_TDSK") <= 999999
0 <= R("R_LYSDC") <= 999999
0 <= R("R_LYSt3pp") <= 999999
0 <= R("R_LYSabcpp") <= 999999
0 <= R("R_LYSt2pp") <= 999999
0 <= R("R_CADVtpp") <= 999999
0 <= R("R_LYXI") <= 999999
0 <= R("R_LYXt2pp") <= 999999
0 <= R("R_MACPD") <= 999999
0 <= R("R_KAS14") <= 999999
0 <= R("R_KAS15") <= 999999
0 <= R("R_3OAS60") <= 999999
0 <= R("R_3OAS141") <= 999999
0 <= R("R_3OAS121") <= 999999
0 <= R("R_3OAS120") <= 999999
0 <= R("R_3OAS140") <= 999999
0 <= R("R_3OAS181") <= 999999
0 <= R("R_3OAS80") <= 999999
0 <= R("R_MALDt2_2pp") <= 999999
0 <= R("R_MALt2_2pp") <= 999999
0 <= R("R_MALt2_3pp") <= 999999
0 <= R("R_MALTtexi") <= 999999
0 <= R("R_MLTG5") <= 999999
0 <= R("R_MLTG4") <= 999999
0 <= R("R_AMALT4") <= 999999
0 <= R("R_MALTHXtexi") <= 999999
0 <= R("R_MALTHXabcpp") <= 999999
0 <= R("R_MALTabcpp") <= 999999
0 <= R("R_MLTG3") <= 999999
0 <= R("R_AMALT3") <= 999999
0 <= R("R_MALTPTtexi") <= 999999
0 <= R("R_MALTPTabcpp") <= 999999
0 <= R("R_MLTG1") <= 999999
0 <= R("R_AMALT1") <= 999999
0 <= R("R_MALTTRtexi") <= 999999
0 <= R("R_MALTTRabcpp") <= 999999
0 <= R("R_MLTG2") <= 999999
0 <= R("R_AMALT2") <= 999999
0 <= R("R_MALTTTRtexi") <= 999999
0 <= R("R_MALTTTRabcpp") <= 999999
0 <= R("R_MAN1PT2") <= 999999
0 <= R("R_MN6PP") <= 999999
0 <= R("R_MANPGH") <= 999999
0 <= R("R_MNNH") <= 999999
0 <= R("R_HEX4") <= 999999
0 <= R("R_GALS3") <= 999999
0 <= R("R_MELIBt3ipp") <= 999999
0 <= R("R_MELIBt2pp") <= 999999
0 <= R("R_MCPST") <= 999999
0 <= R("R_METDabcpp") <= 999999
0 <= R("R_METOX2s") <= 999999
0 <= R("R_METOX1s") <= 999999
0 <= R("R_METAT") <= 999999
0 <= R("R_METabcpp") <= 999999
0 <= R("R_METSOX2abcpp") <= 999999
0 <= R("R_METSOX1abcpp") <= 999999
0 <= R("R_FMETTRS") <= 999999
0 <= R("R_MG2tpp") <= 999999
0 <= R("R_MG2uabcpp") <= 999999
0 <= R("R_MI1PP") <= 999999
0 <= R("R_MINOHPtexi") <= 999999
0 <= R("R_PHYTSpp") <= 999999
0 <= R("R_TMDS") <= 999999
0 <= R("R_MOHMT") <= 999999
0 <= R("R_MMCD") <= 999999
0 <= R("R_HCYSMT2") <= 999999
0 <= R("R_MMETt2pp") <= 999999
0 <= R("R_MN2t3pp") <= 999999
0 <= R("R_MN2tpp") <= 999999
0 <= R("R_MNt2pp") <= 999999
0 <= R("R_MOBDabcpp") <= 999999
0 <= R("R_DMSOR1") <= 999999
0 <= R("R_DMSOR1pp") <= 999999
0 <= R("R_FRD2") <= 999999
0 <= R("R_ASPO4") <= 999999
0 <= R("R_DHORD5") <= 999999
0 <= R("R_DSBAO2") <= 999999
0 <= R("R_G3PD6") <= 999999
0 <= R("R_GLYCTO3") <= 999999
0 <= R("R_FDH5pp") <= 999999
0 <= R("R_HYD2pp") <= 999999
0 <= R("R_L_LACD3") <= 999999
0 <= R("R_MDH3") <= 999999
0 <= R("R_MSO3abcpp") <= 999999
0 <= R("R_LGTHL") <= 999999
0 <= R("R_GLYOX3") <= 999999
0 <= R("R_MLTGY3pp") <= 999999
0 <= R("R_MLDEP1pp") <= 999999
0 <= R("R_MDDEP2pp") <= 999999
0 <= R("R_MLTGY2pp") <= 999999
0 <= R("R_MLDCP2Bpp") <= 999999
0 <= R("R_MLTGY1pp") <= 999999
0 <= R("R_MLDCP1Bpp") <= 999999
0 <= R("R_MLTGY4pp") <= 999999
0 <= R("R_MDDEP1pp") <= 999999
0 <= R("R_MDDEP4pp") <= 999999
0 <= R("R_MDDCP5pp") <= 999999
0 <= R("R_MDDCP4pp") <= 999999
0 <= R("R_MCTP2App") <= 999999
0 <= R("R_MCTP1Bpp") <= 999999
0 <= R("R_MCTP1App") <= 999999
0 <= R("R_MLDCP2App") <= 999999
0 <= R("R_MDDCP3pp") <= 999999
0 <= R("R_MLDCP3App") <= 999999
0 <= R("R_MLDEP2pp") <= 999999
0 <= R("R_MLDCP1App") <= 999999
0 <= R("R_MDDCP1pp") <= 999999
0 <= R("R_MDDEP3pp") <= 999999
0 <= R("R_MDDCP2pp") <= 999999
0 <= R("R_AGPAT140") <= 999999
0 <= R("R_G3PAT140") <= 999999
0 <= R("R_FA140ACPHi") <= 999999
0 <= R("R_EDTXS2") <= 999999
0 <= R("R_EDTXS4") <= 999999
0 <= R("R_3OAS160") <= 999999
0 <= R("R_NAt3pp") <= 999999
0 <= R("R_NAt3_2pp") <= 999999
0 <= R("R_NAt3_1p5pp") <= 999999
0 <= R("R_ACt4pp") <= 999999
0 <= R("R_ALAt4pp") <= 999999
0 <= R("R_GLUt4pp") <= 999999
0 <= R("R_GLYCLTt4pp") <= 999999
0 <= R("R_GLYt4pp") <= 999999
0 <= R("R_INOSTt4pp") <= 999999
0 <= R("R_NACtpp") <= 999999
0 <= R("R_IPMD") <= 999999
0 <= R("R_HACD1i") <= 999999
0 <= R("R_HACD4i") <= 999999
0 <= R("R_HACD5i") <= 999999
0 <= R("R_HACD2i") <= 999999
0 <= R("R_HACD7i") <= 999999
0 <= R("R_HACD3i") <= 999999
0 <= R("R_HACD8i") <= 999999
0 <= R("R_HACD6i") <= 999999
0 <= R("R_PGCD") <= 999999
0 <= R("R_NADK") <= 999999
0 <= R("R_DHPPD") <= 999999
0 <= R("R_DHCIND") <= 999999
0 <= R("R_AKGDH") <= 999999
0 <= R("R_SHCHD2") <= 999999
0 <= R("R_GALCTLO") <= 999999
0 <= R("R_GLYCDx") <= 999999
0 <= R("R_NADN") <= 999999
0 <= R("R_NADDP") <= 999999
0 <= R("R_ABUTD") <= 999999
0 <= R("R_ALDD2x") <= 999999
0 <= R("R_AMPMS2") <= 999999
0 <= R("R_BETALDHx") <= 999999
0 <= R("R_ALDD4") <= 999999
0 <= R("R_GCALDD") <= 999999
0 <= R("R_P5CD") <= 999999
0 <= R("R_OMMBLHX3") <= 999999
0 <= R("R_OMPHHX3") <= 999999
0 <= R("R_OPHHX3") <= 999999
0 <= R("R_HXAND") <= 999999
0 <= R("R_IMPD") <= 999999
0 <= R("R_LCADi") <= 999999
0 <= R("R_MALDDH") <= 999999
0 <= R("R_ME1") <= 999999
0 <= R("R_DKMPPD3") <= 999999
0 <= R("R_HISTD") <= 999999
0 <= R("R_MECDPDH2") <= 999999
0 <= R("R_DOGULNR") <= 999999
0 <= R("R_DKGLCNR2x") <= 999999
0 <= R("R_2DGLCNRx") <= 999999
0 <= R("R_2DGULRx") <= 999999
0 <= R("R_NADH9") <= 999999
0 <= R("R_ALR4x") <= 999999
0 <= R("R_BSORx") <= 999999
0 <= R("R_EAR40x") <= 999999
0 <= R("R_FADRx") <= 999999
0 <= R("R_FMNRx") <= 999999
0 <= R("R_GLYCLTDx") <= 999999
0 <= R("R_DMPPS") <= 999999
0 <= R("R_IPDPS") <= 999999
0 <= R("R_NADH18pp") <= 999999
0 <= R("R_HPYRRx") <= 999999
0 <= R("R_MTHFR2") <= 999999
0 <= R("R_NADH10") <= 999999
0 <= R("R_NADH17pp") <= 999999
0 <= R("R_LALDO2x") <= 999999
0 <= R("R_GND") <= 999999
0 <= R("R_NADPPPS") <= 999999
0 <= R("R_ALDD2y") <= 999999
0 <= R("R_BETALDHy") <= 999999
0 <= R("R_ME2") <= 999999
0 <= R("R_THD2pp") <= 999999
0 <= R("R_DHDPRy") <= 999999
0 <= R("R_DKGLCNR1") <= 999999
0 <= R("R_DKGLCNR2y") <= 999999
0 <= R("R_2DGLCNRy") <= 999999
0 <= R("R_2DGULRy") <= 999999
0 <= R("R_DPR") <= 999999
0 <= R("R_NADPHQR4") <= 999999
0 <= R("R_3OAR121") <= 999999
0 <= R("R_3OAR141") <= 999999
0 <= R("R_3OAR161") <= 999999
0 <= R("R_3OAR181") <= 999999
0 <= R("R_APRAUR") <= 999999
0 <= R("R_IDOND2") <= 999999
0 <= R("R_BSORy") <= 999999
0 <= R("R_EAR40y") <= 999999
0 <= R("R_TDPDRR") <= 999999
0 <= R("R_DXPRIi") <= 999999
0 <= R("R_FADRx2") <= 999999
0 <= R("R_FLDR") <= 999999
0 <= R("R_FMNRx2") <= 999999
0 <= R("R_GOFUCR") <= 999999
0 <= R("R_GLUSy") <= 999999
0 <= R("R_G5SD") <= 999999
0 <= R("R_GLUTRR") <= 999999
0 <= R("R_GLYCLTDy") <= 999999
0 <= R("R_P5CR") <= 999999
0 <= R("R_GMPR") <= 999999
0 <= R("R_HPYRRy") <= 999999
0 <= R("R_NADPHQR3") <= 999999
0 <= R("R_ALR2") <= 999999
0 <= R("R_NADTRHD") <= 999999
0 <= R("R_NNAM") <= 999999
0 <= R("R_ASNS2") <= 999999
0 <= R("R_NADS1") <= 999999
0 <= R("R_GLNS") <= 999999
0 <= R("R_NI2abcpp") <= 999999
0 <= R("R_NI2t3pp") <= 999999
0 <= R("R_NI2tpp") <= 999999
0 <= R("R_NI2uabcpp") <= 999999
0 <= R("R_NNDMBRT") <= 999999
0 <= R("R_NMNDA") <= 999999
0 <= R("R_NMNN") <= 999999
0 <= R("R_NMNAT") <= 999999
0 <= R("R_NMNPtpp") <= 999999
0 <= R("R_NMNt7pp") <= 999999
0 <= R("R_NTRIR2x") <= 999999
0 <= R("R_NTRIR4pp") <= 999999
0 <= R("R_NO3R2pp") <= 999999
0 <= R("R_NO3R2bpp") <= 999999
0 <= R("R_NO3t7pp") <= 999999
0 <= R("R_NHFRBO") <= 999999
0 <= R("R_O16A4COLIPAtex") <= 999999
0 <= R("R_O16A4Lpp") <= 999999
0 <= R("R_O16AUNDtpp") <= 999999
0 <= R("R_O16AP1pp") <= 999999
0 <= R("R_O16AP2pp") <= 999999
0 <= R("R_O16AP3pp") <= 999999
0 <= R("R_OMMBLHX") <= 999999
0 <= R("R_OMPHHX") <= 999999
0 <= R("R_OPHHX") <= 999999
0 <= R("R_CYTBD2pp") <= 999999
0 <= R("R_ASPO6") <= 999999
0 <= R("R_DHCINDO") <= 999999
0 <= R("R_HPPPNDO") <= 999999
0 <= R("R_FDMO4") <= 999999
0 <= R("R_FDMO3") <= 999999
0 <= R("R_AACTOOR") <= 999999
0 <= R("R_DKMPPD") <= 999999
0 <= R("R_GGPTRCO") <= 999999
0 <= R("R_CPPPGO") <= 999999
0 <= R("R_FDMO") <= 999999
0 <= R("R_FDMO2") <= 999999
0 <= R("R_3HCINNMH") <= 999999
0 <= R("R_3HPPPNH") <= 999999
0 <= R("R_CINNDO") <= 999999
0 <= R("R_MTRPOX") <= 999999
0 <= R("R_H2SO") <= 999999
0 <= R("R_QMO3") <= 999999
0 <= R("R_NODOx") <= 999999
0 <= R("R_NODOy") <= 999999
0 <= R("R_42A12BOOXpp") <= 999999
0 <= R("R_CU1Opp") <= 999999
0 <= R("R_FEROpp") <= 999999
0 <= R("R_SPODM") <= 999999
0 <= R("R_SPODMpp") <= 999999
0 <= R("R_PPCK") <= 999999
0 <= R("R_CS") <= 999999
0 <= R("R_OAADC") <= 999999
0 <= R("R_FA80ACPHi") <= 999999
0 <= R("R_3OAS100") <= 999999
0 <= R("R_FACOAE80") <= 999999
0 <= R("R_AGPAT180") <= 999999
0 <= R("R_G3PAT180") <= 999999
0 <= R("R_2AGPEAT180") <= 999999
0 <= R("R_2AGPGAT180") <= 999999
0 <= R("R_AACPS6") <= 999999
0 <= R("R_OCDCAtexi") <= 999999
0 <= R("R_FACOAL180t2pp") <= 999999
0 <= R("R_2AGPEAT181") <= 999999
0 <= R("R_2AGPGAT181") <= 999999
0 <= R("R_AACPS5") <= 999999
0 <= R("R_OCDCEAtexi") <= 999999
0 <= R("R_FACOAL181t2pp") <= 999999
0 <= R("R_AACPS9") <= 999999
0 <= R("R_FACOAL80t2pp") <= 999999
0 <= R("R_HBZOPT") <= 999999
0 <= R("R_DHNAOT4") <= 999999
0 <= R("R_AGPAT181") <= 999999
0 <= R("R_G3PAT181") <= 999999
0 <= R("R_ODECOAI") <= 999999
0 <= R("R_FACOAE181") <= 999999
0 <= R("R_OP4ENH") <= 999999
0 <= R("R_ORNDC") <= 999999
0 <= R("R_ORNabcpp") <= 999999
0 <= R("R_OMPDC") <= 999999
0 <= R("R_OROTt2_2pp") <= 999999
0 <= R("R_DM_OXAM") <= 999999
0 <= R("R_PAPA120") <= 999999
0 <= R("R_PA120abcpp") <= 999999
0 <= R("R_DASYN120") <= 999999
0 <= R("R_PLIPA1A120pp") <= 999999
0 <= R("R_PLIPA2A120pp") <= 999999
0 <= R("R_PAPA120pp") <= 999999
0 <= R("R_PAPA140") <= 999999
0 <= R("R_PA140abcpp") <= 999999
0 <= R("R_DASYN140") <= 999999
0 <= R("R_PAPA140pp") <= 999999
0 <= R("R_PLIPA1A140pp") <= 999999
0 <= R("R_PLIPA2A140pp") <= 999999
0 <= R("R_PAPA141") <= 999999
0 <= R("R_PA141abcpp") <= 999999
0 <= R("R_DASYN141") <= 999999
0 <= R("R_PAPA141pp") <= 999999
0 <= R("R_PLIPA1A141pp") <= 999999
0 <= R("R_PLIPA2A141pp") <= 999999
0 <= R("R_PAPA160") <= 999999
0 <= R("R_PA160abcpp") <= 999999
0 <= R("R_DASYN160") <= 999999
0 <= R("R_PLIPA1A160pp") <= 999999
0 <= R("R_PLIPA2A160pp") <= 999999
0 <= R("R_PAPA160pp") <= 999999
0 <= R("R_PAPA161") <= 999999
0 <= R("R_PA161abcpp") <= 999999
0 <= R("R_DASYN161") <= 999999
0 <= R("R_PLIPA1A161pp") <= 999999
0 <= R("R_PLIPA2A161pp") <= 999999
0 <= R("R_PAPA161pp") <= 999999
0 <= R("R_PAPA180") <= 999999
0 <= R("R_PA180abcpp") <= 999999
0 <= R("R_DASYN180") <= 999999
0 <= R("R_PLIPA1A180pp") <= 999999
0 <= R("R_PLIPA2A180pp") <= 999999
0 <= R("R_PAPA180pp") <= 999999
0 <= R("R_PAPA181") <= 999999
0 <= R("R_PA181abcpp") <= 999999
0 <= R("R_DASYN181") <= 999999
0 <= R("R_PLIPA1A181pp") <= 999999
0 <= R("R_PLIPA2A181pp") <= 999999
0 <= R("R_PAPA181pp") <= 999999
0 <= R("R_ALDD19x") <= 999999
0 <= R("R_PACCOAL") <= 999999
0 <= R("R_AGPAT160") <= 999999
0 <= R("R_G3PAT160") <= 999999
0 <= R("R_FA160ACPHi") <= 999999
0 <= R("R_3OAS180") <= 999999
0 <= R("R_PTPATi") <= 999999
0 <= R("R_PANTS") <= 999999
0 <= R("R_BPNT") <= 999999
0 <= R("R_PAPSR2") <= 999999
0 <= R("R_PDXPP") <= 999999
0 <= R("R_PDX5PO2") <= 999999
0 <= R("R_PDX5POi") <= 999999
0 <= R("R_PE120abcpp") <= 999999
0 <= R("R_PLIPA2E120pp") <= 999999
0 <= R("R_PLIPA1E120pp") <= 999999
0 <= R("R_PE140abcpp") <= 999999
0 <= R("R_PLIPA2E140pp") <= 999999
0 <= R("R_PLIPA1E140pp") <= 999999
0 <= R("R_PE141abcpp") <= 999999
0 <= R("R_PLIPA2E141pp") <= 999999
0 <= R("R_PLIPA1E141pp") <= 999999
0 <= R("R_PE160abcpp") <= 999999
0 <= R("R_ALPATE160pp") <= 999999
0 <= R("R_PLIPA2E160pp") <= 999999
0 <= R("R_PLIPA1E160pp") <= 999999
0 <= R("R_CFAS160E") <= 999999
0 <= R("R_PE161abcpp") <= 999999
0 <= R("R_PLIPA2E161pp") <= 999999
0 <= R("R_PLIPA1E161pp") <= 999999
0 <= R("R_PETNT161pp") <= 999999
0 <= R("R_PE180abcpp") <= 999999
0 <= R("R_PLIPA2E180pp") <= 999999
0 <= R("R_PLIPA1E180pp") <= 999999
0 <= R("R_CFAS180E") <= 999999
0 <= R("R_PE181abcpp") <= 999999
0 <= R("R_PLIPA2E181pp") <= 999999
0 <= R("R_PLIPA1E181pp") <= 999999
0 <= R("R_PETNT181pp") <= 999999
0 <= R("R_PEAMNOpp") <= 999999
0 <= R("R_ACGAptspp") <= 999999
0 <= R("R_ACMANAptspp") <= 999999
0 <= R("R_ACMUMptspp") <= 999999
0 <= R("R_ASCBptspp") <= 999999
0 <= R("R_DHAPT") <= 999999
0 <= R("R_FRUptspp") <= 999999
0 <= R("R_FRUpts2pp") <= 999999
0 <= R("R_GALTptspp") <= 999999
0 <= R("R_GAMptspp") <= 999999
0 <= R("R_GLCptspp") <= 999999
0 <= R("R_KDOPS") <= 999999
0 <= R("R_PPC") <= 999999
0 <= R("R_DDPA") <= 999999
0 <= R("R_PYK") <= 999999
0 <= R("R_MALTptspp") <= 999999
0 <= R("R_MANGLYCptspp") <= 999999
0 <= R("R_MANptspp") <= 999999
0 <= R("R_MNLptspp") <= 999999
0 <= R("R_LPLIPAL2ATE120") <= 999999
0 <= R("R_LPLIPAL2ATG120") <= 999999
0 <= R("R_PG120abcpp") <= 999999
0 <= R("R_PLIPA2G120pp") <= 999999
0 <= R("R_PLIPA1G120pp") <= 999999
0 <= R("R_LPLIPAL2ATE140") <= 999999
0 <= R("R_LPLIPAL2ATG140") <= 999999
0 <= R("R_PG140abcpp") <= 999999
0 <= R("R_PLIPA2G140pp") <= 999999
0 <= R("R_PLIPA1G140pp") <= 999999
0 <= R("R_LPLIPAL2ATE141") <= 999999
0 <= R("R_LPLIPAL2ATG141") <= 999999
0 <= R("R_PG141abcpp") <= 999999
0 <= R("R_PLIPA2G141pp") <= 999999
0 <= R("R_PLIPA1G141pp") <= 999999
0 <= R("R_LPLIPAL2ATE160") <= 999999
0 <= R("R_LPLIPAL2ATG160") <= 999999
0 <= R("R_PG160abcpp") <= 999999
0 <= R("R_ALPATG160pp") <= 999999
0 <= R("R_PLIPA2G160pp") <= 999999
0 <= R("R_PLIPA1G160pp") <= 999999
0 <= R("R_LPLIPAL2ATE161") <= 999999
0 <= R("R_LPLIPAL2ATG161") <= 999999
0 <= R("R_CFAS160G") <= 999999
0 <= R("R_PG161abcpp") <= 999999
0 <= R("R_PLIPA2G161pp") <= 999999
0 <= R("R_PLIPA1G161pp") <= 999999
0 <= R("R_LPLIPAL2ATE180") <= 999999
0 <= R("R_LPLIPAL2ATG180") <= 999999
0 <= R("R_PG180abcpp") <= 999999
0 <= R("R_PLIPA2G180pp") <= 999999
0 <= R("R_PLIPA1G180pp") <= 999999
0 <= R("R_LPLIPAL2ATE181") <= 999999
0 <= R("R_LPLIPAL2ATG181") <= 999999
0 <= R("R_CFAS180G") <= 999999
0 <= R("R_PG181abcpp") <= 999999
0 <= R("R_PLIPA2G181pp") <= 999999
0 <= R("R_PLIPA1G181pp") <= 999999
0 <= R("R_PGPP120") <= 999999
0 <= R("R_PGP120abcpp") <= 999999
0 <= R("R_PGPP120pp") <= 999999
0 <= R("R_PGPP140") <= 999999
0 <= R("R_PGP140abcpp") <= 999999
0 <= R("R_PGPP140pp") <= 999999
0 <= R("R_PGPP141") <= 999999
0 <= R("R_PGP141abcpp") <= 999999
0 <= R("R_PGPP141pp") <= 999999
0 <= R("R_PGPP160") <= 999999
0 <= R("R_PGP160abcpp") <= 999999
0 <= R("R_PGPP160pp") <= 999999
0 <= R("R_PGPP161") <= 999999
0 <= R("R_PGP161abcpp") <= 999999
0 <= R("R_PGPP161pp") <= 999999
0 <= R("R_PGPP180") <= 999999
0 <= R("R_PGP180abcpp") <= 999999
0 <= R("R_PGPP180pp") <= 999999
0 <= R("R_PGPP181") <= 999999
0 <= R("R_PGP181abcpp") <= 999999
0 <= R("R_PGPP181pp") <= 999999
0 <= R("R_PHEMEabcpp") <= 999999
0 <= R("R_HEMEOS") <= 999999
0 <= R("R_PHEMEtiex") <= 999999
0 <= R("R_HEPT3") <= 999999
0 <= R("R_THRS") <= 999999
0 <= R("R_MOAT3C") <= 999999
0 <= R("R_4HTHRS") <= 999999
0 <= R("R_PDX5PS") <= 999999
0 <= R("R_GLCP2") <= 999999
0 <= R("R_GLYC3Pt6pp") <= 999999
0 <= R("R_GLCP") <= 999999
0 <= R("R_OXAMTC") <= 999999
0 <= R("R_F6Pt6_2pp") <= 999999
0 <= R("R_G6Pt6_2pp") <= 999999
0 <= R("R_GAM6Pt6_2pp") <= 999999
0 <= R("R_MAN6Pt6_2pp") <= 999999
0 <= R("R_PIuabcpp") <= 999999
0 <= R("R_FACOAE160") <= 999999
0 <= R("R_PNTK") <= 999999
0 <= R("R_PNTOt4pp") <= 999999
0 <= R("R_ACCOAL") <= 999999
0 <= R("R_ALDD3y") <= 999999
0 <= R("R_PPAt4pp") <= 999999
0 <= R("R_HMBS") <= 999999
0 <= R("R_MCITS") <= 999999
0 <= R("R_PTA2") <= 999999
0 <= R("R_PPGPPDP") <= 999999
0 <= R("R_PPNDH") <= 999999
0 <= R("R_PPND") <= 999999
0 <= R("R_PPA") <= 999999
0 <= R("R_FCLT") <= 999999
0 <= R("R_PPPGO3") <= 999999
0 <= R("R_PPPGO") <= 999999
0 <= R("R_PPA2") <= 999999
0 <= R("R_PPPNDO") <= 999999
0 <= R("R_PPTHpp") <= 999999
0 <= R("R_PRAIi") <= 999999
0 <= R("R_PRAMPC") <= 999999
0 <= R("R_PRATPP") <= 999999
0 <= R("R_IG3PS") <= 999999
0 <= R("R_AMPTASEPG") <= 999999
0 <= R("R_PROGLYabcpp") <= 999999
0 <= R("R_PROD2") <= 999999
0 <= R("R_PROabcpp") <= 999999
0 <= R("R_PROt4pp") <= 999999
0 <= R("R_ADPT") <= 999999
0 <= R("R_ANPRT") <= 999999
0 <= R("R_ATPPRT") <= 999999
0 <= R("R_GUAPRT") <= 999999
0 <= R("R_GLUPRT") <= 999999
0 <= R("R_HXPRT") <= 999999
0 <= R("R_NAMNPP") <= 999999
0 <= R("R_PSD120") <= 999999
0 <= R("R_PSD140") <= 999999
0 <= R("R_PSD141") <= 999999
0 <= R("R_PSD160") <= 999999
0 <= R("R_PSD161") <= 999999
0 <= R("R_PSD180") <= 999999
0 <= R("R_PSD181") <= 999999
0 <= R("R_PSCLYSt2pp") <= 999999
0 <= R("R_PSP_L") <= 999999
0 <= R("R_PSP_Lpp") <= 999999
0 <= R("R_PTRCTA") <= 999999
0 <= R("R_SPMS") <= 999999
0 <= R("R_GGPTRCS") <= 999999
0 <= R("R_PTRCabcpp") <= 999999
0 <= R("R_PTRCt2pp") <= 999999
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
0 <= R("R_DSBAO1") <= 999999
0 <= R("R_G3PD5") <= 999999
0 <= R("R_GLYCTO2") <= 999999
0 <= R("R_GLCDpp") <= 999999
0 <= R("R_FDH4pp") <= 999999
0 <= R("R_HYD1pp") <= 999999
0 <= R("R_LDH_D2") <= 999999
0 <= R("R_L_LACD2") <= 999999
0 <= R("R_MDH2") <= 999999
0 <= R("R_NADH5") <= 999999
0 <= R("R_NADH16pp") <= 999999
0 <= R("R_NADPHQR2") <= 999999
0 <= R("R_POX") <= 999999
0 <= R("R_NO3R1pp") <= 999999
0 <= R("R_NO3R1bpp") <= 999999
0 <= R("R_CYTBDpp") <= 999999
0 <= R("R_CYTBO3_4pp") <= 999999
0 <= R("R_QMO2") <= 999999
0 <= R("R_NTRIR3pp") <= 999999
0 <= R("R_NNDPR") <= 999999
0 <= R("R_R15BPK") <= 999999
0 <= R("R_R1PK") <= 999999
0 <= R("R_R5PP") <= 999999
0 <= R("R_R5PPpp") <= 999999
0 <= R("R_O16AT") <= 999999
0 <= R("R_ARBTNR3") <= 999999
0 <= R("R_CPGNR3") <= 999999
0 <= R("R_FE3HOXR3") <= 999999
0 <= R("R_FECRMR3") <= 999999
0 <= R("R_FEENTERR3") <= 999999
0 <= R("R_FEOXAMR3") <= 999999
0 <= R("R_RBK_L1") <= 999999
0 <= R("R_ADOCBLS") <= 999999
0 <= R("R_RHCCE") <= 999999
0 <= R("R_RBK") <= 999999
0 <= R("R_RIBabcpp") <= 999999
0 <= R("R_RBFK") <= 999999
0 <= R("R_FLVRx") <= 999999
0 <= R("R_FLVR") <= 999999
0 <= R("R_RMK") <= 999999
0 <= R("R_RMNtpp") <= 999999
0 <= R("R_DB4PS") <= 999999
0 <= R("R_S7PI") <= 999999
0 <= R("R_SARCOX") <= 999999
0 <= R("R_SBTptspp") <= 999999
0 <= R("R_NPHS") <= 999999
0 <= R("R_BTS4") <= 999999
0 <= R("R_SHCHF") <= 999999
0 <= R("R_SELNPS") <= 999999
0 <= R("R_ENTCS") <= 999999
0 <= R("R_SERD_D") <= 999999
0 <= R("R_DSERt2pp") <= 999999
0 <= R("R_SERD_L") <= 999999
0 <= R("R_TRPS1") <= 999999
0 <= R("R_PSSA120") <= 999999
0 <= R("R_PSSA161") <= 999999
0 <= R("R_PSSA160") <= 999999
0 <= R("R_PSSA181") <= 999999
0 <= R("R_PSSA180") <= 999999
0 <= R("R_PSSA141") <= 999999
0 <= R("R_PSSA140") <= 999999
0 <= R("R_TRPS2") <= 999999
0 <= R("R_SERt4pp") <= 999999
0 <= R("R_SELCYSS") <= 999999
0 <= R("R_SFGTHi") <= 999999
0 <= R("R_SHKK") <= 999999
0 <= R("R_SKMt2pp") <= 999999
0 <= R("R_SDPDS") <= 999999
0 <= R("R_SULRi") <= 999999
0 <= R("R_SADT2") <= 999999
0 <= R("R_SULabcpp") <= 999999
0 <= R("R_SPMDAT1") <= 999999
0 <= R("R_SPMDAT2") <= 999999
0 <= R("R_GSPMDS") <= 999999
0 <= R("R_SPMDabcpp") <= 999999
0 <= R("R_SHCHCS2") <= 999999
0 <= R("R_FACOAE180") <= 999999
0 <= R("R_FFSD") <= 999999
0 <= R("R_SADH") <= 999999
0 <= R("R_SUCBZL") <= 999999
0 <= R("R_CITt7pp") <= 999999
0 <= R("R_SUCCt3pp") <= 999999
0 <= R("R_PPCSCT") <= 999999
0 <= R("R_SUCDi") <= 999999
0 <= R("R_MMM2") <= 999999
0 <= R("R_AST") <= 999999
0 <= R("R_HSST") <= 999999
0 <= R("R_SUCCt2_2pp") <= 999999
0 <= R("R_SUCCt2_3pp") <= 999999
0 <= R("R_SGDS") <= 999999
0 <= R("R_SGSAD") <= 999999
0 <= R("R_SHSL1") <= 999999
0 <= R("R_SOTA") <= 999999
0 <= R("R_SUCptspp") <= 999999
0 <= R("R_SSALx") <= 999999
0 <= R("R_SSALy") <= 999999
0 <= R("R_FDMO6") <= 999999
0 <= R("R_SULFACabcpp") <= 999999
0 <= R("R_EAR181x") <= 999999
0 <= R("R_EAR181y") <= 999999
0 <= R("R_EAR121x") <= 999999
0 <= R("R_EAR121y") <= 999999
0 <= R("R_EAR141x") <= 999999
0 <= R("R_EAR141y") <= 999999
0 <= R("R_EAR161x") <= 999999
0 <= R("R_EAR161y") <= 999999
0 <= R("R_PFK_2") <= 999999
0 <= R("R_TARTD") <= 999999
0 <= R("R_TAUDO") <= 999999
0 <= R("R_TAURabcpp") <= 999999
0 <= R("R_FACOAE140") <= 999999
0 <= R("R_EAR120x") <= 999999
0 <= R("R_EAR120y") <= 999999
0 <= R("R_AGPAT141") <= 999999
0 <= R("R_G3PAT141") <= 999999
0 <= R("R_FA141ACPHi") <= 999999
0 <= R("R_3OAS161") <= 999999
0 <= R("R_EAR100x") <= 999999
0 <= R("R_EAR100y") <= 999999
0 <= R("R_TDECOAI") <= 999999
0 <= R("R_FACOAE141") <= 999999
0 <= R("R_THDPS") <= 999999
0 <= R("R_EAR60x") <= 999999
0 <= R("R_EAR60y") <= 999999
0 <= R("R_GLYCL") <= 999999
0 <= R("R_TMK") <= 999999
0 <= R("R_TMPK") <= 999999
0 <= R("R_THMabcpp") <= 999999
0 <= R("R_TDP") <= 999999
0 <= R("R_OXGDC2") <= 999999
0 <= R("R_THRAi") <= 999999
0 <= R("R_THRD_L") <= 999999
0 <= R("R_THRt2pp") <= 999999
0 <= R("R_THRD") <= 999999
0 <= R("R_THRabcpp") <= 999999
0 <= R("R_THRt4pp") <= 999999
0 <= R("R_PTHRpp") <= 999999
0 <= R("R_THYMt3pp") <= 999999
0 <= R("R_TMDK1") <= 999999
0 <= R("R_THMDt2pp") <= 999999
0 <= R("R_TMAOR2") <= 999999
0 <= R("R_TMAOR1") <= 999999
0 <= R("R_TMAOR2pp") <= 999999
0 <= R("R_TMAOR1pp") <= 999999
0 <= R("R_EAR140x") <= 999999
0 <= R("R_EAR140y") <= 999999
0 <= R("R_EAR80x") <= 999999
0 <= R("R_EAR80y") <= 999999
0 <= R("R_EAR180x") <= 999999
0 <= R("R_EAR180y") <= 999999
0 <= R("R_EAR160x") <= 999999
0 <= R("R_EAR160y") <= 999999
0 <= R("R_TRDR") <= 999999
0 <= R("R_RNDR1") <= 999999
0 <= R("R_RNDR3") <= 999999
0 <= R("R_DSBDR") <= 999999
0 <= R("R_RNDR2") <= 999999
0 <= R("R_THIORDXi") <= 999999
0 <= R("R_METSOXR2") <= 999999
0 <= R("R_METSOXR1") <= 999999
0 <= R("R_PAPSR") <= 999999
0 <= R("R_TRE6PH") <= 999999
0 <= R("R_TRE6PP") <= 999999
0 <= R("R_TREH") <= 999999
0 <= R("R_TREHpp") <= 999999
0 <= R("R_TREptspp") <= 999999
0 <= R("R_ALATRS") <= 999999
0 <= R("R_ARGTRS") <= 999999
0 <= R("R_ASNTRS") <= 999999
0 <= R("R_ASPTRS") <= 999999
0 <= R("R_CYSTRS") <= 999999
0 <= R("R_GLNTRS") <= 999999
0 <= R("R_GLUTRS") <= 999999
0 <= R("R_GLYTRS") <= 999999
0 <= R("R_HISTRS") <= 999999
0 <= R("R_ILETRS") <= 999999
0 <= R("R_LEUTRS") <= 999999
0 <= R("R_LYSTRS") <= 999999
0 <= R("R_METTRS") <= 999999
0 <= R("R_PHETRS") <= 999999
0 <= R("R_PROTRS") <= 999999
0 <= R("R_SERTRS2") <= 999999
0 <= R("R_SERTRS") <= 999999
0 <= R("R_THRTRS") <= 999999
0 <= R("R_TRPTRS") <= 999999
0 <= R("R_CYANST") <= 999999
0 <= R("R_CYANSTpp") <= 999999
0 <= R("R_TSULabcpp") <= 999999
0 <= R("R_2AGPEAT140") <= 999999
0 <= R("R_2AGPGAT140") <= 999999
0 <= R("R_AACPS1") <= 999999
0 <= R("R_TTDCAtexi") <= 999999
0 <= R("R_FACOAL140t2pp") <= 999999
0 <= R("R_2AGPEAT141") <= 999999
0 <= R("R_2AGPGAT141") <= 999999
0 <= R("R_AACPS2") <= 999999
0 <= R("R_TTDCEAtexi") <= 999999
0 <= R("R_FACOAL141t2pp") <= 999999
0 <= R("R_TUNGSabcpp") <= 999999
0 <= R("R_TYROXDApp") <= 999999
0 <= R("R_THZPSN") <= 999999
0 <= R("R_TYRTRS") <= 999999
0 <= R("R_TYRPpp") <= 999999
0 <= R("R_USHD") <= 999999
0 <= R("R_LPADSS") <= 999999
0 <= R("R_UHGADA") <= 999999
0 <= R("R_U23GAAT") <= 999999
0 <= R("R_MPTG") <= 999999
0 <= R("R_MPTG2") <= 999999
0 <= R("R_UAPGR") <= 999999
0 <= R("R_UAGCVT") <= 999999
0 <= R("R_UACGAMPpp") <= 999999
0 <= R("R_UACMAMO") <= 999999
0 <= R("R_UAGPT3") <= 999999
0 <= R("R_UAMAGS") <= 999999
0 <= R("R_UAAGDS") <= 999999
0 <= R("R_UAMAS") <= 999999
0 <= R("R_UM3PL") <= 999999
0 <= R("R_UM4PL") <= 999999
0 <= R("R_UDCPDP") <= 999999
0 <= R("R_UDCPDPpp") <= 999999
0 <= R("R_ACGAMT") <= 999999
0 <= R("R_UDCPPtppi") <= 999999
0 <= R("R_UACGALPpp") <= 999999
0 <= R("R_RNDR4b") <= 999999
0 <= R("R_RNDR4") <= 999999
0 <= R("R_UDPGALM") <= 999999
0 <= R("R_O16GALFT") <= 999999
0 <= R("R_UDPGALPpp") <= 999999
0 <= R("R_O16GLCT1") <= 999999
0 <= R("R_TRE6PS") <= 999999
0 <= R("R_GLCTR2") <= 999999
0 <= R("R_O16GLCT2") <= 999999
0 <= R("R_GLCTR3") <= 999999
0 <= R("R_GALT1") <= 999999
0 <= R("R_GLCTR1") <= 999999
0 <= R("R_UDPGD") <= 999999
0 <= R("R_UDPGDC") <= 999999
0 <= R("R_UGLCURPpp") <= 999999
0 <= R("R_UDPGPpp") <= 999999
0 <= R("R_UPLA4FNT") <= 999999
0 <= R("R_ULA4NFT") <= 999999
0 <= R("R_PAPPT3") <= 999999
0 <= R("R_UGMDDS") <= 999999
0 <= R("R_UPLA4FNF") <= 999999
0 <= R("R_ULA4Ntppi") <= 999999
0 <= R("R_LA4NTpp") <= 999999
0 <= R("R_UM4PCP") <= 999999
0 <= R("R_NTD2") <= 999999
0 <= R("R_NTD2pp") <= 999999
0 <= R("R_ACMAMUT") <= 999999
0 <= R("R_AADDGT") <= 999999
0 <= R("R_ECAtpp") <= 999999
0 <= R("R_ECAP1pp") <= 999999
0 <= R("R_ECAP2pp") <= 999999
0 <= R("R_ECAP3pp") <= 999999
0 <= R("R_UPP3MT") <= 999999
0 <= R("R_UPPDC1") <= 999999
0 <= R("R_UPPRT") <= 999999
0 <= R("R_URAt2pp") <= 999999
0 <= R("R_URIC") <= 999999
0 <= R("R_UGLYCH") <= 999999
0 <= R("R_URDGLYCD") <= 999999
0 <= R("R_URIK2") <= 999999
0 <= R("R_URIH") <= 999999
0 <= R("R_URIt2pp") <= 999999
0 <= R("R_RNTR4c") <= 999999
0 <= R("R_NTPP8") <= 999999
0 <= R("R_CTPS2") <= 999999
0 <= R("R_UAGDP") <= 999999
0 <= R("R_GALUi") <= 999999
0 <= R("R_VALTRS") <= 999999
0 <= R("R_VALabcpp") <= 999999
0 <= R("R_XAND") <= 999999
0 <= R("R_XPPT") <= 999999
0 <= R("R_XANt2pp") <= 999999
0 <= R("R_NTD10") <= 999999
0 <= R("R_GMPS2") <= 999999
0 <= R("R_NTD10pp") <= 999999
0 <= R("R_NTP12") <= 999999
0 <= R("R_NTPP11") <= 999999
0 <= R("R_XTSNH") <= 999999
0 <= R("R_X5PL3E") <= 999999
0 <= R("R_XYLabcpp") <= 999999
0 <= R("R_XYLt2pp") <= 999999
0 <= R("R_XYLK") <= 999999
0 <= R("R_XYLK2") <= 999999
0 <= R("R_XYLUt2pp") <= 999999
0 <= R("R_Ec_biomass_iAF1260_core_59p81M") <= 999999
0 <= R("R_ZN2abcpp") <= 999999
0 <= R("R_ZN2t3pp") <= 999999
0 <= R("R_ZN2tpp") <= 999999
0 <= R("R_ZNabcpp") <= 999999


\ Some model modifications follow

-999999 <= R("M12ppdRb_Transp") <= 0
-999999 <= R("M12ppdSb_Transp") <= 0
-999999 <= R("M14glucanb_Transp") <= 0
-999999 <= R("M15dapb_Transp") <= 0
-999999 <= R("M23campb_Transp") <= 0
-999999 <= R("M23ccmpb_Transp") <= 0
-999999 <= R("M23cgmpb_Transp") <= 0
-999999 <= R("M23cumpb_Transp") <= 0
-999999 <= R("M23dappab_Transp") <= 0
-999999 <= R("M26dapMb_Transp") <= 0
-999999 <= R("M2ddglcnb_Transp") <= 0
-999999 <= R("M34dhpacb_Transp") <= 0
-999999 <= R("M3ampb_Transp") <= 0
-999999 <= R("M3cmpb_Transp") <= 0
-999999 <= R("M3gmpb_Transp") <= 0
-999999 <= R("M3hcinnmb_Transp") <= 0
-999999 <= R("M3hpppnb_Transp") <= 0
-999999 <= R("M3umpb_Transp") <= 0
-999999 <= R("M4abutb_Transp") <= 0
-999999 <= R("M4hbab_Transp") <= 0
-999999 <= R("M4hoxpacdb_Transp") <= 0
-999999 <= R("M5dglcnb_Transp") <= 0
-999999 <= R("M5dribb_Transp") <= 0
-999999 <= R("Maacaldb_Transp") <= 0
-999999 <= R("Macacb_Transp") <= 0
-999999 <= R("Macaldb_Transp") <= 0
-999999 <= R("Macb_Transp") <= 0
-999999 <= R("Macgal1pb_Transp") <= 0
-999999 <= R("Macgalb_Transp") <= 0
-999999 <= R("Macgam1pb_Transp") <= 0
-999999 <= R("Macgamb_Transp") <= 0
-999999 <= R("Macmanab_Transp") <= 0
-999999 <= R("Macmumb_Transp") <= 0
-999999 <= R("Macnamb_Transp") <= 0
-999999 <= R("Macolipab_Transp") <= 0
-999999 <= R("Macserb_Transp") <= 0
-999999 <= R("Madeb_Transp") <= 0
-999999 <= R("Madnb_Transp") <= 0
-999999 <= R("Madocblb_Transp") <= 0
-999999 <= R("Magb_Transp") <= 0
-999999 <= R("Magmb_Transp") <= 0
-999999 <= R("Makgb_Transp") <= 0
-999999 <= R("Malaalab_Transp") <= 0
-999999 <= R("MalaBb_Transp") <= 0
-999999 <= R("MalaDb_Transp") <= 0
-999999 <= R("MalaLb_Transp") <= 0
-999999 <= R("MallDb_Transp") <= 0
-999999 <= R("Malltnb_Transp") <= 0
-999999 <= R("Mampb_Transp") <= 0
-999999 <= R("Manhgmb_Transp") <= 0
-999999 <= R("MarabLb_Transp") <= 0
-999999 <= R("Marbtnb_Transp") <= 0
-999999 <= R("Marbtnfe3b_Transp") <= 0
-999999 <= R("MargLb_Transp") <= 0
-999999 <= R("MascbLb_Transp") <= 0
-999999 <= R("MasnLb_Transp") <= 0
-999999 <= R("Maso3b_Transp") <= 0
-999999 <= R("MaspLb_Transp") <= 0
-999999 <= R("Mbutb_Transp") <= 0
-999999 <= R("Mbutso3b_Transp") <= 0
-999999 <= R("Mca2b_Transp") <= 0
-999999 <= R("Mcbib_Transp") <= 0
-999999 <= R("Mcbl1b_Transp") <= 0
-999999 <= R("Mcd2b_Transp") <= 0
-999999 <= R("Mcglyb_Transp") <= 0
-999999 <= R("Mcholb_Transp") <= 0
-999999 <= R("Mcitb_Transp") <= 0
-999999 <= R("Mclb_Transp") <= 0
-999999 <= R("Mcmpb_Transp") <= 0
-999999 <= R("Mco2b_Transp") <= 0
-999999 <= R("Mcobalt2b_Transp") <= 0
-999999 <= R("Mcolipab_Transp") <= 0
-999999 <= R("Mcpgnb_Transp") <= 0
-999999 <= R("Mcpgnunb_Transp") <= 0
-999999 <= R("Mcrnb_Transp") <= 0
-999999 <= R("Mcsnb_Transp") <= 0
-999999 <= R("Mcu2b_Transp") <= 0
-999999 <= R("Mcub_Transp") <= 0
-999999 <= R("Mcyanb_Transp") <= 0
-999999 <= R("Mcyntb_Transp") <= 0
-999999 <= R("McysDb_Transp") <= 0
-999999 <= R("McysLb_Transp") <= 0
-999999 <= R("Mcytdb_Transp") <= 0
-999999 <= R("Mdad2b_Transp") <= 0
-999999 <= R("Mdampb_Transp") <= 0
-999999 <= R("Mdcab_Transp") <= 0
-999999 <= R("Mdcmpb_Transp") <= 0
-999999 <= R("Mdcytb_Transp") <= 0
-999999 <= R("Mddcab_Transp") <= 0
-999999 <= R("Mdgmpb_Transp") <= 0
-999999 <= R("Mdgsnb_Transp") <= 0
-999999 <= R("Mdhab_Transp") <= 0
-999999 <= R("Mdimpb_Transp") <= 0
-999999 <= R("Mdinb_Transp") <= 0
-999999 <= R("Mdmsb_Transp") <= 0
-999999 <= R("Mdmsob_Transp") <= 0
-999999 <= R("Mdopab_Transp") <= 0
-999999 <= R("Mdtmpb_Transp") <= 0
-999999 <= R("Mdumpb_Transp") <= 0
-999999 <= R("Mdurib_Transp") <= 0
-999999 <= R("Meca4colipab_Transp") <= 0
-999999 <= R("Menlipab_Transp") <= 0
-999999 <= R("Menterb_Transp") <= 0
-999999 <= R("Methab_Transp") <= 0
-999999 <= R("Methso3b_Transp") <= 0
-999999 <= R("Metohb_Transp") <= 0
-999999 <= R("Mf6pb_Transp") <= 0
-999999 <= R("Mfaldb_Transp") <= 0
-999999 <= R("Mfe2b_Transp") <= 0
-999999 <= R("Mfe3b_Transp") <= 0
-999999 <= R("Mfe3dcitb_Transp") <= 0
-999999 <= R("Mfe3dhbzsb_Transp") <= 0
-999999 <= R("Mfe3hoxb_Transp") <= 0
-999999 <= R("Mfe3hoxunb_Transp") <= 0
-999999 <= R("Mfecrmb_Transp") <= 0
-999999 <= R("Mfecrmunb_Transp") <= 0
-999999 <= R("Mfeenterb_Transp") <= 0
-999999 <= R("Mfeoxamb_Transp") <= 0
-999999 <= R("Mfeoxamunb_Transp") <= 0
-999999 <= R("Mforb_Transp") <= 0
-999999 <= R("Mfrub_Transp") <= 0
-999999 <= R("Mfrulysb_Transp") <= 0
-999999 <= R("Mfruurb_Transp") <= 0
-999999 <= R("MfucLb_Transp") <= 0
-999999 <= R("Mfumb_Transp") <= 0
-999999 <= R("Mg1pb_Transp") <= 0
-999999 <= R("Mg3pcb_Transp") <= 0
-999999 <= R("Mg3peb_Transp") <= 0
-999999 <= R("Mg3pgb_Transp") <= 0
-999999 <= R("Mg3pib_Transp") <= 0
-999999 <= R("Mg3psb_Transp") <= 0
-999999 <= R("Mg6pb_Transp") <= 0
-999999 <= R("Mgal1pb_Transp") <= 0
-999999 <= R("Mgalb_Transp") <= 0
-999999 <= R("MgalbDb_Transp") <= 0
-999999 <= R("MgalctDb_Transp") <= 0
-999999 <= R("MgalctnDb_Transp") <= 0
-999999 <= R("MgalctnLb_Transp") <= 0
-999999 <= R("Mgaltb_Transp") <= 0
-999999 <= R("Mgalurb_Transp") <= 0
-999999 <= R("Mgam6pb_Transp") <= 0
-999999 <= R("Mgamb_Transp") <= 0
-999999 <= R("Mgbbtnb_Transp") <= 0
-999999 <= R("Mgdpb_Transp") <= 0
-999999 <= R("MglcDb_Transp") <= 0
-999999 <= R("Mglcnb_Transp") <= 0
-999999 <= R("Mglcrb_Transp") <= 0
-999999 <= R("Mglcur1pb_Transp") <= 0
-999999 <= R("Mglcurb_Transp") <= 0
-999999 <= R("MglnLb_Transp") <= 0
-999999 <= R("MgluLb_Transp") <= 0
-999999 <= R("Mglyaldb_Transp") <= 0
-999999 <= R("Mglyb_Transp") <= 0
-999999 <= R("Mglybb_Transp") <= 0
-999999 <= R("Mglyc2pb_Transp") <= 0
-999999 <= R("Mglyc3pb_Transp") <= 0
-999999 <= R("Mglycb_Transp") <= 0
-999999 <= R("Mglycltb_Transp") <= 0
-999999 <= R("MglycRb_Transp") <= 0
-999999 <= R("Mgmpb_Transp") <= 0
-999999 <= R("Mgsnb_Transp") <= 0
-999999 <= R("Mgthoxb_Transp") <= 0
-999999 <= R("Mgthrdb_Transp") <= 0
-999999 <= R("Mgtpb_Transp") <= 0
-999999 <= R("Mguab_Transp") <= 0
-999999 <= R("Mh2b_Transp") <= 0
-999999 <= R("Mh2o2b_Transp") <= 0
-999999 <= R("Mh2ob_Transp") <= 0
-999999 <= R("Mh2sb_Transp") <= 0
-999999 <= R("Mhacolipab_Transp") <= 0
-999999 <= R("Mhalipab_Transp") <= 0
-999999 <= R("Mhb_Transp") <= 0
-999999 <= R("Mhdcab_Transp") <= 0
-999999 <= R("Mhdceab_Transp") <= 0
-999999 <= R("Mhg2b_Transp") <= 0
-999999 <= R("MhisLb_Transp") <= 0
-999999 <= R("Mhmfurnb_Transp") <= 0
-999999 <= R("MhomLb_Transp") <= 0
-999999 <= R("Mhxab_Transp") <= 0
-999999 <= R("Mhxanb_Transp") <= 0
-999999 <= R("MidonLb_Transp") <= 0
-999999 <= R("MileLb_Transp") <= 0
-999999 <= R("Mimpb_Transp") <= 0
-999999 <= R("Mindoleb_Transp") <= 0
-999999 <= R("Minostb_Transp") <= 0
-999999 <= R("Minsb_Transp") <= 0
-999999 <= R("Misetacb_Transp") <= 0
-999999 <= R("Mkb_Transp") <= 0
-999999 <= R("Mkdo2lipid4b_Transp") <= 0
-999999 <= R("MlacDb_Transp") <= 0
-999999 <= R("MlacLb_Transp") <= 0
-999999 <= R("MLalaDgluMdapb_Transp") <= 0
-999999 <= R("MLalaDgluMdapDalab_Transp") <= 0
-999999 <= R("Mlctsb_Transp") <= 0
-999999 <= R("MleuLb_Transp") <= 0
-999999 <= R("Mlipab_Transp") <= 0
-999999 <= R("Mlipacoldb_Transp") <= 0
-999999 <= R("MlysLb_Transp") <= 0
-999999 <= R("MlyxLb_Transp") <= 0
-999999 <= R("MmalDb_Transp") <= 0
-999999 <= R("MmalLb_Transp") <= 0
-999999 <= R("Mmaltb_Transp") <= 0
-999999 <= R("Mmalthxb_Transp") <= 0
-999999 <= R("Mmaltptb_Transp") <= 0
-999999 <= R("Mmalttrb_Transp") <= 0
-999999 <= R("Mmaltttrb_Transp") <= 0
-999999 <= R("Mman6pb_Transp") <= 0
-999999 <= R("Mmanb_Transp") <= 0
-999999 <= R("Mmanglycb_Transp") <= 0
-999999 <= R("Mmelibb_Transp") <= 0
-999999 <= R("MmetDb_Transp") <= 0
-999999 <= R("MmetLb_Transp") <= 0
-999999 <= R("MmetsoxRLb_Transp") <= 0
-999999 <= R("MmetsoxSLb_Transp") <= 0
-999999 <= R("Mmg2b_Transp") <= 0
-999999 <= R("Mminohpb_Transp") <= 0
-999999 <= R("Mmmetb_Transp") <= 0
-999999 <= R("Mmn2b_Transp") <= 0
-999999 <= R("Mmnlb_Transp") <= 0
-999999 <= R("Mmobdb_Transp") <= 0
-999999 <= R("Mmso3b_Transp") <= 0
-999999 <= R("Mn2ob_Transp") <= 0
-999999 <= R("Mna1b_Transp") <= 0
-999999 <= R("Mnacb_Transp") <= 0
-999999 <= R("Mnh4b_Transp") <= 0
-999999 <= R("Mni2b_Transp") <= 0
-999999 <= R("Mnmnb_Transp") <= 0
-999999 <= R("Mno2b_Transp") <= 0
-999999 <= R("Mno3b_Transp") <= 0
-999999 <= R("Mnob_Transp") <= 0
-999999 <= R("Mo16a4colipab_Transp") <= 0
-999999 <= R("Mo2b_Transp") <= 0
-999999 <= R("Mo2sb_Transp") <= 0
-999999 <= R("Mocdcab_Transp") <= 0
-999999 <= R("Mocdceab_Transp") <= 0
-999999 <= R("Moctab_Transp") <= 0
-999999 <= R("Mornb_Transp") <= 0
-999999 <= R("Morotb_Transp") <= 0
-999999 <= R("Moxamb_Transp") <= 0
-999999 <= R("Mpacaldb_Transp") <= 0
-999999 <= R("Mpeamnb_Transp") <= 0
-999999 <= R("MpheLb_Transp") <= 0
-999999 <= R("Mphemeb_Transp") <= 0
-999999 <= R("Mpib_Transp") <= 0
-999999 <= R("MpntoRb_Transp") <= 0
-999999 <= R("Mppab_Transp") <= 0
-999999 <= R("Mppalb_Transp") <= 0
-999999 <= R("Mpppnb_Transp") <= 0
-999999 <= R("Mpptb_Transp") <= 0
-999999 <= R("Mproglyb_Transp") <= 0
-999999 <= R("MproLb_Transp") <= 0
-999999 <= R("Mpsclysb_Transp") <= 0
-999999 <= R("MpserLb_Transp") <= 0
-999999 <= R("Mptrcb_Transp") <= 0
-999999 <= R("Mpyrb_Transp") <= 0
-999999 <= R("Mr5pb_Transp") <= 0
-999999 <= R("MribDb_Transp") <= 0
-999999 <= R("Mrmnb_Transp") <= 0
-999999 <= R("MsbtDb_Transp") <= 0
-999999 <= R("MserDb_Transp") <= 0
-999999 <= R("MserLb_Transp") <= 0
-999999 <= R("Mskmb_Transp") <= 0
-999999 <= R("Mso2b_Transp") <= 0
-999999 <= R("Mso3b_Transp") <= 0
-999999 <= R("Mso4b_Transp") <= 0
-999999 <= R("Mspmdb_Transp") <= 0
-999999 <= R("Msuccb_Transp") <= 0
-999999 <= R("Msucrb_Transp") <= 0
-999999 <= R("Msulfacb_Transp") <= 0
-999999 <= R("MtartrLb_Transp") <= 0
-999999 <= R("Mtaurb_Transp") <= 0
-999999 <= R("Mtcyntb_Transp") <= 0
-999999 <= R("Mthmb_Transp") <= 0
-999999 <= R("MthrLb_Transp") <= 0
-999999 <= R("Mthrpb_Transp") <= 0
-999999 <= R("Mthymb_Transp") <= 0
-999999 <= R("Mthymdb_Transp") <= 0
-999999 <= R("Mtmab_Transp") <= 0
-999999 <= R("Mtmaob_Transp") <= 0
-999999 <= R("Mtreb_Transp") <= 0
-999999 <= R("MtrpLb_Transp") <= 0
-999999 <= R("Mtsulb_Transp") <= 0
-999999 <= R("Mttdcab_Transp") <= 0
-999999 <= R("Mttdceab_Transp") <= 0
-999999 <= R("Mtungsb_Transp") <= 0
-999999 <= R("Mtymb_Transp") <= 0
-999999 <= R("MtyrLb_Transp") <= 0
-999999 <= R("Mtyrpb_Transp") <= 0
-999999 <= R("Muacgamb_Transp") <= 0
-999999 <= R("Mudpacgalb_Transp") <= 0
-999999 <= R("Mudpgalb_Transp") <= 0
-999999 <= R("Mudpgb_Transp") <= 0
-999999 <= R("Mudpglcurb_Transp") <= 0
-999999 <= R("Mumpb_Transp") <= 0
-999999 <= R("Murab_Transp") <= 0
-999999 <= R("Mureab_Transp") <= 0
-999999 <= R("Murib_Transp") <= 0
-999999 <= R("MvalLb_Transp") <= 0
-999999 <= R("Mxanb_Transp") <= 0
-999999 <= R("Mxmpb_Transp") <= 0
-999999 <= R("Mxtsnb_Transp") <= 0
-999999 <= R("MxylDb_Transp") <= 0
-999999 <= R("MxyluLb_Transp") <= 0
-999999 <= R("Mzn2b_Transp") <= 0


\ Some model modifications follow

-999999 <= R("Maacaldb_Transp") <= 999999
-999999 <= R("M5dribb_Transp") <= 999999
-999999 <= R("Moxamb_Transp") <= 999999
-999999 <= R("Mhmfurnb_Transp") <= 999999
-999999 <= R("M4hbab_Transp") <= 999999


\ Some model modifications follow

0 <= R("Mna1b_Transp") <= 20
0 <= R("Mkb_Transp") <= 20
0 <= R("Mca2b_Transp") <= 20
0 <= R("Mcu2b_Transp") <= 20
0 <= R("Mmg2b_Transp") <= 20
0 <= R("Mzn2b_Transp") <= 20
0 <= R("Mmobdb_Transp") <= 20
0 <= R("Mfe2b_Transp") <= 20
0 <= R("Mfe3b_Transp") <= 20
0 <= R("Mcobalt2b_Transp") <= 20
0 <= R("Mmn2b_Transp") <= 20
0 <= R("Mclb_Transp") <= 20
