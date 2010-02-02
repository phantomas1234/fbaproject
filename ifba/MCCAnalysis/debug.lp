\* Problem: Unknown *\

Maximize
 obj: + Mvaltrnac

Subject To
 R("R_AICART"): + Mthfc + Mfpricac - Maicarc - M10fthfc = 0
 R("R_AICART_Rev"): - Mthfc - Mfpricac + Maicarc + M10fthfc = 0
 R("R_FMETTRS"): + Mthfc - Mmettrnac + Mhc + Mfmettrnac - M10fthfc = 0
 R("R_FTHFD"): + Mthfc + Mhc - Mh2oc + Mforc - M10fthfc = 0
 R("R_GARFT"): + Mthfc + Mhc - Mgarc + Mfgamc - M10fthfc = 0
 R("R_GARFT_Rev"): - Mthfc - Mhc + Mgarc - Mfgamc + M10fthfc = 0
 R("R_MTHFC"): - Mmethfc + Mhc - Mh2oc + M10fthfc = 0
 R("R_MTHFC_Rev"): + Mmethfc - Mhc + Mh2oc - M10fthfc = 0
 R("R_ULA4NFT"): - MudpLa4nc + MudpLa4fnc + Mthfc + Mhc - M10fthfc = 0
 R("R_12DGR120tipp"): - M12dgr120p + M12dgr120c = 0
 R("R_DAGK120"): + Mpa120c + Mhc - Matpc + Madpc - M12dgr120c = 0
 R("R_PAPA120"): + Mpic - Mpa120c - Mh2oc + M12dgr120c = 0
 R("R_PAPA120pp"): + Mpip - Mpa120p - Mh2op + M12dgr120p = 0
 R("R_12DGR140tipp"): - M12dgr140p + M12dgr140c = 0
 R("R_DAGK140"): + Mpa140c + Mhc - Matpc + Madpc - M12dgr140c = 0
 R("R_PAPA140"): + Mpic - Mpa140c - Mh2oc + M12dgr140c = 0
 R("R_PAPA140pp"): + Mpip - Mpa140p - Mh2op + M12dgr140p = 0
 R("R_12DGR141tipp"): - M12dgr141p + M12dgr141c = 0
 R("R_DAGK141"): + Mpa141c + Mhc - Matpc + Madpc - M12dgr141c = 0
 R("R_PAPA141"): + Mpic - Mpa141c - Mh2oc + M12dgr141c = 0
 R("R_PAPA141pp"): + Mpip - Mpa141p - Mh2op + M12dgr141p = 0
 R("R_12DGR160tipp"): - M12dgr160p + M12dgr160c = 0
 R("R_DAGK160"): + Mpa160c + Mhc - Matpc + Madpc - M12dgr160c = 0
 R("R_PAPA160"): + Mpic - Mpa160c - Mh2oc + M12dgr160c = 0
 R("R_PAPA160pp"): + Mpip - Mpa160p - Mh2op + M12dgr160p = 0
 R("R_12DGR161tipp"): - M12dgr161p + M12dgr161c = 0
 R("R_DAGK161"): + Mpa161c + Mhc - Matpc + Madpc - M12dgr161c = 0
 R("R_PAPA161"): + Mpic - Mpa161c - Mh2oc + M12dgr161c = 0
 R("R_PAPA161pp"): + Mpip - Mpa161p - Mh2op + M12dgr161p = 0
 R("R_PETNT161pp"): - Mpe161p - Mlipap + Menlipap + M12dgr161p = 0
 R("R_12DGR180tipp"): - M12dgr180p + M12dgr180c = 0
 R("R_DAGK180"): + Mpa180c + Mhc - Matpc + Madpc - M12dgr180c = 0
 R("R_PAPA180"): + Mpic - Mpa180c - Mh2oc + M12dgr180c = 0
 R("R_PAPA180pp"): + Mpip - Mpa180p - Mh2op + M12dgr180p = 0
 R("R_12DGR181tipp"): - M12dgr181p + M12dgr181c = 0
 R("R_DAGK181"): + Mpa181c + Mhc - Matpc + Madpc - M12dgr181c = 0
 R("R_PAPA181"): + Mpic - Mpa181c - Mh2oc + M12dgr181c = 0
 R("R_PAPA181pp"): + Mpip - Mpa181p - Mh2op + M12dgr181p = 0
 R("R_PETNT181pp"): - Mpe181p - Mlipap + Menlipap + M12dgr181p = 0
 R("M12ppdRb_Transp"): 0 M10fthfc = 0
 R("R_EX_12ppd_R_e_"): - M12ppdRe = 0
 R("R_EX_12ppd_R_e__Rev"): + M12ppdRe = 0
 R("R_12PPDRtpp"): - M12ppdRp + M12ppdRc = 0
 R("R_12PPDRtpp_Rev"): + M12ppdRp - M12ppdRc = 0
 R("R_ALR4x"): - Mnadhc + Mnadc - Mhc - Macetolc + M12ppdRc = 0
 R("R_LCARR"): - Mnadhc + Mnadc - MlaldDc - Mhc + M12ppdRc = 0
 R("R_LCARR_Rev"): + Mnadhc - Mnadc + MlaldDc + Mhc - M12ppdRc = 0
 R("R_12PPDRtex"): + M12ppdRp - M12ppdRe = 0
 R("R_12PPDRtex_Rev"): - M12ppdRp + M12ppdRe = 0
 R("M12ppdSb_Transp"): 0 M10fthfc = 0
 R("R_EX_12ppd_S_e_"): - M12ppdSe = 0
 R("R_EX_12ppd_S_e__Rev"): + M12ppdSe = 0
 R("R_12PPDStpp"): - M12ppdSp + M12ppdSc = 0
 R("R_12PPDStpp_Rev"): + M12ppdSp - M12ppdSc = 0
 R("R_LCARS"): - Mnadhc + Mnadc - MlaldLc - Mhc + M12ppdSc = 0
 R("R_LCARS_Rev"): + Mnadhc - Mnadc + MlaldLc + Mhc - M12ppdSc = 0
 R("R_12PPDStex"): + M12ppdSp - M12ppdSe = 0
 R("R_12PPDStex_Rev"): - M12ppdSp + M12ppdSe = 0
 R("R_GAPD"): - Mpic + Mnadhc - Mnadc + Mhc - Mg3pc + M13dpgc = 0
 R("R_GAPD_Rev"): + Mpic - Mnadhc + Mnadc - Mhc + Mg3pc - M13dpgc = 0
 R("R_PGK"): - Matpc + Madpc - M3pgc + M13dpgc = 0
 R("R_PGK_Rev"): + Matpc - Madpc + M3pgc - M13dpgc = 0
 R("M14glucanb_Transp"): 0 M10fthfc = 0
 R("R_EX_14glucan_e_"): - M14glucane = 0
 R("R_EX_14glucan_e__Rev"): + M14glucane = 0
 R("R_14GLUCANabcpp"): + Mpic + Mhc - Mh2oc - Matpc + Madpc - M14glucanp
 + M14glucanc = 0
 R("R_AAMYL"): + Mmalthxc - M14glucanc = 0
 R("R_14GLUCANtexi"): + M14glucanp - M14glucane = 0
 R("R_AAMYLpp"): + Mmalthxp - M14glucanp = 0
 R("M15dapb_Transp"): 0 M10fthfc = 0
 R("R_EX_15dap_e_"): - M15dape = 0
 R("R_EX_15dap_e__Rev"): + M15dape = 0
 R("R_CADVtpp"): - MlysLp + MlysLc - Mhp + Mhc + M15dapp - M15dapc = 0
 R("R_LYSDC"): - MlysLc - Mhc + Mco2c + M15dapc = 0
 R("R_DAPtex"): + M15dapp - M15dape = 0
 R("R_DAPtex_Rev"): - M15dapp + M15dape = 0
 R("R_LPLIPAL1E120pp"): + Mhp - Mh2op + Mg3pep + Mddcap - M1agpe120p = 0
 R("R_PLIPA2E120pp"): - Mpe120p + Mhp - Mh2op + Mddcap + M1agpe120p = 0
 R("R_LPLIPAL1E140pp"): + Mttdcap + Mhp - Mh2op + Mg3pep - M1agpe140p
 = 0
 R("R_PLIPA2E140pp"): + Mttdcap - Mpe140p + Mhp - Mh2op + M1agpe140p = 0
 R("R_LPLIPAL1E141pp"): + Mttdceap + Mhp - Mh2op + Mg3pep - M1agpe141p
 = 0
 R("R_PLIPA2E141pp"): + Mttdceap - Mpe141p + Mhp - Mh2op + M1agpe141p
 = 0
 R("R_LPLIPAL1E160pp"): + Mhp + Mhdcap - Mh2op + Mg3pep - M1agpe160p = 0
 R("R_PLIPA2E160pp"): - Mpe160p + Mhp + Mhdcap - Mh2op + M1agpe160p = 0
 R("R_LPLIPAL1E161pp"): + Mhp + Mhdceap - Mh2op + Mg3pep - M1agpe161p
 = 0
 R("R_PLIPA2E161pp"): - Mpe161p + Mhp + Mhdceap - Mh2op + M1agpe161p = 0
 R("R_LPLIPAL1E180pp"): + Mocdcap + Mhp - Mh2op + Mg3pep - M1agpe180p
 = 0
 R("R_PLIPA2E180pp"): - Mpe180p + Mocdcap + Mhp - Mh2op + M1agpe180p = 0
 R("R_LPLIPAL1E181pp"): + Mocdceap + Mhp - Mh2op + Mg3pep - M1agpe181p
 = 0
 R("R_PLIPA2E181pp"): - Mpe181p + Mocdceap + Mhp - Mh2op + M1agpe181p
 = 0
 R("R_LPLIPAL1G120pp"): + Mhp - Mh2op + Mg3pgp + Mddcap - M1agpg120p = 0
 R("R_PLIPA2G120pp"): - Mpg120p + Mhp - Mh2op + Mddcap + M1agpg120p = 0
 R("R_LPLIPAL1G140pp"): + Mttdcap + Mhp - Mh2op + Mg3pgp - M1agpg140p
 = 0
 R("R_PLIPA2G140pp"): + Mttdcap - Mpg140p + Mhp - Mh2op + M1agpg140p = 0
 R("R_LPLIPAL1G141pp"): + Mttdceap + Mhp - Mh2op + Mg3pgp - M1agpg141p
 = 0
 R("R_PLIPA2G141pp"): + Mttdceap - Mpg141p + Mhp - Mh2op + M1agpg141p
 = 0
 R("R_LPLIPAL1G160pp"): + Mhp + Mhdcap - Mh2op + Mg3pgp - M1agpg160p = 0
 R("R_PLIPA2G160pp"): - Mpg160p + Mhp + Mhdcap - Mh2op + M1agpg160p = 0
 R("R_LPLIPAL1G161pp"): + Mhp + Mhdceap - Mh2op + Mg3pgp - M1agpg161p
 = 0
 R("R_PLIPA2G161pp"): - Mpg161p + Mhp + Mhdceap - Mh2op + M1agpg161p = 0
 R("R_LPLIPAL1G180pp"): + Mocdcap + Mhp - Mh2op + Mg3pgp - M1agpg180p
 = 0
 R("R_PLIPA2G180pp"): - Mpg180p + Mocdcap + Mhp - Mh2op + M1agpg180p = 0
 R("R_LPLIPAL1G181pp"): + Mocdceap + Mhp - Mh2op + Mg3pgp - M1agpg181p
 = 0
 R("R_PLIPA2G181pp"): - Mpg181p + Mocdceap + Mhp - Mh2op + M1agpg181p
 = 0
 R("R_AGPAT120"): + Mpa120c - MddcaACPc + MACPc - M1ddecg3pc = 0
 R("R_G3PAT120"): - Mglyc3pc - MddcaACPc + MACPc + M1ddecg3pc = 0
 R("R_LPLIPAL1A120pp"): + Mhp - Mh2op + Mglyc3pp + Mddcap - M1ddecg3pp
 = 0
 R("R_PLIPA2A120pp"): - Mpa120p + Mhp - Mh2op + Mddcap + M1ddecg3pp = 0
 R("R_AGPAT161"): + Mpa161c - MhdeACPc + MACPc - M1hdec9eg3pc = 0
 R("R_G3PAT161"): - MhdeACPc - Mglyc3pc + MACPc + M1hdec9eg3pc = 0
 R("R_LPLIPAL1A161pp"): + Mhp + Mhdceap - Mh2op + Mglyc3pp
 - M1hdec9eg3pp = 0
 R("R_PLIPA2A161pp"): - Mpa161p + Mhp + Mhdceap - Mh2op + M1hdec9eg3pp
 = 0
 R("R_AGPAT160"): - MpalmACPc + Mpa160c + MACPc - M1hdecg3pc = 0
 R("R_G3PAT160"): - MpalmACPc - Mglyc3pc + MACPc + M1hdecg3pc = 0
 R("R_LPLIPAL1A160pp"): + Mhp + Mhdcap - Mh2op + Mglyc3pp - M1hdecg3pp
 = 0
 R("R_PLIPA2A160pp"): - Mpa160p + Mhp + Mhdcap - Mh2op + M1hdecg3pp = 0
 R("R_AGPAT181"): + Mpa181c - MocteACPc + MACPc - M1odec11eg3pc = 0
 R("R_G3PAT181"): - MocteACPc - Mglyc3pc + MACPc + M1odec11eg3pc = 0
 R("R_LPLIPAL1A181pp"): + Mocdceap + Mhp - Mh2op + Mglyc3pp
 - M1odec11eg3pp = 0
 R("R_PLIPA2A181pp"): - Mpa181p + Mocdceap + Mhp - Mh2op + M1odec11eg3pp
 = 0
 R("R_AGPAT180"): + Mpa180c - MocdcaACPc + MACPc - M1odecg3pc = 0
 R("R_G3PAT180"): - MocdcaACPc - Mglyc3pc + MACPc + M1odecg3pc = 0
 R("R_LPLIPAL1A180pp"): + Mocdcap + Mhp - Mh2op + Mglyc3pp - M1odecg3pp
 = 0
 R("R_PLIPA2A180pp"): - Mpa180p + Mocdcap + Mhp - Mh2op + M1odecg3pp = 0
 R("R_G5SADs"): + Mhc + Mh2oc - Mglu5sac + M1pyr5cc = 0
 R("R_P5CD"): + Mnadhc - Mnadc + Mhc - 2 Mh2oc + MgluLc - M1pyr5cc = 0
 R("R_P5CR"): + MproLc - Mnadphc + Mnadpc - 2 Mhc - M1pyr5cc = 0
 R("R_PROD2"): - MproLc + Mhc + Mfadh2c - Mfadc + M1pyr5cc = 0
 R("R_AGPAT141"): - MtdeACPc + Mpa141c + MACPc - M1tdec7eg3pc = 0
 R("R_G3PAT141"): - MtdeACPc - Mglyc3pc + MACPc + M1tdec7eg3pc = 0
 R("R_LPLIPAL1A141pp"): + Mttdceap + Mhp - Mh2op + Mglyc3pp
 - M1tdec7eg3pp = 0
 R("R_PLIPA2A141pp"): + Mttdceap - Mpa141p + Mhp - Mh2op + M1tdec7eg3pp
 = 0
 R("R_AGPAT140"): + Mpa140c - MmyrsACPc + MACPc - M1tdecg3pc = 0
 R("R_G3PAT140"): - MmyrsACPc - Mglyc3pc + MACPc + M1tdecg3pc = 0
 R("R_LPLIPAL1A140pp"): + Mttdcap + Mhp - Mh2op + Mglyc3pp - M1tdecg3pp
 = 0
 R("R_PLIPA2A140pp"): + Mttdcap - Mpa140p + Mhp - Mh2op + M1tdecg3pp = 0
 R("M23campb_Transp"): 0 M10fthfc = 0
 R("R_EX_23camp_e_"): - M23campe = 0
 R("R_EX_23camp_e__Rev"): + M23campe = 0
 R("R_23CAMPtex"): + M23campp - M23campe = 0
 R("R_23CAMPtex_Rev"): - M23campp + M23campe = 0
 R("R_23PDE7pp"): + Mhp - Mh2op + M3ampp - M23campp = 0
 R("M23ccmpb_Transp"): 0 M10fthfc = 0
 R("R_EX_23ccmp_e_"): - M23ccmpe = 0
 R("R_EX_23ccmp_e__Rev"): + M23ccmpe = 0
 R("R_23CCMPtex"): + M23ccmpp - M23ccmpe = 0
 R("R_23CCMPtex_Rev"): - M23ccmpp + M23ccmpe = 0
 R("R_23PDE4pp"): + Mhp - Mh2op + M3cmpp - M23ccmpp = 0
 R("M23cgmpb_Transp"): 0 M10fthfc = 0
 R("R_EX_23cgmp_e_"): - M23cgmpe = 0
 R("R_EX_23cgmp_e__Rev"): + M23cgmpe = 0
 R("R_23CGMPtex"): + M23cgmpp - M23cgmpe = 0
 R("R_23CGMPtex_Rev"): - M23cgmpp + M23cgmpe = 0
 R("R_23PDE9pp"): + Mhp - Mh2op + M3gmpp - M23cgmpp = 0
 R("M23cumpb_Transp"): 0 M10fthfc = 0
 R("R_EX_23cump_e_"): - M23cumpe = 0
 R("R_EX_23cump_e__Rev"): + M23cumpe = 0
 R("R_23CUMPtex"): + M23cumpp - M23cumpe = 0
 R("R_23CUMPtex_Rev"): - M23cumpp + M23cumpe = 0
 R("R_23PDE2pp"): + Mhp - Mh2op + M3umpp - M23cumpp = 0
 R("M23dappab_Transp"): 0 M10fthfc = 0
 R("R_EX_23dappa_e_"): - M23dappae = 0
 R("R_EX_23dappa_e__Rev"): + M23dappae = 0
 R("R_23DAPPAt2pp"): - Mhp + Mhc - M23dappap + M23dappac = 0
 R("R_DAPAL"): + Mpyrc + 2 Mnh4c - Mh2oc - M23dappac = 0
 R("R_23DAPPAtex"): + M23dappap - M23dappae = 0
 R("R_23DAPPAtex_Rev"): - M23dappap + M23dappae = 0
 R("R_DHBD"): + Mnadhc - Mnadc + Mhc + M23dhbc - M23ddhbc = 0
 R("R_DHBD_Rev"): - Mnadhc + Mnadc - Mhc - M23dhbc + M23ddhbc = 0
 R("R_ICHORT"): + Mpyrc - Michorc - Mh2oc + M23ddhbc = 0
 R("R_DHBS"): + Mppic - Mhc - Matpc - M23dhbc + M23dhbac = 0
 R("R_ENTCS"): - 3 Mserampc + 9 Mhc + Menterc + 6 Mampc - 3 M23dhbac = 0
 R("R_DHBSH"): + MserLc - Mh2oc - M23dhbzsc + M23dhbc = 0
 R("R_ENTERES"): + 3 Mhc - 3 Mh2oc - Menterc + 3 M23dhbzsc = 0
 R("R_ENTERES2"): + 3 Mhc - 3 Mh2oc - Mfeenterc + Mfe3c + 3 M23dhbzsc
 = 0
 R("R_FE3DHBZR"): - Mfe3dhbzsc + Mfe3c + M23dhbzsc = 0
 R("R_DHDPRy"): + Mthdpc - Mnadphc + Mnadpc - Mhc - M23dhdpc = 0
 R("R_DHDPS"): - Mpyrc + Mhc + 2 Mh2oc - Maspsac + M23dhdpc = 0
 R("R_DHAD1"): + Mh2oc + M3mobc - M23dhmbc = 0
 R("R_KARA1"): + Mnadphc - Mnadpc + Mhc + MalacSc - M23dhmbc = 0
 R("R_KARA1_Rev"): - Mnadphc + Mnadpc - Mhc - MalacSc + M23dhmbc = 0
 R("R_DHAD2"): + Mh2oc + M3mopc - M23dhmpc = 0
 R("R_KARA2"): - Mnadphc + Mnadpc - Mhc - M2ahbutc + M23dhmpc = 0
 R("R_KARA2_Rev"): + Mnadphc - Mnadpc + Mhc + M2ahbutc - M23dhmpc = 0
 R("R_DOGULNR"): - Mnadhc + Mnadc - Mhc + M3dhgulnc - M23dogulnc = 0
 R("R_ADSL2r"): + Mfumc + Maicarc - M25aicsc = 0
 R("R_ADSL2r_Rev"): - Mfumc - Maicarc + M25aicsc = 0
 R("R_PRASCSi"): + Mpic + Mhc - Matpc - MaspLc + Madpc - M5aizcc
 + M25aicsc = 0
 R("R_DKGLCNR1"): - Mnadphc + Mnadpc - Mhc + M2dhgulnc - M25dkglcnc = 0
 R("R_DKGLCNR2x"): - Mnadhc + Mnadc - Mhc + M5dglcnc - M25dkglcnc = 0
 R("R_DKGLCNR2y"): - Mnadphc + Mnadpc - Mhc + M5dglcnc - M25dkglcnc = 0
 R("R_DHPPDA2"): + Mnh4c - Mhc - Mh2oc + M5apruc - M25drappc = 0
 R("R_GTPCII2"): + Mppic + 2 Mhc - 3 Mh2oc - Mgtpc + Mforc + M25drappc
 = 0
 R("R_DAPE"): + M26dapMc - M26dapLLc = 0
 R("R_DAPE_Rev"): - M26dapMc + M26dapLLc = 0
 R("R_SDPDS"): + Msuccc - Msl26dac - Mh2oc + M26dapLLc = 0
 R("M26dapMb_Transp"): 0 M10fthfc = 0
 R("R_EX_26dap_M_e_"): - M26dapMe = 0
 R("R_EX_26dap_M_e__Rev"): + M26dapMe = 0
 R("R_DAPabcpp"): + Mpic + Mhc - Mh2oc - Matpc + Madpc - M26dapMp
 + M26dapMc = 0
 R("R_DAPDC"): + MlysLc - Mhc + Mco2c - M26dapMc = 0
 R("R_LADGMDH"): - MLalaDgluMdapc + MLalaDgluc - Mh2oc + M26dapMc = 0
 R("R_UAAGDS"): + Mugmdc - Muamagc + Mpic + Mhc - Matpc + Madpc
 - M26dapMc = 0
 R("R_26DAHtex"): + M26dapMp - M26dapMe = 0
 R("R_26DAHtex_Rev"): - M26dapMp + M26dapMe = 0
 R("R_2AGPE120tipp"): - M2agpe120p + M2agpe120c = 0
 R("R_2AGPEAT120"): + Mppic + Mpe120c - Mddcac - Matpc + Mampc
 - M2agpe120c = 0
 R("R_LPLIPAL2ATE120"): - Mpg120c + Mg3pec + Mapg120c - M2agpe120c = 0
 R("R_LPLIPAL2E120"): + Mhc - Mh2oc + Mg3pec + Mddcac - M2agpe120c = 0
 R("R_PLIPA1E120pp"): - Mpe120p + Mhp - Mh2op + Mddcap + M2agpe120p = 0
 R("R_2AGPE140tipp"): - M2agpe140p + M2agpe140c = 0
 R("R_2AGPEAT140"): - Mttdcac + Mppic + Mpe140c - Matpc + Mampc
 - M2agpe140c = 0
 R("R_LPLIPAL2ATE140"): - Mpg140c + Mg3pec + Mapg140c - M2agpe140c = 0
 R("R_LPLIPAL2E140"): + Mttdcac + Mhc - Mh2oc + Mg3pec - M2agpe140c = 0
 R("R_PLIPA1E140pp"): + Mttdcap - Mpe140p + Mhp - Mh2op + M2agpe140p = 0
 R("R_2AGPE141tipp"): - M2agpe141p + M2agpe141c = 0
 R("R_2AGPEAT141"): - Mttdceac + Mppic + Mpe141c - Matpc + Mampc
 - M2agpe141c = 0
 R("R_LPLIPAL2ATE141"): - Mpg141c + Mg3pec + Mapg141c - M2agpe141c = 0
 R("R_LPLIPAL2E141"): + Mttdceac + Mhc - Mh2oc + Mg3pec - M2agpe141c = 0
 R("R_PLIPA1E141pp"): + Mttdceap - Mpe141p + Mhp - Mh2op + M2agpe141p
 = 0
 R("R_2AGPE160tipp"): - M2agpe160p + M2agpe160c = 0
 R("R_2AGPEAT160"): + Mppic + Mpe160c - Mhdcac - Matpc + Mampc
 - M2agpe160c = 0
 R("R_LPLIPAL2ATE160"): - Mpg160c + Mg3pec + Mapg160c - M2agpe160c = 0
 R("R_LPLIPAL2E160"): + Mhdcac + Mhc - Mh2oc + Mg3pec - M2agpe160c = 0
 R("R_ALPATE160pp"): - Mpe160p + Mlppp - Malppp + M2agpe160p = 0
 R("R_PLIPA1E160pp"): - Mpe160p + Mhp + Mhdcap - Mh2op + M2agpe160p = 0
 R("R_2AGPE161tipp"): - M2agpe161p + M2agpe161c = 0
 R("R_2AGPEAT161"): + Mppic + Mpe161c - Mhdceac - Matpc + Mampc
 - M2agpe161c = 0
 R("R_LPLIPAL2ATE161"): - Mpg161c + Mg3pec + Mapg161c - M2agpe161c = 0
 R("R_LPLIPAL2E161"): + Mhdceac + Mhc - Mh2oc + Mg3pec - M2agpe161c = 0
 R("R_PLIPA1E161pp"): - Mpe161p + Mhp + Mhdceap - Mh2op + M2agpe161p = 0
 R("R_2AGPE180tipp"): - M2agpe180p + M2agpe180c = 0
 R("R_2AGPEAT180"): + Mppic + Mpe180c - Mocdcac - Matpc + Mampc
 - M2agpe180c = 0
 R("R_LPLIPAL2ATE180"): - Mpg180c + Mg3pec + Mapg180c - M2agpe180c = 0
 R("R_LPLIPAL2E180"): + Mocdcac + Mhc - Mh2oc + Mg3pec - M2agpe180c = 0
 R("R_PLIPA1E180pp"): - Mpe180p + Mocdcap + Mhp - Mh2op + M2agpe180p = 0
 R("R_2AGPE181tipp"): - M2agpe181p + M2agpe181c = 0
 R("R_2AGPEAT181"): + Mppic + Mpe181c - Mocdceac - Matpc + Mampc
 - M2agpe181c = 0
 R("R_LPLIPAL2ATE181"): - Mpg181c + Mg3pec + Mapg181c - M2agpe181c = 0
 R("R_LPLIPAL2E181"): + Mocdceac + Mhc - Mh2oc + Mg3pec - M2agpe181c = 0
 R("R_PLIPA1E181pp"): - Mpe181p + Mocdceap + Mhp - Mh2op + M2agpe181p
 = 0
 R("R_2AGPG120tipp"): - M2agpg120p + M2agpg120c = 0
 R("R_2AGPGAT120"): + Mppic + Mpg120c - Mddcac - Matpc + Mampc
 - M2agpg120c = 0
 R("R_LPLIPAL2ATG120"): - Mpg120c + Mg3pgc + Mapg120c - M2agpg120c = 0
 R("R_LPLIPAL2G120"): + Mhc - Mh2oc + Mg3pgc + Mddcac - M2agpg120c = 0
 R("R_PLIPA1G120pp"): - Mpg120p + Mhp - Mh2op + Mddcap + M2agpg120p = 0
 R("R_2AGPG140tipp"): - M2agpg140p + M2agpg140c = 0
 R("R_2AGPGAT140"): - Mttdcac + Mppic + Mpg140c - Matpc + Mampc
 - M2agpg140c = 0
 R("R_LPLIPAL2ATG140"): - Mpg140c + Mg3pgc + Mapg140c - M2agpg140c = 0
 R("R_LPLIPAL2G140"): + Mttdcac + Mhc - Mh2oc + Mg3pgc - M2agpg140c = 0
 R("R_PLIPA1G140pp"): + Mttdcap - Mpg140p + Mhp - Mh2op + M2agpg140p = 0
 R("R_2AGPG141tipp"): - M2agpg141p + M2agpg141c = 0
 R("R_2AGPGAT141"): - Mttdceac + Mppic + Mpg141c - Matpc + Mampc
 - M2agpg141c = 0
 R("R_LPLIPAL2ATG141"): - Mpg141c + Mg3pgc + Mapg141c - M2agpg141c = 0
 R("R_LPLIPAL2G141"): + Mttdceac + Mhc - Mh2oc + Mg3pgc - M2agpg141c = 0
 R("R_PLIPA1G141pp"): + Mttdceap - Mpg141p + Mhp - Mh2op + M2agpg141p
 = 0
 R("R_2AGPG160tipp"): - M2agpg160p + M2agpg160c = 0
 R("R_2AGPGAT160"): + Mppic + Mpg160c - Mhdcac - Matpc + Mampc
 - M2agpg160c = 0
 R("R_LPLIPAL2ATG160"): - Mpg160c + Mg3pgc + Mapg160c - M2agpg160c = 0
 R("R_LPLIPAL2G160"): + Mhdcac + Mhc - Mh2oc + Mg3pgc - M2agpg160c = 0
 R("R_ALPATG160pp"): - Mpg160p + Mlppp - Malppp + M2agpg160p = 0
 R("R_PLIPA1G160pp"): - Mpg160p + Mhp + Mhdcap - Mh2op + M2agpg160p = 0
 R("R_2AGPG161tipp"): - M2agpg161p + M2agpg161c = 0
 R("R_2AGPGAT161"): + Mppic + Mpg161c - Mhdceac - Matpc + Mampc
 - M2agpg161c = 0
 R("R_LPLIPAL2ATG161"): - Mpg161c + Mg3pgc + Mapg161c - M2agpg161c = 0
 R("R_LPLIPAL2G161"): + Mhdceac + Mhc - Mh2oc + Mg3pgc - M2agpg161c = 0
 R("R_PLIPA1G161pp"): - Mpg161p + Mhp + Mhdceap - Mh2op + M2agpg161p = 0
 R("R_2AGPG180tipp"): - M2agpg180p + M2agpg180c = 0
 R("R_2AGPGAT180"): + Mppic + Mpg180c - Mocdcac - Matpc + Mampc
 - M2agpg180c = 0
 R("R_LPLIPAL2ATG180"): - Mpg180c + Mg3pgc + Mapg180c - M2agpg180c = 0
 R("R_LPLIPAL2G180"): + Mocdcac + Mhc - Mh2oc + Mg3pgc - M2agpg180c = 0
 R("R_PLIPA1G180pp"): - Mpg180p + Mocdcap + Mhp - Mh2op + M2agpg180p = 0
 R("R_2AGPG181tipp"): - M2agpg181p + M2agpg181c = 0
 R("R_2AGPGAT181"): + Mppic + Mpg181c - Mocdceac - Matpc + Mampc
 - M2agpg181c = 0
 R("R_LPLIPAL2ATG181"): - Mpg181c + Mg3pgc + Mapg181c - M2agpg181c = 0
 R("R_LPLIPAL2G181"): + Mocdceac + Mhc - Mh2oc + Mg3pgc - M2agpg181c = 0
 R("R_PLIPA1G181pp"): - Mpg181p + Mocdceap + Mhp - Mh2op + M2agpg181p
 = 0
 R("R_ACHBS"): - Mpyrc - Mhc + Mco2c - M2obutc + M2ahbutc = 0
 R("R_DSERDHr"): - MserDc + Mnadphc - Mnadpc + Mhc + M2amsac = 0
 R("R_DSERDHr_Rev"): + MserDc - Mnadphc + Mnadpc - Mhc - M2amsac = 0
 R("R_LSERDHr"): - MserLc + Mnadphc - Mnadpc + Mhc + M2amsac = 0
 R("R_LSERDHr_Rev"): + MserLc - Mnadphc + Mnadpc - Mhc - M2amsac = 0
 R("R_AOBUTDs"): - Mhc + Mco2c + Maactc - M2aobutc = 0
 R("R_ATHRDHr"): + Mnadphc - Mnadpc + Mhc - MathrLc + M2aobutc = 0
 R("R_ATHRDHr_Rev"): - Mnadphc + Mnadpc - Mhc + MathrLc - M2aobutc = 0
 R("R_GLYAT"): - Mglyc + Mcoac - Maccoac + M2aobutc = 0
 R("R_GLYAT_Rev"): + Mglyc - Mcoac + Maccoac - M2aobutc = 0
 R("R_THRD"): - MthrLc + Mnadhc - Mnadc + Mhc + M2aobutc = 0
 R("R_IGPS"): - Mhc + Mh2oc + Mco2c + M3ig3pc - M2cpr5pc = 0
 R("R_PRAIi"): - Mpranc + M2cpr5pc = 0
 R("R_DDPA"): + Mpic - Mpepc - Mh2oc - Me4pc + M2dda7pc = 0
 R("R_DHQS"): + Mpic + M3dhqc - M2dda7pc = 0
 R("R_2AGPA120tipp"): - M2ddecg3pp + M2ddecg3pc = 0
 R("R_LPLIPAL2A120"): + 2 Mhc - Mh2oc + Mglyc3pc + Mddcac - M2ddecg3pc
 = 0
 R("R_PLIPA1A120pp"): - Mpa120p - Mh2op + Mddcap + M2ddecg3pp = 0
 R("R_DDGLK"): + Mhc - Matpc + Madpc - M2ddglcnc + M2ddg6pc = 0
 R("R_EDA"): + Mpyrc + Mg3pc - M2ddg6pc = 0
 R("R_EDD"): + Mh2oc - M6pgcc + M2ddg6pc = 0
 R("M2ddglcnb_Transp"): 0 M10fthfc = 0
 R("R_EX_2ddglcn_e_"): - M2ddglcne = 0
 R("R_EX_2ddglcn_e__Rev"): + M2ddglcne = 0
 R("R_ALTRH"): + Mh2oc - Maltrnc + M2ddglcnc = 0
 R("R_DDGLCNt2rpp"): - Mhp + Mhc - M2ddglcnp + M2ddglcnc = 0
 R("R_DDGLCNt2rpp_Rev"): + Mhp - Mhc + M2ddglcnp - M2ddglcnc = 0
 R("R_MNNH"): - Mmanac + Mh2oc + M2ddglcnc = 0
 R("R_DDGLCNtex"): + M2ddglcnp - M2ddglcne = 0
 R("R_DDGLCNtex_Rev"): - M2ddglcnp + M2ddglcne = 0
 R("R_DDGALK"): + Mhc - Matpc + Madpc - M2dh3dgalc + M2dh3dgal6pc = 0
 R("R_DDPGALA"): + Mpyrc + Mg3pc - M2dh3dgal6pc = 0
 R("R_DDPGALA_Rev"): - Mpyrc - Mg3pc + M2dh3dgal6pc = 0
 R("R_GALCTND"): + Mh2oc - MgalctnDc + M2dh3dgalc = 0
 R("R_2DGLCNRx"): - Mnadhc + Mnadc - Mhc + Mglcnc - M2dhglcnc = 0
 R("R_2DGLCNRy"): - Mnadphc + Mnadpc - Mhc + Mglcnc - M2dhglcnc = 0
 R("R_2DGULRx"): - Mnadhc + Mnadc + MidonLc - Mhc - M2dhgulnc = 0
 R("R_2DGULRy"): - Mnadphc + Mnadpc + MidonLc - Mhc - M2dhgulnc = 0
 R("R_DPR"): + MpantRc - Mnadphc + Mnadpc - Mhc - M2dhpc = 0
 R("R_MOHMT"): + Mthfc - Mmlthfc - Mh2oc - M3mobc + M2dhpc = 0
 R("R_DMSOR2"): + Mh2oc - Mdmsoc + Mdmsc - M2dmmql8c + M2dmmq8c = 0
 R("R_DMSOR2pp"): + Mh2op + Mdmsp - Mdmsop - M2dmmql8c + M2dmmq8c = 0
 R("R_FRD3"): + Msuccc - Mfumc - M2dmmql8c + M2dmmq8c = 0
 R("R_G3PD7"): - Mglyc3pc + Mdhapc + M2dmmql8c - M2dmmq8c = 0
 R("R_GLYCTO4"): - Mglycltc + Mglxc + M2dmmql8c - M2dmmq8c = 0
 R("R_HYD3pp"): + 2 Mhp - 2 Mhc - Mh2c + M2dmmql8c - M2dmmq8c = 0
 R("R_NADH18pp"): - Mnadhc + Mnadc + 3 Mhp - 4 Mhc + M2dmmql8c
 - M2dmmq8c = 0
 R("R_NADH9"): - Mnadhc + Mnadc - Mhc + M2dmmql8c - M2dmmq8c = 0
 R("R_NADPHQR4"): - Mnadphc + Mnadpc - Mhc + M2dmmql8c - M2dmmq8c = 0
 R("R_TMAOR2"): - Mtmaoc + Mtmac - Mhc + Mh2oc - M2dmmql8c + M2dmmq8c
 = 0
 R("R_TMAOR2pp"): + Mtmap - Mtmaop - Mhp + Mh2op - M2dmmql8c + M2dmmq8c
 = 0
 R("R_AMMQLT8"): + Mmql8c + Mhc - Mametc + Mahcysc - M2dmmql8c = 0
 R("R_DHNAOT4"): + Mppic - Moctdpc - Mhc - Mdhnac + Mco2c + M2dmmql8c
 = 0
 R("R_DURIPP"): + Murac - Mpic - Mduric + M2dr1pc = 0
 R("R_DURIPP_Rev"): - Murac + Mpic + Mduric - M2dr1pc = 0
 R("R_PPM2"): + M2dr5pc - M2dr1pc = 0
 R("R_PPM2_Rev"): - M2dr5pc + M2dr1pc = 0
 R("R_PUNP2"): - Mpic - Mdad2c + Madec + M2dr1pc = 0
 R("R_PUNP2_Rev"): + Mpic + Mdad2c - Madec - M2dr1pc = 0
 R("R_PUNP4"): - Mpic + Mguac - Mdgsnc + M2dr1pc = 0
 R("R_PUNP4_Rev"): + Mpic - Mguac + Mdgsnc - M2dr1pc = 0
 R("R_PUNP6"): - Mpic + Mhxanc - Mdinc + M2dr1pc = 0
 R("R_PUNP6_Rev"): + Mpic - Mhxanc + Mdinc - M2dr1pc = 0
 R("R_TMDPP"): - Mthymdc + Mthymc - Mpic + M2dr1pc = 0
 R("R_TMDPP_Rev"): + Mthymdc - Mthymc + Mpic - M2dr1pc = 0
 R("R_DRPA"): + Mg3pc + Macaldc - M2dr5pc = 0
 R("R_GLCRAL"): + Mpyrc - M5dh4dglcc + M2h3oppanc = 0
 R("R_GLXCL"): - Mhc - 2 Mglxc + Mco2c + M2h3oppanc = 0
 R("R_HPYRI"): - Mhpyrc + M2h3oppanc = 0
 R("R_HPYRI_Rev"): + Mhpyrc - M2h3oppanc = 0
 R("R_TRSARr"): - Mnadhc + Mnadc - Mhc + MglycRc - M2h3oppanc = 0
 R("R_TRSARr_Rev"): + Mnadhc - Mnadc + Mhc - MglycRc + M2h3oppanc = 0
 R("R_2AGPA161tipp"): - M2hdec9eg3pp + M2hdec9eg3pc = 0
 R("R_LPLIPAL2A161"): + Mhdceac + 2 Mhc - Mh2oc + Mglyc3pc
 - M2hdec9eg3pc = 0
 R("R_PLIPA1A161pp"): - Mpa161p + Mhdceap - Mh2op + M2hdec9eg3pp = 0
 R("R_2AGPA160tipp"): - M2hdecg3pp + M2hdecg3pc = 0
 R("R_LPLIPAL2A160"): + Mhdcac + 2 Mhc - Mh2oc + Mglyc3pc - M2hdecg3pc
 = 0
 R("R_PLIPA1A160pp"): - Mpa160p + Mhdcap - Mh2op + M2hdecg3pp = 0
 R("R_IPPMIa"): + Mh2oc - M3c2hmpc + M2ippmc = 0
 R("R_IPPMIa_Rev"): - Mh2oc + M3c2hmpc - M2ippmc = 0
 R("R_IPPMIb"): - Mh2oc + M3c3hmpc - M2ippmc = 0
 R("R_IPPMIb_Rev"): + Mh2oc - M3c3hmpc + M2ippmc = 0
 R("R_DKMPPD"): + Mpic - Mo2c + 2 Mhc - Mh2oc + Mforc - Mdkmppc + M2kmbc
 = 0
 R("R_DKMPPD3"): + Mpic + 2 Mnadhc - 2 Mnadc + 4 Mhc - 3 Mh2oc + Mforc
 - Mdkmppc + M2kmbc = 0
 R("R_UNK3"): + MmetLc - MgluLc + Makgc - M2kmbc = 0
 R("R_2MAHMP"): + Mpic + Mhc - Mh2oc + M4ampmc - M2mahmpc = 0
 R("R_PMPK"): - Matpc + Madpc - M4ampmc + M2mahmpc = 0
 R("R_TMPPP"): + Mthmmpc + Mppic - Mhc - M4mpetzc - M2mahmpc = 0
 R("R_MCITD"): + Mh2oc - M2mcitc + M2mcacnc = 0
 R("R_MICITD"): + Mmicitc - Mh2oc - M2mcacnc = 0
 R("R_MCITS"): - Mppcoac - Moaac + Mhc - Mh2oc + Mcoac + M2mcitc = 0
 R("R_DXPRIi"): - Mnadphc + Mnadpc - Mhc - Mdxyl5pc + M2me4pc = 0
 R("R_MEPCT"): + Mppic - Mhc - Mctpc + M4c2mec - M2me4pc = 0
 R("R_MECDPDH2"): - Mnadhc + Mnadc + Mh2oc + Mh2mb4pc - M2mecdpc = 0
 R("R_MECDPS"): + Mcmpc - M2p4c2mec + M2mecdpc = 0
 R("R_OBTFL"): + Mppcoac + Mforc - Mcoac - M2obutc = 0
 R("R_THRD_L"): - MthrLc + Mnh4c + M2obutc = 0
 R("R_2AGPA181tipp"): - M2odec11eg3pp + M2odec11eg3pc = 0
 R("R_LPLIPAL2A181"): + Mocdceac + 2 Mhc - Mh2oc + Mglyc3pc
 - M2odec11eg3pc = 0
 R("R_PLIPA1A181pp"): - Mpa181p + Mocdceap - Mh2op + M2odec11eg3pp = 0
 R("R_2AGPA180tipp"): - M2odecg3pp + M2odecg3pc = 0
 R("R_LPLIPAL2A180"): + Mocdcac + 2 Mhc - Mh2oc + Mglyc3pc - M2odecg3pc
 = 0
 R("R_PLIPA1A180pp"): - Mpa180p + Mocdcap - Mh2op + M2odecg3pp = 0
 R("R_OHPHM"): + Mhc - Mametc + Mahcysc + M2omphc - M2ohphc = 0
 R("R_OPHHX"): - 0.5 Mo2c - M2ophc + M2ohphc = 0
 R("R_OPHHX3"): + 2 Mpic + Mnadhc - Mnadc + 3 Mhc - 3 Mh2oc - 2 Matpc
 + 2 Madpc - M2ophc + M2ohphc = 0
 R("R_OMBZLM"): + Mhc - Mametc + Mahcysc + M2ommblc - M2ombzlc = 0
 R("R_OMPHHX"): - 0.5 Mo2c - M2omphc + M2ombzlc = 0
 R("R_OMPHHX3"): + 2 Mpic + Mnadhc - Mnadc + 3 Mhc - 3 Mh2oc - 2 Matpc
 + 2 Madpc - M2omphc + M2ombzlc = 0
 R("R_DMQMT"): + Mq8h2c + Mhc - Mametc + Mahcysc - M2omhmblc = 0
 R("R_OMMBLHX"): - 0.5 Mo2c - M2ommblc + M2omhmblc = 0
 R("R_OMMBLHX3"): + 2 Mpic + Mnadhc - Mnadc + 3 Mhc - 3 Mh2oc - 2 Matpc
 + 2 Madpc - M2ommblc + M2omhmblc = 0
 R("R_OPHBDC"): - Mhc + Mco2c - M3ophbc + M2ophc = 0
 R("R_CDPMEK"): + Mhc - Matpc + Madpc - M4c2mec + M2p4c2mec = 0
 R("R_ENO"): + Mpepc + Mh2oc - M2pgc = 0
 R("R_ENO_Rev"): - Mpepc - Mh2oc + M2pgc = 0
 R("R_GLYCK2"): + Mhc - MglycRc - Matpc + Madpc + M2pgc = 0
 R("R_PGM"): + M3pgc - M2pgc = 0
 R("R_PGM_Rev"): - M3pgc + M2pgc = 0
 R("R_PGLYCP"): + Mpic - Mh2oc + Mglycltc - M2pglycc = 0
 R("R_SHCHCS2"): + Mthmppc - Mssaltppc + Mpyrc - Michorc + M2shchcc = 0
 R("R_SUCBZS"): + Msucbzc + Mh2oc - M2shchcc = 0
 R("R_2AGPA141tipp"): - M2tdec7eg3pp + M2tdec7eg3pc = 0
 R("R_LPLIPAL2A141"): + Mttdceac + 2 Mhc - Mh2oc + Mglyc3pc
 - M2tdec7eg3pc = 0
 R("R_PLIPA1A141pp"): + Mttdceap - Mpa141p - Mh2op + M2tdec7eg3pp = 0
 R("R_2AGPA140tipp"): - M2tdecg3pp + M2tdecg3pc = 0
 R("R_LPLIPAL2A140"): + Mttdcac + 2 Mhc - Mh2oc + Mglyc3pc - M2tdecg3pc
 = 0
 R("R_PLIPA1A140pp"): + Mttdcap - Mpa140p - Mh2op + M2tdecg3pp = 0
 R("R_TPRDCOAS"): - Mdpcoac - Matpc + Madec + M2tpr3dpcoac = 0
 R("M34dhpacb_Transp"): 0 M10fthfc = 0
 R("R_EX_34dhpac_e_"): - M34dhpace = 0
 R("R_EX_34dhpac_e__Rev"): + M34dhpace = 0
 R("R_34dhpactex"): + M34dhpacp - M34dhpace = 0
 R("R_34dhpactex_Rev"): - M34dhpacp + M34dhpace = 0
 R("R_42A12BOOXpp"): - Mo2p + Mnh4p - Mh2op + Mh2o2p - Mdopap
 + M34dhpacp = 0
 R("R_PPND"): - Mpphnc + Mnadhc - Mnadc + Mco2c + M34hppc = 0
 R("R_TYRTA"): - MtyrLc + MgluLc - Makgc + M34hppc = 0
 R("R_TYRTA_Rev"): + MtyrLc - MgluLc + Makgc - M34hppc = 0
 R("M3ampb_Transp"): 0 M10fthfc = 0
 R("R_EX_3amp_e_"): - M3ampe = 0
 R("R_EX_3amp_e__Rev"): + M3ampe = 0
 R("R_3AMPtex"): + M3ampp - M3ampe = 0
 R("R_3AMPtex_Rev"): - M3ampp + M3ampe = 0
 R("R_3NTD7pp"): + Mpip - Mh2op + Madnp - M3ampp = 0
 R("R_IPMD"): + Mnadhc - Mnadc + Mhc + M3c4mopc - M3c2hmpc = 0
 R("R_IPPS"): + Mhc - Mh2oc + Mcoac - Maccoac - M3mobc + M3c3hmpc = 0
 R("R_OMCDC"): - Mhc + Mco2c + M4mopc - M3c4mopc = 0
 R("M3cmpb_Transp"): 0 M10fthfc = 0
 R("R_EX_3cmp_e_"): - M3cmpe = 0
 R("R_EX_3cmp_e__Rev"): + M3cmpe = 0
 R("R_3CMPtex"): + M3cmpp - M3cmpe = 0
 R("R_3CMPtex_Rev"): - M3cmpp + M3cmpe = 0
 R("R_3NTD4pp"): + Mpip - Mh2op + Mcytdp - M3cmpp = 0
 R("R_3KGK"): + Mhc - Matpc + Madpc + M3dhgulnpc - M3dhgulnc = 0
 R("R_ASCBPL"): + Mhc - Mh2oc - Mascb6pc + M3dhgulnpc = 0
 R("R_KG6PDC"): + Mxu5pLc - Mhc + Mco2c - M3dhgulnpc = 0
 R("R_DHQTi"): + Mh2oc + M3dhskc - M3dhqc = 0
 R("R_SHK3Dr"): + Mskmc - Mnadphc + Mnadpc - Mhc - M3dhskc = 0
 R("R_SHK3Dr_Rev"): - Mskmc + Mnadphc - Mnadpc + Mhc + M3dhskc = 0
 R("M3gmpb_Transp"): 0 M10fthfc = 0
 R("R_EX_3gmp_e_"): - M3gmpe = 0
 R("R_EX_3gmp_e__Rev"): + M3gmpe = 0
 R("R_3GMPtex"): + M3gmpp - M3gmpe = 0
 R("R_3GMPtex_Rev"): - M3gmpp + M3gmpe = 0
 R("R_3NTD9pp"): + Mpip - Mh2op + Mgsnp - M3gmpp = 0
 R("R_3HAD40"): + Mh2oc + Mbut2eACPc - M3haACPc = 0
 R("R_3OAR40"): - Mnadphc + Mnadpc - Mhc - MactACPc + M3haACPc = 0
 R("R_3OAR40_Rev"): + Mnadphc - Mnadpc + Mhc + MactACPc - M3haACPc = 0
 R("R_ECOAH1"): + Mh2oc + Mb2coac - M3hbcoac = 0
 R("R_ECOAH1_Rev"): - Mh2oc - Mb2coac + M3hbcoac = 0
 R("R_HACD1i"): + Mnadhc - Mnadc + Mhc + Maacoac - M3hbcoac = 0
 R("R_3HAD121"): + Mt3c5ddeceACPc + Mh2oc - M3hcddec5eACPc = 0
 R("R_3OAR121"): - Mnadphc + Mnadpc - Mhc - M3ocddec5eACPc
 + M3hcddec5eACPc = 0
 R("M3hcinnmb_Transp"): 0 M10fthfc = 0
 R("R_EX_3hcinnm_e_"): - M3hcinnme = 0
 R("R_EX_3hcinnm_e__Rev"): + M3hcinnme = 0
 R("R_3HCINNMH"): - Mo2c - Mnadhc + Mnadc - Mhc + Mh2oc + Mdhcinnmc
 - M3hcinnmc = 0
 R("R_HCINNMt2rpp"): - Mhp + Mhc - M3hcinnmp + M3hcinnmc = 0
 R("R_HCINNMt2rpp_Rev"): + Mhp - Mhc + M3hcinnmp - M3hcinnmc = 0
 R("R_HCINNMtex"): + M3hcinnmp - M3hcinnme = 0
 R("R_HCINNMtex_Rev"): - M3hcinnmp + M3hcinnme = 0
 R("R_3HAD141"): + Mt3c7mrseACPc + Mh2oc - M3hcmrs7eACPc = 0
 R("R_3OAR141"): - Mnadphc + Mnadpc - Mhc - M3ocmrs7eACPc
 + M3hcmrs7eACPc = 0
 R("R_3HAD161"): + Mt3c9palmeACPc + Mh2oc - M3hcpalm9eACPc = 0
 R("R_3OAR161"): - Mnadphc + Mnadpc - Mhc - M3ocpalm9eACPc
 + M3hcpalm9eACPc = 0
 R("R_3HAD181"): + Mt3c11vaceACPc + Mh2oc - M3hcvac11eACPc = 0
 R("R_3OAR181"): - Mnadphc + Mnadpc - Mhc - M3ocvac11eACPc
 + M3hcvac11eACPc = 0
 R("R_ECOAH4"): + Mh2oc + Mdc2coac - M3hdcoac = 0
 R("R_ECOAH4_Rev"): - Mh2oc - Mdc2coac + M3hdcoac = 0
 R("R_HACD4i"): + Mnadhc - Mnadc + Mhc + M3odcoac - M3hdcoac = 0
 R("R_ECOAH5"): + Mh2oc + Mdd2coac - M3hddcoac = 0
 R("R_ECOAH5_Rev"): - Mh2oc - Mdd2coac + M3hddcoac = 0
 R("R_HACD5i"): + Mnadhc - Mnadc + Mhc + M3oddcoac - M3hddcoac = 0
 R("R_3HAD120"): + Mtddec2eACPc + Mh2oc - M3hddecACPc = 0
 R("R_3OAR120"): - Mnadphc + Mnadpc - Mhc - M3oddecACPc + M3hddecACPc
 = 0
 R("R_3OAR120_Rev"): + Mnadphc - Mnadpc + Mhc + M3oddecACPc
 - M3hddecACPc = 0
 R("R_3HAD100"): + Mtdec2eACPc + Mh2oc - M3hdecACPc = 0
 R("R_3OAR100"): - Mnadphc + Mnadpc - Mhc - M3odecACPc + M3hdecACPc = 0
 R("R_3OAR100_Rev"): + Mnadphc - Mnadpc + Mhc + M3odecACPc - M3hdecACPc
 = 0
 R("R_ECOAH2"): + Mhx2coac + Mh2oc - M3hhcoac = 0
 R("R_ECOAH2_Rev"): - Mhx2coac - Mh2oc + M3hhcoac = 0
 R("R_HACD2i"): + Mnadhc - Mnadc + Mhc + M3ohcoac - M3hhcoac = 0
 R("R_ECOAH7"): + Mhdd2coac + Mh2oc - M3hhdcoac = 0
 R("R_ECOAH7_Rev"): - Mhdd2coac - Mh2oc + M3hhdcoac = 0
 R("R_HACD7i"): + Mnadhc - Mnadc + Mhc + M3ohdcoac - M3hhdcoac = 0
 R("R_3HAD60"): + Mthex2eACPc + Mh2oc - M3hhexACPc = 0
 R("R_3OAR60"): - Mnadphc + Mnadpc - Mhc - M3ohexACPc + M3hhexACPc = 0
 R("R_3OAR60_Rev"): + Mnadphc - Mnadpc + Mhc + M3ohexACPc - M3hhexACPc
 = 0
 R("R_3HAD140"): + Mtmrs2eACPc + Mh2oc - M3hmrsACPc = 0
 R("R_3OAR140"): - Mnadphc + Mnadpc - Mhc - M3omrsACPc + M3hmrsACPc = 0
 R("R_3OAR140_Rev"): + Mnadphc - Mnadpc + Mhc + M3omrsACPc - M3hmrsACPc
 = 0
 R("R_U23GAAT"): - Mu3hgac + Mu23gac + Mhc + MACPc - M3hmrsACPc = 0
 R("R_UAGAAT"): - Muacgamc + Mu3agac + MACPc - M3hmrsACPc = 0
 R("R_UAGAAT_Rev"): + Muacgamc - Mu3agac - MACPc + M3hmrsACPc = 0
 R("R_ECOAH3"): + Moc2coac + Mh2oc - M3hocoac = 0
 R("R_ECOAH3_Rev"): - Moc2coac - Mh2oc + M3hocoac = 0
 R("R_HACD3i"): + Mnadhc - Mnadc + Mhc + M3oocoac - M3hocoac = 0
 R("R_3HAD180"): + Mtoctd2eACPc + Mh2oc - M3hoctaACPc = 0
 R("R_3OAR180"): - Mnadphc + Mnadpc - Mhc - M3ooctdACPc + M3hoctaACPc
 = 0
 R("R_3OAR180_Rev"): + Mnadphc - Mnadpc + Mhc + M3ooctdACPc
 - M3hoctaACPc = 0
 R("R_3HAD80"): + Mtoct2eACPc + Mh2oc - M3hoctACPc = 0
 R("R_3OAR80"): - Mnadphc + Mnadpc - Mhc - M3ooctACPc + M3hoctACPc = 0
 R("R_3OAR80_Rev"): + Mnadphc - Mnadpc + Mhc + M3ooctACPc - M3hoctACPc
 = 0
 R("R_ECOAH8"): + Mod2coac + Mh2oc - M3hodcoac = 0
 R("R_ECOAH8_Rev"): - Mod2coac - Mh2oc + M3hodcoac = 0
 R("R_HACD8i"): + Mnadhc - Mnadc + Mhc + M3oodcoac - M3hodcoac = 0
 R("R_3HAD160"): + Mtpalm2eACPc + Mh2oc - M3hpalmACPc = 0
 R("R_3OAR160"): - Mnadphc + Mnadpc - Mhc - M3opalmACPc + M3hpalmACPc
 = 0
 R("R_3OAR160_Rev"): + Mnadphc - Mnadpc + Mhc + M3opalmACPc
 - M3hpalmACPc = 0
 R("M3hpppnb_Transp"): 0 M10fthfc = 0
 R("R_EX_3hpppn_e_"): - M3hpppne = 0
 R("R_EX_3hpppn_e__Rev"): + M3hpppne = 0
 R("R_3HPPPNH"): - Mo2c - Mnadhc + Mnadc - Mhc + Mh2oc + Mdhpppnc
 - M3hpppnc = 0
 R("R_HPPPNt2rpp"): - Mhp + Mhc - M3hpppnp + M3hpppnc = 0
 R("R_HPPPNt2rpp_Rev"): + Mhp - Mhc + M3hpppnp - M3hpppnc = 0
 R("R_HPPPNtex"): + M3hpppnp - M3hpppne = 0
 R("R_HPPPNtex_Rev"): - M3hpppnp + M3hpppne = 0
 R("R_ECOAH6"): + Mtd2coac + Mh2oc - M3htdcoac = 0
 R("R_ECOAH6_Rev"): - Mtd2coac - Mh2oc + M3htdcoac = 0
 R("R_HACD6i"): + Mnadhc - Mnadc + Mhc + M3otdcoac - M3htdcoac = 0
 R("R_TRPS1"): + MtrpLc - MserLc + Mh2oc + Mg3pc - M3ig3pc = 0
 R("R_TRPS3"): + Mindolec + Mg3pc - M3ig3pc = 0
 R("R_VALTA"): - MvalLc + MgluLc - Makgc + M3mobc = 0
 R("R_VALTA_Rev"): + MvalLc - MgluLc + Makgc - M3mobc = 0
 R("R_VPAMT"): + MvalLc + Mpyrc - MalaLc - M3mobc = 0
 R("R_ILETA"): - MileLc + MgluLc - Makgc + M3mopc = 0
 R("R_ILETA_Rev"): + MileLc - MgluLc + Makgc - M3mopc = 0
 R("R_3OAS121"): - MmalACPc - Mhc + Mco2c - Mcdec3eACPc + MACPc
 + M3ocddec5eACPc = 0
 R("R_3OAS141"): - MmalACPc - Mhc + Mco2c - Mcddec5eACPc + MACPc
 + M3ocmrs7eACPc = 0
 R("R_3OAS161"): - MtdeACPc - MmalACPc - Mhc + Mco2c + MACPc
 + M3ocpalm9eACPc = 0
 R("R_3OAS181"): - MmalACPc - MhdeACPc - Mhc + Mco2c + MACPc
 + M3ocvac11eACPc = 0
 R("R_KAT4"): + Moccoac - Mcoac + Maccoac - M3odcoac = 0
 R("R_KAT5"): + Mdcacoac - Mcoac + Maccoac - M3oddcoac = 0
 R("R_3OAS120"): - MmalACPc - Mhc - MdcaACPc + Mco2c + MACPc
 + M3oddecACPc = 0
 R("R_3OAS100"): - MocACPc - MmalACPc - Mhc + Mco2c + MACPc + M3odecACPc
 = 0
 R("R_KAT2"): - Mcoac + Mbtcoac + Maccoac - M3ohcoac = 0
 R("R_KAT7"): + Mtdcoac - Mcoac + Maccoac - M3ohdcoac = 0
 R("R_3OAS60"): - MmalACPc - Mhc + Mco2c - MbutACPc + MACPc + M3ohexACPc
 = 0
 R("R_3OAS140"): - MmalACPc - Mhc - MddcaACPc + Mco2c + MACPc
 + M3omrsACPc = 0
 R("R_KAT3"): + Mhxcoac - Mcoac + Maccoac - M3oocoac = 0
 R("R_3OAS80"): - MmalACPc - MhexACPc - Mhc + Mco2c + MACPc + M3ooctACPc
 = 0
 R("R_3OAS180"): - MpalmACPc - MmalACPc - Mhc + Mco2c + MACPc
 + M3ooctdACPc = 0
 R("R_KAT8"): + Mpmtcoac - Mcoac + Maccoac - M3oodcoac = 0
 R("R_3OAS160"): - MmyrsACPc - MmalACPc - Mhc + Mco2c + MACPc
 + M3opalmACPc = 0
 R("R_HBZOPT"): + Mppic - Moctdpc - M4hbzc + M3ophbc = 0
 R("R_KAT6"): + Mddcacoac - Mcoac + Maccoac - M3otdcoac = 0
 R("R_GLYCK"): + Mhc - MglycRc - Matpc + Madpc + M3pgc = 0
 R("R_PGCD"): + Mnadhc - Mnadc + Mhc + M3phpc - M3pgc = 0
 R("R_PSERT"): + MpserLc - MgluLc + Makgc - M3phpc = 0
 R("R_CHORS"): + Mpic + Mchorc - M3psmec = 0
 R("R_PSCVT"): - Mskm5pc + Mpic - Mpepc + M3psmec = 0
 R("R_PSCVT_Rev"): + Mskm5pc - Mpic + Mpepc - M3psmec = 0
 R("R_CYSSADS"): + Mso2c - 2 Mhc + MalaLc - M3salac = 0
 R("M3umpb_Transp"): 0 M10fthfc = 0
 R("R_EX_3ump_e_"): - M3umpe = 0
 R("R_EX_3ump_e__Rev"): + M3umpe = 0
 R("R_3UMPtex"): + M3umpp - M3umpe = 0
 R("R_3UMPtex_Rev"): - M3umpp + M3umpe = 0
 R("R_3NTD2pp"): + Murip + Mpip - Mh2op - M3umpp = 0
 R("M4abutb_Transp"): 0 M10fthfc = 0
 R("R_EX_4abut_e_"): - M4abute = 0
 R("R_EX_4abut_e__Rev"): + M4abute = 0
 R("R_ABTA"): + Msucsalc + MgluLc - Makgc - M4abutc = 0
 R("R_ABUTD"): + Mnadhc - Mnadc + 2 Mhc - Mh2oc - M4abutnc + M4abutc = 0
 R("R_ABUTt2pp"): - Mhp + Mhc - M4abutp + M4abutc = 0
 R("R_GGGABAH"): - Mh2oc + MgluLc - Mgg4abutc + M4abutc = 0
 R("R_GLUABUTt7pp"): - MgluLp + MgluLc + M4abutp - M4abutc = 0
 R("R_GLUABUTt7pp_Rev"): + MgluLp - MgluLc - M4abutp + M4abutc = 0
 R("R_GLUDC"): - Mhc - MgluLc + Mco2c + M4abutc = 0
 R("R_ABUTtex"): + M4abutp - M4abute = 0
 R("R_ABUTtex_Rev"): - M4abutp + M4abute = 0
 R("R_PTRCTA"): - Mptrcc + MgluLc - Makgc + M4abutnc = 0
 R("R_ADCL"): + Mpyrc + Mhc - M4adchoc + M4abzc = 0
 R("R_DHPS2"): + Mppic + Mdhptc - M6hmhptppc - M4abzc = 0
 R("R_ADCS"): + MgluLc - MglnLc - Mchorc + M4adchoc = 0
 R("R_HMPK1"): + Mhc - Matpc + Madpc + M4ampmc - M4ahmmpc = 0
 R("R_AMPMS2"): + Mnadhc - Mnadc + 3 Mhc - Mh2oc + 2 Mforc - Mairc
 + M4ampmc = 0
 R("R_HOPNTAL"): + Mpyrc + Macaldc - M4h2opntnc = 0
 R("R_OP4ENH"): - Mop4enc - Mh2oc + M4h2opntnc = 0
 R("M4hbab_Transp"): 0 M10fthfc = 0
 R("R_DM_4HBA"): - M4hbac = 0
 R("R_THZPSN"): - MtyrLc + Mppic + Mhc + Mh2oc - Mdxyl5pc - McysLc
 + Mco2c - Matpc + Mampc + MalaLc + M4mpetzc + M4hbac = 0
 R("R_CHRPL"): + Mpyrc - Mchorc + M4hbzc = 0
 R("M4hoxpacdb_Transp"): 0 M10fthfc = 0
 R("R_EX_4hoxpacd_e_"): - M4hoxpacde = 0
 R("R_EX_4hoxpacd_e__Rev"): + M4hoxpacde = 0
 R("R_4HOXPACDtex"): + M4hoxpacdp - M4hoxpacde = 0
 R("R_4HOXPACDtex_Rev"): - M4hoxpacdp + M4hoxpacde = 0
 R("R_TYROXDApp"): - Mtymp - Mo2p + Mnh4p - Mh2op + Mh2o2p + M4hoxpacdp
 = 0
 R("R_4HTHRS"): + Mpic - Mphthrc - Mh2oc + M4hthrc = 0
 R("R_HETZK"): + Mhc - Matpc + Madpc + M4mpetzc - M4mhetzc = 0
 R("R_LEUTAi"): + MleuLc - MgluLc + Makgc - M4mopc = 0
 R("R_ASAD"): - Mpic + Mnadphc - Mnadpc + Mhc - Maspsac + M4paspc = 0
 R("R_ASAD_Rev"): + Mpic - Mnadphc + Mnadpc - Mhc + Maspsac - M4paspc
 = 0
 R("R_ASPK"): - Matpc - MaspLc + Madpc + M4paspc = 0
 R("R_ASPK_Rev"): + Matpc + MaspLc - Madpc - M4paspc = 0
 R("R_E4PD"): + Mnadhc - Mnadc + 2 Mhc - Mh2oc - Me4pc + M4perc = 0
 R("R_E4PD_Rev"): - Mnadhc + Mnadc - 2 Mhc + Mh2oc + Me4pc - M4perc = 0
 R("R_PERD"): + Mohpbc + Mnadhc - Mnadc + Mhc - M4perc = 0
 R("R_PERD_Rev"): - Mohpbc - Mnadhc + Mnadc - Mhc + M4perc = 0
 R("R_PNTK"): - MpntoRc + Mhc - Matpc + Madpc + M4ppanc = 0
 R("R_PPNCL2"): + Mppic + Mhc - McysLc - Mctpc + Mcmpc + M4ppcysc
 - M4ppanc = 0
 R("R_PPCDC"): + Mpan4pc - Mhc + Mco2c - M4ppcysc = 0
 R("R_PMDPHT"): + Mpic - Mh2oc - M5aprbuc + M4r5auc = 0
 R("R_RBFSa"): + Mpic + 2 Mh2oc + Mdmlzc - Mdb4pc - M4r5auc = 0
 R("R_RBFSb"): + Mribflvc - 2 Mdmlzc + M4r5auc = 0
 R("R_AIRC3"): + M5caizc - M5aizcc = 0
 R("R_AIRC3_Rev"): - M5caizc + M5aizcc = 0
 R("R_G1SAT"): - Mglu1sac + M5aopc = 0
 R("R_G1SAT_Rev"): + Mglu1sac - M5aopc = 0
 R("R_PPBNGS"): + Mppbngc + Mhc + 2 Mh2oc - 2 M5aopc = 0
 R("R_APRAUR"): - Mnadphc + Mnadpc - Mhc - M5apruc + M5aprbuc = 0
 R("R_AIRC2"): + Mpic - Mhco3c + Mhc - Matpc - Mairc + Madpc + M5caizc
 = 0
 R("M5dglcnb_Transp"): 0 M10fthfc = 0
 R("R_EX_5dglcn_e_"): - M5dglcne = 0
 R("R_EX_5dglcn_e__Rev"): + M5dglcne = 0
 R("R_5DGLCNR"): - Mnadphc + Mnadpc - Mhc + Mglcnc - M5dglcnc = 0
 R("R_5DGLCNR_Rev"): + Mnadphc - Mnadpc + Mhc - Mglcnc + M5dglcnc = 0
 R("R_5DGLCNt2rpp"): - Mhp + Mhc - M5dglcnp + M5dglcnc = 0
 R("R_5DGLCNt2rpp_Rev"): + Mhp - Mhc + M5dglcnp - M5dglcnc = 0
 R("R_IDOND"): - Mnadhc + Mnadc + MidonLc - Mhc - M5dglcnc = 0
 R("R_IDOND2"): - Mnadphc + Mnadpc + MidonLc - Mhc - M5dglcnc = 0
 R("R_IDOND_Rev"): + Mnadhc - Mnadc - MidonLc + Mhc + M5dglcnc = 0
 R("R_5DGLCNtex"): + M5dglcnp - M5dglcne = 0
 R("R_5DGLCNtex_Rev"): - M5dglcnp + M5dglcne = 0
 R("R_GALCTD"): + Mh2oc - MgalctDc + M5dh4dglcc = 0
 R("R_GLCRD"): + Mh2oc - Mglcrc + M5dh4dglcc = 0
 R("M5dribb_Transp"): 0 M10fthfc = 0
 R("R_DM_5DRIB"): - M5dribc = 0
 R("R_5DOAN"): - Mh2oc - Mdad5c + Madec + M5dribc = 0
 R("R_MTRI"): + M5mdru1pc - M5mdr1pc = 0
 R("R_MTRI_Rev"): - M5mdru1pc + M5mdr1pc = 0
 R("R_MTRK"): + Mhc - Matpc + Madpc - M5mtrc + M5mdr1pc = 0
 R("R_MDRPD"): + Mh2oc + Mdkmppc - M5mdru1pc = 0
 R("R_MTAN"): - Mh2oc + Madec + M5mtrc - M5mtac = 0
 R("R_SPMS"): + Mspmdc - Mptrcc + Mhc - Mametamc + M5mtac = 0
 R("R_METS"): + Mthfc + MmetLc - MhcysLc + Mhc - M5mthfc = 0
 R("R_MTHFR2"): - Mnadhc + Mnadc - Mmlthfc - 2 Mhc + M5mthfc = 0
 R("R_NNDMBRT"): - Mnicrntc + Mnacc + Mhc - Mdmbzidc + M5prdmbzc = 0
 R("R_RZ5PP"): + Mrdmbzic + Mpic - Mh2oc - M5prdmbzc = 0
 R("R_DHNPA2"): + Mgcaldc - Mdhnptc + M6hmhptc = 0
 R("R_HPPK2"): + Mhc - Matpc + Mampc + M6hmhptppc - M6hmhptc = 0
 R("R_GND"): + Mru5pDc + Mnadphc - Mnadpc + Mco2c - M6pgcc = 0
 R("R_GNK"): + Mhc - Mglcnc - Matpc + Madpc + M6pgcc = 0
 R("R_PGL"): + Mhc - Mh2oc - M6pglc + M6pgcc = 0
 R("R_G6PDH2r"): + Mnadphc - Mnadpc + Mhc - Mg6pc + M6pglc = 0
 R("R_G6PDH2r_Rev"): - Mnadphc + Mnadpc - Mhc + Mg6pc - M6pglc = 0
 R("R_AMAOTr"): + Mdannc + Mamobc - Mametc - M8aonnc = 0
 R("R_AMAOTr_Rev"): - Mdannc - Mamobc + Mametc + M8aonnc = 0
 R("R_AOXSr"): - Mpmcoac - Mhc + Mcoac + Mco2c - MalaLc + M8aonnc = 0
 R("R_AOXSr_Rev"): + Mpmcoac + Mhc - Mcoac - Mco2c + MalaLc - M8aonnc
 = 0
 R("Maacaldb_Transp"): 0 M10fthfc = 0
 R("R_DM_AACALD"): - Maacaldc = 0
 R("R_TAUDO"): - Mtaurc + Msuccc + Mso3c - Mo2c + Mhc + Mco2c - Makgc
 + Maacaldc = 0
 R("R_ACACCT"): - Maccoac + Macc - Macacc + Maacoac = 0
 R("R_ACACT1r"): + Mcoac - 2 Maccoac + Maacoac = 0
 R("R_ACACT1r_Rev"): - Mcoac + 2 Maccoac - Maacoac = 0
 R("R_KAT1"): - Mcoac + 2 Maccoac - Maacoac = 0
 R("R_AACTOOR"): - Mo2c + Mnh4c + Mmthgxlc - Mh2oc + Mh2o2c - Maactc = 0
 R("Macacb_Transp"): 0 M10fthfc = 0
 R("R_EX_acac_e_"): - Macace = 0
 R("R_EX_acac_e__Rev"): + Macace = 0
 R("R_ACACt2pp"): - Mhp + Mhc - Macacp + Macacc = 0
 R("R_ACACt2pp_Rev"): + Mhp - Mhc + Macacp - Macacc = 0
 R("R_ACACtex"): + Macacp - Macace = 0
 R("R_ACACtex_Rev"): - Macacp + Macace = 0
 R("R_ACOATA"): + Mcoac - MACPc - Maccoac + MacACPc = 0
 R("R_ACOATA_Rev"): - Mcoac + MACPc + Maccoac - MacACPc = 0
 R("R_KAS14"): - MmalACPc - Mhc + Mco2c + MactACPc + MACPc - MacACPc = 0
 R("R_MACPD"): - MmalACPc - Mhc + Mco2c + MacACPc = 0
 R("Macaldb_Transp"): 0 M10fthfc = 0
 R("R_EX_acald_e_"): - Macalde = 0
 R("R_EX_acald_e__Rev"): + Macalde = 0
 R("R_ACALD"): + Mnadhc - Mnadc + Mhc - Mcoac + Maccoac - Macaldc = 0
 R("R_ACALD_Rev"): - Mnadhc + Mnadc - Mhc + Mcoac - Maccoac + Macaldc
 = 0
 R("R_ACALDtpp"): - Macaldp + Macaldc = 0
 R("R_ACALDtpp_Rev"): + Macaldp - Macaldc = 0
 R("R_ALCD2x"): + Mnadhc - Mnadc + Mhc - Metohc + Macaldc = 0
 R("R_ALCD2x_Rev"): - Mnadhc + Mnadc - Mhc + Metohc - Macaldc = 0
 R("R_ALDD2x"): + Mnadhc - Mnadc + 2 Mhc - Mh2oc + Macc - Macaldc = 0
 R("R_ALDD2y"): + Mnadphc - Mnadpc + 2 Mhc - Mh2oc + Macc - Macaldc = 0
 R("R_ETHAAL"): + Mnh4c - Methac + Macaldc = 0
 R("R_FDMO3"): + Mso3c - Mo2c + Mhc + Mh2oc - Mfmnh2c + Mfmnc - Methso3c
 + Macaldc = 0
 R("R_THRA2i"): + Mglyc - MathrLc + Macaldc = 0
 R("R_THRAi"): - MthrLc + Mglyc + Macaldc = 0
 R("R_ACALDtex"): + Macaldp - Macalde = 0
 R("R_ACALDtex_Rev"): - Macaldp + Macalde = 0
 R("R_ACANTHAT"): + Mcoac - Manthc - Maccoac + Macanthc = 0
 R("Macb_Transp"): 0 M10fthfc = 0
 R("R_EX_ac_e_"): - Mace = 0
 R("R_EX_ac_e__Rev"): + Mace = 0
 R("R_ACKr"): - Matpc + Madpc + Mactpc - Macc = 0
 R("R_ACKr_Rev"): + Matpc - Madpc - Mactpc + Macc = 0
 R("R_ACODA"): + Mornc - Mh2oc - Macornc + Macc = 0
 R("R_ACS"): + Mppic - Mcoac - Matpc + Mampc + Maccoac - Macc = 0
 R("R_ACt2rpp"): - Mhp + Mhc - Macp + Macc = 0
 R("R_ACt2rpp_Rev"): + Mhp - Mhc + Macp - Macc = 0
 R("R_ACt4pp"): - Mna1p + Mna1c - Macp + Macc = 0
 R("R_AGDC"): - Mh2oc + Mgam6pc - Macgam6pc + Macc = 0
 R("R_BUTCT"): - Mbutc + Mbtcoac - Maccoac + Macc = 0
 R("R_CITL"): + Moaac - Mcitc + Macc = 0
 R("R_CYSS"): + Mhc - Mh2sc + McysLc - Macserc + Macc = 0
 R("R_HXCT"): + Mhxcoac - Mhxac - Maccoac + Macc = 0
 R("R_NACODA"): - Mh2oc + Mglu5sac - Macg5sac + Macc = 0
 R("R_POX"): + Mq8h2c - Mq8c - Mpyrc - Mh2oc + Mco2c + Macc = 0
 R("R_UHGADA"): + Mu3hgac - Mu3agac - Mh2oc + Macc = 0
 R("R_ACCOAC"): + Mpic + Mmalcoac - Mhco3c + Mhc - Matpc + Madpc
 - Maccoac = 0
 R("R_ACGS"): + Mhc - MgluLc + Mcoac + Macgluc - Maccoac = 0
 R("R_CS"): - Moaac + Mhc - Mh2oc + Mcoac + Mcitc - Maccoac = 0
 R("R_G1PACT"): + Mhc - Mgam1pc + Mcoac + Macgam1pc - Maccoac = 0
 R("R_GLCATr"): - MglcDc + Mcoac + MacglcDc - Maccoac = 0
 R("R_GLCATr_Rev"): + MglcDc - Mcoac - MacglcDc + Maccoac = 0
 R("R_KAS15"): - MmalACPc - Mhc + Mcoac + Mco2c + MactACPc - Maccoac = 0
 R("R_MALS"): + MmalLc + Mhc - Mh2oc - Mglxc + Mcoac - Maccoac = 0
 R("R_MALTATr"): - Mmaltc + Mcoac + Macmaltc - Maccoac = 0
 R("R_MALTATr_Rev"): + Mmaltc - Mcoac - Macmaltc + Maccoac = 0
 R("R_O16AT"): - Mragundc + Mcoac + Maragundc - Maccoac = 0
 R("R_PDH"): - Mpyrc + Mnadhc - Mnadc - Mcoac + Mco2c + Maccoac = 0
 R("R_PFL"): - Mpyrc + Mforc - Mcoac + Maccoac = 0
 R("R_PTAr"): - Mpic + Mcoac + Mactpc - Maccoac = 0
 R("R_PTAr_Rev"): + Mpic - Mcoac - Mactpc + Maccoac = 0
 R("R_SERAT"): - MserLc + Mcoac + Macserc - Maccoac = 0
 R("R_SERAT_Rev"): + MserLc - Mcoac - Macserc + Maccoac = 0
 R("R_SPMDAT1"): - Mspmdc + MN1aspmdc + Mhc + Mcoac - Maccoac = 0
 R("R_SPMDAT2"): - Mspmdc + Mn8aspmdc + Mhc + Mcoac - Maccoac = 0
 R("R_TDPADGAT"): + Mhc - Mdtdp4addgc + Mdtdp4aaddgc + Mcoac - Maccoac
 = 0
 R("R_ACtex"): + Macp - Mace = 0
 R("R_ACtex_Rev"): - Macp + Mace = 0
 R("R_ALR2"): - Mnadphc + Mnadpc - Mmthgxlc - Mhc + Macetolc = 0
 R("R_ACGK"): - Matpc + Madpc - Macgluc + Macg5pc = 0
 R("R_AGPR"): - Mpic + Mnadphc - Mnadpc + Mhc - Macg5sac + Macg5pc = 0
 R("R_AGPR_Rev"): + Mpic - Mnadphc + Mnadpc - Mhc + Macg5sac - Macg5pc
 = 0
 R("R_ACOTA"): + MgluLc - Makgc - Macornc + Macg5sac = 0
 R("R_ACOTA_Rev"): - MgluLc + Makgc + Macornc - Macg5sac = 0
 R("Macgal1pb_Transp"): 0 M10fthfc = 0
 R("R_EX_acgal1p_e_"): - Macgal1pe = 0
 R("R_EX_acgal1p_e__Rev"): + Macgal1pe = 0
 R("R_ACGAL1Ptex"): + Macgal1pp - Macgal1pe = 0
 R("R_ACGAL1Ptex_Rev"): - Macgal1pp + Macgal1pe = 0
 R("R_ACGAL1PPpp"): + Mpip - Mh2op + Macgalp - Macgal1pp = 0
 R("R_UACGALPpp"): + Mumpp - Mudpacgalp + 2 Mhp - Mh2op + Macgal1pp = 0
 R("Macgalb_Transp"): 0 M10fthfc = 0
 R("R_EX_acgal_e_"): - Macgale = 0
 R("R_EX_acgal_e__Rev"): + Macgale = 0
 R("R_ACGALtex"): + Macgalp - Macgale = 0
 R("R_ACGALtex_Rev"): - Macgalp + Macgale = 0
 R("Macgam1pb_Transp"): 0 M10fthfc = 0
 R("R_EX_acgam1p_e_"): - Macgam1pe = 0
 R("R_EX_acgam1p_e__Rev"): + Macgam1pe = 0
 R("R_UAGDP"): - Mutpc + Muacgamc + Mppic - Mhc - Macgam1pc = 0
 R("R_ACGAM1Ptex"): + Macgam1pp - Macgam1pe = 0
 R("R_ACGAM1Ptex_Rev"): - Macgam1pp + Macgam1pe = 0
 R("R_ACGAM1PPpp"): + Mpip - Mh2op + Macgamp - Macgam1pp = 0
 R("R_UACGAMPpp"): + Mumpp - Muacgamp + 2 Mhp - Mh2op + Macgam1pp = 0
 R("R_ACGAMK"): + Mhc - Matpc + Madpc - Macgamc + Macgam6pc = 0
 R("R_ACGAptspp"): + Mpyrc - Mpepc - Macgamp + Macgam6pc = 0
 R("R_ACM6PH"): + MlacDc - Mh2oc - Macmum6pc + Macgam6pc = 0
 R("R_AMANAPEr"): - Macmanapc + Macgam6pc = 0
 R("R_AMANAPEr_Rev"): + Macmanapc - Macgam6pc = 0
 R("Macgamb_Transp"): 0 M10fthfc = 0
 R("R_EX_acgam_e_"): - Macgame = 0
 R("R_EX_acgam_e__Rev"): + Macgame = 0
 R("R_AGM3PH"): - Mh2oc + Manhm3pc - Manhgm3pc + Macgamc = 0
 R("R_AGM4PH"): - Mh2oc + Manhm4pc - Manhgm4pc + Macgamc = 0
 R("R_AGMH"): - Mh2oc + Manhmc - Manhgmc + Macgamc = 0
 R("R_ACGAtex"): + Macgamp - Macgame = 0
 R("R_ACGAtex_Rev"): - Macgamp + Macgame = 0
 R("Macmanab_Transp"): 0 M10fthfc = 0
 R("R_EX_acmana_e_"): - Macmanae = 0
 R("R_EX_acmana_e__Rev"): + Macmanae = 0
 R("R_ACNML"): + Mpyrc - Macnamc + Macmanac = 0
 R("R_AMANK"): + Mhc - Matpc + Madpc + Macmanapc - Macmanac = 0
 R("R_ACMANAtex"): + Macmanap - Macmanae = 0
 R("R_ACMANAtex_Rev"): - Macmanap + Macmanae = 0
 R("R_ACMANAptspp"): + Mpyrc - Mpepc + Macmanapc - Macmanap = 0
 R("R_ACMUMptspp"): + Mpyrc - Mpepc - Macmump + Macmum6pc = 0
 R("R_ANHMK"): + Mhc - Mh2oc - Matpc - Manhmc + Madpc + Macmum6pc = 0
 R("Macmumb_Transp"): 0 M10fthfc = 0
 R("R_EX_acmum_e_"): - Macmume = 0
 R("R_EX_acmum_e__Rev"): + Macmume = 0
 R("R_ACMUMtex"): + Macmump - Macmume = 0
 R("Macnamb_Transp"): 0 M10fthfc = 0
 R("R_EX_acnam_e_"): - Macname = 0
 R("R_EX_acnam_e__Rev"): + Macname = 0
 R("R_ACNAMt2pp"): - Mhp + Mhc - Macnamp + Macnamc = 0
 R("R_ACNAMtex"): + Macnamp - Macname = 0
 R("R_ACNAMtex_Rev"): - Macnamp + Macname = 0
 R("Macolipab_Transp"): 0 M10fthfc = 0
 R("R_EX_acolipa_e_"): - Macolipae = 0
 R("R_EX_acolipa_e__Rev"): + Macolipae = 0
 R("R_ACOLIPAtex"): - Macolipap + Macolipae = 0
 R("R_LA4NTpp"): - MuLa4np + Mudcppp - Mcolipap + Macolipap = 0
 R("R_ACONIs"): - MaconTc + MaconCc = 0
 R("R_ACONIs_Rev"): + MaconTc - MaconCc = 0
 R("R_ACONTa"): + Mh2oc - Mcitc + MaconCc = 0
 R("R_ACONTa_Rev"): - Mh2oc + Mcitc - MaconCc = 0
 R("R_ACONTb"): + Micitc - Mh2oc - MaconCc = 0
 R("R_ACONTb_Rev"): - Micitc + Mh2oc + MaconCc = 0
 R("R_ACONMT"): - Mametc + Mahcysc - MaconTc + Maconmc = 0
 R("R_AACPS1"): - Mttdcac + Mppic + MmyrsACPc - Matpc + Mampc - MACPc
 = 0
 R("R_AACPS2"): - Mttdceac + MtdeACPc + Mppic - Matpc + Mampc - MACPc
 = 0
 R("R_AACPS3"): + Mppic + MpalmACPc - Mhdcac - Matpc + Mampc - MACPc = 0
 R("R_AACPS4"): + Mppic + MhdeACPc - Mhdceac - Matpc + Mampc - MACPc = 0
 R("R_AACPS5"): + Mppic + MocteACPc - Mocdceac - Matpc + Mampc - MACPc
 = 0
 R("R_AACPS6"): + Mppic - Mocdcac + MocdcaACPc - Matpc + Mampc - MACPc
 = 0
 R("R_AACPS7"): + Mppic - Mddcac + MddcaACPc - Matpc + Mampc - MACPc = 0
 R("R_AACPS8"): + Mppic - Mdcac + MdcaACPc - Matpc + Mampc - MACPc = 0
 R("R_AACPS9"): + Mppic - Moctac + MocACPc - Matpc + Mampc - MACPc = 0
 R("R_ACPS1"): + Mpapc + Mhc - Mcoac - MapoACPc + MACPc = 0
 R("R_EDTXS1"): + Mkdo2lipid4Lc - Mkdo2lipid4c - MddcaACPc + MACPc = 0
 R("R_EDTXS2"): - MmyrsACPc + Mlipac - Mkdo2lipid4Lc + MACPc = 0
 R("R_EDTXS3"): + Mkdo2lipid4pc - Mkdo2lipid4c - MhdeACPc + MACPc = 0
 R("R_EDTXS4"): - MmyrsACPc + Mlipacoldc - Mkdo2lipid4pc + MACPc = 0
 R("R_FA100ACPHi"): + Mhc - Mh2oc + Mdcac - MdcaACPc + MACPc = 0
 R("R_FA120ACPHi"): + Mhc - Mh2oc + Mddcac - MddcaACPc + MACPc = 0
 R("R_FA140ACPHi"): + Mttdcac - MmyrsACPc + Mhc - Mh2oc + MACPc = 0
 R("R_FA141ACPHi"): + Mttdceac - MtdeACPc + Mhc - Mh2oc + MACPc = 0
 R("R_FA160ACPHi"): - MpalmACPc + Mhdcac + Mhc - Mh2oc + MACPc = 0
 R("R_FA161ACPHi"): - MhdeACPc + Mhdceac + Mhc - Mh2oc + MACPc = 0
 R("R_FA80ACPHi"): + Moctac - MocACPc + Mhc - Mh2oc + MACPc = 0
 R("R_MCOATA"): - Mmalcoac + MmalACPc + Mcoac - MACPc = 0
 R("R_MCOATA_Rev"): + Mmalcoac - MmalACPc - Mcoac + MACPc = 0
 R("Macserb_Transp"): 0 M10fthfc = 0
 R("R_EX_acser_e_"): - Macsere = 0
 R("R_EX_acser_e__Rev"): + Macsere = 0
 R("R_ACSERtpp"): + Macserp - Macserc = 0
 R("R_ACSERtex"): + Macserp - Macsere = 0
 R("R_ACSERtex_Rev"): - Macserp + Macsere = 0
 R("Madeb_Transp"): 0 M10fthfc = 0
 R("R_EX_ade_e_"): - Madee = 0
 R("R_EX_ade_e__Rev"): + Madee = 0
 R("R_ADD"): + Mnh4c + Mhxanc - Mhc - Mh2oc - Madec = 0
 R("R_ADEt2rpp"): - Mhp + Mhc - Madep + Madec = 0
 R("R_ADEt2rpp_Rev"): + Mhp - Mhc + Madep - Madec = 0
 R("R_ADNUC"): + MribDc - Mh2oc - Madnc + Madec = 0
 R("R_ADPT"): - Mprppc + Mppic + Mampc - Madec = 0
 R("R_AHCYSNS"): + Mrhcysc - Mh2oc - Mahcysc + Madec = 0
 R("R_AMPN"): + Mr5pc - Mh2oc - Mampc + Madec = 0
 R("R_PUNP1"): + Mr1pc - Mpic - Madnc + Madec = 0
 R("R_PUNP1_Rev"): - Mr1pc + Mpic + Madnc - Madec = 0
 R("R_ADEtex"): + Madep - Madee = 0
 R("R_ADEtex_Rev"): - Madep + Madee = 0
 R("Madnb_Transp"): 0 M10fthfc = 0
 R("R_EX_adn_e_"): - Madne = 0
 R("R_EX_adn_e__Rev"): + Madne = 0
 R("R_ADA"): + Mnh4c + Minsc - Mhc - Mh2oc - Madnc = 0
 R("R_ADNK1"): + Mhc - Matpc + Mampc + Madpc - Madnc = 0
 R("R_ADNt2pp"): - Mhp + Mhc - Madnp + Madnc = 0
 R("R_ADNt2rpp"): - Mhp + Mhc - Madnp + Madnc = 0
 R("R_ADNt2rpp_Rev"): + Mhp - Mhc + Madnp - Madnc = 0
 R("R_NTD7"): + Mpic - Mh2oc - Mampc + Madnc = 0
 R("R_ADNtex"): + Madnp - Madne = 0
 R("R_ADNtex_Rev"): - Madnp + Madne = 0
 R("R_NTD7pp"): + Mpip - Mh2op - Mampp + Madnp = 0
 R("R_ADOCBIK"): + Mhc - Matpc + Madpc + Madocbipc - Madocbic = 0
 R("R_CBIAT"): + Mpppic - Mhc - Mcbic - Matpc + Madocbic = 0
 R("R_CBIAT_Rev"): - Mpppic + Mhc + Mcbic + Matpc - Madocbic = 0
 R("R_ACBIPGT"): + Mppic - Mhc - Mgtpc + Magdpcbic - Madocbipc = 0
 R("Madocblb_Transp"): 0 M10fthfc = 0
 R("R_EX_adocbl_e_"): - Madocble = 0
 R("R_EX_adocbl_e__Rev"): + Madocble = 0
 R("R_ADOCBLabcpp"): + Mpic + Mhc - Mh2oc - Matpc + Madpc - Madocblp
 + Madocblc = 0
 R("R_ADOCBLS"): - Mrdmbzic + Mhc + Mgmpc - Magdpcbic + Madocblc = 0
 R("R_CBLAT"): + Mpppic - Mhc - Mcbl1c - Matpc + Madocblc = 0
 R("R_CBLAT_Rev"): - Mpppic + Mhc + Mcbl1c + Matpc - Madocblc = 0
 R("R_ADOCBLtonex"): - Mhp + Mhc + Madocblp - Madocble = 0
 R("R_3PEPTabcpp"): + Mpic - MLalaDgluMdapp + MLalaDgluMdapc + Mhc
 - Mh2oc - Matpc + Madpc = 0
 R("R_4PEPTabcpp"): + Mpic - MLalaDgluMdapDalap + MLalaDgluMdapDalac
 + Mhc - Mh2oc - Matpc + Madpc = 0
 R("R_ACCOAL"): + Mppcoac - Mppac + Mpic - Mcoac - Matpc + Madpc = 0
 R("R_ADK1"): - Matpc - Mampc + 2 Madpc = 0
 R("R_ADK1_Rev"): + Matpc + Mampc - 2 Madpc = 0
 R("R_ADK3"): - Mgtpc + Mgdpc - Mampc + Madpc = 0
 R("R_ADK3_Rev"): + Mgtpc - Mgdpc + Mampc - Madpc = 0
 R("R_ADK4"): - Mitpc + Midpc - Mampc + Madpc = 0
 R("R_ADK4_Rev"): + Mitpc - Midpc + Mampc - Madpc = 0
 R("R_ADSK"): + Mpapsc + Mhc - Matpc - Mapsc + Madpc = 0
 R("R_ALAabcpp"): + Mpic + Mhc - Mh2oc - Matpc - MalaLp + MalaLc + Madpc
 = 0
 R("R_ALAALAabcpp"): + Mpic + Mhc - Mh2oc - Matpc - Malaalap + Malaalac
 + Madpc = 0
 R("R_ALAALAr"): + Mpic + Mhc - Matpc - 2 MalaDc + Malaalac + Madpc = 0
 R("R_ALAALAr_Rev"): - Mpic - Mhc + Matpc + 2 MalaDc - Malaalac - Madpc
 = 0
 R("R_ALLabcpp"): + Mpic + Mhc - Mh2oc - Matpc - MallDp + MallDc + Madpc
 = 0
 R("R_ALLK"): + Mhc - Matpc - MallDc + Mall6pc + Madpc = 0
 R("R_AP4AH"): + 2 Mhc - Mh2oc - Map4ac + 2 Madpc = 0
 R("R_AP5AH"): + 2 Mhc - Mh2oc + Matpc - Map5ac + Madpc = 0
 R("R_ARBabcpp"): + Mpic + Mhc - Mh2oc - Matpc - MarabLp + MarabLc
 + Madpc = 0
 R("R_ARBTNabcpp"): + Mpic + Mhc - Mh2oc - Matpc - Marbtnfe3p
 + Marbtnfe3c + Madpc = 0
 R("R_ARGabcpp"): + Mpic + Mhc - Mh2oc - Matpc - MargLp + MargLc + Madpc
 = 0
 R("R_ASNabcpp"): + Mpic + Mhc - Mh2oc - Matpc - MasnLp + MasnLc + Madpc
 = 0
 R("R_ASO3t8pp"): + Mpic + Mhc - Mh2oc - Matpc + Maso3p - Maso3c + Madpc
 = 0
 R("R_ASPabcpp"): + Mpic + Mhc - Mh2oc - Matpc - MaspLp + MaspLc + Madpc
 = 0
 R("R_ATPM"): + Mpic + Mhc - Mh2oc - Matpc + Madpc = 0
 R("R_ATPS4rpp"): - Mpic - 4 Mhp + 3 Mhc + Mh2oc + Matpc - Madpc = 0
 R("R_ATPS4rpp_Rev"): + Mpic + 4 Mhp - 3 Mhc - Mh2oc - Matpc + Madpc = 0
 R("R_BUTSO3abcpp"): + Mpic + Mhc - Mh2oc - Mbutso3p + Mbutso3c - Matpc
 + Madpc = 0
 R("R_CBIuabcpp"): + Mpic + Mhc - Mh2oc - Mcbip + Mcbic - Matpc + Madpc
 = 0
 R("R_CBL1abcpp"): + Mpic + Mhc - Mh2oc - Mcbl1p + Mcbl1c - Matpc
 + Madpc = 0
 R("R_CBMKr"): - Mnh4c + 2 Mhc - Mco2c + Mcbpc - Matpc + Madpc = 0
 R("R_CBMKr_Rev"): + Mnh4c - 2 Mhc + Mco2c - Mcbpc + Matpc - Madpc = 0
 R("R_CBPS"): + Mpic - Mhco3c + 2 Mhc - Mh2oc + MgluLc - MglnLc + Mcbpc
 - 2 Matpc + 2 Madpc = 0
 R("R_CD2abcpp"): + Mpic + Mhc - Mh2oc + Mcd2p - Mcd2c - Matpc + Madpc
 = 0
 R("R_CGLYabcpp"): + Mpic + Mhc - Mh2oc - Mcglyp + Mcglyc - Matpc
 + Madpc = 0
 R("R_CHLabcpp"): + Mpic + Mhc - Mh2oc - Mcholp + Mcholc - Matpc + Madpc
 = 0
 R("R_COBALT2abcpp"): + Mpic + Mhc - Mh2oc + Mcobalt2p - Mcobalt2c
 - Matpc + Madpc = 0
 R("R_COLIPAabcpp"): + Mpic + Mhc - Mh2oc + Mcolipap - Mcolipac - Matpc
 + Madpc = 0
 R("R_CPGNabcpp"): + Mpic + Mhc - Mh2oc - Mcpgnp + Mcpgnc - Matpc
 + Madpc = 0
 R("R_CRNabcpp"): + Mpic + Mhc - Mh2oc - Mcrnp + Mcrnc - Matpc + Madpc
 = 0
 R("R_CRNCAL2"): + Mpic + Mcrncoac - Mcrnc - Mcoac - Matpc + Madpc = 0
 R("R_CRNDabcpp"): + Mpic + Mhc - Mh2oc - McrnDp + McrnDc - Matpc
 + Madpc = 0
 R("R_CRNDCAL2"): + Mpic + McrnDcoac - McrnDc - Mcoac - Matpc + Madpc
 = 0
 R("R_CTBTabcpp"): + Mpic + Mhc - Mh2oc - Mctbtp + Mctbtc - Matpc
 + Madpc = 0
 R("R_CTBTCAL2"): + Mpic + Mctbtcoac - Mctbtc - Mcoac - Matpc + Madpc
 = 0
 R("R_CTPS2"): - Mutpc + Mpic + 2 Mhc - Mh2oc + MgluLc - MglnLc + Mctpc
 - Matpc + Madpc = 0
 R("R_CU1abcpp"): + Mpic + Mhc - Mh2oc + Mcup - Mcuc - Matpc + Madpc = 0
 R("R_CU2abcpp"): + Mpic + Mhc - Mh2oc + Mcu2p - Mcu2c - Matpc + Madpc
 = 0
 R("R_CYSabc2pp"): + Mpic + Mhc - Mh2oc + McysLp - McysLc - Matpc
 + Madpc = 0
 R("R_CYSabcpp"): + Mpic + Mhc - Mh2oc - McysLp + McysLc - Matpc + Madpc
 = 0
 R("R_CYSDabcpp"): + Mpic + Mhc - Mh2oc - McysDp + McysDc - Matpc
 + Madpc = 0
 R("R_CYTK1"): - Mcmpc + Mcdpc - Matpc + Madpc = 0
 R("R_CYTK1_Rev"): + Mcmpc - Mcdpc + Matpc - Madpc = 0
 R("R_CYTK2"): - Mdcmpc + Mdcdpc - Matpc + Madpc = 0
 R("R_CYTK2_Rev"): + Mdcmpc - Mdcdpc + Matpc - Madpc = 0
 R("R_DADK"): - Mdampc + Mdadpc - Matpc + Madpc = 0
 R("R_DADK_Rev"): + Mdampc - Mdadpc + Matpc - Madpc = 0
 R("R_DBTS"): + Mpic + 3 Mhc + Mdtbtc - Mdannc - Mco2c - Matpc + Madpc
 = 0
 R("R_DGK1"): - Mdgmpc + Mdgdpc - Matpc + Madpc = 0
 R("R_DGK1_Rev"): + Mdgmpc - Mdgdpc + Matpc - Madpc = 0
 R("R_DHFS"): + Mpic + Mhc - MgluLc - Mdhptc + Mdhfc - Matpc + Madpc = 0
 R("R_DPCOAK"): + Mhc - Mdpcoac + Mcoac - Matpc + Madpc = 0
 R("R_DTMPK"): - Mdtmpc + Mdtdpc - Matpc + Madpc = 0
 R("R_DTMPK_Rev"): + Mdtmpc - Mdtdpc + Matpc - Madpc = 0
 R("R_DURIK1"): + Mhc - Mduric + Mdumpc - Matpc + Madpc = 0
 R("R_DXYLK"): + Mhc - Mdxylc + Mdxyl5pc - Matpc + Madpc = 0
 R("R_ETHSO3abcpp"): + Mpic + Mhc - Mh2oc - Methso3p + Methso3c - Matpc
 + Madpc = 0
 R("R_FCLK"): + Mhc - MfclLc + Mfc1pc - Matpc + Madpc = 0
 R("R_FE2abcpp"): + Mpic + Mhc - Mh2oc - Mfe2p + Mfe2c - Matpc + Madpc
 = 0
 R("R_FE3abcpp"): + Mpic + Mhc - Mh2oc - Mfe3p + Mfe3c - Matpc + Madpc
 = 0
 R("R_FE3DCITabcpp"): + Mpic + Mhc - Mh2oc - Mfe3dcitp + Mfe3c + 2 Mcitc
 - Matpc + Madpc = 0
 R("R_FE3DHBZSabcpp"): + Mpic + Mhc - Mh2oc - Mfe3dhbzsp + Mfe3dhbzsc
 - Matpc + Madpc = 0
 R("R_FE3HOXabcpp"): + Mpic + Mhc - Mh2oc - Mfe3hoxp + Mfe3hoxc - Matpc
 + Madpc = 0
 R("R_FECRMabcpp"): + Mpic + Mhc - Mh2oc - Mfecrmp + Mfecrmc - Matpc
 + Madpc = 0
 R("R_FEENTERabcpp"): + Mpic + Mhc - Mh2oc - Mfeenterp + Mfeenterc
 - Matpc + Madpc = 0
 R("R_FEOXAMabcpp"): + Mpic + Mhc - Mh2oc - Mfeoxamp + Mfeoxamc - Matpc
 + Madpc = 0
 R("R_FRUK"): + Mhc + Mfdpc - Mf1pc - Matpc + Madpc = 0
 R("R_FRULYSK"): + Mhc + Mfrulyspc - Mfrulysc - Matpc + Madpc = 0
 R("R_G3PCabcpp"): + Mpic + Mhc - Mh2oc - Mg3pcp + Mg3pcc - Matpc
 + Madpc = 0
 R("R_G3PEabcpp"): + Mpic + Mhc - Mh2oc - Mg3pep + Mg3pec - Matpc
 + Madpc = 0
 R("R_G3PGabcpp"): + Mpic + Mhc - Mh2oc - Mg3pgp + Mg3pgc - Matpc
 + Madpc = 0
 R("R_G3PIabcpp"): + Mpic + Mhc - Mh2oc - Mg3pip + Mg3pic - Matpc
 + Madpc = 0
 R("R_G3PSabcpp"): + Mpic + Mhc - Mh2oc - Mg3psp + Mg3psc - Matpc
 + Madpc = 0
 R("R_GALabcpp"): + Mpic + Mhc - Mh2oc - Mgalp + Mgalc - Matpc + Madpc
 = 0
 R("R_GALKr"): + Mhc - Mgalc + Mgal1pc - Matpc + Madpc = 0
 R("R_GALKr_Rev"): - Mhc + Mgalc - Mgal1pc + Matpc - Madpc = 0
 R("R_GART"): + Mpic + Mhc - Mgarc - Mforc + Mfgamc - Matpc + Madpc = 0
 R("R_GGPTRCS"): - Mptrcc + Mpic + Mhc - MgluLc + Mggptrcc - Matpc
 + Madpc = 0
 R("R_GK1"): - Mgmpc + Mgdpc - Matpc + Madpc = 0
 R("R_GK1_Rev"): + Mgmpc - Mgdpc + Matpc - Madpc = 0
 R("R_GLCabcpp"): + Mpic + Mhc - Mh2oc - MglcDp + MglcDc - Matpc + Madpc
 = 0
 R("R_GLCS1"): + Mhc + Mglycogenc - Madpglcc + Madpc = 0
 R("R_GLNabcpp"): + Mpic + Mhc - Mh2oc - MglnLp + MglnLc - Matpc + Madpc
 = 0
 R("R_GLNS"): + Mpic - Mnh4c + Mhc - MgluLc + MglnLc - Matpc + Madpc = 0
 R("R_GLU5K"): - MgluLc + Mglu5pc - Matpc + Madpc = 0
 R("R_GLUabcpp"): + Mpic + Mhc - Mh2oc - MgluLp + MgluLc - Matpc + Madpc
 = 0
 R("R_GLUCYS"): + Mpic + Mhc - MgluLc + Mglucysc - McysLc - Matpc
 + Madpc = 0
 R("R_GLYBabcpp"): + Mpic + Mhc - Mh2oc - Mglybp + Mglybc - Matpc
 + Madpc = 0
 R("R_GLYC3Pabcpp"): + Mpic + Mhc - Mh2oc - Mglyc3pp + Mglyc3pc - Matpc
 + Madpc = 0
 R("R_GLYK"): + Mhc - Mglycc + Mglyc3pc - Matpc + Madpc = 0
 R("R_GMHEPK"): + Mhc - Mgmhep7pc + Mgmhep17bpc - Matpc + Madpc = 0
 R("R_GSNK"): + Mhc - Mgsnc + Mgmpc - Matpc + Madpc = 0
 R("R_GSPMDS"): - Mspmdc + Mpic + Mhc + Mgtspmdc - Mgthrdc - Matpc
 + Madpc = 0
 R("R_GTHRDabc2pp"): + Mpic + Mhc - Mh2oc + Mgthrdp - Mgthrdc - Matpc
 + Madpc = 0
 R("R_GTHRDabcpp"): + Mpic + Mhc - Mh2oc - Mgthrdp + Mgthrdc - Matpc
 + Madpc = 0
 R("R_GTHS"): + Mpic + Mhc + Mgthrdc - Mglyc - Mglucysc - Matpc + Madpc
 = 0
 R("R_HEPK1"): + Mphhlipac - Mhhlipac + Mhc - Matpc + Madpc = 0
 R("R_HEPK2"): + Mphphhlipac - Mhphhlipac + Mhc - Matpc + Madpc = 0
 R("R_HEPT1"): - Mlipac + Mhlipac + Mhc - MadphepLDc + Madpc = 0
 R("R_HEPT2"): - Mhlipac + Mhhlipac + Mhc - MadphepLDc + Madpc = 0
 R("R_HEPT3"): - Mphhlipac + Mhphhlipac + Mhc - MadphepLDc + Madpc = 0
 R("R_HEPT4"): + Mhc - Mgggagicolipac + Mcolipac - MadphepLDc + Madpc
 = 0
 R("R_HEX1"): + Mhc - MglcDc + Mg6pc - Matpc + Madpc = 0
 R("R_HEX4"): - Mmanc + Mman6pc + Mhc - Matpc + Madpc = 0
 R("R_HEX7"): + Mhc - Mfruc + Mf6pc - Matpc + Madpc = 0
 R("R_HG2abcpp"): + Mpic + Mhg2p - Mhg2c + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_HISabcpp"): + Mpic - MhisLp + MhisLc + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_HSK"): + Mphomc - MhomLc + Mhc - Matpc + Madpc = 0
 R("R_ILEabcpp"): + Mpic - MileLp + MileLc + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_INSK"): - Minsc + Mimpc + Mhc - Matpc + Madpc = 0
 R("R_ISETACabcpp"): + Mpic - Misetacp + Misetacc + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_K2L4Aabcpp"): + Mpic + Mkdo2lipid4p - Mkdo2lipid4c + Mhc - Mh2oc
 - Matpc + Madpc = 0
 R("R_Kabcpp"): + Mpic - Mkp + Mkc + Mhc - Mh2oc - Matpc + Madpc = 0
 R("R_LEUabcpp"): + Mpic - MleuLp + MleuLc + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_LIPAabcpp"): + Mpic + Mlipap - Mlipac + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_LIPACabcpp"): + Mpic + Mlipacoldp - Mlipacoldc + Mhc - Mh2oc
 - Matpc + Madpc = 0
 R("R_LYSabcpp"): + Mpic - MlysLp + MlysLc + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_MALTabcpp"): + Mpic - Mmaltp + Mmaltc + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_MALTHXabcpp"): + Mpic - Mmalthxp + Mmalthxc + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_MALTPTabcpp"): + Mpic - Mmaltptp + Mmaltptc + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_MALTTRabcpp"): + Mpic - Mmalttrp + Mmalttrc + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_MALTTTRabcpp"): + Mpic - Mmaltttrp + Mmaltttrc + Mhc - Mh2oc
 - Matpc + Madpc = 0
 R("R_METabcpp"): + Mpic - MmetLp + MmetLc + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_METDabcpp"): + Mpic - MmetDp + MmetDc + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_METSOX1abcpp"): + Mpic - MmetsoxSLp + MmetsoxSLc + Mhc - Mh2oc
 - Matpc + Madpc = 0
 R("R_METSOX2abcpp"): + Mpic - MmetsoxRLp + MmetsoxRLc + Mhc - Mh2oc
 - Matpc + Madpc = 0
 R("R_MG2uabcpp"): + Mpic - Mmg2p + Mmg2c + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_MOBDabcpp"): + Mpic - Mmobdp + Mmobdc + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_MSO3abcpp"): + Mpic - Mmso3p + Mmso3c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_NADK"): + Mnadpc - Mnadc + Mhc - Matpc + Madpc = 0
 R("R_NAMNPP"): - Mprppc + Mppic + Mpic + Mnicrntc - Mnacc - Mh2oc
 - Matpc + Madpc = 0
 R("R_NDPK1"): + Mgtpc - Mgdpc - Matpc + Madpc = 0
 R("R_NDPK1_Rev"): - Mgtpc + Mgdpc + Matpc - Madpc = 0
 R("R_NDPK2"): + Mutpc - Mudpc - Matpc + Madpc = 0
 R("R_NDPK2_Rev"): - Mutpc + Mudpc + Matpc - Madpc = 0
 R("R_NDPK3"): + Mctpc - Mcdpc - Matpc + Madpc = 0
 R("R_NDPK3_Rev"): - Mctpc + Mcdpc + Matpc - Madpc = 0
 R("R_NDPK4"): + Mdttpc - Mdtdpc - Matpc + Madpc = 0
 R("R_NDPK4_Rev"): - Mdttpc + Mdtdpc + Matpc - Madpc = 0
 R("R_NDPK5"): + Mdgtpc - Mdgdpc - Matpc + Madpc = 0
 R("R_NDPK5_Rev"): - Mdgtpc + Mdgdpc + Matpc - Madpc = 0
 R("R_NDPK6"): + Mdutpc - Mdudpc - Matpc + Madpc = 0
 R("R_NDPK6_Rev"): - Mdutpc + Mdudpc + Matpc - Madpc = 0
 R("R_NDPK7"): + Mdctpc - Mdcdpc - Matpc + Madpc = 0
 R("R_NDPK7_Rev"): - Mdctpc + Mdcdpc + Matpc - Madpc = 0
 R("R_NDPK8"): + Mdatpc - Mdadpc - Matpc + Madpc = 0
 R("R_NDPK8_Rev"): - Mdatpc + Mdadpc + Matpc - Madpc = 0
 R("R_NI2abcpp"): + Mpic + Mni2p - Mni2c + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_NI2uabcpp"): + Mpic - Mni2p + Mni2c + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_NTP1"): + Mpic + Mhc - Mh2oc - Matpc + Madpc = 0
 R("R_ORNabcpp"): + Mpic - Mornp + Mornc + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_PA120abcpp"): + Mpic + Mpa120p - Mpa120c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PA140abcpp"): + Mpic + Mpa140p - Mpa140c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PA141abcpp"): + Mpic + Mpa141p - Mpa141c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PA160abcpp"): + Mpic + Mpa160p - Mpa160c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PA161abcpp"): + Mpic + Mpa161p - Mpa161c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PA180abcpp"): + Mpic + Mpa180p - Mpa180c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PA181abcpp"): + Mpic + Mpa181p - Mpa181c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PE120abcpp"): + Mpic + Mpe120p - Mpe120c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PE140abcpp"): + Mpic + Mpe140p - Mpe140c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PE141abcpp"): + Mpic + Mpe141p - Mpe141c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PE160abcpp"): + Mpic + Mpe160p - Mpe160c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PE161abcpp"): + Mpic + Mpe161p - Mpe161c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PE180abcpp"): + Mpic + Mpe180p - Mpe180c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PE181abcpp"): + Mpic + Mpe181p - Mpe181c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PFK"): + Mhc + Mfdpc - Mf6pc - Matpc + Madpc = 0
 R("R_PFK_2"): + MtagdpDc - Mtag6pDc + Mhc - Matpc + Madpc = 0
 R("R_PG120abcpp"): + Mpic + Mpg120p - Mpg120c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PG140abcpp"): + Mpic + Mpg140p - Mpg140c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PG141abcpp"): + Mpic + Mpg141p - Mpg141c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PG160abcpp"): + Mpic + Mpg160p - Mpg160c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PG161abcpp"): + Mpic + Mpg161p - Mpg161c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PG180abcpp"): + Mpic + Mpg180p - Mpg180c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PG181abcpp"): + Mpic + Mpg181p - Mpg181c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PGP120abcpp"): + Mpic + Mpgp120p - Mpgp120c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PGP140abcpp"): + Mpic + Mpgp140p - Mpgp140c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PGP141abcpp"): + Mpic + Mpgp141p - Mpgp141c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PGP160abcpp"): + Mpic + Mpgp160p - Mpgp160c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PGP161abcpp"): + Mpic + Mpgp161p - Mpgp161c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PGP180abcpp"): + Mpic + Mpgp180p - Mpgp180c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PGP181abcpp"): + Mpic + Mpgp181p - Mpgp181c + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PHEMEabcpp"): + Mpic + Mphemep - Mphemec + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PIuabcpp"): - Mpip + 2 Mpic + Mhc - Mh2oc - Matpc + Madpc = 0
 R("R_PPAKr"): - Mppapc + Mppac + Matpc - Madpc = 0
 R("R_PPAKr_Rev"): + Mppapc - Mppac - Matpc + Madpc = 0
 R("R_PPCK"): + Mpepc - Moaac + Mco2c - Matpc + Madpc = 0
 R("R_PPK2r"): + Mpppic - Mppic - Matpc + Madpc = 0
 R("R_PPK2r_Rev"): - Mpppic + Mppic + Matpc - Madpc = 0
 R("R_PPKr"): + Mppic - Mpic - Matpc + Madpc = 0
 R("R_PPKr_Rev"): - Mppic + Mpic + Matpc - Madpc = 0
 R("R_PRAGSr"): - Mpramc + Mpic + Mhc - Mglyc + Mgarc - Matpc + Madpc
 = 0
 R("R_PRAGSr_Rev"): + Mpramc - Mpic - Mhc + Mglyc - Mgarc + Matpc
 - Madpc = 0
 R("R_PRAIS"): + Mpic + 2 Mhc - Mfpramc - Matpc + Mairc + Madpc = 0
 R("R_PRFGS"): + Mpic + Mhc - Mh2oc + MgluLc - MglnLc + Mfpramc - Mfgamc
 - Matpc + Madpc = 0
 R("R_PROabcpp"): - MproLp + MproLc + Mpic + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_PROGLYabcpp"): - Mproglyp + Mproglyc + Mpic + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PTRCabcpp"): - Mptrcp + Mptrcc + Mpic + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_PYDAMK"): - Mpydamc + Mpyam5pc + Mhc - Matpc + Madpc = 0
 R("R_PYDXK"): - Mpydxc + Mpydx5pc + Mhc - Matpc + Madpc = 0
 R("R_PYDXNK"): - Mpydxnc + Mpdx5pc + Mhc - Matpc + Madpc = 0
 R("R_PYK"): + Mpyrc - Mpepc - Mhc + Matpc - Madpc = 0
 R("R_R15BPK"): - Mr15bpc + Mprppc - Matpc + Madpc = 0
 R("R_R1PK"): - Mr1pc + Mr15bpc + Mhc - Matpc + Madpc = 0
 R("R_RBFK"): - Mribflvc + Mhc + Mfmnc - Matpc + Madpc = 0
 R("R_RBK"): - MribDc + Mr5pc + Mhc - Matpc + Madpc = 0
 R("R_RBK_L1"): + Mru5pLc - MrblLc + Mhc - Matpc + Madpc = 0
 R("R_RIBabcpp"): - MribDp + MribDc + Mpic + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_RMK"): - Mrmlc + Mrml1pc + Mhc - Matpc + Madpc = 0
 R("R_RNDR1"): - Mtrdrdc + Mtrdoxc + Mh2oc + Mdadpc - Madpc = 0
 R("R_RNDR1b"): + Mh2oc - Mgrxrdc + Mgrxoxc + Mdadpc - Madpc = 0
 R("R_SHKK"): - Mskmc + Mskm5pc + Mhc - Matpc + Madpc = 0
 R("R_SPMDabcpp"): - Mspmdp + Mspmdc + Mpic + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_SUCOAS"): + Msuccoac - Msuccc + Mpic - Mcoac - Matpc + Madpc = 0
 R("R_SUCOAS_Rev"): - Msuccoac + Msuccc - Mpic + Mcoac + Matpc - Madpc
 = 0
 R("R_SULabcpp"): - Mso4p + Mso4c + Mpic + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_SULFACabcpp"): - Msulfacp + Msulfacc + Mpic + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_TAURabcpp"): - Mtaurp + Mtaurc + Mpic + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_TDSK"): - MlipidAdsc + MlipidAc + Mhc - Matpc + Madpc = 0
 R("R_THMabcpp"): - Mthmp + Mthmc + Mpic + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_THRabcpp"): - MthrLp + MthrLc + Mpic + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_TMDK1"): - Mthymdc + Mhc + Mdtmpc - Matpc + Madpc = 0
 R("R_TMK"): + Mthmmpc - Mthmc + Mhc - Matpc + Madpc = 0
 R("R_TMPK"): + Mthmppc - Mthmmpc - Matpc + Madpc = 0
 R("R_TSULabcpp"): - Mtsulp + Mtsulc + Mpic + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_TUNGSabcpp"): - Mtungsp + Mtungsc + Mpic + Mhc - Mh2oc - Matpc
 + Madpc = 0
 R("R_UAMAGS"): + Muamagc - Muamac + Mpic + Mhc - MgluDc - Matpc + Madpc
 = 0
 R("R_UAMAS"): - Muamrc + Muamac + Mpic + Mhc - Matpc - MalaLc + Madpc
 = 0
 R("R_UGMDDS"): - Mugmdc + Mugmdac + Mpic + Mhc - Matpc - Malaalac
 + Madpc = 0
 R("R_UM3PL"): + Mugmdc - Muamrc + Mpic - MLalaDgluMdapc + Mhc - Matpc
 + Madpc = 0
 R("R_UM4PL"): + Mum4pc - Muamrc + Mpic - MLalaDgluMdapDalac + Mhc
 - Matpc + Madpc = 0
 R("R_UMPK"): - Mumpc + Mudpc - Matpc + Madpc = 0
 R("R_UMPK_Rev"): + Mumpc - Mudpc + Matpc - Madpc = 0
 R("R_URIDK2r"): - Mdumpc + Mdudpc - Matpc + Madpc = 0
 R("R_URIDK2r_Rev"): + Mdumpc - Mdudpc + Matpc - Madpc = 0
 R("R_VALabcpp"): - MvalLp + MvalLc + Mpic + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_XYLabcpp"): - MxylDp + MxylDc + Mpic + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_XYLK"): - MxyluDc + Mxu5pDc + Mhc - Matpc + Madpc = 0
 R("R_XYLK2"): - MxyluLc + Mxu5pLc + Mhc - Matpc + Madpc = 0
 R("R_ZN2abcpp"): + Mzn2p - Mzn2c + Mpic + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_ZNabcpp"): - Mzn2p + Mzn2c + Mpic + Mhc - Mh2oc - Matpc + Madpc
 = 0
 R("R_GLGC"): + Mppic - Mhc - Mg1pc - Matpc + Madpglcc = 0
 R("R_AGMHE"): + MadphepLDc - MadphepDDc = 0
 R("R_GMHEPAT"): + Mppic - Mhc - Mgmhep1pc - Matpc + MadphepDDc = 0
 R("R_ADPRDP"): + Mr5pc + 2 Mhc - Mh2oc + Mampc - Madpribc = 0
 R("R_NADN"): + Mncamc - Mnadc + Mhc - Mh2oc + Madpribc = 0
 R("Magb_Transp"): 0 M10fthfc = 0
 R("R_EX_ag_e_"): - Mage = 0
 R("R_EX_ag_e__Rev"): + Mage = 0
 R("R_AGt3"): - Mhe + Mhc + Mage - Magc = 0
 R("Magmb_Transp"): 0 M10fthfc = 0
 R("R_EX_agm_e_"): - Magme = 0
 R("R_EX_agm_e__Rev"): + Magme = 0
 R("R_AGMT"): + Mureac + Mptrcc - Mh2oc - Magmc = 0
 R("R_ARGAGMt7pp"): - MargLp + MargLc + Magmp - Magmc = 0
 R("R_ARGAGMt7pp_Rev"): + MargLp - MargLc - Magmp + Magmc = 0
 R("R_ARGDC"): - Mhc + Mco2c - MargLc + Magmc = 0
 R("R_AGMtex"): + Magmp - Magme = 0
 R("R_AGMtex_Rev"): - Magmp + Magme = 0
 R("R_ARGDCpp"): - Mhp + Mco2p - MargLp + Magmp = 0
 R("R_CFAS160E"): - Mpe161c + 2 Mhc + Mcpe160c - 2 Mametc + 2 Mahcysc
 = 0
 R("R_CFAS160G"): - Mpg161c + 2 Mhc + Mcpg160c - 2 Mametc + 2 Mahcysc
 = 0
 R("R_CFAS180E"): - Mpe181c + 2 Mhc + Mcpe180c - 2 Mametc + 2 Mahcysc
 = 0
 R("R_CFAS180G"): - Mpg181c + 2 Mhc + Mcpg180c - 2 Mametc + 2 Mahcysc
 = 0
 R("R_HCYSMT"): + MmetLc - MhcysLc + Mhc - Mametc + Mahcysc = 0
 R("R_UPP3MT"): - Muppg3c + Mhc + Mdsclc - 2 Mametc + 2 Mahcysc = 0
 R("R_DHPTPE"): + Mdhmptpc - Mahdtc = 0
 R("R_DHPTPE_Rev"): - Mdhmptpc + Mahdtc = 0
 R("R_DNTPPA"): + Mppic + Mhc - Mh2oc + Mdhpmpc - Mahdtc = 0
 R("R_GTPCI"): + Mhc - Mh2oc - Mgtpc + Mforc + Mahdtc = 0
 R("R_IG3PS"): - Mprlpc + Mhc + MgluLc - MglnLc + Meig3pc + Maicarc = 0
 R("Makgb_Transp"): 0 M10fthfc = 0
 R("R_EX_akg_e_"): - Makge = 0
 R("R_EX_akg_e__Rev"): + Makge = 0
 R("R_AKGDH"): + Msuccoac + Mnadhc - Mnadc - Mcoac + Mco2c - Makgc = 0
 R("R_AKGt2rpp"): - Mhp + Mhc - Makgp + Makgc = 0
 R("R_AKGt2rpp_Rev"): + Mhp - Mhc + Makgp - Makgc = 0
 R("R_ALATA_L"): + Mpyrc + MgluLc - MalaLc - Makgc = 0
 R("R_ALATA_L_Rev"): - Mpyrc - MgluLc + MalaLc + Makgc = 0
 R("R_ASPTA"): + Moaac + MgluLc - MaspLc - Makgc = 0
 R("R_ASPTA_Rev"): - Moaac - MgluLc + MaspLc + Makgc = 0
 R("R_GLUDy"): + Mnh4c + Mnadphc - Mnadpc + Mhc - Mh2oc - MgluLc + Makgc
 = 0
 R("R_GLUDy_Rev"): - Mnh4c - Mnadphc + Mnadpc - Mhc + Mh2oc + MgluLc
 - Makgc = 0
 R("R_GLUSy"): - Mnadphc + Mnadpc - Mhc + 2 MgluLc - MglnLc - Makgc = 0
 R("R_HSTPT"): - Mimacpc + Mhispc - MgluLc + Makgc = 0
 R("R_ICDHyr"): + Mnadphc - Mnadpc - Micitc + Mco2c + Makgc = 0
 R("R_ICDHyr_Rev"): - Mnadphc + Mnadpc + Micitc - Mco2c - Makgc = 0
 R("R_OHPBAT"): + Mphthrc - Mohpbc - MgluLc + Makgc = 0
 R("R_OHPBAT_Rev"): - Mphthrc + Mohpbc + MgluLc - Makgc = 0
 R("R_OXGDC2"): - Mthmppc + Mssaltppc - Mhc + Mco2c - Makgc = 0
 R("R_PHETA1"): + Mphpyrc - MpheLc + MgluLc - Makgc = 0
 R("R_PHETA1_Rev"): - Mphpyrc + MpheLc - MgluLc + Makgc = 0
 R("R_SDPTA"): + Msl2a6oc - Msl26dac + MgluLc - Makgc = 0
 R("R_SDPTA_Rev"): - Msl2a6oc + Msl26dac - MgluLc + Makgc = 0
 R("R_SOTA"): - Msucornc + Msucgsac + MgluLc - Makgc = 0
 R("R_TDPAGTA"): - MgluLc - Mdtdp4d6dgc + Mdtdp4addgc + Makgc = 0
 R("R_UDPKAAT"): - MudpLa4oc + MudpLa4nc - MgluLc + Makgc = 0
 R("R_UDPKAAT_Rev"): + MudpLa4oc - MudpLa4nc + MgluLc - Makgc = 0
 R("R_AKGtex"): + Makgp - Makge = 0
 R("R_AKGtex_Rev"): - Makgp + Makge = 0
 R("Malaalab_Transp"): 0 M10fthfc = 0
 R("R_EX_alaala_e_"): - Malaalae = 0
 R("R_EX_alaala_e__Rev"): + Malaalae = 0
 R("R_ALAALAD"): - Mh2oc + 2 MalaDc - Malaalac = 0
 R("R_ALAALAtex"): + Malaalap - Malaalae = 0
 R("R_ALAALAtex_Rev"): - Malaalap + Malaalae = 0
 R("R_MCTP1Bpp"): + Mmurein5px3pp - Mmurein5p5pp + Malaalap = 0
 R("R_MLDCP1App"): - Mmurein5px4pp + Mmurein3px4pp - Mh2op + Malaalap
 = 0
 R("R_MLDCP2App"): - Mmurein5p5pp + Mmurein5p3pp - Mh2op + Malaalap = 0
 R("R_MLDCP3App"): - Mmurein5px3pp + Mmurein3px3pp - Mh2op + Malaalap
 = 0
 R("MalaBb_Transp"): 0 M10fthfc = 0
 R("R_EX_ala_B_e_"): - MalaBe = 0
 R("R_EX_ala_B_e__Rev"): + MalaBe = 0
 R("R_ASP1DC"): - Mhc + Mco2c - MaspLc + MalaBc = 0
 R("R_BALAt2pp"): - Mhp + Mhc - MalaBp + MalaBc = 0
 R("R_PANTS"): + Mppic + MpntoRc - MpantRc + Mhc - Matpc + Mampc
 - MalaBc = 0
 R("R_BALAtex"): + MalaBp - MalaBe = 0
 R("R_BALAtex_Rev"): - MalaBp + MalaBe = 0
 R("R_ACLS"): - 2 Mpyrc - Mhc + Mco2c + MalacSc = 0
 R("MalaDb_Transp"): 0 M10fthfc = 0
 R("R_EX_ala_D_e_"): - MalaDe = 0
 R("R_EX_ala_D_e__Rev"): + MalaDe = 0
 R("R_4PCP"): - MLalaDgluMdapDalac + MLalaDgluMdapc - Mh2oc + MalaDc = 0
 R("R_AGM4PCP"): - Mh2oc - Manhgm4pc + Manhgm3pc + MalaDc = 0
 R("R_ALAR"): - MalaLc + MalaDc = 0
 R("R_ALAR_Rev"): + MalaLc - MalaDc = 0
 R("R_ALATA_D2"): + Mpyrc - Mpydx5pc + Mpyam5pc - MalaDc = 0
 R("R_AM4PCP"): - Mh2oc - Manhm4pc + Manhm3pc + MalaDc = 0
 R("R_DAAD"): + Mpyrc + Mnh4c - Mh2oc + Mfadh2c - Mfadc - MalaDc = 0
 R("R_DALAt2pp"): - Mhp + Mhc - MalaDp + MalaDc = 0
 R("R_UM4PCP"): - Mum4pc + Mugmdc - Mh2oc + MalaDc = 0
 R("R_DALAtex"): + MalaDp - MalaDe = 0
 R("R_DALAtex_Rev"): - MalaDp + MalaDe = 0
 R("R_4PCPpp"): + MLalaDgluMdapp - MLalaDgluMdapDalap - Mh2op + MalaDp
 = 0
 R("R_AGM4PCPpp"): - Mh2op - Manhgm4pp + Manhgm3pp + MalaDp = 0
 R("R_MCTP1App"): + Mmurein5px4pp - Mmurein5p5pp + MalaDp = 0
 R("R_MCTP2App"): + Mmurein5px4px4pp - Mmurein5p5p5pp + 2 MalaDp = 0
 R("R_MDDCP1pp"): - Mmurein5px4pp + Mmurein4px4pp - Mh2op + MalaDp = 0
 R("R_MDDCP2pp"): - Mmurein5px4px4pp + Mmurein4px4px4pp - Mh2op + MalaDp
 = 0
 R("R_MDDCP3pp"): - Mmurein5p5pp + Mmurein5p4pp - Mh2op + MalaDp = 0
 R("R_MDDCP4pp"): - Mmurein5p4pp + Mmurein4p4pp - Mh2op + MalaDp = 0
 R("R_MDDCP5pp"): - Mmurein5p3pp + Mmurein4p3pp - Mh2op + MalaDp = 0
 R("R_MLDCP1Bpp"): - Mmurein4p4pp + Mmurein4p3pp - Mh2op + MalaDp = 0
 R("R_MLDCP2Bpp"): - Mmurein4p3pp + Mmurein3p3pp - Mh2op + MalaDp = 0
 R("MalaLb_Transp"): 0 M10fthfc = 0
 R("R_EX_ala_L_e_"): - MalaLe = 0
 R("R_EX_ala_L_e__Rev"): + MalaLe = 0
 R("R_ALAt2pp"): - Mhp + Mhc - MalaLp + MalaLc = 0
 R("R_ALAt4pp"): - Mna1p + Mna1c - MalaLp + MalaLc = 0
 R("R_ALATA_L2"): + Mpyrc - Mpydx5pc + Mpyam5pc - MalaLc = 0
 R("R_ALATRS"): - Mtrnaalac + Mppic - Matpc + Mampc + Malatrnac - MalaLc
 = 0
 R("R_LALGP"): - MLalaLgluc - Mh2oc + MgluLc + MalaLc = 0
 R("R_ALAtex"): + MalaLp - MalaLe = 0
 R("R_ALAtex_Rev"): - MalaLp + MalaLe = 0
 R("R_ALLPI"): + Mallul6pc - Mall6pc = 0
 R("R_ALLPI_Rev"): - Mallul6pc + Mall6pc = 0
 R("MallDb_Transp"): 0 M10fthfc = 0
 R("R_EX_all_D_e_"): - MallDe = 0
 R("R_EX_all_D_e__Rev"): + MallDe = 0
 R("R_ALLtex"): + MallDp - MallDe = 0
 R("R_ALLtex_Rev"): - MallDp + MallDe = 0
 R("Malltnb_Transp"): 0 M10fthfc = 0
 R("R_EX_alltn_e_"): - Malltne = 0
 R("R_EX_alltn_e__Rev"): + Malltne = 0
 R("R_ALLTN"): + Mhc - Mh2oc + Mallttc - Malltnc = 0
 R("R_ALLTNt2rpp"): - Mhp + Mhc - Malltnp + Malltnc = 0
 R("R_ALLTNt2rpp_Rev"): + Mhp - Mhc + Malltnp - Malltnc = 0
 R("R_URIC"): - Muratec - Mo2c - 2 Mh2oc + Mh2o2c + Mco2c + Malltnc = 0
 R("R_ALLTNtex"): + Malltnp - Malltne = 0
 R("R_ALLTNtex_Rev"): - Malltnp + Malltne = 0
 R("R_ALLTAMH"): + Murdglycc + 2 Mnh4c - 2 Mhc - 2 Mh2oc + Mco2c
 - Mallttc = 0
 R("R_ALLULPE"): + Mf6pc - Mallul6pc = 0
 R("R_ALLULPE_Rev"): - Mf6pc + Mallul6pc = 0
 R("R_TAGURr"): + Mtagurc + Mnadhc - Mnadc + Mhc - Maltrnc = 0
 R("R_TAGURr_Rev"): - Mtagurc - Mnadhc + Mnadc - Mhc + Maltrnc = 0
 R("R_ADMDC"): - Mhc + Mco2c - Mametc + Mametamc = 0
 R("R_BTS4"): - Msc + MmetLc + Mhc - Mdtbtc + Mdad5c + Mbtnc - Mametc
 = 0
 R("R_CPPPGO2"): + Mpppg9c + 2 MmetLc + 2 Mdad5c - Mcpppg3c + 2 Mco2c
 - 2 Mametc = 0
 R("R_METAT"): + Mppic + Mpic - MmetLc - Mh2oc - Matpc + Mametc = 0
 R("Mampb_Transp"): 0 M10fthfc = 0
 R("R_EX_amp_e_"): - Mampe = 0
 R("R_EX_amp_e__Rev"): + Mampe = 0
 R("R_ADSL1r"): + Mfumc - Mdcampc + Mampc = 0
 R("R_ADSL1r_Rev"): - Mfumc + Mdcampc - Mampc = 0
 R("R_ARGSS"): + Mppic + Mhc - McitrLc - Matpc - MaspLc + Margsucc
 + Mampc = 0
 R("R_ARGTRS"): - Mtrnaargc + Mppic - Matpc + Margtrnac - MargLc + Mampc
 = 0
 R("R_ASNS1"): + Mppic + Mhc - Mh2oc + MgluLc - MglnLc - Matpc - MaspLc
 + MasnLc + Mampc = 0
 R("R_ASNS2"): + Mppic - Mnh4c + Mhc - Matpc - MaspLc + MasnLc + Mampc
 = 0
 R("R_ASNTRS"): - Mtrnaasnc + Mppic - Matpc + Masntrnac - MasnLc + Mampc
 = 0
 R("R_ASPTRS"): - Mtrnaaspc + Mppic - Matpc + Masptrnac - MaspLc + Mampc
 = 0
 R("R_BPNT"): + Mpic - Mpapc - Mh2oc + Mampc = 0
 R("R_CYSTRS"): - Mtrnacysc + Mppic + Mcystrnac - McysLc - Matpc + Mampc
 = 0
 R("R_FACOAL100t2pp"): + Mppic - Mhp + Mhc - Mdcap + Mdcacoac - Mcoac
 - Matpc + Mampc = 0
 R("R_FACOAL120t2pp"): + Mppic - Mhp + Mhc - Mddcap + Mddcacoac - Mcoac
 - Matpc + Mampc = 0
 R("R_FACOAL140t2pp"): - Mttdcap + Mtdcoac + Mppic - Mhp + Mhc - Mcoac
 - Matpc + Mampc = 0
 R("R_FACOAL141t2pp"): - Mttdceap + Mtdecoac + Mppic - Mhp + Mhc - Mcoac
 - Matpc + Mampc = 0
 R("R_FACOAL160t2pp"): + Mppic + Mpmtcoac - Mhp - Mhdcap + Mhc - Mcoac
 - Matpc + Mampc = 0
 R("R_FACOAL161t2pp"): + Mppic - Mhp + Mhdcoac - Mhdceap + Mhc - Mcoac
 - Matpc + Mampc = 0
 R("R_FACOAL180t2pp"): + Mstcoac + Mppic - Mocdcap - Mhp + Mhc - Mcoac
 - Matpc + Mampc = 0
 R("R_FACOAL181t2pp"): + Mppic + Modecoac - Mocdceap - Mhp + Mhc - Mcoac
 - Matpc + Mampc = 0
 R("R_FACOAL60t2pp"): + Mppic + Mhxcoac - Mhxap - Mhp + Mhc - Mcoac
 - Matpc + Mampc = 0
 R("R_FACOAL80t2pp"): + Mppic - Moctap + Moccoac - Mhp + Mhc - Mcoac
 - Matpc + Mampc = 0
 R("R_GDPDPK"): + Mppgppc + Mhc - Mgdpc - Matpc + Mampc = 0
 R("R_GLNTRS"): - Mtrnaglnc + Mppic + Mglntrnac - MglnLc - Matpc + Mampc
 = 0
 R("R_GLUTRS"): - Mtrnagluc + Mppic + Mglutrnac - MgluLc - Matpc + Mampc
 = 0
 R("R_GLYTRS"): - Mtrnaglyc + Mppic + Mglytrnac - Mglyc - Matpc + Mampc
 = 0
 R("R_GMPS2"): - Mxmpc + Mppic + 2 Mhc - Mh2oc + Mgmpc + MgluLc - MglnLc
 - Matpc + Mampc = 0
 R("R_GTPDPK"): + Mhc - Mgtpc + Mgdptpc - Matpc + Mampc = 0
 R("R_HISTRS"): - Mtrnahisc + Mppic + Mhistrnac - MhisLc - Matpc + Mampc
 = 0
 R("R_ILETRS"): - Mtrnailec + Mppic + Miletrnac - MileLc - Matpc + Mampc
 = 0
 R("R_LEUTRS"): - Mtrnaleuc + Mppic + Mleutrnac - MleuLc - Matpc + Mampc
 = 0
 R("R_LYSTRS"): - Mtrnalysc + Mppic + Mlystrnac - MlysLc - Matpc + Mampc
 = 0
 R("R_METTRS"): - Mtrnametc + Mppic + Mmettrnac - MmetLc - Matpc + Mampc
 = 0
 R("R_NADDP"): + Mnmnc - Mnadc + 2 Mhc - Mh2oc + Mampc = 0
 R("R_NADS1"): + Mppic - Mnh4c + Mnadc + Mhc - Mdnadc - Matpc + Mampc
 = 0
 R("R_NTPP6"): + Mppic + Mhc - Mh2oc - Matpc + Mampc = 0
 R("R_PACCOAL"): + Mppic + Mphaccoac - Mpacc - Mcoac - Matpc + Mampc = 0
 R("R_PDE1"): + Mhc - Mh2oc - Mcampc + Mampc = 0
 R("R_PHETRS"): - Mtrnaphec + Mppic + Mphetrnac - MpheLc - Matpc + Mampc
 = 0
 R("R_PPS"): - Mpyrc + Mpic + Mpepc + 2 Mhc - Mh2oc - Matpc + Mampc = 0
 R("R_PROTRS"): - Mtrnaproc + Mprotrnac - MproLc + Mppic - Matpc + Mampc
 = 0
 R("R_PRPPS"): - Mr5pc + Mprppc + Mhc - Matpc + Mampc = 0
 R("R_PRPPS_Rev"): + Mr5pc - Mprppc - Mhc + Matpc - Mampc = 0
 R("R_SELNPS"): + Mselnpc - Mselnc + Mpic - Mh2oc - Matpc + Mampc = 0
 R("R_SERTRS"): - Mtrnaserc + Msertrnac - MserLc + Mppic - Matpc + Mampc
 = 0
 R("R_SERTRS2"): - Mtrnasecysc + Msertrnasecoooc - MserLc + Mppic
 - Matpc + Mampc = 0
 R("R_SUCBZL"): - Msucbzc + Msbzcoac + Mppic - Mcoac - Matpc + Mampc = 0
 R("R_THRTRS"): - Mtrnathrc + Mthrtrnac - MthrLc + Mppic - Matpc + Mampc
 = 0
 R("R_TRPTRS"): + Mtrptrnac - MtrpLc - Mtrnatrpc + Mppic - Matpc + Mampc
 = 0
 R("R_TYRTRS"): + Mtyrtrnac - MtyrLc - Mtrnatyrc + Mppic - Matpc + Mampc
 = 0
 R("R_VALTRS"): + Mvaltrnac - MvalLc - Mtrnavalc + Mppic - Matpc + Mampc
 = 0
 R("R_AMPtex"): + Mampp - Mampe = 0
 R("R_AMPtex_Rev"): - Mampp + Mampe = 0
 R("R_AGM3PA"): + MLalaDgluMdapc - Mh2oc + Manhgmc - Manhgm3pc = 0
 R("R_AGM3Pt2pp"): - Mhp + Mhc - Manhgm3pp + Manhgm3pc = 0
 R("R_AGM3PApp"): + MLalaDgluMdapp - Mh2op + Manhgmp - Manhgm3pp = 0
 R("R_MLTGY2pp"): - Mmurein4p3pp + Manhgm4pp + Manhgm3pp = 0
 R("R_MLTGY3pp"): - Mmurein3p3pp + 2 Manhgm3pp = 0
 R("R_AGM4PA"): + MLalaDgluMdapDalac - Mh2oc + Manhgmc - Manhgm4pc = 0
 R("R_AGM4Pt2pp"): - Mhp + Mhc - Manhgm4pp + Manhgm4pc = 0
 R("R_AGM4PApp"): + MLalaDgluMdapDalap - Mh2op + Manhgmp - Manhgm4pp = 0
 R("R_MLTGY1pp"): - Mmurein4p4pp + 2 Manhgm4pp = 0
 R("R_MLTGY4pp"): + Mmurein4px4pp - Mmurein4px4p4pp + Manhgm4pp = 0
 R("Manhgmb_Transp"): 0 M10fthfc = 0
 R("R_EX_anhgm_e_"): - Manhgme = 0
 R("R_EX_anhgm_e__Rev"): + Manhgme = 0
 R("R_AGMt2pp"): - Mhp + Mhc - Manhgmp + Manhgmc = 0
 R("R_ANHGMtex"): + Manhgmp - Manhgme = 0
 R("R_ANHGMtex_Rev"): - Manhgmp + Manhgme = 0
 R("R_AM3PA"): + MLalaDgluMdapc - Mh2oc + Manhmc - Manhm3pc = 0
 R("R_AM4PA"): + MLalaDgluMdapDalac - Mh2oc + Manhmc - Manhm4pc = 0
 R("R_ANPRT"): - Mprppc + Mpranc + Mppic - Manthc = 0
 R("R_ANS"): + Mpyrc + Mhc + MgluLc - MglnLc - Mchorc + Manthc = 0
 R("R_SADT2"): - Mso4c + Mppic + Mpic - Mh2oc - Mgtpc + Mgdpc - Matpc
 + Mapsc = 0
 R("R_A5PISO"): - Mru5pDc + Mara5pc = 0
 R("R_A5PISO_Rev"): + Mru5pDc - Mara5pc = 0
 R("R_KDOPS"): + Mpic - Mpepc + Mkdo8pc - Mh2oc - Mara5pc = 0
 R("MarabLb_Transp"): 0 M10fthfc = 0
 R("R_EX_arab_L_e_"): - MarabLe = 0
 R("R_EX_arab_L_e__Rev"): + MarabLe = 0
 R("R_ARAI"): + MrblLc - MarabLc = 0
 R("R_ARAI_Rev"): - MrblLc + MarabLc = 0
 R("R_ARBt2rpp"): - Mhp + Mhc - MarabLp + MarabLc = 0
 R("R_ARBt2rpp_Rev"): + Mhp - Mhc + MarabLp - MarabLc = 0
 R("R_ARBt3ipp"): - Mhp + Mhc + MarabLp - MarabLc = 0
 R("R_ARBtex"): + MarabLp - MarabLe = 0
 R("R_ARBtex_Rev"): - MarabLp + MarabLe = 0
 R("R_O16GLCT1"): - Mudpgc + Mudpc + Mhc + Mgaragundc - Maragundc = 0
 R("R_AB6PGH"): + Mhqnc - Mh2oc + Mg6pc - Marbt6pc = 0
 R("Marbtnb_Transp"): 0 M10fthfc = 0
 R("R_EX_arbtn_e_"): - Marbtne = 0
 R("R_EX_arbtn_e__Rev"): + Marbtne = 0
 R("R_ARBTNR1"): + 2 Mhc + 2 Mfe2c - Mfadh2c + Mfadc - 2 Marbtnfe3c
 + 2 Marbtnc = 0
 R("R_ARBTNR2"): + 2 Mhc - Mfmnh2c + Mfmnc + 2 Mfe2c - 2 Marbtnfe3c
 + 2 Marbtnc = 0
 R("R_ARBTNR3"): + Mribflvc - Mrbflvrdc + 2 Mhc + 2 Mfe2c - 2 Marbtnfe3c
 + 2 Marbtnc = 0
 R("R_ARBTNtpp"): - Mhp + Mhc + Marbtnp - Marbtnc = 0
 R("R_ARBTNexs"): - Mfe3e + Marbtnfe3e - Marbtne = 0
 R("R_ARBTNtex"): - Mhp + Mhc - Marbtnp + Marbtne = 0
 R("Marbtnfe3b_Transp"): 0 M10fthfc = 0
 R("R_EX_arbtn_fe3_e_"): - Marbtnfe3e = 0
 R("R_EX_arbtn_fe3_e__Rev"): + Marbtnfe3e = 0
 R("R_ARBTNtonex"): - Mhp + Mhc + Marbtnfe3p - Marbtnfe3e = 0
 R("MargLb_Transp"): 0 M10fthfc = 0
 R("R_EX_arg_L_e_"): - MargLe = 0
 R("R_EX_arg_L_e__Rev"): + MargLe = 0
 R("R_ARGORNt7pp"): + Mornp - Mornc - MargLp + MargLc = 0
 R("R_ARGORNt7pp_Rev"): - Mornp + Mornc + MargLp - MargLc = 0
 R("R_ARGSL"): + Mfumc - Margsucc + MargLc = 0
 R("R_ARGSL_Rev"): - Mfumc + Margsucc - MargLc = 0
 R("R_ARGt3pp"): - Mhp + Mhc + MargLp - MargLc = 0
 R("R_AST"): - Msuccoac + Msucargc + Mhc + Mcoac - MargLc = 0
 R("R_ARGtex"): + MargLp - MargLe = 0
 R("R_ARGtex_Rev"): - MargLp + MargLe = 0
 R("R_ASCBptspp"): + Mpyrc - Mpepc - MascbLp + Mascb6pc = 0
 R("MascbLb_Transp"): 0 M10fthfc = 0
 R("R_EX_ascb_L_e_"): - MascbLe = 0
 R("R_EX_ascb_L_e__Rev"): + MascbLe = 0
 R("R_ASCBtex"): + MascbLp - MascbLe = 0
 R("R_ASCBtex_Rev"): - MascbLp + MascbLe = 0
 R("MasnLb_Transp"): 0 M10fthfc = 0
 R("R_EX_asn_L_e_"): - MasnLe = 0
 R("R_EX_asn_L_e__Rev"): + MasnLe = 0
 R("R_ASNN"): + Mnh4c - Mh2oc + MaspLc - MasnLc = 0
 R("R_ASNt2rpp"): - Mhp + Mhc - MasnLp + MasnLc = 0
 R("R_ASNt2rpp_Rev"): + Mhp - Mhc + MasnLp - MasnLc = 0
 R("R_ASNtex"): + MasnLp - MasnLe = 0
 R("R_ASNtex_Rev"): - MasnLp + MasnLe = 0
 R("R_ASNNpp"): + Mnh4p - Mh2op + MaspLp - MasnLp = 0
 R("Maso3b_Transp"): 0 M10fthfc = 0
 R("R_EX_aso3_e_"): - Maso3e = 0
 R("R_EX_aso3_e__Rev"): + Maso3e = 0
 R("R_ASR"): + Mh2oc - 2 Mgthrdc + Mgthoxc - Maso4c + Maso3c = 0
 R("R_ASO3tex"): + Maso3p - Maso3e = 0
 R("R_ASO3tex_Rev"): - Maso3p + Maso3e = 0
 R("MaspLb_Transp"): 0 M10fthfc = 0
 R("R_EX_asp_L_e_"): - MaspLe = 0
 R("R_EX_asp_L_e__Rev"): + MaspLe = 0
 R("R_ADSS"): + Mpic - Mimpc + 2 Mhc - Mgtpc + Mgdpc + Mdcampc - MaspLc
 = 0
 R("R_ASPCT"): + Mpic + Mhc - Mcbpc + Mcbaspc - MaspLc = 0
 R("R_ASPO3"): + Mq8h2c - Mq8c + Miaspc + Mhc - MaspLc = 0
 R("R_ASPO4"): - Mmqn8c + Mmql8c + Miaspc + Mhc - MaspLc = 0
 R("R_ASPO5"): + Msuccc + Miaspc + Mhc - Mfumc - MaspLc = 0
 R("R_ASPO6"): - Mo2c + Miaspc + Mhc + Mh2o2c - MaspLc = 0
 R("R_ASPT"): + Mnh4c + Mfumc - MaspLc = 0
 R("R_ASPt2_2pp"): - 2 Mhp + 2 Mhc - MaspLp + MaspLc = 0
 R("R_ASPt2_3pp"): - 3 Mhp + 3 Mhc - MaspLp + MaspLc = 0
 R("R_ASPt2pp"): - Mhp + Mhc - MaspLp + MaspLc = 0
 R("R_ASPtex"): + MaspLp - MaspLe = 0
 R("R_ASPtex_Rev"): - MaspLp + MaspLe = 0
 R("R_HSDy"): + Mnadphc - Mnadpc - MhomLc + Mhc + Maspsac = 0
 R("R_HSDy_Rev"): - Mnadphc + Mnadpc + MhomLc - Mhc - Maspsac = 0
 R("R_ADNCYC"): + Mppic + Mcampc - Matpc = 0
 R("R_ATPHs"): + Mnh4c + Mitpc - Mhc - Mh2oc - Matpc = 0
 R("R_ATPPRT"): - Mprppc + Mprbatpc + Mppic - Matpc = 0
 R("R_FMNAT"): + Mppic - Mhc - Mfmnc + Mfadc - Matpc = 0
 R("R_NMNAT"): + Mppic - Mnmnc + Mnadc - Mhc - Matpc = 0
 R("R_NNATr"): + Mppic - Mnicrntc - Mhc + Mdnadc - Matpc = 0
 R("R_NNATr_Rev"): - Mppic + Mnicrntc + Mhc - Mdnadc + Matpc = 0
 R("R_PTPATi"): + Mppic - Mpan4pc - Mhc + Mdpcoac - Matpc = 0
 R("R_RNTR1c"): + Mh2oc - Mfldrdc + Mfldoxc + Mdatpc - Matpc = 0
 R("R_SERASr"): - MserLc + Mserampc + Mppic - Mhc - Matpc = 0
 R("R_SERASr_Rev"): + MserLc - Mserampc - Mppic + Mhc + Matpc = 0
 R("R_ACOAD1f"): + Mfadh2c - Mfadc - Mbtcoac + Mb2coac = 0
 R("R_ACOAD1f_Rev"): - Mfadh2c + Mfadc + Mbtcoac - Mb2coac = 0
 R("R_CRNBTCT"): + Mgbbtnc + Mcrncoac - Mcrnc - Mbbtcoac = 0
 R("R_CRNBTCT_Rev"): - Mgbbtnc - Mcrncoac + Mcrnc + Mbbtcoac = 0
 R("R_BETALDHx"): + Mnadhc - Mnadc + 2 Mhc - Mh2oc + Mglybc - Mbetaldc
 = 0
 R("R_BETALDHy"): + Mnadphc - Mnadpc + 2 Mhc - Mh2oc + Mglybc - Mbetaldc
 = 0
 R("R_GLBRAN2"): - Mglycogenc + Mbglycogenc = 0
 R("R_GLCP2"): - Mpic + Mg1pc - Mbglycogenc = 0
 R("R_GLDBRAN2"): + Mglycogenc - Mbglycogenc = 0
 R("R_ALDD4"): + Mnadhc - Mnadc + 2 Mhc - Mh2oc + Mbutc - Mbtalc = 0
 R("R_FDMO4"): + Mso3c - Mo2c + Mhc + Mh2oc - Mfmnh2c + Mfmnc - Mbutso3c
 + Mbtalc = 0
 R("R_BSORx"): - Mnadhc + Mnadc - Mhc + Mh2oc - Mbtnsoc + Mbtnc = 0
 R("R_BSORy"): - Mnadphc + Mnadpc - Mhc + Mh2oc - Mbtnsoc + Mbtnc = 0
 R("R_EAR40x"): - Mnadhc + Mnadc - Mhc + MbutACPc - Mbut2eACPc = 0
 R("R_EAR40y"): - Mnadphc + Mnadpc - Mhc + MbutACPc - Mbut2eACPc = 0
 R("Mbutb_Transp"): 0 M10fthfc = 0
 R("R_EX_but_e_"): - Mbute = 0
 R("R_EX_but_e__Rev"): + Mbute = 0
 R("R_BUTt2rpp"): - Mhp + Mhc - Mbutp + Mbutc = 0
 R("R_BUTt2rpp_Rev"): + Mhp - Mhc + Mbutp - Mbutc = 0
 R("R_BUTtex"): + Mbutp - Mbute = 0
 R("R_BUTtex_Rev"): - Mbutp + Mbute = 0
 R("Mbutso3b_Transp"): 0 M10fthfc = 0
 R("R_EX_butso3_e_"): - Mbutso3e = 0
 R("R_EX_butso3_e__Rev"): + Mbutso3e = 0
 R("R_BUTSO3tex"): + Mbutso3p - Mbutso3e = 0
 R("R_BUTSO3tex_Rev"): - Mbutso3p + Mbutso3e = 0
 R("Mca2b_Transp"): 0 M10fthfc = 0
 R("R_EX_ca2_e_"): - Mca2e = 0
 R("R_EX_ca2_e__Rev"): + Mca2e = 0
 R("R_CA2t3pp"): - Mhp + Mhc + Mca2p - Mca2c = 0
 R("R_CAt6pp"): - Mna1p + Mna1c + Mca2p - Mca2c = 0
 R("R_CAt6pp_Rev"): + Mna1p - Mna1c - Mca2p + Mca2c = 0
 R("R_CA2tex"): + Mca2p - Mca2e = 0
 R("R_CA2tex_Rev"): - Mca2p + Mca2e = 0
 R("R_DHORTS"): + Mhc - Mh2oc - MdhorSc + Mcbaspc = 0
 R("R_DHORTS_Rev"): - Mhc + Mh2oc + MdhorSc - Mcbaspc = 0
 R("Mcbib_Transp"): 0 M10fthfc = 0
 R("R_EX_cbi_e_"): - Mcbie = 0
 R("R_EX_cbi_e__Rev"): + Mcbie = 0
 R("R_CBItonex"): - Mhp + Mhc + Mcbip - Mcbie = 0
 R("Mcbl1b_Transp"): 0 M10fthfc = 0
 R("R_EX_cbl1_e_"): - Mcbl1e = 0
 R("R_EX_cbl1_e__Rev"): + Mcbl1e = 0
 R("R_CBL1tonex"): - Mhp + Mhc + Mcbl1p - Mcbl1e = 0
 R("R_OCBT"): + Mpic - Mornc + Mhc + McitrLc - Mcbpc = 0
 R("R_OCBT_Rev"): - Mpic + Mornc - Mhc - McitrLc + Mcbpc = 0
 R("R_OXAMTC"): - Mpic - Moxurc + Moxamc + Mcbpc = 0
 R("Mcd2b_Transp"): 0 M10fthfc = 0
 R("R_EX_cd2_e_"): - Mcd2e = 0
 R("R_EX_cd2_e__Rev"): + Mcd2e = 0
 R("R_CD2t3pp"): - Mhp + Mhc + Mcd2p - Mcd2c = 0
 R("R_CD2tpp"): - Mcd2p + Mcd2c = 0
 R("R_CD2tex"): + Mcd2p - Mcd2e = 0
 R("R_CD2tex_Rev"): - Mcd2p + Mcd2e = 0
 R("R_EAR121x"): - Mt3c5ddeceACPc - Mnadhc + Mnadc - Mhc + Mcddec5eACPc
 = 0
 R("R_EAR121y"): - Mt3c5ddeceACPc - Mnadphc + Mnadpc - Mhc
 + Mcddec5eACPc = 0
 R("R_T2DECAI"): - Mtdec2eACPc + Mcdec3eACPc = 0
 R("R_T2DECAI_Rev"): + Mtdec2eACPc - Mcdec3eACPc = 0
 R("R_NTP5"): + Mpic + Mhc - Mh2oc - Mctpc + Mcdpc = 0
 R("R_RNDR3"): - Mtrdrdc + Mtrdoxc + Mh2oc + Mdcdpc - Mcdpc = 0
 R("R_RNDR3b"): + Mh2oc - Mgrxrdc + Mgrxoxc + Mdcdpc - Mcdpc = 0
 R("R_CDAPPA120"): + Mpa120c + 2 Mhc - Mh2oc + Mcmpc - Mcdpdddecgc = 0
 R("R_DASYN120"): + Mppic - Mpa120c - Mhc - Mctpc + Mcdpdddecgc = 0
 R("R_PGSA120"): + Mpgp120c + Mhc - Mglyc3pc + Mcmpc - Mcdpdddecgc = 0
 R("R_PSSA120"): - MserLc + Mps120c + Mhc + Mcmpc - Mcdpdddecgc = 0
 R("R_CDAPPA161"): + Mpa161c + 2 Mhc - Mh2oc + Mcmpc - Mcdpdhdec9egc = 0
 R("R_DASYN161"): + Mppic - Mpa161c - Mhc - Mctpc + Mcdpdhdec9egc = 0
 R("R_PGSA161"): + Mpgp161c + Mhc - Mglyc3pc + Mcmpc - Mcdpdhdec9egc = 0
 R("R_PSSA161"): - MserLc + Mps161c + Mhc + Mcmpc - Mcdpdhdec9egc = 0
 R("R_CDAPPA160"): + Mpa160c + 2 Mhc - Mh2oc + Mcmpc - Mcdpdhdecgc = 0
 R("R_DASYN160"): + Mppic - Mpa160c - Mhc - Mctpc + Mcdpdhdecgc = 0
 R("R_PGSA160"): + Mpgp160c + Mhc - Mglyc3pc + Mcmpc - Mcdpdhdecgc = 0
 R("R_PSSA160"): - MserLc + Mps160c + Mhc + Mcmpc - Mcdpdhdecgc = 0
 R("R_CDAPPA181"): + Mpa181c + 2 Mhc - Mh2oc + Mcmpc - Mcdpdodec11egc
 = 0
 R("R_DASYN181"): + Mppic - Mpa181c - Mhc - Mctpc + Mcdpdodec11egc = 0
 R("R_PGSA181"): + Mpgp181c + Mhc - Mglyc3pc + Mcmpc - Mcdpdodec11egc
 = 0
 R("R_PSSA181"): - MserLc + Mps181c + Mhc + Mcmpc - Mcdpdodec11egc = 0
 R("R_CDAPPA180"): + Mpa180c + 2 Mhc - Mh2oc + Mcmpc - Mcdpdodecgc = 0
 R("R_DASYN180"): + Mppic - Mpa180c - Mhc - Mctpc + Mcdpdodecgc = 0
 R("R_PGSA180"): + Mpgp180c + Mhc - Mglyc3pc + Mcmpc - Mcdpdodecgc = 0
 R("R_PSSA180"): - MserLc + Mps180c + Mhc + Mcmpc - Mcdpdodecgc = 0
 R("R_CDAPPA141"): + Mpa141c + 2 Mhc - Mh2oc + Mcmpc - Mcdpdtdec7egc = 0
 R("R_DASYN141"): + Mppic - Mpa141c - Mhc - Mctpc + Mcdpdtdec7egc = 0
 R("R_PGSA141"): + Mpgp141c + Mhc - Mglyc3pc + Mcmpc - Mcdpdtdec7egc = 0
 R("R_PSSA141"): - MserLc + Mps141c + Mhc + Mcmpc - Mcdpdtdec7egc = 0
 R("R_CDAPPA140"): + Mpa140c + 2 Mhc - Mh2oc + Mcmpc - Mcdpdtdecgc = 0
 R("R_DASYN140"): + Mppic - Mpa140c - Mhc - Mctpc + Mcdpdtdecgc = 0
 R("R_PGSA140"): + Mpgp140c + Mhc - Mglyc3pc + Mcmpc - Mcdpdtdecgc = 0
 R("R_PSSA140"): - MserLc + Mps140c + Mhc + Mcmpc - Mcdpdtdecgc = 0
 R("R_DHPPD"): + Mnadhc - Mnadc + Mhc + Mdhpppnc - Mcechdddc = 0
 R("R_PPPNDO"): - Mpppnc - Mo2c - Mnadhc + Mnadc - Mhc + Mcechdddc = 0
 R("R_CINNDO"): - Mo2c - Mnadhc + Mnadc - Mhc - Mcinnmc + Mcenchdddc = 0
 R("R_DHCIND"): + Mnadhc - Mnadc + Mhc + Mdhcinnmc - Mcenchdddc = 0
 R("Mcglyb_Transp"): 0 M10fthfc = 0
 R("R_EX_cgly_e_"): - Mcglye = 0
 R("R_EX_cgly_e__Rev"): + Mcglye = 0
 R("R_AMPTASECG"): - Mh2oc + Mglyc + McysLc - Mcglyc = 0
 R("R_CGLYtex"): + Mcglyp - Mcglye = 0
 R("R_CGLYtex_Rev"): - Mcglyp + Mcglye = 0
 R("R_GTHRDHpp"): - Mh2op - Mgthrdp + MgluLp + Mcglyp = 0
 R("Mcholb_Transp"): 0 M10fthfc = 0
 R("R_EX_chol_e_"): - Mchole = 0
 R("R_EX_chol_e__Rev"): + Mchole = 0
 R("R_CHLt2pp"): - Mhp + Mhc - Mcholp + Mcholc = 0
 R("R_GPDDA1"): + Mhc - Mh2oc + Mglyc3pc - Mg3pcc + Mcholc = 0
 R("R_CHLtex"): + Mcholp - Mchole = 0
 R("R_CHLtex_Rev"): - Mcholp + Mchole = 0
 R("R_GPDDA1pp"): + Mhp - Mh2op + Mglyc3pp - Mg3pcp + Mcholp = 0
 R("R_CHORM"): + Mpphnc - Mchorc = 0
 R("R_ICHORS"): + Michorc - Mchorc = 0
 R("R_ICHORSi"): + Michorc - Mchorc = 0
 R("R_ICHORS_Rev"): - Michorc + Mchorc = 0
 R("Mcitb_Transp"): 0 M10fthfc = 0
 R("R_EX_cit_e_"): - Mcite = 0
 R("R_EX_cit_e__Rev"): + Mcite = 0
 R("R_CITt7pp"): + Msuccp - Msuccc - Mcitp + Mcitc = 0
 R("R_CITtex"): + Mcitp - Mcite = 0
 R("R_CITtex_Rev"): - Mcitp + Mcite = 0
 R("R_KDOCT2"): + Mppic - Mkdoc - Mctpc + Mckdoc = 0
 R("R_MOAT"): - MlipidAc + Mkdolipid4c + Mhc + Mcmpc - Mckdoc = 0
 R("R_MOAT2"): - Mkdolipid4c + Mkdo2lipid4c + Mhc + Mcmpc - Mckdoc = 0
 R("R_MOAT3C"): - Mphphhlipac + Mkphphhlipac + Mhc + Mcmpc - Mckdoc = 0
 R("Mclb_Transp"): 0 M10fthfc = 0
 R("R_EX_cl_e_"): - Mcle = 0
 R("R_EX_cl_e__Rev"): + Mcle = 0
 R("R_CLt3_2pp"): + Mhp - Mhc - 2 Mclp + 2 Mclc = 0
 R("R_CLtex"): + Mclp - Mcle = 0
 R("R_CLtex_Rev"): - Mclp + Mcle = 0
 R("R_CLPNH120pp"): + Mpg120p + Mpa120p + Mhp - Mh2op - Mclpn120p = 0
 R("R_CLPNS120pp"): - 2 Mpg120p + Mglycp + Mclpn120p = 0
 R("R_CLPNS120pp_Rev"): + 2 Mpg120p - Mglycp - Mclpn120p = 0
 R("R_CLPNH140pp"): + Mpg140p + Mpa140p + Mhp - Mh2op - Mclpn140p = 0
 R("R_CLPNS140pp"): - 2 Mpg140p + Mglycp + Mclpn140p = 0
 R("R_CLPNS140pp_Rev"): + 2 Mpg140p - Mglycp - Mclpn140p = 0
 R("R_CLPNH141pp"): + Mpg141p + Mpa141p + Mhp - Mh2op - Mclpn141p = 0
 R("R_CLPNS141pp"): - 2 Mpg141p + Mglycp + Mclpn141p = 0
 R("R_CLPNS141pp_Rev"): + 2 Mpg141p - Mglycp - Mclpn141p = 0
 R("R_CLPNH160pp"): + Mpg160p + Mpa160p + Mhp - Mh2op - Mclpn160p = 0
 R("R_CLPNS160pp"): - 2 Mpg160p + Mglycp + Mclpn160p = 0
 R("R_CLPNS160pp_Rev"): + 2 Mpg160p - Mglycp - Mclpn160p = 0
 R("R_CLPNH161pp"): + Mpg161p + Mpa161p + Mhp - Mh2op - Mclpn161p = 0
 R("R_CLPNS161pp"): - 2 Mpg161p + Mglycp + Mclpn161p = 0
 R("R_CLPNS161pp_Rev"): + 2 Mpg161p - Mglycp - Mclpn161p = 0
 R("R_CLPNH180pp"): + Mpg180p + Mpa180p + Mhp - Mh2op - Mclpn180p = 0
 R("R_CLPNS180pp"): - 2 Mpg180p + Mglycp + Mclpn180p = 0
 R("R_CLPNS180pp_Rev"): + 2 Mpg180p - Mglycp - Mclpn180p = 0
 R("R_CLPNH181pp"): + Mpg181p + Mpa181p + Mhp - Mh2op - Mclpn181p = 0
 R("R_CLPNS181pp"): - 2 Mpg181p + Mglycp + Mclpn181p = 0
 R("R_CLPNS181pp_Rev"): + 2 Mpg181p - Mglycp - Mclpn181p = 0
 R("Mcmpb_Transp"): 0 M10fthfc = 0
 R("R_EX_cmp_e_"): - Mcmpe = 0
 R("R_EX_cmp_e__Rev"): + Mcmpe = 0
 R("R_CMPN"): + Mr5pc - Mh2oc + Mcsnc - Mcmpc = 0
 R("R_CYTDK2"): + Mhc - Mgtpc + Mgdpc - Mcytdc + Mcmpc = 0
 R("R_NTD4"): + Mpic - Mh2oc + Mcytdc - Mcmpc = 0
 R("R_NTPP4"): + Mppic + Mhc - Mh2oc - Mctpc + Mcmpc = 0
 R("R_CMPtex"): + Mcmpp - Mcmpe = 0
 R("R_CMPtex_Rev"): - Mcmpp + Mcmpe = 0
 R("R_NTD4pp"): + Mpip - Mh2op + Mcytdp - Mcmpp = 0
 R("Mco2b_Transp"): 0 M10fthfc = 0
 R("R_EX_co2_e_"): - Mco2e = 0
 R("R_EX_co2_e__Rev"): + Mco2e = 0
 R("R_CO2tpp"): - Mco2p + Mco2c = 0
 R("R_CO2tpp_Rev"): + Mco2p - Mco2c = 0
 R("R_CPPPGO"): + Mpppg9c - Mo2c - 2 Mhc + 2 Mh2oc - Mcpppg3c + 2 Mco2c
 = 0
 R("R_CYNTAH"): + Mnh4c - Mhco3c - 3 Mhc - Mcyntc + 2 Mco2c = 0
 R("R_DXPS"): - Mpyrc - Mhc - Mg3pc + Mdxyl5pc + Mco2c = 0
 R("R_FDH4pp"): + Mq8h2c - Mq8c + Mhp - 2 Mhc - Mforp + Mco2c = 0
 R("R_FDH5pp"): - Mmqn8c + Mmql8c + Mhp - 2 Mhc - Mforp + Mco2c = 0
 R("R_FHL"): - Mhc + Mh2c - Mforc + Mco2c = 0
 R("R_GLYCL"): - Mthfc + Mnh4c + Mnadhc - Mnadc + Mmlthfc - Mglyc
 + Mco2c = 0
 R("R_HCO3E"): + Mhco3c + Mhc - Mh2oc - Mco2c = 0
 R("R_HCO3E_Rev"): - Mhco3c - Mhc + Mh2oc + Mco2c = 0
 R("R_MALDDH"): + Mpyrc + Mnadhc - Mnadc - MmalDc + Mco2c = 0
 R("R_ME1"): + Mpyrc + Mnadhc - Mnadc - MmalLc + Mco2c = 0
 R("R_ME2"): + Mpyrc + Mnadphc - Mnadpc - MmalLc + Mco2c = 0
 R("R_MMCD"): + Mppcoac - MmmcoaSc - Mhc + Mco2c = 0
 R("R_NNDPR"): - Mqulnc - Mprppc + Mppic + Mnicrntc - 2 Mhc + Mco2c = 0
 R("R_OAADC"): + Mpyrc - Moaac - Mhc + Mco2c = 0
 R("R_OMPDC"): + Mumpc - Morot5pc - Mhc + Mco2c = 0
 R("R_ORNDC"): + Mptrcc - Mornc - Mhc + Mco2c = 0
 R("R_PDX5PS"): + Mpic - Mphthrc + Mpdx5pc + Mnadhc - Mnadc + Mhc
 + 2 Mh2oc - Mdxyl5pc + Mco2c = 0
 R("R_PPC"): + Mpic - Mpepc + Moaac + Mhc - Mh2oc - Mco2c = 0
 R("R_PPNDH"): - Mpphnc + Mphpyrc - Mhc + Mh2oc + Mco2c = 0
 R("R_PSD120"): - Mps120c + Mpe120c - Mhc + Mco2c = 0
 R("R_PSD140"): - Mps140c + Mpe140c - Mhc + Mco2c = 0
 R("R_PSD141"): - Mps141c + Mpe141c - Mhc + Mco2c = 0
 R("R_PSD160"): - Mps160c + Mpe160c - Mhc + Mco2c = 0
 R("R_PSD161"): - Mps161c + Mpe161c - Mhc + Mco2c = 0
 R("R_PSD180"): - Mps180c + Mpe180c - Mhc + Mco2c = 0
 R("R_PSD181"): - Mps181c + Mpe181c - Mhc + Mco2c = 0
 R("R_SADH"): + Msucornc - Msucargc + 2 Mnh4c - 2 Mhc - 2 Mh2oc + Mco2c
 = 0
 R("R_UDPGDC"): + MudpLa4oc - Mudpglcurc + Mnadhc - Mnadc + Mco2c = 0
 R("R_UGLYCH"): - Murdglycc + 2 Mnh4c - 2 Mhc - Mh2oc + Mglxc + Mco2c
 = 0
 R("R_UPPDC1"): - Muppg3c - 4 Mhc + Mcpppg3c + 4 Mco2c = 0
 R("R_CO2tex"): + Mco2p - Mco2e = 0
 R("R_CO2tex_Rev"): - Mco2p + Mco2e = 0
 R("R_FACOAE100"): + Mhc - Mh2oc - Mdcacoac + Mdcac + Mcoac = 0
 R("R_FACOAE120"): + Mhc - Mh2oc - Mddcacoac + Mddcac + Mcoac = 0
 R("R_FACOAE140"): + Mttdcac - Mtdcoac + Mhc - Mh2oc + Mcoac = 0
 R("R_FACOAE141"): + Mttdceac - Mtdecoac + Mhc - Mh2oc + Mcoac = 0
 R("R_FACOAE160"): - Mpmtcoac + Mhdcac + Mhc - Mh2oc + Mcoac = 0
 R("R_FACOAE161"): - Mhdcoac + Mhdceac + Mhc - Mh2oc + Mcoac = 0
 R("R_FACOAE180"): - Mstcoac + Mocdcac + Mhc - Mh2oc + Mcoac = 0
 R("R_FACOAE181"): - Modecoac + Mocdceac + Mhc - Mh2oc + Mcoac = 0
 R("R_FACOAE60"): - Mhxcoac + Mhxac + Mhc - Mh2oc + Mcoac = 0
 R("R_FACOAE80"): + Moctac - Moccoac + Mhc - Mh2oc + Mcoac = 0
 R("R_HSST"): + Msuchmsc - Msuccoac - MhomLc + Mcoac = 0
 R("R_NPHS"): - Msbzcoac + Mdhnac + Mcoac = 0
 R("R_PTA2"): - Mppcoac + Mppapc - Mpic + Mcoac = 0
 R("R_THDPS"): - Mthdpc - Msuccoac + Msl2a6oc - Mh2oc + Mcoac = 0
 R("Mcobalt2b_Transp"): 0 M10fthfc = 0
 R("R_EX_cobalt2_e_"): - Mcobalt2e = 0
 R("R_EX_cobalt2_e__Rev"): + Mcobalt2e = 0
 R("R_COBALT2t3pp"): - Mhp + Mhc + Mcobalt2p - Mcobalt2c = 0
 R("R_COBALT2tpp"): - Mcobalt2p + Mcobalt2c = 0
 R("R_COBALT2tex"): + Mcobalt2p - Mcobalt2e = 0
 R("R_COBALT2tex_Rev"): - Mcobalt2p + Mcobalt2e = 0
 R("Mcolipab_Transp"): 0 M10fthfc = 0
 R("R_EX_colipa_e_"): - Mcolipae = 0
 R("R_EX_colipa_e__Rev"): + Mcolipae = 0
 R("R_COLIPAtex"): - Mcolipap + Mcolipae = 0
 R("R_LIPAHT2ex"): - Mhe - Mhdcae + Mhacolipae + Mh2oe - Mcolipae = 0
 R("R_ECA4OALpp"): + Mudcpdpp + Mhp - Meca4undp + Meca4colipap
 - Mcolipap = 0
 R("R_O16A4Lpp"): + Mudcpdpp - Mo16a4undp + Mo16a4colipap + Mhp
 - Mcolipap = 0
 R("Mcpgnb_Transp"): 0 M10fthfc = 0
 R("R_EX_cpgn_e_"): - Mcpgne = 0
 R("R_EX_cpgn_e__Rev"): + Mcpgne = 0
 R("R_CPGNR1"): + 2 Mhc + 2 Mfe2c - Mfadh2c + Mfadc + 2 Mcpgnunc
 - 2 Mcpgnc = 0
 R("R_CPGNR2"): + 2 Mhc - Mfmnh2c + Mfmnc + 2 Mfe2c + 2 Mcpgnunc
 - 2 Mcpgnc = 0
 R("R_CPGNR3"): + Mribflvc - Mrbflvrdc + 2 Mhc + 2 Mfe2c + 2 Mcpgnunc
 - 2 Mcpgnc = 0
 R("R_CPGNexs"): - Mfe3e - Mcpgnune + Mcpgne = 0
 R("R_CPGNtonex"): - Mhp + Mhc + Mcpgnp - Mcpgne = 0
 R("Mcpgnunb_Transp"): 0 M10fthfc = 0
 R("R_EX_cpgn_un_e_"): - Mcpgnune = 0
 R("R_EX_cpgn_un_e__Rev"): + Mcpgnune = 0
 R("R_CPGNUtpp"): - Mhp + Mhc + Mcpgnunp - Mcpgnunc = 0
 R("R_CPGNUtex"): - Mhp + Mhc - Mcpgnunp + Mcpgnune = 0
 R("Mcrnb_Transp"): 0 M10fthfc = 0
 R("R_EX_crn_e_"): - Mcrne = 0
 R("R_EX_crn_e__Rev"): + Mcrne = 0
 R("R_CRNCBCT"): - Mctbtcoac + Mctbtc + Mcrncoac - Mcrnc = 0
 R("R_CRNCBCT_Rev"): + Mctbtcoac - Mctbtc - Mcrncoac + Mcrnc = 0
 R("R_CRNt2rpp"): - Mhp + Mhc - Mcrnp + Mcrnc = 0
 R("R_CRNt2rpp_Rev"): + Mhp - Mhc + Mcrnp - Mcrnc = 0
 R("R_CRNt7pp"): + Mgbbtnp - Mgbbtnc - Mcrnp + Mcrnc = 0
 R("R_CRNt8pp"): - Mcrnp + McrnDp - McrnDc + Mcrnc = 0
 R("R_CRNCAR"): + McrnDcoac - Mcrncoac = 0
 R("R_CRNCAR_Rev"): - McrnDcoac + Mcrncoac = 0
 R("R_CRNCDH"): + Mh2oc + Mctbtcoac - Mcrncoac = 0
 R("R_CRNCDH_Rev"): - Mh2oc - Mctbtcoac + Mcrncoac = 0
 R("R_CRNDt2rpp"): - Mhp + Mhc - McrnDp + McrnDc = 0
 R("R_CRNDt2rpp_Rev"): + Mhp - Mhc + McrnDp - McrnDc = 0
 R("R_CRNtex"): + Mcrnp - Mcrne = 0
 R("R_CRNtex_Rev"): - Mcrnp + Mcrne = 0
 R("Mcsnb_Transp"): 0 M10fthfc = 0
 R("R_EX_csn_e_"): - Mcsne = 0
 R("R_EX_csn_e__Rev"): + Mcsne = 0
 R("R_CSND"): + Murac + Mnh4c - Mhc - Mh2oc - Mcsnc = 0
 R("R_CSNt2pp"): - Mhp + Mhc - Mcsnp + Mcsnc = 0
 R("R_CYTDH"): + MribDc - Mh2oc - Mcytdc + Mcsnc = 0
 R("R_CSNtex"): + Mcsnp - Mcsne = 0
 R("R_CSNtex_Rev"): - Mcsnp + Mcsne = 0
 R("R_CTBTt2rpp"): - Mhp + Mhc - Mctbtp + Mctbtc = 0
 R("R_CTBTt2rpp_Rev"): + Mhp - Mhc + Mctbtp - Mctbtc = 0
 R("R_RNTR3c"): + Mh2oc - Mfldrdc + Mfldoxc + Mdctpc - Mctpc = 0
 R("Mcu2b_Transp"): 0 M10fthfc = 0
 R("R_EX_cu2_e_"): - Mcu2e = 0
 R("R_EX_cu2_e__Rev"): + Mcu2e = 0
 R("R_CU2tpp"): - Mcu2p + Mcu2c = 0
 R("R_CU2tex"): + Mcu2p - Mcu2e = 0
 R("R_CU2tex_Rev"): - Mcu2p + Mcu2e = 0
 R("R_CU1Opp"): - Mo2p - 4 Mhp + 2 Mh2op - 4 Mcup + 4 Mcu2p = 0
 R("Mcub_Transp"): 0 M10fthfc = 0
 R("R_EX_cu_e_"): - Mcue = 0
 R("R_EX_cu_e__Rev"): + Mcue = 0
 R("R_CUt3"): - Mhe + Mhc + Mcue - Mcuc = 0
 R("R_CUtex"): + Mcup - Mcue = 0
 R("R_CUtex_Rev"): - Mcup + Mcue = 0
 R("Mcyanb_Transp"): 0 M10fthfc = 0
 R("R_EX_cyan_e_"): - Mcyane = 0
 R("R_EX_cyan_e__Rev"): + Mcyane = 0
 R("R_CYANST"): - Mtsulc + Mtcyntc + Mso3c + Mhc - Mcyanc = 0
 R("R_MCPST"): + Mtcyntc + Mpyrc - Mmercppyrc + Mhc - Mcyanc = 0
 R("R_CYANtex"): + Mcyanp - Mcyane = 0
 R("R_CYANtex_Rev"): - Mcyanp + Mcyane = 0
 R("R_CYANSTpp"): - Mtsulp + Mtcyntp + Mso3p + Mhp - Mcyanp = 0
 R("Mcyntb_Transp"): 0 M10fthfc = 0
 R("R_EX_cynt_e_"): - Mcynte = 0
 R("R_EX_cynt_e__Rev"): + Mcynte = 0
 R("R_CYNTt2pp"): - Mhp + Mhc - Mcyntp + Mcyntc = 0
 R("R_CYNTtex"): + Mcyntp - Mcynte = 0
 R("R_CYNTtex_Rev"): - Mcyntp + Mcynte = 0
 R("McysDb_Transp"): 0 M10fthfc = 0
 R("R_EX_cys_D_e_"): - McysDe = 0
 R("R_EX_cys_D_e__Rev"): + McysDe = 0
 R("R_CYSDDS"): + Mpyrc + Mnh4c + Mh2sc - Mh2oc - McysDc = 0
 R("R_CYSDtex"): + McysDp - McysDe = 0
 R("R_CYSDtex_Rev"): - McysDp + McysDe = 0
 R("McysLb_Transp"): 0 M10fthfc = 0
 R("R_EX_cys_L_e_"): - McysLe = 0
 R("R_EX_cys_L_e__Rev"): + McysLe = 0
 R("R_CYSDS"): + Mpyrc + Mnh4c + Mh2sc - Mh2oc - McysLc = 0
 R("R_CYStpp"): + McysLp - McysLc = 0
 R("R_SHSL1"): - Msuchmsc + Msuccc + Mhc + McystLc - McysLc = 0
 R("R_CYStex"): + McysLp - McysLe = 0
 R("R_CYStex_Rev"): - McysLp + McysLe = 0
 R("R_CYSTL"): + Mpyrc + Mnh4c + MhcysLc - Mh2oc - McystLc = 0
 R("Mcytdb_Transp"): 0 M10fthfc = 0
 R("R_EX_cytd_e_"): - Mcytde = 0
 R("R_EX_cytd_e__Rev"): + Mcytde = 0
 R("R_CYTD"): + Muric + Mnh4c - Mhc - Mh2oc - Mcytdc = 0
 R("R_CYTDt2pp"): - Mhp + Mhc - Mcytdp + Mcytdc = 0
 R("R_CYTDt2rpp"): - Mhp + Mhc - Mcytdp + Mcytdc = 0
 R("R_CYTDt2rpp_Rev"): + Mhp - Mhc + Mcytdp - Mcytdc = 0
 R("R_CYTDtex"): + Mcytdp - Mcytde = 0
 R("R_CYTDtex_Rev"): - Mcytdp + Mcytde = 0
 R("Mdad2b_Transp"): 0 M10fthfc = 0
 R("R_EX_dad_2_e_"): - Mdad2e = 0
 R("R_EX_dad_2_e__Rev"): + Mdad2e = 0
 R("R_DADA"): + Mnh4c - Mhc - Mh2oc + Mdinc - Mdad2c = 0
 R("R_DADNt2pp"): - Mhp + Mhc - Mdad2p + Mdad2c = 0
 R("R_NTD6"): + Mpic - Mh2oc - Mdampc + Mdad2c = 0
 R("R_DADNtex"): + Mdad2p - Mdad2e = 0
 R("R_DADNtex_Rev"): - Mdad2p + Mdad2e = 0
 R("R_NTD6pp"): + Mpip - Mh2op - Mdampp + Mdad2p = 0
 R("Mdampb_Transp"): 0 M10fthfc = 0
 R("R_EX_damp_e_"): - Mdampe = 0
 R("R_EX_damp_e__Rev"): + Mdampe = 0
 R("R_NTPP5"): + Mppic + Mhc - Mh2oc - Mdatpc + Mdampc = 0
 R("R_DAMPtex"): + Mdampp - Mdampe = 0
 R("R_DAMPtex_Rev"): - Mdampp + Mdampe = 0
 R("R_DATPHs"): + Mnh4c - Mhc - Mh2oc + Mditpc - Mdatpc = 0
 R("R_DB4PS"): - Mru5pDc + Mhc + Mforc + Mdb4pc = 0
 R("R_ACOAD4f"): + Mfadh2c - Mfadc - Mdcacoac + Mdc2coac = 0
 R("R_ACOAD4f_Rev"): - Mfadh2c + Mfadc + Mdcacoac - Mdc2coac = 0
 R("R_EAR100x"): - Mtdec2eACPc - Mnadhc + Mnadc - Mhc + MdcaACPc = 0
 R("R_EAR100y"): - Mtdec2eACPc - Mnadphc + Mnadpc - Mhc + MdcaACPc = 0
 R("Mdcab_Transp"): 0 M10fthfc = 0
 R("R_EX_dca_e_"): - Mdcae = 0
 R("R_EX_dca_e__Rev"): + Mdcae = 0
 R("R_DCAtex"): + Mdcap - Mdcae = 0
 R("R_DCAtex_Rev"): - Mdcap + Mdcae = 0
 R("Mdcmpb_Transp"): 0 M10fthfc = 0
 R("R_EX_dcmp_e_"): - Mdcmpe = 0
 R("R_EX_dcmp_e__Rev"): + Mdcmpe = 0
 R("R_NTD3"): + Mpic - Mh2oc + Mdcytc - Mdcmpc = 0
 R("R_NTPP3"): + Mppic + Mhc - Mh2oc - Mdctpc + Mdcmpc = 0
 R("R_DCMPtex"): + Mdcmpp - Mdcmpe = 0
 R("R_DCMPtex_Rev"): - Mdcmpp + Mdcmpe = 0
 R("R_NTD3pp"): + Mpip - Mh2op + Mdcytp - Mdcmpp = 0
 R("R_DCTPD"): + Mnh4c - Mhc - Mh2oc + Mdutpc - Mdctpc = 0
 R("Mdcytb_Transp"): 0 M10fthfc = 0
 R("R_EX_dcyt_e_"): - Mdcyte = 0
 R("R_EX_dcyt_e__Rev"): + Mdcyte = 0
 R("R_DCYTD"): + Mnh4c - Mhc - Mh2oc + Mduric - Mdcytc = 0
 R("R_DCYTt2pp"): - Mhp + Mhc - Mdcytp + Mdcytc = 0
 R("R_DCYTtex"): + Mdcytp - Mdcyte = 0
 R("R_DCYTtex_Rev"): - Mdcytp + Mdcyte = 0
 R("R_ACOAD5f"): + Mfadh2c - Mfadc - Mddcacoac + Mdd2coac = 0
 R("R_ACOAD5f_Rev"): - Mfadh2c + Mfadc + Mddcacoac - Mdd2coac = 0
 R("R_EAR120x"): - Mtddec2eACPc - Mnadhc + Mnadc - Mhc + MddcaACPc = 0
 R("R_EAR120y"): - Mtddec2eACPc - Mnadphc + Mnadpc - Mhc + MddcaACPc = 0
 R("Mddcab_Transp"): 0 M10fthfc = 0
 R("R_EX_ddca_e_"): - Mddcae = 0
 R("R_EX_ddca_e__Rev"): + Mddcae = 0
 R("R_DDCAtexi"): + Mddcap - Mddcae = 0
 R("R_RNDR2"): - Mtrdrdc + Mtrdoxc + Mh2oc - Mgdpc + Mdgdpc = 0
 R("R_RNDR2b"): + Mh2oc - Mgrxrdc + Mgrxoxc - Mgdpc + Mdgdpc = 0
 R("Mdgmpb_Transp"): 0 M10fthfc = 0
 R("R_EX_dgmp_e_"): - Mdgmpe = 0
 R("R_EX_dgmp_e__Rev"): + Mdgmpe = 0
 R("R_NTD8"): + Mpic - Mh2oc + Mdgsnc - Mdgmpc = 0
 R("R_NTPP1"): + Mppic + Mhc - Mh2oc - Mdgtpc + Mdgmpc = 0
 R("R_DGMPtex"): + Mdgmpp - Mdgmpe = 0
 R("R_DGMPtex_Rev"): - Mdgmpp + Mdgmpe = 0
 R("R_NTD8pp"): + Mpip - Mh2op + Mdgsnp - Mdgmpp = 0
 R("Mdgsnb_Transp"): 0 M10fthfc = 0
 R("R_EX_dgsn_e_"): - Mdgsne = 0
 R("R_EX_dgsn_e__Rev"): + Mdgsne = 0
 R("R_DGSNt2pp"): - Mhp + Mhc - Mdgsnp + Mdgsnc = 0
 R("R_NTPTP1"): + Mpppic - Mh2oc - Mdgtpc + Mdgsnc = 0
 R("R_DGSNtex"): + Mdgsnp - Mdgsne = 0
 R("R_DGSNtex_Rev"): - Mdgsnp + Mdgsne = 0
 R("R_RNTR2c"): + Mh2oc - Mgtpc - Mfldrdc + Mfldoxc + Mdgtpc = 0
 R("Mdhab_Transp"): 0 M10fthfc = 0
 R("R_EX_dha_e_"): - Mdhae = 0
 R("R_EX_dha_e__Rev"): + Mdhae = 0
 R("R_DHAPT"): + Mpyrc - Mpepc + Mdhapc - Mdhac = 0
 R("R_DHAtpp"): - Mdhap + Mdhac = 0
 R("R_DHAtpp_Rev"): + Mdhap - Mdhac = 0
 R("R_F6PA"): + Mg3pc - Mf6pc + Mdhac = 0
 R("R_F6PA_Rev"): - Mg3pc + Mf6pc - Mdhac = 0
 R("R_GLYCDx"): + Mnadhc - Mnadc + Mhc - Mglycc + Mdhac = 0
 R("R_DHAtex"): + Mdhap - Mdhae = 0
 R("R_DHAtex_Rev"): - Mdhap + Mdhae = 0
 R("R_FBA"): + Mg3pc - Mfdpc + Mdhapc = 0
 R("R_FBA_Rev"): - Mg3pc + Mfdpc - Mdhapc = 0
 R("R_FCLPA"): + MlaldLc - Mfc1pc + Mdhapc = 0
 R("R_FCLPA_Rev"): - MlaldLc + Mfc1pc - Mdhapc = 0
 R("R_G3PD2"): + Mnadphc - Mnadpc + Mhc - Mglyc3pc + Mdhapc = 0
 R("R_G3PD2_Rev"): - Mnadphc + Mnadpc - Mhc + Mglyc3pc - Mdhapc = 0
 R("R_G3PD5"): + Mq8h2c - Mq8c - Mglyc3pc + Mdhapc = 0
 R("R_G3PD6"): - Mmqn8c + Mmql8c - Mglyc3pc + Mdhapc = 0
 R("R_MGSA"): + Mpic + Mmthgxlc - Mdhapc = 0
 R("R_QULNS"): + Mqulnc + Mpic - Miaspc + 2 Mh2oc - Mdhapc = 0
 R("R_RMPA"): - Mrml1pc + MlaldLc + Mdhapc = 0
 R("R_RMPA_Rev"): + Mrml1pc - MlaldLc - Mdhapc = 0
 R("R_TGBPA"): - MtagdpDc + Mg3pc + Mdhapc = 0
 R("R_TGBPA_Rev"): + MtagdpDc - Mg3pc - Mdhapc = 0
 R("R_TPI"): + Mg3pc - Mdhapc = 0
 R("R_TPI_Rev"): - Mg3pc + Mdhapc = 0
 R("R_DHCINDO"): - Mo2c + Mhkntdc + Mhc - Mdhcinnmc = 0
 R("R_DHFR"): + Mthfc - Mnadphc + Mnadpc - Mhc - Mdhfc = 0
 R("R_DHFR_Rev"): - Mthfc + Mnadphc - Mnadpc + Mhc + Mdhfc = 0
 R("R_TMDS"): - Mmlthfc - Mdumpc + Mdtmpc + Mdhfc = 0
 R("R_DNMPPA"): + Mpic - Mh2oc - Mdhpmpc + Mdhnptc = 0
 R("R_DHORD2"): + Mq8h2c - Mq8c + Morotc - MdhorSc = 0
 R("R_DHORD5"): + Morotc - Mmqn8c + Mmql8c - MdhorSc = 0
 R("R_HPPPNDO"): - Mo2c + Mhknddc + Mhc - Mdhpppnc = 0
 R("R_DHPTDCs"): + Mhmfurnc + Mh2oc - Mdhptdc = 0
 R("R_RHCCE"): - Mrhcysc + MhcysLc + Mdhptdc = 0
 R("R_NTP11"): + Mpic + Mhc - Mh2oc - Mditpc + Mdidpc = 0
 R("Mdimpb_Transp"): 0 M10fthfc = 0
 R("R_EX_dimp_e_"): - Mdimpe = 0
 R("R_EX_dimp_e__Rev"): + Mdimpe = 0
 R("R_NTD12"): + Mpic - Mh2oc + Mdinc - Mdimpc = 0
 R("R_NTPP10"): + Mppic + Mhc - Mh2oc - Mditpc + Mdimpc = 0
 R("R_DIMPtex"): + Mdimpp - Mdimpe = 0
 R("R_DIMPtex_Rev"): - Mdimpp + Mdimpe = 0
 R("R_NTD12pp"): + Mpip - Mh2op + Mdinp - Mdimpp = 0
 R("Mdinb_Transp"): 0 M10fthfc = 0
 R("R_EX_din_e_"): - Mdine = 0
 R("R_EX_din_e__Rev"): + Mdine = 0
 R("R_DINSt2pp"): - Mhp + Mhc - Mdinp + Mdinc = 0
 R("R_DINStex"): + Mdinp - Mdine = 0
 R("R_DINStex_Rev"): - Mdinp + Mdine = 0
 R("R_DMATT"): + Mppic - Mipdpc + Mgrdpc - Mdmppc = 0
 R("R_DMPPS"): - Mnadhc + Mnadc - Mhc + Mh2oc - Mh2mb4pc + Mdmppc = 0
 R("R_IPDDI"): - Mipdpc + Mdmppc = 0
 R("R_IPDDI_Rev"): + Mipdpc - Mdmppc = 0
 R("Mdmsb_Transp"): 0 M10fthfc = 0
 R("R_EX_dms_e_"): - Mdmse = 0
 R("R_EX_dms_e__Rev"): + Mdmse = 0
 R("R_DMSOR1"): + Mmqn8c - Mmql8c + Mh2oc - Mdmsoc + Mdmsc = 0
 R("R_DMStex"): + Mdmsp - Mdmse = 0
 R("R_DMStex_Rev"): - Mdmsp + Mdmse = 0
 R("Mdmsob_Transp"): 0 M10fthfc = 0
 R("R_EX_dmso_e_"): - Mdmsoe = 0
 R("R_EX_dmso_e__Rev"): + Mdmsoe = 0
 R("R_DMSOtex"): + Mdmsop - Mdmsoe = 0
 R("R_DMSOtex_Rev"): - Mdmsop + Mdmsoe = 0
 R("R_DMSOR1pp"): + Mmqn8c - Mmql8c + Mh2op + Mdmsp - Mdmsop = 0
 R("Mdopab_Transp"): 0 M10fthfc = 0
 R("R_EX_dopa_e_"): - Mdopae = 0
 R("R_EX_dopa_e__Rev"): + Mdopae = 0
 R("R_DOPAtex"): + Mdopap - Mdopae = 0
 R("R_DOPAtex_Rev"): - Mdopap + Mdopae = 0
 R("R_DSBAO1"): + Mq8h2c - Mq8c - Mdsbardp + Mdsbaoxp = 0
 R("R_DSBAO2"): - Mmqn8c + Mmql8c - Mdsbardp + Mdsbaoxp = 0
 R("R_DSBCGT"): - 2 Mgthrdp + Mgthoxp + Mdsbcrdp - Mdsbcoxp = 0
 R("R_TDSR1"): - Mdsbdrdc + Mdsbdoxc + Mdsbcrdp - Mdsbcoxp = 0
 R("R_DSBDR"): - Mtrdrdc + Mtrdoxc + Mdsbdrdc - Mdsbdoxc = 0
 R("R_TDSR2"): + Mdsbgrdp - Mdsbgoxp - Mdsbdrdc + Mdsbdoxc = 0
 R("R_DSBGGT"): - 2 Mgthrdp + Mgthoxp + Mdsbgrdp - Mdsbgoxp = 0
 R("R_SHCHD2"): + Msclc + Mnadhc - Mnadc + Mhc - Mdsclc = 0
 R("R_AADDGT"): + Munagamufc - Munagamuc + Mhc + Mdtdpc - Mdtdp4aaddgc
 = 0
 R("R_TDPDRE"): + Mdtdp4d6dmc - Mdtdp4d6dgc = 0
 R("R_TDPGDH"): + Mh2oc - Mdtdpgluc + Mdtdp4d6dgc = 0
 R("R_TDPDRR"): - Mnadphc + Mnadpc - Mhc + Mdtdprmnc - Mdtdp4d6dmc = 0
 R("R_RHAT1"): - Mkphphhlipac + Micolipac + Mhc - Mdtdprmnc + Mdtdpc = 0
 R("R_G1PTT"): + Mppic - Mhc - Mg1pc - Mdttpc + Mdtdpgluc = 0
 R("Mdtmpb_Transp"): 0 M10fthfc = 0
 R("R_EX_dtmp_e_"): - Mdtmpe = 0
 R("R_EX_dtmp_e__Rev"): + Mdtmpe = 0
 R("R_NTD5"): + Mthymdc + Mpic - Mh2oc - Mdtmpc = 0
 R("R_NTPP7"): + Mppic + Mhc - Mh2oc - Mdttpc + Mdtmpc = 0
 R("R_DTMPtex"): + Mdtmpp - Mdtmpe = 0
 R("R_DTMPtex_Rev"): - Mdtmpp + Mdtmpe = 0
 R("R_NTD5pp"): + Mthymdp + Mpip - Mh2op - Mdtmpp = 0
 R("R_RNDR4"): - Mudpc - Mtrdrdc + Mtrdoxc + Mh2oc + Mdudpc = 0
 R("R_RNDR4b"): - Mudpc + Mh2oc - Mgrxrdc + Mgrxoxc + Mdudpc = 0
 R("Mdumpb_Transp"): 0 M10fthfc = 0
 R("R_EX_dump_e_"): - Mdumpe = 0
 R("R_EX_dump_e__Rev"): + Mdumpe = 0
 R("R_DUTPDP"): + Mppic + Mhc - Mh2oc - Mdutpc + Mdumpc = 0
 R("R_NTD1"): + Mpic - Mh2oc + Mduric - Mdumpc = 0
 R("R_DUMPtex"): + Mdumpp - Mdumpe = 0
 R("R_DUMPtex_Rev"): - Mdumpp + Mdumpe = 0
 R("R_NTD1pp"): + Mpip - Mh2op + Mdurip - Mdumpp = 0
 R("Mdurib_Transp"): 0 M10fthfc = 0
 R("R_EX_duri_e_"): - Mdurie = 0
 R("R_EX_duri_e__Rev"): + Mdurie = 0
 R("R_DURIt2pp"): - Mhp + Mhc - Mdurip + Mduric = 0
 R("R_DURItex"): + Mdurip - Mdurie = 0
 R("R_DURItex_Rev"): - Mdurip + Mdurie = 0
 R("R_RNTR4c"): - Mutpc + Mh2oc - Mfldrdc + Mfldoxc + Mdutpc = 0
 R("R_TALA"): - Ms7pc - Mg3pc + Mf6pc + Me4pc = 0
 R("R_TALA_Rev"): + Ms7pc + Mg3pc - Mf6pc - Me4pc = 0
 R("R_TKT2"): - Mxu5pDc + Mg3pc + Mf6pc - Me4pc = 0
 R("R_TKT2_Rev"): + Mxu5pDc - Mg3pc - Mf6pc + Me4pc = 0
 R("R_ECAP1pp"): - 2 Munagamufp + Mudcpdpp + Mhp + Meca2undp = 0
 R("R_ECAP2pp"): - Munagamufp + Mudcpdpp + Mhp + Meca3undp - Meca2undp
 = 0
 R("R_ECAP3pp"): - Munagamufp + Mudcpdpp + Mhp + Meca4undp - Meca3undp
 = 0
 R("Meca4colipab_Transp"): 0 M10fthfc = 0
 R("R_EX_eca4colipa_e_"): - Meca4colipae = 0
 R("R_EX_eca4colipa_e__Rev"): + Meca4colipae = 0
 R("R_ECA4COLIPAtex"): - Meca4colipap + Meca4colipae = 0
 R("R_IGPDH"): + Mimacpc + Mh2oc - Meig3pc = 0
 R("Menlipab_Transp"): 0 M10fthfc = 0
 R("R_EX_enlipa_e_"): - Menlipae = 0
 R("R_EX_enlipa_e__Rev"): + Menlipae = 0
 R("R_ENLIPAtex"): - Menlipap + Menlipae = 0
 R("Menterb_Transp"): 0 M10fthfc = 0
 R("R_EX_enter_e_"): - Mentere = 0
 R("R_EX_enter_e__Rev"): + Mentere = 0
 R("R_FEENTERR1"): + 2 Mhc - 2 Mfeenterc + 2 Mfe2c - Mfadh2c + Mfadc
 + 2 Menterc = 0
 R("R_FEENTERR2"): + 2 Mhc - Mfmnh2c + Mfmnc - 2 Mfeenterc + 2 Mfe2c
 + 2 Menterc = 0
 R("R_FEENTERR3"): + Mribflvc - Mrbflvrdc + 2 Mhc - 2 Mfeenterc
 + 2 Mfe2c + 2 Menterc = 0
 R("R_FEENTERtpp"): - Mhp + Mhc + Menterp - Menterc = 0
 R("R_FEENTERexs"): + Mfeentere - Mfe3e - Mentere = 0
 R("R_FEENTERtex"): - Mhp + Mhc - Menterp + Mentere = 0
 R("Methab_Transp"): 0 M10fthfc = 0
 R("R_EX_etha_e_"): - Methae = 0
 R("R_EX_etha_e__Rev"): + Methae = 0
 R("R_ETHAt2pp"): - Mhp + Mhc - Methap + Methac = 0
 R("R_GPDDA2"): + Mhc - Mh2oc + Mglyc3pc - Mg3pec + Methac = 0
 R("R_ETHAtex"): + Methap - Methae = 0
 R("R_ETHAtex_Rev"): - Methap + Methae = 0
 R("R_GPDDA2pp"): + Mhp - Mh2op + Mglyc3pp - Mg3pep + Methap = 0
 R("Methso3b_Transp"): 0 M10fthfc = 0
 R("R_EX_ethso3_e_"): - Methso3e = 0
 R("R_EX_ethso3_e__Rev"): + Methso3e = 0
 R("R_ETHSO3tex"): + Methso3p - Methso3e = 0
 R("R_ETHSO3tex_Rev"): - Methso3p + Methso3e = 0
 R("Metohb_Transp"): 0 M10fthfc = 0
 R("R_EX_etoh_e_"): - Metohe = 0
 R("R_EX_etoh_e__Rev"): + Metohe = 0
 R("R_ETOHt2rpp"): - Mhp + Mhc - Metohp + Metohc = 0
 R("R_ETOHt2rpp_Rev"): + Mhp - Mhc + Metohp - Metohc = 0
 R("R_ETOHtex"): + Metohp - Metohe = 0
 R("R_ETOHtex_Rev"): - Metohp + Metohe = 0
 R("R_FRUptspp"): + Mpyrc - Mpepc - Mfrup + Mf1pc = 0
 R("Mf6pb_Transp"): 0 M10fthfc = 0
 R("R_EX_f6p_e_"): - Mf6pe = 0
 R("R_EX_f6p_e__Rev"): + Mf6pe = 0
 R("R_F6PP"): + Mpic - Mh2oc + Mfruc - Mf6pc = 0
 R("R_F6Pt6_2pp"): + 2 Mpip - 2 Mpic - Mf6pp + Mf6pc = 0
 R("R_FBP"): + Mpic - Mh2oc - Mfdpc + Mf6pc = 0
 R("R_FRUpts2pp"): + Mpyrc - Mpepc - Mfrup + Mf6pc = 0
 R("R_G6PDA"): + Mnh4c - Mh2oc - Mgam6pc + Mf6pc = 0
 R("R_GF6PTA"): + MgluLc - MglnLc + Mgam6pc - Mf6pc = 0
 R("R_M1PD"): + Mnadhc - Mnadc - Mmnl1pc + Mhc + Mf6pc = 0
 R("R_M1PD_Rev"): - Mnadhc + Mnadc + Mmnl1pc - Mhc - Mf6pc = 0
 R("R_MAN6PI"): - Mman6pc + Mf6pc = 0
 R("R_MAN6PI_Rev"): + Mman6pc - Mf6pc = 0
 R("R_PGI"): - Mg6pc + Mf6pc = 0
 R("R_PGI_Rev"): + Mg6pc - Mf6pc = 0
 R("R_SBTPD"): - Msbt6pc + Mnadhc - Mnadc + Mhc + Mf6pc = 0
 R("R_SBTPD_Rev"): + Msbt6pc - Mnadhc + Mnadc - Mhc - Mf6pc = 0
 R("R_F6Ptex"): + Mf6pp - Mf6pe = 0
 R("R_F6Ptex_Rev"): - Mf6pp + Mf6pe = 0
 R("R_ACOAD2f"): - Mhxcoac + Mhx2coac + Mfadh2c - Mfadc = 0
 R("R_ACOAD2f_Rev"): + Mhxcoac - Mhx2coac - Mfadh2c + Mfadc = 0
 R("R_ACOAD3f"): - Moccoac + Moc2coac + Mfadh2c - Mfadc = 0
 R("R_ACOAD3f_Rev"): + Moccoac - Moc2coac - Mfadh2c + Mfadc = 0
 R("R_ACOAD6f"): - Mtdcoac + Mtd2coac + Mfadh2c - Mfadc = 0
 R("R_ACOAD6f_Rev"): + Mtdcoac - Mtd2coac - Mfadh2c + Mfadc = 0
 R("R_ACOAD7f"): - Mpmtcoac + Mhdd2coac + Mfadh2c - Mfadc = 0
 R("R_ACOAD7f_Rev"): + Mpmtcoac - Mhdd2coac - Mfadh2c + Mfadc = 0
 R("R_ACOAD8f"): - Mstcoac + Mod2coac + Mfadh2c - Mfadc = 0
 R("R_ACOAD8f_Rev"): + Mstcoac - Mod2coac - Mfadh2c + Mfadc = 0
 R("R_FADRx"): - Mnadhc + Mnadc - Mhc + Mfadh2c - Mfadc = 0
 R("R_FADRx2"): - Mnadphc + Mnadpc - Mhc + Mfadh2c - Mfadc = 0
 R("R_FE3HOXR1"): + 2 Mhc + 2 Mfe3hoxunc - 2 Mfe3hoxc + 2 Mfe2c
 - Mfadh2c + Mfadc = 0
 R("R_FE3Ri"): + 2 Mhc - 2 Mfe3c + 2 Mfe2c - Mfadh2c + Mfadc = 0
 R("R_FECRMR1"): + 2 Mhc + 2 Mfecrmunc - 2 Mfecrmc + 2 Mfe2c - Mfadh2c
 + Mfadc = 0
 R("R_FEOXAMR1"): + 2 Mhc + 2 Mfeoxamunc - 2 Mfeoxamc + 2 Mfe2c
 - Mfadh2c + Mfadc = 0
 R("Mfaldb_Transp"): 0 M10fthfc = 0
 R("R_EX_fald_e_"): - Mfalde = 0
 R("R_EX_fald_e__Rev"): + Mfalde = 0
 R("R_FALDtpp"): - Mfaldp + Mfaldc = 0
 R("R_FALDtpp_Rev"): + Mfaldp - Mfaldc = 0
 R("R_FALGTHLs"): + Mhmgthc - Mgthrdc - Mfaldc = 0
 R("R_FALGTHLs_Rev"): - Mhmgthc + Mgthrdc + Mfaldc = 0
 R("R_FDMO2"): + Mso3c - Mo2c - Mmso3c + Mhc + Mh2oc - Mfmnh2c + Mfmnc
 + Mfaldc = 0
 R("R_MTRPOX"): + MtrpLc - Mo2c - MNmtrpc - Mh2oc + Mh2o2c + Mfaldc = 0
 R("R_SARCOX"): - Msarcsc - Mo2c - Mh2oc + Mh2o2c + Mglyc + Mfaldc = 0
 R("R_FALDtex"): + Mfaldp - Mfalde = 0
 R("R_FALDtex_Rev"): - Mfaldp + Mfalde = 0
 R("R_FCI"): - MfucLc + MfclLc = 0
 R("R_FCI_Rev"): + MfucLc - MfclLc = 0
 R("Mfe2b_Transp"): 0 M10fthfc = 0
 R("R_EX_fe2_e_"): - Mfe2e = 0
 R("R_EX_fe2_e__Rev"): + Mfe2e = 0
 R("R_FCLT"): - Mppp9c + Mphemec + 2 Mhc - Mfe2c = 0
 R("R_FE2t2pp"): - Mhp + Mhc - Mfe2p + Mfe2c = 0
 R("R_FE2t3pp"): - Mhp + Mhc + Mfe2p - Mfe2c = 0
 R("R_FE2tpp"): - Mfe2p + Mfe2c = 0
 R("R_FE3HOXR2"): + 2 Mhc - Mfmnh2c + Mfmnc + 2 Mfe3hoxunc - 2 Mfe3hoxc
 + 2 Mfe2c = 0
 R("R_FE3HOXR3"): + Mribflvc - Mrbflvrdc + 2 Mhc + 2 Mfe3hoxunc
 - 2 Mfe3hoxc + 2 Mfe2c = 0
 R("R_FECRMR2"): + 2 Mhc - Mfmnh2c + Mfmnc + 2 Mfecrmunc - 2 Mfecrmc
 + 2 Mfe2c = 0
 R("R_FECRMR3"): + Mribflvc - Mrbflvrdc + 2 Mhc + 2 Mfecrmunc
 - 2 Mfecrmc + 2 Mfe2c = 0
 R("R_FEOXAMR2"): + 2 Mhc - Mfmnh2c + Mfmnc + 2 Mfeoxamunc - 2 Mfeoxamc
 + 2 Mfe2c = 0
 R("R_FEOXAMR3"): + Mribflvc - Mrbflvrdc + 2 Mhc + 2 Mfeoxamunc
 - 2 Mfeoxamc + 2 Mfe2c = 0
 R("R_SHCHF"): + Mshemec - Msclc + 3 Mhc - Mfe2c = 0
 R("R_FE2tex"): + Mfe2p - Mfe2e = 0
 R("R_FE2tex_Rev"): - Mfe2p + Mfe2e = 0
 R("R_FEROpp"): - Mo2p - 4 Mhp + 2 Mh2op + 4 Mfe3p - 4 Mfe2p = 0
 R("Mfe3b_Transp"): 0 M10fthfc = 0
 R("R_EX_fe3_e_"): - Mfe3e = 0
 R("R_EX_fe3_e__Rev"): + Mfe3e = 0
 R("Mfe3dcitb_Transp"): 0 M10fthfc = 0
 R("R_EX_fe3dcit_e_"): - Mfe3dcite = 0
 R("R_EX_fe3dcit_e__Rev"): + Mfe3dcite = 0
 R("R_FE3DCITtonex"): - Mhp + Mhc + Mfe3dcitp - Mfe3dcite = 0
 R("Mfe3dhbzsb_Transp"): 0 M10fthfc = 0
 R("R_EX_fe3dhbzs_e_"): - Mfe3dhbzse = 0
 R("R_EX_fe3dhbzs_e__Rev"): + Mfe3dhbzse = 0
 R("R_FE3DHBZStonex"): - Mhp + Mhc + Mfe3dhbzsp - Mfe3dhbzse = 0
 R("R_FE3HOXexs"): - Mfe3hoxune + Mfe3hoxe - Mfe3e = 0
 R("R_FE3tex"): + Mfe3p - Mfe3e = 0
 R("R_FE3tex_Rev"): - Mfe3p + Mfe3e = 0
 R("R_FECRMexs"): - Mfecrmune + Mfecrme - Mfe3e = 0
 R("R_FEOXAMexs"): - Mfeoxamune + Mfeoxame - Mfe3e = 0
 R("Mfe3hoxb_Transp"): 0 M10fthfc = 0
 R("R_EX_fe3hox_e_"): - Mfe3hoxe = 0
 R("R_EX_fe3hox_e__Rev"): + Mfe3hoxe = 0
 R("R_FE3HOXtonex"): - Mhp + Mhc + Mfe3hoxp - Mfe3hoxe = 0
 R("Mfe3hoxunb_Transp"): 0 M10fthfc = 0
 R("R_EX_fe3hox_un_e_"): - Mfe3hoxune = 0
 R("R_EX_fe3hox_un_e__Rev"): + Mfe3hoxune = 0
 R("R_FE3HOXUtex"): - Mhp + Mhc + Mfe3hoxunp - Mfe3hoxunc = 0
 R("R_FE3HOXUtpp"): - Mhp + Mhc - Mfe3hoxunp + Mfe3hoxune = 0
 R("Mfecrmb_Transp"): 0 M10fthfc = 0
 R("R_EX_fecrm_e_"): - Mfecrme = 0
 R("R_EX_fecrm_e__Rev"): + Mfecrme = 0
 R("R_FECRMtonex"): - Mhp + Mhc + Mfecrmp - Mfecrme = 0
 R("Mfecrmunb_Transp"): 0 M10fthfc = 0
 R("R_EX_fecrm_un_e_"): - Mfecrmune = 0
 R("R_EX_fecrm_un_e__Rev"): + Mfecrmune = 0
 R("R_FECRMUtpp"): - Mhp + Mhc + Mfecrmunp - Mfecrmunc = 0
 R("R_FECRMUtex"): - Mhp + Mhc - Mfecrmunp + Mfecrmune = 0
 R("Mfeenterb_Transp"): 0 M10fthfc = 0
 R("R_EX_feenter_e_"): - Mfeentere = 0
 R("R_EX_feenter_e__Rev"): + Mfeentere = 0
 R("R_FEENTERtonex"): - Mhp + Mhc + Mfeenterp - Mfeentere = 0
 R("Mfeoxamb_Transp"): 0 M10fthfc = 0
 R("R_EX_feoxam_e_"): - Mfeoxame = 0
 R("R_EX_feoxam_e__Rev"): + Mfeoxame = 0
 R("R_FEOXAMtonex"): - Mhp + Mhc + Mfeoxamp - Mfeoxame = 0
 R("Mfeoxamunb_Transp"): 0 M10fthfc = 0
 R("R_EX_feoxam_un_e_"): - Mfeoxamune = 0
 R("R_EX_feoxam_un_e__Rev"): + Mfeoxamune = 0
 R("R_FEOXAMUtpp"): - Mhp + Mhc + Mfeoxamunp - Mfeoxamunc = 0
 R("R_FEOXAMUtex"): - Mhp + Mhc - Mfeoxamunp + Mfeoxamune = 0
 R("R_FLDR"): - Mnadphc + Mnadpc - Mhc + Mfldrdc - Mfldoxc = 0
 R("R_FDMO"): + Mso3c - Mo2c - Misetacc + Mhc + Mh2oc + Mgcaldc
 - Mfmnh2c + Mfmnc = 0
 R("R_FDMO6"): - Msulfacc + Mso3c - Mo2c + Mhc + Mh2oc + Mglxc - Mfmnh2c
 + Mfmnc = 0
 R("R_FMNRx"): - Mnadhc + Mnadc - Mhc + Mfmnh2c - Mfmnc = 0
 R("R_FMNRx2"): - Mnadphc + Mnadpc - Mhc + Mfmnh2c - Mfmnc = 0
 R("Mforb_Transp"): 0 M10fthfc = 0
 R("R_EX_for_e_"): - Mfore = 0
 R("R_EX_for_e__Rev"): + Mfore = 0
 R("R_FORt2pp"): - Mhp + Mhc - Mforp + Mforc = 0
 R("R_FORtppi"): + Mforp - Mforc = 0
 R("R_SFGTHi"): - MSfglutthc + Mhc - Mh2oc + Mgthrdc + Mforc = 0
 R("R_UPLA4FNF"): + MuLa4nc - MuLa4fnc - Mh2oc + Mforc = 0
 R("R_FORtex"): + Mforp - Mfore = 0
 R("R_FORtex_Rev"): - Mforp + Mfore = 0
 R("R_IMPC"): - Mimpc - Mh2oc + Mfpricac = 0
 R("R_IMPC_Rev"): + Mimpc + Mh2oc - Mfpricac = 0
 R("R_GRTT"): + Mppic - Mipdpc - Mgrdpc + Mfrdpc = 0
 R("R_HEMEOS"): + Mppic - Mphemec + MhemeOc - Mh2oc - Mfrdpc = 0
 R("R_OCTDPS"): + 5 Mppic + Moctdpc - 5 Mipdpc - Mfrdpc = 0
 R("R_UDCPDPS"): + Mudcpdpc + 8 Mppic - 8 Mipdpc - Mfrdpc = 0
 R("Mfrub_Transp"): 0 M10fthfc = 0
 R("R_EX_fru_e_"): - Mfrue = 0
 R("R_EX_fru_e__Rev"): + Mfrue = 0
 R("R_FFSD"): - Msuc6pc - Mh2oc + Mg6pc + Mfruc = 0
 R("R_XYLI2"): - MglcDc + Mfruc = 0
 R("R_XYLI2_Rev"): + MglcDc - Mfruc = 0
 R("R_FRUtex"): + Mfrup - Mfrue = 0
 R("R_FRUtex_Rev"): - Mfrup + Mfrue = 0
 R("Mfrulysb_Transp"): 0 M10fthfc = 0
 R("R_EX_frulys_e_"): - Mfrulyse = 0
 R("R_EX_frulys_e__Rev"): + Mfrulyse = 0
 R("R_FRULYSE"): - Mpsclysc + Mfrulysc = 0
 R("R_FRULYSE_Rev"): + Mpsclysc - Mfrulysc = 0
 R("R_FRULYSt2pp"): - Mhp + Mhc - Mfrulysp + Mfrulysc = 0
 R("R_FRULYStex"): + Mfrulysp - Mfrulyse = 0
 R("R_FRULYStex_Rev"): - Mfrulysp + Mfrulyse = 0
 R("R_FRULYSDG"): + MlysLc - Mh2oc + Mg6pc - Mfrulyspc = 0
 R("R_FRULYSDG_Rev"): - MlysLc + Mh2oc - Mg6pc + Mfrulyspc = 0
 R("Mfruurb_Transp"): 0 M10fthfc = 0
 R("R_EX_fruur_e_"): - Mfruure = 0
 R("R_EX_fruur_e__Rev"): + Mfruure = 0
 R("R_FRUURt2rpp"): - Mhp + Mhc - Mfruurp + Mfruurc = 0
 R("R_FRUURt2rpp_Rev"): + Mhp - Mhc + Mfruurp - Mfruurc = 0
 R("R_GUI1"): - Mglcurc + Mfruurc = 0
 R("R_GUI1_Rev"): + Mglcurc - Mfruurc = 0
 R("R_MANAO"): + Mnadhc - Mnadc - Mmanac + Mhc + Mfruurc = 0
 R("R_MANAO_Rev"): - Mnadhc + Mnadc + Mmanac - Mhc - Mfruurc = 0
 R("R_FRUURtex"): + Mfruurp - Mfruure = 0
 R("R_FRUURtex_Rev"): - Mfruurp + Mfruure = 0
 R("MfucLb_Transp"): 0 M10fthfc = 0
 R("R_EX_fuc_L_e_"): - MfucLe = 0
 R("R_EX_fuc_L_e__Rev"): + MfucLe = 0
 R("R_FUCtpp"): - Mhp + Mhc - MfucLp + MfucLc = 0
 R("R_FUCtpp_Rev"): + Mhp - Mhc + MfucLp - MfucLc = 0
 R("R_FUCtex"): + MfucLp - MfucLe = 0
 R("R_FUCtex_Rev"): - MfucLp + MfucLe = 0
 R("Mfumb_Transp"): 0 M10fthfc = 0
 R("R_EX_fum_e_"): - Mfume = 0
 R("R_EX_fum_e__Rev"): + Mfume = 0
 R("R_FRD2"): + Msuccc + Mmqn8c - Mmql8c - Mfumc = 0
 R("R_FUM"): + MmalLc - Mh2oc - Mfumc = 0
 R("R_FUM_Rev"): - MmalLc + Mh2oc + Mfumc = 0
 R("R_FUMt2_2pp"): - 2 Mhp + 2 Mhc - Mfump + Mfumc = 0
 R("R_FUMt2_3pp"): - 3 Mhp + 3 Mhc - Mfump + Mfumc = 0
 R("R_HKNTDH"): + Mop4enc - Mhkntdc + Mhc - Mh2oc + Mfumc = 0
 R("R_PPPGO3"): + 3 Msuccc - Mpppg9c + Mppp9c - 3 Mfumc = 0
 R("R_SUCDi"): - Msuccc + Mq8h2c - Mq8c + Mfumc = 0
 R("R_SUCFUMtpp"): + Msuccp - Msuccc - Mfump + Mfumc = 0
 R("R_SUCFUMtpp_Rev"): - Msuccp + Msuccc + Mfump - Mfumc = 0
 R("R_FUMtex"): + Mfump - Mfume = 0
 R("R_FUMtex_Rev"): - Mfump + Mfume = 0
 R("Mg1pb_Transp"): 0 M10fthfc = 0
 R("R_EX_g1p_e_"): - Mg1pe = 0
 R("R_EX_g1p_e__Rev"): + Mg1pe = 0
 R("R_GALUi"): - Mutpc + Mudpgc + Mppic - Mhc - Mg1pc = 0
 R("R_GLCP"): - Mpic - Mglycogenc + Mg1pc = 0
 R("R_MLTP1"): - Mpic + Mmaltttrc - Mmaltptc + Mg1pc = 0
 R("R_MLTP1_Rev"): + Mpic - Mmaltttrc + Mmaltptc - Mg1pc = 0
 R("R_MLTP2"): - Mpic + Mmaltptc - Mmalthxc + Mg1pc = 0
 R("R_MLTP2_Rev"): + Mpic - Mmaltptc + Mmalthxc - Mg1pc = 0
 R("R_MLTP3"): - Mpic + Mmalthxc - Mmalthpc + Mg1pc = 0
 R("R_MLTP3_Rev"): + Mpic - Mmalthxc + Mmalthpc - Mg1pc = 0
 R("R_PGMT"): + Mg6pc - Mg1pc = 0
 R("R_PGMT_Rev"): - Mg6pc + Mg1pc = 0
 R("R_UGLT"): - Mudpgc + Mudpgalc - Mgal1pc + Mg1pc = 0
 R("R_UGLT_Rev"): + Mudpgc - Mudpgalc + Mgal1pc - Mg1pc = 0
 R("R_G1Ptex"): + Mg1pp - Mg1pe = 0
 R("R_G1Ptex_Rev"): - Mg1pp + Mg1pe = 0
 R("R_G1PPpp"): + Mpip - Mh2op + MglcDp - Mg1pp = 0
 R("R_UDPGPpp"): + Mumpp - Mudpgp + 2 Mhp - Mh2op + Mg1pp = 0
 R("R_TKT1"): - Mxu5pDc + Ms7pc - Mr5pc + Mg3pc = 0
 R("R_TKT1_Rev"): + Mxu5pDc - Ms7pc + Mr5pc - Mg3pc = 0
 R("Mg3pcb_Transp"): 0 M10fthfc = 0
 R("R_EX_g3pc_e_"): - Mg3pce = 0
 R("R_EX_g3pc_e__Rev"): + Mg3pce = 0
 R("R_G3PCtex"): + Mg3pcp - Mg3pce = 0
 R("R_G3PCtex_Rev"): - Mg3pcp + Mg3pce = 0
 R("Mg3peb_Transp"): 0 M10fthfc = 0
 R("R_EX_g3pe_e_"): - Mg3pee = 0
 R("R_EX_g3pe_e__Rev"): + Mg3pee = 0
 R("R_G3PEtex"): + Mg3pep - Mg3pee = 0
 R("R_G3PEtex_Rev"): - Mg3pep + Mg3pee = 0
 R("Mg3pgb_Transp"): 0 M10fthfc = 0
 R("R_EX_g3pg_e_"): - Mg3pge = 0
 R("R_EX_g3pg_e__Rev"): + Mg3pge = 0
 R("R_GPDDA4"): + Mhc - Mh2oc + Mglycc + Mglyc3pc - Mg3pgc = 0
 R("R_G3PGtex"): + Mg3pgp - Mg3pge = 0
 R("R_G3PGtex_Rev"): - Mg3pgp + Mg3pge = 0
 R("R_GPDDA4pp"): + Mhp - Mh2op + Mglycp + Mglyc3pp - Mg3pgp = 0
 R("Mg3pib_Transp"): 0 M10fthfc = 0
 R("R_EX_g3pi_e_"): - Mg3pie = 0
 R("R_EX_g3pi_e__Rev"): + Mg3pie = 0
 R("R_GPDDA5"): + Minostc + Mhc - Mh2oc + Mglyc3pc - Mg3pic = 0
 R("R_G3PItex"): + Mg3pip - Mg3pie = 0
 R("R_G3PItex_Rev"): - Mg3pip + Mg3pie = 0
 R("R_GPDDA5pp"): + Minostp + Mhp - Mh2op + Mglyc3pp - Mg3pip = 0
 R("Mg3psb_Transp"): 0 M10fthfc = 0
 R("R_EX_g3ps_e_"): - Mg3pse = 0
 R("R_EX_g3ps_e__Rev"): + Mg3pse = 0
 R("R_GPDDA3"): + MserLc + Mhc - Mh2oc + Mglyc3pc - Mg3psc = 0
 R("R_G3PStex"): + Mg3psp - Mg3pse = 0
 R("R_G3PStex_Rev"): - Mg3psp + Mg3pse = 0
 R("R_GPDDA3pp"): + MserLp + Mhp - Mh2op + Mglyc3pp - Mg3psp = 0
 R("Mg6pb_Transp"): 0 M10fthfc = 0
 R("R_EX_g6p_e_"): - Mg6pe = 0
 R("R_EX_g6p_e__Rev"): + Mg6pe = 0
 R("R_G6PP"): + Mpic - Mh2oc + MglcDc - Mg6pc = 0
 R("R_G6Pt6_2pp"): + 2 Mpip - 2 Mpic - Mg6pp + Mg6pc = 0
 R("R_GLCptspp"): + Mpyrc - Mpepc - MglcDp + Mg6pc = 0
 R("R_TRE6PH"): - Mtre6pc - Mh2oc + MglcDc + Mg6pc = 0
 R("R_TRE6PS"): - Mudpgc + Mudpc + Mtre6pc + Mhc - Mg6pc = 0
 R("R_G6Ptex"): + Mg6pp - Mg6pe = 0
 R("R_G6Ptex_Rev"): - Mg6pp + Mg6pe = 0
 R("R_GALT1"): - Mudpgc + Mudpc + Mhc - Mgicolipac + Mgagicolipac = 0
 R("R_GLCTR2"): - Mudpgc + Mudpc + Mhc + Mggagicolipac - Mgagicolipac
 = 0
 R("Mgal1pb_Transp"): 0 M10fthfc = 0
 R("R_EX_gal1p_e_"): - Mgal1pe = 0
 R("R_EX_gal1p_e__Rev"): + Mgal1pe = 0
 R("R_GAL1Ptex"): + Mgal1pp - Mgal1pe = 0
 R("R_GAL1Ptex_Rev"): - Mgal1pp + Mgal1pe = 0
 R("R_GAL1PPpp"): + Mpip - Mh2op + Mgalp - Mgal1pp = 0
 R("R_UDPGALPpp"): + Mumpp - Mudpgalp + 2 Mhp - Mh2op + Mgal1pp = 0
 R("Mgalb_Transp"): 0 M10fthfc = 0
 R("R_EX_gal_e_"): - Mgale = 0
 R("R_EX_gal_e__Rev"): + Mgale = 0
 R("MgalbDb_Transp"): 0 M10fthfc = 0
 R("R_EX_gal_bD_e_"): - MgalbDe = 0
 R("R_EX_gal_bD_e__Rev"): + MgalbDe = 0
 R("R_GALBDtex"): + MgalbDp - MgalbDe = 0
 R("R_GALBDtex_Rev"): - MgalbDp + MgalbDe = 0
 R("R_GALM2pp"): + Mgalp - MgalbDp = 0
 R("R_GALS3"): - Mmelibc - Mh2oc + MglcDc + Mgalc = 0
 R("R_GALt2pp"): - Mhp + Mhc - Mgalp + Mgalc = 0
 R("R_LACZ"): - Mlctsc - Mh2oc + MglcDc + Mgalc = 0
 R("MgalctDb_Transp"): 0 M10fthfc = 0
 R("R_EX_galct_D_e_"): - MgalctDe = 0
 R("R_EX_galct_D_e__Rev"): + MgalctDe = 0
 R("R_GALCTt2rpp"): - Mhp + Mhc - MgalctDp + MgalctDc = 0
 R("R_GALCTt2rpp_Rev"): + Mhp - Mhc + MgalctDp - MgalctDc = 0
 R("R_GALCTtex"): + MgalctDp - MgalctDe = 0
 R("R_GALCTtex_Rev"): - MgalctDp + MgalctDe = 0
 R("MgalctnDb_Transp"): 0 M10fthfc = 0
 R("R_EX_galctn_D_e_"): - MgalctnDe = 0
 R("R_EX_galctn_D_e__Rev"): + MgalctnDe = 0
 R("R_GALCTNt2pp"): - Mhp + Mhc - MgalctnDp + MgalctnDc = 0
 R("R_GALCTNtex"): + MgalctnDp - MgalctnDe = 0
 R("R_GALCTNtex_Rev"): - MgalctnDp + MgalctnDe = 0
 R("MgalctnLb_Transp"): 0 M10fthfc = 0
 R("R_EX_galctn_L_e_"): - MgalctnLe = 0
 R("R_EX_galctn_L_e__Rev"): + MgalctnLe = 0
 R("R_GALCTLO"): + Mtagurc + Mnadhc - Mnadc + Mhc - MgalctnLc = 0
 R("R_GALCTNLt2pp"): - Mhp + Mhc - MgalctnLp + MgalctnLc = 0
 R("R_GALCTNLtex"): + MgalctnLp - MgalctnLe = 0
 R("R_GALCTNLtex_Rev"): - MgalctnLp + MgalctnLe = 0
 R("R_GALtex"): + Mgalp - Mgale = 0
 R("R_GALtex_Rev"): - Mgalp + Mgale = 0
 R("R_LACZpp"): - Mlctsp - Mh2op + MglcDp + Mgalp = 0
 R("R_GALTptspp"): + Mpyrc - Mpepc - Mgaltp + Mgalt1pc = 0
 R("R_GLTPD"): + Mtag6pDc + Mnadhc - Mnadc + Mhc - Mgalt1pc = 0
 R("R_GLTPD_Rev"): - Mtag6pDc - Mnadhc + Mnadc - Mhc + Mgalt1pc = 0
 R("Mgaltb_Transp"): 0 M10fthfc = 0
 R("R_EX_galt_e_"): - Mgalte = 0
 R("R_EX_galt_e__Rev"): + Mgalte = 0
 R("R_GALTtex"): + Mgaltp - Mgalte = 0
 R("R_GALTtex_Rev"): - Mgaltp + Mgalte = 0
 R("Mgalurb_Transp"): 0 M10fthfc = 0
 R("R_EX_galur_e_"): - Mgalure = 0
 R("R_EX_galur_e__Rev"): + Mgalure = 0
 R("R_GALURt2rpp"): - Mhp + Mhc - Mgalurp + Mgalurc = 0
 R("R_GALURt2rpp_Rev"): + Mhp - Mhc + Mgalurp - Mgalurc = 0
 R("R_GUI2"): + Mtagurc - Mgalurc = 0
 R("R_GUI2_Rev"): - Mtagurc + Mgalurc = 0
 R("R_GALURtex"): + Mgalurp - Mgalure = 0
 R("R_GALURtex_Rev"): - Mgalurp + Mgalure = 0
 R("R_PGAMT"): + Mgam6pc - Mgam1pc = 0
 R("R_PGAMT_Rev"): - Mgam6pc + Mgam1pc = 0
 R("Mgam6pb_Transp"): 0 M10fthfc = 0
 R("R_EX_gam6p_e_"): - Mgam6pe = 0
 R("R_EX_gam6p_e__Rev"): + Mgam6pe = 0
 R("R_GAM6Pt6_2pp"): + 2 Mpip - 2 Mpic - Mgam6pp + Mgam6pc = 0
 R("R_GAMptspp"): + Mpyrc - Mpepc - Mgamp + Mgam6pc = 0
 R("R_GAMAN6Ptex"): + Mgam6pp - Mgam6pe = 0
 R("R_GAMAN6Ptex_Rev"): - Mgam6pp + Mgam6pe = 0
 R("Mgamb_Transp"): 0 M10fthfc = 0
 R("R_EX_gam_e_"): - Mgame = 0
 R("R_EX_gam_e__Rev"): + Mgame = 0
 R("R_GAMtex"): + Mgamp - Mgame = 0
 R("R_GAMtex_Rev"): - Mgamp + Mgame = 0
 R("R_O16GALFT"): - Mudpgalfurc + Mudpc + Mhc + Mgfgaragundc
 - Mgaragundc = 0
 R("Mgbbtnb_Transp"): 0 M10fthfc = 0
 R("R_EX_gbbtn_e_"): - Mgbbtne = 0
 R("R_EX_gbbtn_e__Rev"): + Mgbbtne = 0
 R("R_GBBTNtex"): + Mgbbtnp - Mgbbtne = 0
 R("R_GBBTNtex_Rev"): - Mgbbtnp + Mgbbtne = 0
 R("R_GCALDD"): + Mnadhc - Mnadc + 2 Mhc - Mh2oc + Mglycltc - Mgcaldc
 = 0
 R("Mgdpb_Transp"): 0 M10fthfc = 0
 R("R_EX_gdp_e_"): - Mgdpe = 0
 R("R_EX_gdp_e__Rev"): + Mgdpe = 0
 R("R_GDPMNH"): + Mmanc + Mhc - Mh2oc - Mgdpmannc + Mgdpc = 0
 R("R_GP4GH"): + 2 Mhc - Mh2oc - Mgp4gc + 2 Mgdpc = 0
 R("R_MAN1PT2"): + Mpic - Mman1pc - Mhc + Mgdpmannc - Mgdpc = 0
 R("R_NTP3"): + Mpic + Mhc - Mh2oc - Mgtpc + Mgdpc = 0
 R("R_PPGPPDP"): + Mppic - Mppgppc - Mh2oc + Mgdpc = 0
 R("R_URIK2"): - Muric + Mumpc + Mhc - Mgtpc + Mgdpc = 0
 R("R_GDMANE"): + Mgdpofucc - Mgdpddmanc = 0
 R("R_GMAND"): + Mh2oc - Mgdpmannc + Mgdpddmanc = 0
 R("R_GDPtex"): + Mgdpp - Mgdpe = 0
 R("R_GDPtex_Rev"): - Mgdpp + Mgdpe = 0
 R("R_GOFUCR"): - Mnadphc + Mnadpc - Mhc - Mgdpofucc + Mgdpfucc = 0
 R("R_GDPMNP"): + Mman1pc + 2 Mhc - Mh2oc + Mgmpc - Mgdpmannc = 0
 R("R_NTP3pp"): + Mpip + Mhp - Mh2op - Mgtpp + Mgdpp = 0
 R("R_GTPDPDP"): + Mppgppc + Mpic + Mhc - Mh2oc - Mgdptpc = 0
 R("R_O16GLCT2"): - Mudpgc + Mudpc + Mo16aundc + Mhc - Mgfgaragundc = 0
 R("R_GGGABADr"): + Mnadphc - Mnadpc + 2 Mhc - Mh2oc - Mggbutalc
 + Mgg4abutc = 0
 R("R_GGGABADr_Rev"): - Mnadphc + Mnadpc - 2 Mhc + Mh2oc + Mggbutalc
 - Mgg4abutc = 0
 R("R_GLCTR3"): - Mudpgc + Mudpc + Mhc + Mgggagicolipac - Mggagicolipac
 = 0
 R("R_GGPTRCO"): - Mo2c + Mnh4c - Mh2oc + Mh2o2c - Mggptrcc + Mggbutalc
 = 0
 R("R_GLCTR1"): - Mudpgc + Mudpc - Micolipac + Mhc + Mgicolipac = 0
 R("MglcDb_Transp"): 0 M10fthfc = 0
 R("R_EX_glc_e_"): - MglcDe = 0
 R("R_EX_glc_e__Rev"): + MglcDe = 0
 R("R_AMALT1"): + Mmaltttrc - Mmalttrc - Mmaltc + MglcDc = 0
 R("R_AMALT2"): - Mmaltttrc + Mmaltptc - Mmaltc + MglcDc = 0
 R("R_AMALT3"): - Mmaltptc + Mmalthxc - Mmaltc + MglcDc = 0
 R("R_AMALT4"): - Mmalthxc + Mmalthpc - Mmaltc + MglcDc = 0
 R("R_GLCt2pp"): - Mhp + Mhc - MglcDp + MglcDc = 0
 R("R_MLTG1"): - Mmalttrc + Mmaltc - Mh2oc + MglcDc = 0
 R("R_MLTG2"): - Mmaltttrc + Mmalttrc - Mh2oc + MglcDc = 0
 R("R_MLTG3"): + Mmaltttrc - Mmaltptc - Mh2oc + MglcDc = 0
 R("R_MLTG4"): + Mmaltptc - Mmalthxc - Mh2oc + MglcDc = 0
 R("R_MLTG5"): + Mmalthxc - Mmalthpc - Mh2oc + MglcDc = 0
 R("R_TREH"): - Mtrec - Mh2oc + 2 MglcDc = 0
 R("R_GLCtex"): + MglcDp - MglcDe = 0
 R("R_GLCtexi"): + MglcDp - MglcDe = 0
 R("R_GLCtex_Rev"): - MglcDp + MglcDe = 0
 R("R_GLCDpp"): + Mq8h2c - Mq8c + Mhp - Mh2op + Mglcnp - MglcDp = 0
 R("R_TREHpp"): - Mtrep - Mh2op + 2 MglcDp = 0
 R("Mglcnb_Transp"): 0 M10fthfc = 0
 R("R_EX_glcn_e_"): - Mglcne = 0
 R("R_EX_glcn_e__Rev"): + Mglcne = 0
 R("R_GLCNt2rpp"): - Mhp + Mhc - Mglcnp + Mglcnc = 0
 R("R_GLCNt2rpp_Rev"): + Mhp - Mhc + Mglcnp - Mglcnc = 0
 R("R_GLCNtex"): + Mglcnp - Mglcne = 0
 R("R_GLCNtex_Rev"): - Mglcnp + Mglcne = 0
 R("Mglcrb_Transp"): 0 M10fthfc = 0
 R("R_EX_glcr_e_"): - Mglcre = 0
 R("R_EX_glcr_e__Rev"): + Mglcre = 0
 R("R_GLCRt2rpp"): - Mhp + Mhc - Mglcrp + Mglcrc = 0
 R("R_GLCRt2rpp_Rev"): + Mhp - Mhc + Mglcrp - Mglcrc = 0
 R("R_GLCRtex"): + Mglcrp - Mglcre = 0
 R("R_GLCRtex_Rev"): - Mglcrp + Mglcre = 0
 R("Mglcur1pb_Transp"): 0 M10fthfc = 0
 R("R_EX_glcur1p_e_"): - Mglcur1pe = 0
 R("R_EX_glcur1p_e__Rev"): + Mglcur1pe = 0
 R("R_GLCUR1Ptex"): + Mglcur1pp - Mglcur1pe = 0
 R("R_GLCUR1Ptex_Rev"): - Mglcur1pp + Mglcur1pe = 0
 R("R_GUR1PPpp"): + Mpip - Mh2op + Mglcurp - Mglcur1pp = 0
 R("R_UGLCURPpp"): + Mumpp - Mudpglcurp + 2 Mhp - Mh2op + Mglcur1pp = 0
 R("Mglcurb_Transp"): 0 M10fthfc = 0
 R("R_EX_glcur_e_"): - Mglcure = 0
 R("R_EX_glcur_e__Rev"): + Mglcure = 0
 R("R_GLCURt2rpp"): - Mhp + Mhc - Mglcurp + Mglcurc = 0
 R("R_GLCURt2rpp_Rev"): + Mhp - Mhc + Mglcurp - Mglcurc = 0
 R("R_GLCURtex"): + Mglcurp - Mglcure = 0
 R("R_GLCURtex_Rev"): - Mglcurp + Mglcure = 0
 R("MglnLb_Transp"): 0 M10fthfc = 0
 R("R_EX_gln_L_e_"): - MglnLe = 0
 R("R_EX_gln_L_e__Rev"): + MglnLe = 0
 R("R_GLUN"): + Mnh4c - Mh2oc + MgluLc - MglnLc = 0
 R("R_GLUPRT"): - Mprppc + Mpramc + Mppic - Mh2oc + MgluLc - MglnLc = 0
 R("R_GLNtex"): + MglnLp - MglnLe = 0
 R("R_GLNtex_Rev"): - MglnLp + MglnLe = 0
 R("R_GLUNpp"): + Mnh4p - Mh2op + MgluLp - MglnLp = 0
 R("R_GLUTRR"): + Mtrnagluc - Mnadphc + Mnadpc - Mhc - Mglutrnac
 + Mglu1sac = 0
 R("R_G5SD"): + Mpic - Mnadphc + Mnadpc - Mhc + Mglu5sac - Mglu5pc = 0
 R("R_GLUR"): + MgluLc - MgluDc = 0
 R("R_GLUR_Rev"): - MgluLc + MgluDc = 0
 R("MgluLb_Transp"): 0 M10fthfc = 0
 R("R_EX_glu_L_e_"): - MgluLe = 0
 R("R_EX_glu_L_e__Rev"): + MgluLe = 0
 R("R_GLUt2rpp"): - Mhp + Mhc - MgluLp + MgluLc = 0
 R("R_GLUt2rpp_Rev"): + Mhp - Mhc + MgluLp - MgluLc = 0
 R("R_GLUt4pp"): - Mna1p + Mna1c - MgluLp + MgluLc = 0
 R("R_SGDS"): - Msucgluc + Msuccc - Mh2oc + MgluLc = 0
 R("R_GLUtex"): + MgluLp - MgluLe = 0
 R("R_GLUtex_Rev"): - MgluLp + MgluLe = 0
 R("R_GLYCLTDx"): - Mnadhc + Mnadc - Mhc + Mglycltc - Mglxc = 0
 R("R_GLYCLTDy"): - Mnadphc + Mnadpc - Mhc + Mglycltc - Mglxc = 0
 R("R_GLYCTO2"): + Mq8h2c - Mq8c - Mglycltc + Mglxc = 0
 R("R_GLYCTO3"): - Mmqn8c + Mmql8c - Mglycltc + Mglxc = 0
 R("R_ICL"): + Msuccc - Micitc + Mglxc = 0
 R("Mglyaldb_Transp"): 0 M10fthfc = 0
 R("R_EX_glyald_e_"): - Mglyalde = 0
 R("R_EX_glyald_e__Rev"): + Mglyalde = 0
 R("R_ALCD19"): - Mnadhc + Mnadc - Mhc + Mglycc - Mglyaldc = 0
 R("R_ALCD19_Rev"): + Mnadhc - Mnadc + Mhc - Mglycc + Mglyaldc = 0
 R("R_GLYALDtpp"): - Mglyaldp + Mglyaldc = 0
 R("R_GLYALDtpp_Rev"): + Mglyaldp - Mglyaldc = 0
 R("R_GLYALDtex"): + Mglyaldp - Mglyalde = 0
 R("R_GLYALDtex_Rev"): - Mglyaldp + Mglyalde = 0
 R("Mglyb_Transp"): 0 M10fthfc = 0
 R("R_EX_gly_e_"): - Mglye = 0
 R("R_EX_gly_e__Rev"): + Mglye = 0
 R("Mglybb_Transp"): 0 M10fthfc = 0
 R("R_EX_glyb_e_"): - Mglybe = 0
 R("R_EX_glyb_e__Rev"): + Mglybe = 0
 R("R_GLYBt2pp"): - Mhp + Mhc - Mglybp + Mglybc = 0
 R("R_GLYBtex"): + Mglybp - Mglybe = 0
 R("R_GLYBtex_Rev"): - Mglybp + Mglybe = 0
 R("R_AMPTASEPG"): + MproLc - Mproglyc - Mh2oc + Mglyc = 0
 R("R_GHMT2r"): - Mthfc - MserLc + Mmlthfc + Mh2oc + Mglyc = 0
 R("R_GHMT2r_Rev"): + Mthfc + MserLc - Mmlthfc - Mh2oc - Mglyc = 0
 R("R_GLYt2pp"): - Mhp + Mhc - Mglyp + Mglyc = 0
 R("R_GLYt4pp"): - Mna1p + Mna1c - Mglyp + Mglyc = 0
 R("Mglyc2pb_Transp"): 0 M10fthfc = 0
 R("R_EX_glyc2p_e_"): - Mglyc2pe = 0
 R("R_EX_glyc2p_e__Rev"): + Mglyc2pe = 0
 R("R_GLYC2Ptex"): + Mglyc2pp - Mglyc2pe = 0
 R("R_GLYC2Ptex_Rev"): - Mglyc2pp + Mglyc2pe = 0
 R("R_G2PPpp"): + Mpip - Mh2op + Mglycp - Mglyc2pp = 0
 R("Mglyc3pb_Transp"): 0 M10fthfc = 0
 R("R_EX_glyc3p_e_"): - Mglyc3pe = 0
 R("R_EX_glyc3p_e__Rev"): + Mglyc3pe = 0
 R("R_G3PT"): + Mpic - Mh2oc + Mglycc - Mglyc3pc = 0
 R("R_GLYC3Pt6pp"): + Mpip - Mpic - Mglyc3pp + Mglyc3pc = 0
 R("R_GLYC3Ptex"): + Mglyc3pp - Mglyc3pe = 0
 R("R_GLYC3Ptex_Rev"): - Mglyc3pp + Mglyc3pe = 0
 R("Mglycb_Transp"): 0 M10fthfc = 0
 R("R_EX_glyc_e_"): - Mglyce = 0
 R("R_EX_glyc_e__Rev"): + Mglyce = 0
 R("R_GLYCtpp"): + Mglycp - Mglycc = 0
 R("R_GLYCtpp_Rev"): - Mglycp + Mglycc = 0
 R("R_GLYCtex"): + Mglycp - Mglyce = 0
 R("R_GLYCtex_Rev"): - Mglycp + Mglyce = 0
 R("Mglycltb_Transp"): 0 M10fthfc = 0
 R("R_EX_glyclt_e_"): - Mglyclte = 0
 R("R_EX_glyclt_e__Rev"): + Mglyclte = 0
 R("R_GLYCLTt2rpp"): - Mhp + Mhc - Mglycltp + Mglycltc = 0
 R("R_GLYCLTt2rpp_Rev"): + Mhp - Mhc + Mglycltp - Mglycltc = 0
 R("R_GLYCLTt4pp"): - Mna1p + Mna1c - Mglycltp + Mglycltc = 0
 R("R_GLYCLTtex"): + Mglycltp - Mglyclte = 0
 R("R_GLYCLTtex_Rev"): - Mglycltp + Mglyclte = 0
 R("MglycRb_Transp"): 0 M10fthfc = 0
 R("R_EX_glyc_R_e_"): - MglycRe = 0
 R("R_EX_glyc_R_e__Rev"): + MglycRe = 0
 R("R_GLYCAt2rpp"): - Mhp + Mhc - MglycRp + MglycRc = 0
 R("R_GLYCAt2rpp_Rev"): + Mhp - Mhc + MglycRp - MglycRc = 0
 R("R_HPYRRx"): - Mnadhc + Mnadc - Mhpyrc - Mhc + MglycRc = 0
 R("R_HPYRRy"): - Mnadphc + Mnadpc - Mhpyrc - Mhc + MglycRc = 0
 R("R_MANPGH"): - Mman6pglycc + Mman6pc - Mh2oc + MglycRc = 0
 R("R_GLYCAtex"): + MglycRp - MglycRe = 0
 R("R_GLYCAtex_Rev"): - MglycRp + MglycRe = 0
 R("R_GLYtex"): + Mglyp - Mglye = 0
 R("R_GLYtex_Rev"): - Mglyp + Mglye = 0
 R("R_GMHEPPA"): + Mpic - Mh2oc + Mgmhep1pc - Mgmhep17bpc = 0
 R("R_S7PI"): - Ms7pc + Mgmhep7pc = 0
 R("Mgmpb_Transp"): 0 M10fthfc = 0
 R("R_EX_gmp_e_"): - Mgmpe = 0
 R("R_EX_gmp_e__Rev"): + Mgmpe = 0
 R("R_GMPR"): + Mnh4c - Mnadphc + Mnadpc + Mimpc - 2 Mhc - Mgmpc = 0
 R("R_GUAPRT"): - Mprppc + Mppic - Mguac + Mgmpc = 0
 R("R_NTD9"): + Mpic - Mh2oc + Mgsnc - Mgmpc = 0
 R("R_NTPP2"): + Mppic + Mhc - Mh2oc - Mgtpc + Mgmpc = 0
 R("R_GMPtex"): + Mgmpp - Mgmpe = 0
 R("R_GMPtex_Rev"): - Mgmpp + Mgmpe = 0
 R("R_NTD9pp"): + Mpip - Mh2op + Mgsnp - Mgmpp = 0
 R("R_GRXR"): - 2 Mgthrdc + Mgthoxc + Mgrxrdc - Mgrxoxc = 0
 R("R_PAPSR2"): + Mso3c - Mpapsc + Mpapc + 2 Mhc - Mgrxrdc + Mgrxoxc = 0
 R("Mgsnb_Transp"): 0 M10fthfc = 0
 R("R_EX_gsn_e_"): - Mgsne = 0
 R("R_EX_gsn_e__Rev"): + Mgsne = 0
 R("R_GSNt2pp"): - Mhp + Mhc - Mgsnp + Mgsnc = 0
 R("R_NTPTP2"): + Mpppic - Mh2oc - Mgtpc + Mgsnc = 0
 R("R_PUNP3"): + Mr1pc - Mpic + Mguac - Mgsnc = 0
 R("R_PUNP3_Rev"): - Mr1pc + Mpic - Mguac + Mgsnc = 0
 R("R_GSNtex"): + Mgsnp - Mgsne = 0
 R("R_GSNtex_Rev"): - Mgsnp + Mgsne = 0
 R("Mgthoxb_Transp"): 0 M10fthfc = 0
 R("R_EX_gthox_e_"): - Mgthoxe = 0
 R("R_EX_gthox_e__Rev"): + Mgthoxe = 0
 R("R_GTHOr"): - Mnadphc + Mnadpc - Mhc + 2 Mgthrdc - Mgthoxc = 0
 R("R_GTHOr_Rev"): + Mnadphc - Mnadpc + Mhc - 2 Mgthrdc + Mgthoxc = 0
 R("R_GTHOXtex"): + Mgthoxp - Mgthoxe = 0
 R("R_GTHOXtex_Rev"): - Mgthoxp + Mgthoxe = 0
 R("Mgthrdb_Transp"): 0 M10fthfc = 0
 R("R_EX_gthrd_e_"): - Mgthrde = 0
 R("R_EX_gthrd_e__Rev"): + Mgthrde = 0
 R("R_GLYOX"): - MlgtSc + MlacDc + Mhc - Mh2oc + Mgthrdc = 0
 R("R_GSPMDA"): + Mspmdc - Mh2oc - Mgtspmdc + Mgthrdc = 0
 R("R_LGTHL"): - Mmthgxlc + MlgtSc - Mgthrdc = 0
 R("R_GTHRDtex"): + Mgthrdp - Mgthrde = 0
 R("R_GTHRDtex_Rev"): - Mgthrdp + Mgthrde = 0
 R("Mgtpb_Transp"): 0 M10fthfc = 0
 R("R_EX_gtp_e_"): - Mgtpe = 0
 R("R_EX_gtp_e__Rev"): + Mgtpe = 0
 R("R_GTPHs"): + Mxtpc + Mnh4c - Mhc - Mh2oc - Mgtpc = 0
 R("R_GTPtex"): + Mgtpp - Mgtpe = 0
 R("R_GTPtex_Rev"): - Mgtpp + Mgtpe = 0
 R("Mguab_Transp"): 0 M10fthfc = 0
 R("R_EX_gua_e_"): - Mguae = 0
 R("R_EX_gua_e__Rev"): + Mguae = 0
 R("R_GUAD"): + Mxanc + Mnh4c - Mhc - Mh2oc - Mguac = 0
 R("R_GUAt2pp"): - Mhp + Mhc - Mguap + Mguac = 0
 R("R_GUAtpp"): - Mguap + Mguac = 0
 R("R_GUAtpp_Rev"): + Mguap - Mguac = 0
 R("R_GUAtex"): + Mguap - Mguae = 0
 R("R_GUAtex_Rev"): - Mguap + Mguae = 0
 R("Mh2b_Transp"): 0 M10fthfc = 0
 R("R_EX_h2_e_"): - Mh2e = 0
 R("R_EX_h2_e__Rev"): + Mh2e = 0
 R("R_H2tpp"): - Mh2p + Mh2c = 0
 R("R_H2tpp_Rev"): + Mh2p - Mh2c = 0
 R("R_HYD1pp"): + Mq8h2c - Mq8c + 2 Mhp - 2 Mhc - Mh2c = 0
 R("R_HYD2pp"): - Mmqn8c + Mmql8c + 2 Mhp - 2 Mhc - Mh2c = 0
 R("R_H2tex"): + Mh2p - Mh2e = 0
 R("R_H2tex_Rev"): - Mh2p + Mh2e = 0
 R("R_IPDPS"): - Mnadhc + Mnadc + Mipdpc - Mhc + Mh2oc - Mh2mb4pc = 0
 R("Mh2o2b_Transp"): 0 M10fthfc = 0
 R("R_EX_h2o2_e_"): - Mh2o2e = 0
 R("R_EX_h2o2_e__Rev"): + Mh2o2e = 0
 R("R_CAT"): + Mo2c + 2 Mh2oc - 2 Mh2o2c = 0
 R("R_METOX1s"): + MmetsoxSLc - MmetLc + Mh2oc - Mh2o2c = 0
 R("R_METOX2s"): + MmetsoxRLc - MmetLc + Mh2oc - Mh2o2c = 0
 R("R_PDX5POi"): + Mpydx5pc - Mpdx5pc - Mo2c + Mh2o2c = 0
 R("R_PYAM5PO"): + Mpydx5pc - Mpyam5pc - Mo2c + Mnh4c - Mh2oc + Mh2o2c
 = 0
 R("R_SPODM"): - 2 Mo2sc + Mo2c - 2 Mhc + Mh2o2c = 0
 R("R_THIORDXi"): - Mtrdrdc + Mtrdoxc + 2 Mh2oc - Mh2o2c = 0
 R("R_H2O2tex"): + Mh2o2p - Mh2o2e = 0
 R("R_H2O2tex_Rev"): - Mh2o2p + Mh2o2e = 0
 R("R_PEAMNOpp"): - Mpeamnp + Mpacaldp - Mo2p + Mnh4p - Mh2op + Mh2o2p
 = 0
 R("R_SPODMpp"): - 2 Mo2sp + Mo2p - 2 Mhp + Mh2o2p = 0
 R("Mh2ob_Transp"): 0 M10fthfc = 0
 R("R_EX_h2o_e_"): - Mh2oe = 0
 R("R_EX_h2o_e__Rev"): + Mh2oe = 0
 R("R_ALDD19x"): + Mpacc - Mpacaldc + Mnadhc - Mnadc + 2 Mhc - Mh2oc = 0
 R("R_ALDD3y"): - Mppalc + Mppac + Mnadphc - Mnadpc + 2 Mhc - Mh2oc = 0
 R("R_CYTBD2pp"): - 0.5 Mo2c + Mmqn8c - Mmql8c + 2 Mhp - 2 Mhc + Mh2oc
 = 0
 R("R_CYTBDpp"): - Mq8h2c + Mq8c - 0.5 Mo2c + 2 Mhp - 2 Mhc + Mh2oc = 0
 R("R_CYTBO3_4pp"): - Mq8h2c + Mq8c - 0.5 Mo2c + 4 Mhp - 4 Mhc + Mh2oc
 = 0
 R("R_GLYOX3"): - Mmthgxlc + MlacDc + Mhc - Mh2oc = 0
 R("R_H2Otpp"): - Mh2op + Mh2oc = 0
 R("R_H2Otpp_Rev"): + Mh2op - Mh2oc = 0
 R("R_HISTD"): + 2 Mnadhc - 2 Mnadc - Mhistdc + MhisLc + 3 Mhc - Mh2oc
 = 0
 R("R_HISTP"): + Mpic + Mhistdc - Mhispc - Mh2oc = 0
 R("R_HKNDDH"): + Msuccc + Mop4enc - Mhknddc + Mhc - Mh2oc = 0
 R("R_HMBS"): - 4 Mppbngc + 4 Mnh4c + Mhmbilc - Mh2oc = 0
 R("R_HXAND"): + Mxanc + Mnadhc - Mnadc - Mhxanc + Mhc - Mh2oc = 0
 R("R_HYPOE"): + Mpydamc - Mpyam5pc + Mpic - Mh2oc = 0
 R("R_IMPD"): + Mxmpc + Mnadhc - Mnadc - Mimpc + Mhc - Mh2oc = 0
 R("R_INSH"): + MribDc - Minsc + Mhxanc - Mh2oc = 0
 R("R_KDOPP"): + Mpic + Mkdoc - Mkdo8pc - Mh2oc = 0
 R("R_LCADi"): + Mnadhc - Mnadc - MlaldLc + MlacLc + 2 Mhc - Mh2oc = 0
 R("R_METSOXR1"): - Mtrdrdc + Mtrdoxc - MmetsoxSLc + MmetLc + Mh2oc = 0
 R("R_METSOXR2"): - Mtrdrdc + Mtrdoxc - MmetsoxRLc + MmetLc + Mh2oc = 0
 R("R_MI1PP"): + Mpic - Mmi1pDc + Minostc - Mh2oc = 0
 R("R_MN6PP"): + Mpic + Mmanc - Mman6pc - Mh2oc = 0
 R("R_NADPPPS"): + Mpic - Mnadpc + Mnadc - Mh2oc = 0
 R("R_NHFRBO"): - 2 Mnoc - Mnadhc + Mnadc + Mn2oc - Mhc + Mh2oc = 0
 R("R_NMNDA"): - Mnmnc + Mnicrntc + Mnh4c - Mh2oc = 0
 R("R_NMNN"): + Mr5pc - Mnmnc + Mncamc + Mhc - Mh2oc = 0
 R("R_NMNt7pp"): + Mr5pc - Mnmnp + Mncamc + Mhc - Mh2oc = 0
 R("R_NNAM"): + Mnh4c - Mncamc + Mnacc - Mh2oc = 0
 R("R_NO3R1pp"): - Mq8h2c + Mq8c - Mno3c + Mno2c + 2 Mhp - 2 Mhc + Mh2oc
 = 0
 R("R_NO3R2pp"): - Mno3c + Mno2c + Mmqn8c - Mmql8c + 2 Mhp - 2 Mhc
 + Mh2oc = 0
 R("R_NTD10"): + Mxtsnc - Mxmpc + Mpic - Mh2oc = 0
 R("R_NTD11"): + Mpic + Minsc - Mimpc - Mh2oc = 0
 R("R_NTD2"): + Muric - Mumpc + Mpic - Mh2oc = 0
 R("R_NTP10"): + Mpic - Mitpc + Midpc + Mhc - Mh2oc = 0
 R("R_NTP12"): - Mxtpc + Mxdpc + Mpic + Mhc - Mh2oc = 0
 R("R_NTPP11"): - Mxtpc + Mxmpc + Mppic + Mhc - Mh2oc = 0
 R("R_NTPP8"): - Mutpc + Mumpc + Mppic + Mhc - Mh2oc = 0
 R("R_NTPP9"): + Mppic - Mitpc + Mimpc + Mhc - Mh2oc = 0
 R("R_NTRIR2x"): - Mno2c + Mnh4c - 3 Mnadhc + 3 Mnadc - 5 Mhc + 2 Mh2oc
 = 0
 R("R_PDXPP"): + Mpydxnc + Mpic - Mpdx5pc - Mh2oc = 0
 R("R_PGPP120"): + Mpic - Mpgp120c + Mpg120c - Mh2oc = 0
 R("R_PGPP140"): + Mpic - Mpgp140c + Mpg140c - Mh2oc = 0
 R("R_PGPP141"): + Mpic - Mpgp141c + Mpg141c - Mh2oc = 0
 R("R_PGPP160"): + Mpic - Mpgp160c + Mpg160c - Mh2oc = 0
 R("R_PGPP161"): + Mpic - Mpgp161c + Mpg161c - Mh2oc = 0
 R("R_PGPP180"): + Mpic - Mpgp180c + Mpg180c - Mh2oc = 0
 R("R_PGPP181"): + Mpic - Mpgp181c + Mpg181c - Mh2oc = 0
 R("R_PPA"): - Mppic + 2 Mpic + Mhc - Mh2oc = 0
 R("R_PPA2"): - Mpppic + Mppic + Mpic + Mhc - Mh2oc = 0
 R("R_PPPGO"): - Mpppg9c + Mppp9c - 1.5 Mo2c + 3 Mh2oc = 0
 R("R_PRAMPC"): + Mprfpc - Mprbampc - Mh2oc = 0
 R("R_PRATPP"): - Mprbatpc + Mprbampc + Mppic + Mhc - Mh2oc = 0
 R("R_PSP_L"): + MserLc - MpserLc + Mpic - Mh2oc = 0
 R("R_PYDXPP"): + Mpydxc - Mpydx5pc + Mpic - Mh2oc = 0
 R("R_R5PP"): + MribDc - Mr5pc + Mpic - Mh2oc = 0
 R("R_SGSAD"): - Msucgsac + Msucgluc + Mnadhc - Mnadc + 2 Mhc - Mh2oc
 = 0
 R("R_SSALx"): - Msucsalc + Msuccc + Mnadhc - Mnadc + 2 Mhc - Mh2oc = 0
 R("R_SSALy"): - Msucsalc + Msuccc + Mnadphc - Mnadpc + 2 Mhc - Mh2oc
 = 0
 R("R_SULRi"): - Mso3c - 3 Mnadphc + 3 Mnadpc - 5 Mhc + Mh2sc + 3 Mh2oc
 = 0
 R("R_TARTD"): - MtartrLc + Moaac + Mh2oc = 0
 R("R_TDP"): - Mthmppc + Mthmmpc + Mpic + Mhc - Mh2oc = 0
 R("R_THRS"): + MthrLc + Mpic - Mphomc - Mh2oc = 0
 R("R_TMAOR1"): - Mtmaoc + Mtmac + Mmqn8c - Mmql8c - Mhc + Mh2oc = 0
 R("R_TRE6PP"): + Mtrec - Mtre6pc + Mpic - Mh2oc = 0
 R("R_TRPAS2"): - MtrpLc + Mpyrc + Mnh4c + Mindolec - Mh2oc = 0
 R("R_TRPAS2_Rev"): + MtrpLc - Mpyrc - Mnh4c - Mindolec + Mh2oc = 0
 R("R_TRPS2"): + MtrpLc - MserLc - Mindolec + Mh2oc = 0
 R("R_UACMAMO"): + Muacmamuc - Muacmamc + 2 Mnadhc - 2 Mnadc + 3 Mhc
 - Mh2oc = 0
 R("R_UDCPDP"): + Mudcppc - Mudcpdpc + Mpic + Mhc - Mh2oc = 0
 R("R_UDPGD"): + Mudpglcurc - Mudpgc + 2 Mnadhc - 2 Mnadc + 3 Mhc
 - Mh2oc = 0
 R("R_UPP3S"): + Muppg3c - Mhmbilc + Mh2oc = 0
 R("R_URIH"): - Muric + Murac + MribDc - Mh2oc = 0
 R("R_USHD"): + Mumpc - Mu23gac + MlipidXc + 2 Mhc - Mh2oc = 0
 R("R_XAND"): - Mxanc + Muratec + Mnadhc - Mnadc + Mhc - Mh2oc = 0
 R("R_XTSNH"): - Mxtsnc + Mxanc + MribDc - Mh2oc = 0
 R("R_H2Otex"): + Mh2op - Mh2oe = 0
 R("R_H2Otex_Rev"): - Mh2op + Mh2oe = 0
 R("R_LIPAHTex"): - Mlipae - Mhe - Mhdcae + Mhalipae + Mh2oe = 0
 R("R_MDDEP1pp"): - Mmurein4px4pp + Mmurein4p4pp - Mh2op = 0
 R("R_MDDEP2pp"): + Mmurein4p3pp - Mmurein3px4pp - Mh2op = 0
 R("R_MDDEP3pp"): - Mmurein5px4pp + Mmurein5p4pp - Mh2op = 0
 R("R_MDDEP4pp"): - Mmurein4px4px4pp + Mmurein4px4p4pp - Mh2op = 0
 R("R_MLDEP1pp"): - Mmurein3px3pp + Mmurein3p3pp - Mh2op = 0
 R("R_MLDEP2pp"): - Mmurein5px3pp + Mmurein5p3pp - Mh2op = 0
 R("R_NO3R1bpp"): - Mq8h2c + Mq8c - Mno3p + Mno2p + Mh2op = 0
 R("R_NO3R2bpp"): - Mno3p + Mno2p + Mmqn8c - Mmql8c + Mh2op = 0
 R("R_NTD10pp"): + Mxtsnp - Mxmpp + Mpip - Mh2op = 0
 R("R_NTD11pp"): + Mpip + Minsp - Mimpp - Mh2op = 0
 R("R_NTD2pp"): + Murip - Mumpp + Mpip - Mh2op = 0
 R("R_NTRIR3pp"): - 3 Mq8h2c + 3 Mq8c - Mno2p + Mnh4p - 2 Mhp + 2 Mh2op
 = 0
 R("R_NTRIR4pp"): - Mno2p + Mnh4p + 3 Mmqn8c - 3 Mmql8c - 2 Mhp
 + 2 Mh2op = 0
 R("R_PGPP120pp"): + Mpip - Mpgp120p + Mpg120p - Mh2op = 0
 R("R_PGPP140pp"): + Mpip - Mpgp140p + Mpg140p - Mh2op = 0
 R("R_PGPP141pp"): + Mpip - Mpgp141p + Mpg141p - Mh2op = 0
 R("R_PGPP160pp"): + Mpip - Mpgp160p + Mpg160p - Mh2op = 0
 R("R_PGPP161pp"): + Mpip - Mpgp161p + Mpg161p - Mh2op = 0
 R("R_PGPP180pp"): + Mpip - Mpgp180p + Mpg180p - Mh2op = 0
 R("R_PGPP181pp"): + Mpip - Mpgp181p + Mpg181p - Mh2op = 0
 R("R_PHYTSpp"): + 6 Mpip - Mminohpp + Minostp - 6 Mh2op = 0
 R("R_PPTHpp"): - Mpptp + Mpip + Mh2p - Mh2op = 0
 R("R_PSP_Lpp"): + MserLp - MpserLp + Mpip - Mh2op = 0
 R("R_PTHRpp"): - Mthrpp + MthrLp + Mpip - Mh2op = 0
 R("R_R5PPpp"): + MribDp - Mr5pp + Mpip - Mh2op = 0
 R("R_TMAOR1pp"): + Mtmap - Mtmaop + Mmqn8c - Mmql8c - Mhp + Mh2op = 0
 R("R_TYRPpp"): - Mtyrpp + MtyrLp + Mpip - Mh2op = 0
 R("R_UDCPDPpp"): + Mudcppp - Mudcpdpp + Mpip + Mhp - Mh2op = 0
 R("Mh2sb_Transp"): 0 M10fthfc = 0
 R("R_EX_h2s_e_"): - Mh2se = 0
 R("R_EX_h2s_e__Rev"): + Mh2se = 0
 R("R_H2SO"): + Mso4c - 2 Mo2c + 2 Mhc - Mh2sc = 0
 R("R_H2St1pp"): + Mh2sp - Mh2sc = 0
 R("R_H2Stex"): + Mh2sp - Mh2se = 0
 R("R_H2Stex_Rev"): - Mh2sp + Mh2se = 0
 R("Mhacolipab_Transp"): 0 M10fthfc = 0
 R("R_EX_hacolipa_e_"): - Mhacolipae = 0
 R("R_EX_hacolipa_e__Rev"): + Mhacolipae = 0
 R("Mhalipab_Transp"): 0 M10fthfc = 0
 R("R_EX_halipa_e_"): - Mhalipae = 0
 R("R_EX_halipa_e__Rev"): + Mhalipae = 0
 R("Mhb_Transp"): 0 M10fthfc = 0
 R("R_EX_h_e_"): - Mhe = 0
 R("R_EX_h_e__Rev"): + Mhe = 0
 R("R_ACMAMUT"): + Munagamuc - Munagac + Mudpc - Muacmamuc + Mhc = 0
 R("R_D_LACt2pp"): - MlacDp + MlacDc - Mhp + Mhc = 0
 R("R_D_LACt2pp_Rev"): + MlacDp - MlacDc + Mhp - Mhc = 0
 R("R_DSERt2pp"): - MserDp + MserDc - Mhp + Mhc = 0
 R("R_EAR140x"): - Mtmrs2eACPc - Mnadhc + Mnadc + MmyrsACPc - Mhc = 0
 R("R_EAR140y"): - Mtmrs2eACPc - Mnadphc + Mnadpc + MmyrsACPc - Mhc = 0
 R("R_EAR141x"): + MtdeACPc - Mt3c7mrseACPc - Mnadhc + Mnadc - Mhc = 0
 R("R_EAR141y"): + MtdeACPc - Mt3c7mrseACPc - Mnadphc + Mnadpc - Mhc = 0
 R("R_EAR160x"): - Mtpalm2eACPc + MpalmACPc - Mnadhc + Mnadc - Mhc = 0
 R("R_EAR160y"): - Mtpalm2eACPc + MpalmACPc - Mnadphc + Mnadpc - Mhc = 0
 R("R_EAR161x"): - Mt3c9palmeACPc - Mnadhc + Mnadc + MhdeACPc - Mhc = 0
 R("R_EAR161y"): - Mt3c9palmeACPc - Mnadphc + Mnadpc + MhdeACPc - Mhc
 = 0
 R("R_EAR180x"): - Mtoctd2eACPc + MocdcaACPc - Mnadhc + Mnadc - Mhc = 0
 R("R_EAR180y"): - Mtoctd2eACPc + MocdcaACPc - Mnadphc + Mnadpc - Mhc
 = 0
 R("R_EAR181x"): - Mt3c11vaceACPc + MocteACPc - Mnadhc + Mnadc - Mhc = 0
 R("R_EAR181y"): - Mt3c11vaceACPc + MocteACPc - Mnadphc + Mnadpc - Mhc
 = 0
 R("R_EAR60x"): - Mthex2eACPc - Mnadhc + Mnadc + MhexACPc - Mhc = 0
 R("R_EAR60y"): - Mthex2eACPc - Mnadphc + Mnadpc + MhexACPc - Mhc = 0
 R("R_EAR80x"): - Mtoct2eACPc + MocACPc - Mnadhc + Mnadc - Mhc = 0
 R("R_EAR80y"): - Mtoct2eACPc + MocACPc - Mnadphc + Mnadpc - Mhc = 0
 R("R_FALDH2"): + MSfglutthc + Mnadhc - Mnadc - Mhmgthc + Mhc = 0
 R("R_FALDH2_Rev"): - MSfglutthc - Mnadhc + Mnadc + Mhmgthc - Mhc = 0
 R("R_FLVR"): - Mribflvc + Mrbflvrdc - Mnadphc + Mnadpc - Mhc = 0
 R("R_FLVRx"): - Mribflvc + Mrbflvrdc - Mnadhc + Mnadc - Mhc = 0
 R("R_HCYSMT2"): - Mmmetc + 2 MmetLc - MhcysLc + Mhc = 0
 R("R_HEXt2rpp"): - Mhxap + Mhxac - Mhp + Mhc = 0
 R("R_HEXt2rpp_Rev"): + Mhxap - Mhxac + Mhp - Mhc = 0
 R("R_HG2t3pp"): - Mhp + Mhg2p - Mhg2c + Mhc = 0
 R("R_HISt2rpp"): - Mhp - MhisLp + MhisLc + Mhc = 0
 R("R_HISt2rpp_Rev"): + Mhp + MhisLp - MhisLc - Mhc = 0
 R("R_HOMt2pp"): - Mhp + MhomLp - MhomLc + Mhc = 0
 R("R_IDONt2rpp"): - MidonLp + MidonLc - Mhp + Mhc = 0
 R("R_IDONt2rpp_Rev"): + MidonLp - MidonLc + Mhp - Mhc = 0
 R("R_ILEt2rpp"): - MileLp + MileLc - Mhp + Mhc = 0
 R("R_ILEt2rpp_Rev"): + MileLp - MileLc + Mhp - Mhc = 0
 R("R_INDOLEt2pp"): + Mindolep - Mindolec + Mhp - Mhc = 0
 R("R_INDOLEt2rpp"): - Mindolep + Mindolec - Mhp + Mhc = 0
 R("R_INDOLEt2rpp_Rev"): + Mindolep - Mindolec + Mhp - Mhc = 0
 R("R_INSt2pp"): - Minsp + Minsc - Mhp + Mhc = 0
 R("R_INSt2rpp"): - Minsp + Minsc - Mhp + Mhc = 0
 R("R_INSt2rpp_Rev"): + Minsp - Minsc + Mhp - Mhc = 0
 R("R_Kt2pp"): - Mkp + Mkc - Mhp + Mhc = 0
 R("R_Kt3pp"): + Mkp - Mkc - Mhp + Mhc = 0
 R("R_LALDO2x"): - Mnadhc + Mnadc - Mmthgxlc + MlaldDc - Mhc = 0
 R("R_LCTSt3ipp"): + Mlctsp - Mlctsc - Mhp + Mhc = 0
 R("R_LCTStpp"): - Mlctsp + Mlctsc - Mhp + Mhc = 0
 R("R_LCTStpp_Rev"): + Mlctsp - Mlctsc + Mhp - Mhc = 0
 R("R_LDH_D"): + Mpyrc + Mnadhc - Mnadc - MlacDc + Mhc = 0
 R("R_LDH_D_Rev"): - Mpyrc - Mnadhc + Mnadc + MlacDc - Mhc = 0
 R("R_LEUt2rpp"): - MleuLp + MleuLc - Mhp + Mhc = 0
 R("R_LEUt2rpp_Rev"): + MleuLp - MleuLc + Mhp - Mhc = 0
 R("R_L_LACt2rpp"): - MlacLp + MlacLc - Mhp + Mhc = 0
 R("R_L_LACt2rpp_Rev"): + MlacLp - MlacLc + Mhp - Mhc = 0
 R("R_LPADSS"): + Mudpc - Mu23gac - MlipidXc + MlipidAdsc + Mhc = 0
 R("R_LYSt2pp"): - MlysLp + MlysLc - Mhp + Mhc = 0
 R("R_LYSt3pp"): + MlysLp - MlysLc - Mhp + Mhc = 0
 R("R_LYXt2pp"): - MlyxLp + MlyxLc - Mhp + Mhc = 0
 R("R_MALDt2_2pp"): - MmalDp + MmalDc - 2 Mhp + 2 Mhc = 0
 R("R_MALt2_2pp"): - MmalLp + MmalLc - 2 Mhp + 2 Mhc = 0
 R("R_MALt2_3pp"): - MmalLp + MmalLc - 3 Mhp + 3 Mhc = 0
 R("R_MDH"): + Moaac + Mnadhc - Mnadc - MmalLc + Mhc = 0
 R("R_MDH_Rev"): - Moaac - Mnadhc + Mnadc + MmalLc - Mhc = 0
 R("R_MELIBt2pp"): - Mmelibp + Mmelibc - Mhp + Mhc = 0
 R("R_MELIBt3ipp"): + Mmelibp - Mmelibc - Mhp + Mhc = 0
 R("R_MG2t3_2pp"): - Mmg2p + Mmg2c + 2 Mhp - 2 Mhc = 0
 R("R_MG2t3_2pp_Rev"): + Mmg2p - Mmg2c - 2 Mhp + 2 Mhc = 0
 R("R_MMETt2pp"): - Mmmetp + Mmmetc - Mhp + Mhc = 0
 R("R_MN2t3pp"): + Mmn2p - Mmn2c - Mhp + Mhc = 0
 R("R_MNt2pp"): - Mmn2p + Mmn2c - Mhp + Mhc = 0
 R("R_MPTG"): + 2 Mudcpdpc - 2 Muaagmdac + Mmurein5p5pp + 2 Mhc = 0
 R("R_MPTG2"): + Mudcpdpc - Muaagmdac - Mmurein5p5pp + Mmurein5p5p5pp
 + Mhc = 0
 R("R_NADH10"): - Mnadhc + Mnadc - Mmqn8c + Mmql8c - Mhc = 0
 R("R_NADH16pp"): + Mq8h2c - Mq8c - Mnadhc + Mnadc + 3 Mhp - 4 Mhc = 0
 R("R_NADH17pp"): - Mnadhc + Mnadc - Mmqn8c + Mmql8c + 3 Mhp - 4 Mhc = 0
 R("R_NADH5"): + Mq8h2c - Mq8c - Mnadhc + Mnadc - Mhc = 0
 R("R_NADPHQR2"): + Mq8h2c - Mq8c - Mnadphc + Mnadpc - Mhc = 0
 R("R_NADPHQR3"): - Mnadphc + Mnadpc - Mmqn8c + Mmql8c - Mhc = 0
 R("R_NAt3_1p5pp"): + 2 Mna1p - 2 Mna1c - 3 Mhp + 3 Mhc = 0
 R("R_NAt3_2pp"): + Mna1p - Mna1c - 2 Mhp + 2 Mhc = 0
 R("R_NAt3pp"): + Mna1p - Mna1c - Mhp + Mhc = 0
 R("R_NI2t3pp"): + Mni2p - Mni2c - Mhp + Mhc = 0
 R("R_NO2t2rpp"): - Mno2p + Mno2c - Mhp + Mhc = 0
 R("R_NO2t2rpp_Rev"): + Mno2p - Mno2c + Mhp - Mhc = 0
 R("R_NODOx"): - 2 Mo2c - 2 Mnoc + 2 Mno3c - Mnadhc + Mnadc + Mhc = 0
 R("R_NODOy"): - 2 Mo2c - 2 Mnoc + 2 Mno3c - Mnadphc + Mnadpc + Mhc = 0
 R("R_OROTt2_2pp"): - Morotp + Morotc - 2 Mhp + 2 Mhc = 0
 R("R_PACALDt2rpp"): - Mpacaldp + Mpacaldc - Mhp + Mhc = 0
 R("R_PACALDt2rpp_Rev"): + Mpacaldp - Mpacaldc + Mhp - Mhc = 0
 R("R_PAPSR"): - Mtrdrdc + Mtrdoxc + Mso3c - Mpapsc + Mpapc + 2 Mhc = 0
 R("R_PDX5PO2"): + Mpydx5pc - Mpdx5pc + Mnadhc - Mnadc + Mhc = 0
 R("R_PHEt2rpp"): - MpheLp + MpheLc - Mhp + Mhc = 0
 R("R_PHEt2rpp_Rev"): + MpheLp - MpheLc + Mhp - Mhc = 0
 R("R_PIt2rpp"): - Mpip + Mpic - Mhp + Mhc = 0
 R("R_PIt2rpp_Rev"): + Mpip - Mpic + Mhp - Mhc = 0
 R("R_PPPNt2rpp"): - Mpppnp + Mpppnc - Mhp + Mhc = 0
 R("R_PPPNt2rpp_Rev"): + Mpppnp - Mpppnc + Mhp - Mhc = 0
 R("R_PROt2rpp"): - MproLp + MproLc - Mhp + Mhc = 0
 R("R_PROt2rpp_Rev"): + MproLp - MproLc + Mhp - Mhc = 0
 R("R_PSCLYSt2pp"): - Mpsclysp + Mpsclysc - Mhp + Mhc = 0
 R("R_PTRCt2pp"): - Mptrcp + Mptrcc - Mhp + Mhc = 0
 R("R_PYRt2rpp"): - Mpyrp + Mpyrc - Mhp + Mhc = 0
 R("R_PYRt2rpp_Rev"): + Mpyrp - Mpyrc + Mhp - Mhc = 0
 R("R_QMO2"): - Mq8h2c + Mq8c + 2 Mo2sc - 2 Mo2c + 2 Mhc = 0
 R("R_QMO3"): + 2 Mo2sc - 2 Mo2c + Mmqn8c - Mmql8c + 2 Mhc = 0
 R("R_RMNtpp"): - Mrmnp + Mrmnc - Mhp + Mhc = 0
 R("R_SELCYSS"): - Msertrnasecoooc - Mselnpc + Msectrnac + Mpic + Mhc
 = 0
 R("R_SERt2rpp"): - MserLp + MserLc - Mhp + Mhc = 0
 R("R_SERt2rpp_Rev"): + MserLp - MserLc + Mhp - Mhc = 0
 R("R_SKMt2pp"): - Mskmp + Mskmc - Mhp + Mhc = 0
 R("R_SUCCt2_2pp"): - Msuccp + Msuccc - 2 Mhp + 2 Mhc = 0
 R("R_SUCCt2_3pp"): - Msuccp + Msuccc - 3 Mhp + 3 Mhc = 0
 R("R_SUCCt3pp"): + Msuccp - Msuccc - Mhp + Mhc = 0
 R("R_THD2pp"): + Mnadphc - Mnadpc - Mnadhc + Mnadc - 2 Mhp + 2 Mhc = 0
 R("R_THMDt2pp"): - Mthymdp + Mthymdc - Mhp + Mhc = 0
 R("R_THMDt2rpp"): - Mthymdp + Mthymdc - Mhp + Mhc = 0
 R("R_THMDt2rpp_Rev"): + Mthymdp - Mthymdc + Mhp - Mhc = 0
 R("R_THRt2pp"): + MthrLp - MthrLc - Mhp + Mhc = 0
 R("R_THRt2rpp"): - MthrLp + MthrLc - Mhp + Mhc = 0
 R("R_THRt2rpp_Rev"): + MthrLp - MthrLc + Mhp - Mhc = 0
 R("R_THYMt3pp"): + Mthymp - Mthymc - Mhp + Mhc = 0
 R("R_TRDR"): + Mtrdrdc - Mtrdoxc - Mnadphc + Mnadpc - Mhc = 0
 R("R_TRPt2rpp"): - MtrpLp + MtrpLc - Mhp + Mhc = 0
 R("R_TRPt2rpp_Rev"): + MtrpLp - MtrpLc + Mhp - Mhc = 0
 R("R_TYRt2rpp"): - MtyrLp + MtyrLc - Mhp + Mhc = 0
 R("R_TYRt2rpp_Rev"): + MtyrLp - MtyrLc + Mhp - Mhc = 0
 R("R_UAGPT3"): + Mudpc - Muagmdac - Muacgamc + Muaagmdac + Mhc = 0
 R("R_UAPGR"): + Muamrc - Muaccgc - Mnadphc + Mnadpc - Mhc = 0
 R("R_URAt2pp"): - Murap + Murac - Mhp + Mhc = 0
 R("R_URAt2rpp"): - Murap + Murac - Mhp + Mhc = 0
 R("R_URAt2rpp_Rev"): + Murap - Murac + Mhp - Mhc = 0
 R("R_URDGLYCD"): - Murdglycc + Moxurc + Mnadhc - Mnadc + Mhc = 0
 R("R_URIt2pp"): - Murip + Muric - Mhp + Mhc = 0
 R("R_URIt2rpp"): - Murip + Muric - Mhp + Mhc = 0
 R("R_URIt2rpp_Rev"): + Murip - Muric + Mhp - Mhc = 0
 R("R_VALt2rpp"): - MvalLp + MvalLc - Mhp + Mhc = 0
 R("R_VALt2rpp_Rev"): + MvalLp - MvalLc + Mhp - Mhc = 0
 R("R_XANt2pp"): - Mxanp + Mxanc - Mhp + Mhc = 0
 R("R_XTSNt2rpp"): - Mxtsnp + Mxtsnc - Mhp + Mhc = 0
 R("R_XTSNt2rpp_Rev"): + Mxtsnp - Mxtsnc + Mhp - Mhc = 0
 R("R_XYLt2pp"): - MxylDp + MxylDc - Mhp + Mhc = 0
 R("R_XYLUt2pp"): - MxyluLp + MxyluLc - Mhp + Mhc = 0
 R("R_ZN2t3pp"): + Mzn2p - Mzn2c - Mhp + Mhc = 0
 R("Mhdcab_Transp"): 0 M10fthfc = 0
 R("R_EX_hdca_e_"): - Mhdcae = 0
 R("R_EX_hdca_e__Rev"): + Mhdcae = 0
 R("R_HDCAtexi"): + Mhdcap - Mhdcae = 0
 R("Mhdceab_Transp"): 0 M10fthfc = 0
 R("R_EX_hdcea_e_"): - Mhdceae = 0
 R("R_EX_hdcea_e__Rev"): + Mhdceae = 0
 R("R_HDCEAtexi"): + Mhdceap - Mhdceae = 0
 R("R_HDCOAI"): + Mhdd2coac - Mhdcoac = 0
 R("R_Htex"): + Mhp - Mhe = 0
 R("R_Htex_Rev"): - Mhp + Mhe = 0
 R("Mhg2b_Transp"): 0 M10fthfc = 0
 R("R_EX_hg2_e_"): - Mhg2e = 0
 R("R_EX_hg2_e__Rev"): + Mhg2e = 0
 R("R_HG2tex"): + Mhg2p - Mhg2e = 0
 R("R_HG2tex_Rev"): - Mhg2p + Mhg2e = 0
 R("MhisLb_Transp"): 0 M10fthfc = 0
 R("R_EX_his_L_e_"): - MhisLe = 0
 R("R_EX_his_L_e__Rev"): + MhisLe = 0
 R("R_HIStex"): + MhisLp - MhisLe = 0
 R("R_HIStex_Rev"): - MhisLp + MhisLe = 0
 R("Mhmfurnb_Transp"): 0 M10fthfc = 0
 R("R_DM_HMFURN"): - Mhmfurnc = 0
 R("MhomLb_Transp"): 0 M10fthfc = 0
 R("R_EX_hom_L_e_"): - MhomLe = 0
 R("R_EX_hom_L_e__Rev"): + MhomLe = 0
 R("R_HOMtex"): + MhomLp - MhomLe = 0
 R("R_HOMtex_Rev"): - MhomLp + MhomLe = 0
 R("R_O16AP1pp"): + Mudcpdpp - 2 Mo16aundp + Mo16a2undp + Mhp = 0
 R("R_O16AP2pp"): + Mudcpdpp - Mo16aundp + Mo16a3undp - Mo16a2undp + Mhp
 = 0
 R("R_O16AP3pp"): + Mudcpdpp - Mo16aundp + Mo16a4undp - Mo16a3undp + Mhp
 = 0
 R("Mhxab_Transp"): 0 M10fthfc = 0
 R("R_EX_hxa_e_"): - Mhxae = 0
 R("R_EX_hxa_e__Rev"): + Mhxae = 0
 R("R_HXAtex"): + Mhxap - Mhxae = 0
 R("R_HXAtex_Rev"): - Mhxap + Mhxae = 0
 R("Mhxanb_Transp"): 0 M10fthfc = 0
 R("R_EX_hxan_e_"): - Mhxane = 0
 R("R_EX_hxan_e__Rev"): + Mhxane = 0
 R("R_HXPRT"): - Mprppc + Mppic + Mimpc - Mhxanc = 0
 R("R_HYXNtpp"): - Mhxanp + Mhxanc = 0
 R("R_HYXNtpp_Rev"): + Mhxanp - Mhxanc = 0
 R("R_PUNP5"): + Mr1pc - Mpic - Minsc + Mhxanc = 0
 R("R_PUNP5_Rev"): - Mr1pc + Mpic + Minsc - Mhxanc = 0
 R("R_HYXNtex"): + Mhxanp - Mhxane = 0
 R("R_HYXNtex_Rev"): - Mhxanp + Mhxane = 0
 R("MidonLb_Transp"): 0 M10fthfc = 0
 R("R_EX_idon_L_e_"): - MidonLe = 0
 R("R_EX_idon_L_e__Rev"): + MidonLe = 0
 R("R_IDONtex"): + MidonLp - MidonLe = 0
 R("R_IDONtex_Rev"): - MidonLp + MidonLe = 0
 R("MileLb_Transp"): 0 M10fthfc = 0
 R("R_EX_ile_L_e_"): - MileLe = 0
 R("R_EX_ile_L_e__Rev"): + MileLe = 0
 R("R_ILEtex"): + MileLp - MileLe = 0
 R("R_ILEtex_Rev"): - MileLp + MileLe = 0
 R("Mimpb_Transp"): 0 M10fthfc = 0
 R("R_EX_imp_e_"): - Mimpe = 0
 R("R_EX_imp_e__Rev"): + Mimpe = 0
 R("R_IMPtex"): + Mimpp - Mimpe = 0
 R("R_IMPtex_Rev"): - Mimpp + Mimpe = 0
 R("Mindoleb_Transp"): 0 M10fthfc = 0
 R("R_EX_indole_e_"): - Mindolee = 0
 R("R_EX_indole_e__Rev"): + Mindolee = 0
 R("R_INDOLEtex"): + Mindolep - Mindolee = 0
 R("R_INDOLEtex_Rev"): - Mindolep + Mindolee = 0
 R("Minostb_Transp"): 0 M10fthfc = 0
 R("R_EX_inost_e_"): - Minoste = 0
 R("R_EX_inost_e__Rev"): + Minoste = 0
 R("R_INOSTt4pp"): - Mna1p + Mna1c - Minostp + Minostc = 0
 R("R_INSTtex"): + Minostp - Minoste = 0
 R("R_INSTtex_Rev"): - Minostp + Minoste = 0
 R("Minsb_Transp"): 0 M10fthfc = 0
 R("R_EX_ins_e_"): - Minse = 0
 R("R_EX_ins_e__Rev"): + Minse = 0
 R("R_INStex"): + Minsp - Minse = 0
 R("R_INStex_Rev"): - Minsp + Minse = 0
 R("Misetacb_Transp"): 0 M10fthfc = 0
 R("R_EX_isetac_e_"): - Misetace = 0
 R("R_EX_isetac_e__Rev"): + Misetace = 0
 R("R_ISETACtex"): + Misetacp - Misetace = 0
 R("R_ISETACtex_Rev"): - Misetacp + Misetace = 0
 R("Mkb_Transp"): 0 M10fthfc = 0
 R("R_EX_k_e_"): - Mke = 0
 R("R_EX_k_e__Rev"): + Mke = 0
 R("Mkdo2lipid4b_Transp"): 0 M10fthfc = 0
 R("R_EX_kdo2lipid4_e_"): - Mkdo2lipid4e = 0
 R("R_EX_kdo2lipid4_e__Rev"): + Mkdo2lipid4e = 0
 R("R_K2L4Atex"): - Mkdo2lipid4p + Mkdo2lipid4e = 0
 R("R_Ktex"): + Mkp - Mke = 0
 R("R_Ktex_Rev"): - Mkp + Mke = 0
 R("MlacDb_Transp"): 0 M10fthfc = 0
 R("R_EX_lac_D_e_"): - MlacDe = 0
 R("R_EX_lac_D_e__Rev"): + MlacDe = 0
 R("R_LDH_D2"): + Mq8h2c - Mq8c + Mpyrc - MlacDc = 0
 R("R_D_LACtex"): + MlacDp - MlacDe = 0
 R("R_D_LACtex_Rev"): - MlacDp + MlacDe = 0
 R("MlacLb_Transp"): 0 M10fthfc = 0
 R("R_EX_lac_L_e_"): - MlacLe = 0
 R("R_EX_lac_L_e__Rev"): + MlacLe = 0
 R("R_L_LACD2"): + Mq8h2c - Mq8c + Mpyrc - MlacLc = 0
 R("R_L_LACD3"): + Mpyrc - Mmqn8c + Mmql8c - MlacLc = 0
 R("R_L_LACtex"): + MlacLp - MlacLe = 0
 R("R_L_LACtex_Rev"): - MlacLp + MlacLe = 0
 R("R_ALAGLUE"): + MLalaLgluc - MLalaDgluc = 0
 R("R_ALAGLUE_Rev"): - MLalaLgluc + MLalaDgluc = 0
 R("MLalaDgluMdapb_Transp"): 0 M10fthfc = 0
 R("R_EX_LalaDgluMdap_e_"): - MLalaDgluMdape = 0
 R("R_EX_LalaDgluMdap_e__Rev"): + MLalaDgluMdape = 0
 R("MLalaDgluMdapDalab_Transp"): 0 M10fthfc = 0
 R("R_EX_LalaDgluMdapDala_e_"): - MLalaDgluMdapDalae = 0
 R("R_EX_LalaDgluMdapDala_e__Rev"): + MLalaDgluMdapDalae = 0
 R("R_4PEPTtex"): + MLalaDgluMdapDalap - MLalaDgluMdapDalae = 0
 R("R_4PEPTtex_Rev"): - MLalaDgluMdapDalap + MLalaDgluMdapDalae = 0
 R("R_3PEPTtex"): + MLalaDgluMdapp - MLalaDgluMdape = 0
 R("R_3PEPTtex_Rev"): - MLalaDgluMdapp + MLalaDgluMdape = 0
 R("Mlctsb_Transp"): 0 M10fthfc = 0
 R("R_EX_lcts_e_"): - Mlctse = 0
 R("R_EX_lcts_e__Rev"): + Mlctse = 0
 R("R_LCTStex"): + Mlctsp - Mlctse = 0
 R("R_LCTStex_Rev"): - Mlctsp + Mlctse = 0
 R("MleuLb_Transp"): 0 M10fthfc = 0
 R("R_EX_leu_L_e_"): - MleuLe = 0
 R("R_EX_leu_L_e__Rev"): + MleuLe = 0
 R("R_LEUtex"): + MleuLp - MleuLe = 0
 R("R_LEUtex_Rev"): - MleuLp + MleuLe = 0
 R("Mlipab_Transp"): 0 M10fthfc = 0
 R("R_EX_lipa_e_"): - Mlipae = 0
 R("R_EX_lipa_e__Rev"): + Mlipae = 0
 R("Mlipacoldb_Transp"): 0 M10fthfc = 0
 R("R_EX_lipa_cold_e_"): - Mlipacolde = 0
 R("R_EX_lipa_cold_e__Rev"): + Mlipacolde = 0
 R("R_CLIPAtex"): - Mlipacoldp + Mlipacolde = 0
 R("R_LIPAtex"): - Mlipap + Mlipae = 0
 R("MlysLb_Transp"): 0 M10fthfc = 0
 R("R_EX_lys_L_e_"): - MlysLe = 0
 R("R_EX_lys_L_e__Rev"): + MlysLe = 0
 R("R_LYStex"): + MlysLp - MlysLe = 0
 R("R_LYStex_Rev"): - MlysLp + MlysLe = 0
 R("MlyxLb_Transp"): 0 M10fthfc = 0
 R("R_EX_lyx_L_e_"): - MlyxLe = 0
 R("R_EX_lyx_L_e__Rev"): + MlyxLe = 0
 R("R_LYXI"): + MxyluLc - MlyxLc = 0
 R("R_LYXtex"): + MlyxLp - MlyxLe = 0
 R("R_LYXtex_Rev"): - MlyxLp + MlyxLe = 0
 R("MmalDb_Transp"): 0 M10fthfc = 0
 R("R_EX_mal_D_e_"): - MmalDe = 0
 R("R_EX_mal_D_e__Rev"): + MmalDe = 0
 R("R_MALDtex"): + MmalDp - MmalDe = 0
 R("R_MALDtex_Rev"): - MmalDp + MmalDe = 0
 R("MmalLb_Transp"): 0 M10fthfc = 0
 R("R_EX_mal_L_e_"): - MmalLe = 0
 R("R_EX_mal_L_e__Rev"): + MmalLe = 0
 R("R_MDH2"): + Mq8h2c - Mq8c + Moaac - MmalLc = 0
 R("R_MDH3"): + Moaac - Mmqn8c + Mmql8c - MmalLc = 0
 R("R_MALtex"): + MmalLp - MmalLe = 0
 R("R_MALtex_Rev"): - MmalLp + MmalLe = 0
 R("R_MALTptspp"): + Mpyrc - Mpepc - Mmaltp + Mmalt6pc = 0
 R("Mmaltb_Transp"): 0 M10fthfc = 0
 R("R_EX_malt_e_"): - Mmalte = 0
 R("R_EX_malt_e__Rev"): + Mmalte = 0
 R("R_MALTtexi"): + Mmaltp - Mmalte = 0
 R("Mmalthxb_Transp"): 0 M10fthfc = 0
 R("R_EX_malthx_e_"): - Mmalthxe = 0
 R("R_EX_malthx_e__Rev"): + Mmalthxe = 0
 R("R_MALTHXtexi"): + Mmalthxp - Mmalthxe = 0
 R("Mmaltptb_Transp"): 0 M10fthfc = 0
 R("R_EX_maltpt_e_"): - Mmaltpte = 0
 R("R_EX_maltpt_e__Rev"): + Mmaltpte = 0
 R("R_MALTPTtexi"): + Mmaltptp - Mmaltpte = 0
 R("Mmalttrb_Transp"): 0 M10fthfc = 0
 R("R_EX_malttr_e_"): - Mmalttre = 0
 R("R_EX_malttr_e__Rev"): + Mmalttre = 0
 R("R_MALTTRtexi"): + Mmalttrp - Mmalttre = 0
 R("Mmaltttrb_Transp"): 0 M10fthfc = 0
 R("R_EX_maltttr_e_"): - Mmaltttre = 0
 R("R_EX_maltttr_e__Rev"): + Mmaltttre = 0
 R("R_MALTTTRtexi"): + Mmaltttrp - Mmaltttre = 0
 R("R_PMANM"): + Mman6pc - Mman1pc = 0
 R("R_PMANM_Rev"): - Mman6pc + Mman1pc = 0
 R("Mman6pb_Transp"): 0 M10fthfc = 0
 R("R_EX_man6p_e_"): - Mman6pe = 0
 R("R_EX_man6p_e__Rev"): + Mman6pe = 0
 R("R_MAN6Pt6_2pp"): + 2 Mpip - 2 Mpic - Mman6pp + Mman6pc = 0
 R("R_MANptspp"): + Mpyrc - Mpepc - Mmanp + Mman6pc = 0
 R("R_MAN6Ptex"): + Mman6pp - Mman6pe = 0
 R("R_MAN6Ptex_Rev"): - Mman6pp + Mman6pe = 0
 R("R_MANGLYCptspp"): + Mpyrc - Mpepc - Mmanglycp + Mman6pglycc = 0
 R("Mmanb_Transp"): 0 M10fthfc = 0
 R("R_EX_man_e_"): - Mmane = 0
 R("R_EX_man_e__Rev"): + Mmane = 0
 R("R_MANtex"): + Mmanp - Mmane = 0
 R("R_MANtex_Rev"): - Mmanp + Mmane = 0
 R("Mmanglycb_Transp"): 0 M10fthfc = 0
 R("R_EX_manglyc_e_"): - Mmanglyce = 0
 R("R_EX_manglyc_e__Rev"): + Mmanglyce = 0
 R("R_MANGLYCtex"): + Mmanglycp - Mmanglyce = 0
 R("R_MANGLYCtex_Rev"): - Mmanglycp + Mmanglyce = 0
 R("Mmelibb_Transp"): 0 M10fthfc = 0
 R("R_EX_melib_e_"): - Mmelibe = 0
 R("R_EX_melib_e__Rev"): + Mmelibe = 0
 R("R_MELIBtex"): + Mmelibp - Mmelibe = 0
 R("R_MELIBtex_Rev"): - Mmelibp + Mmelibe = 0
 R("MmetDb_Transp"): 0 M10fthfc = 0
 R("R_EX_met_D_e_"): - MmetDe = 0
 R("R_EX_met_D_e__Rev"): + MmetDe = 0
 R("R_METDtex"): + MmetDp - MmetDe = 0
 R("R_METDtex_Rev"): - MmetDp + MmetDe = 0
 R("R_MTHFD"): + Mnadphc - Mnadpc - Mmlthfc + Mmethfc = 0
 R("R_MTHFD_Rev"): - Mnadphc + Mnadpc + Mmlthfc - Mmethfc = 0
 R("MmetLb_Transp"): 0 M10fthfc = 0
 R("R_EX_met_L_e_"): - MmetLe = 0
 R("R_EX_met_L_e__Rev"): + MmetLe = 0
 R("R_METtex"): + MmetLp - MmetLe = 0
 R("R_METtex_Rev"): - MmetLp + MmetLe = 0
 R("MmetsoxRLb_Transp"): 0 M10fthfc = 0
 R("R_EX_metsox_R_L_e_"): - MmetsoxRLe = 0
 R("R_EX_metsox_R_L_e__Rev"): + MmetsoxRLe = 0
 R("R_METSOX2tex"): + MmetsoxRLp - MmetsoxRLe = 0
 R("R_METSOX2tex_Rev"): - MmetsoxRLp + MmetsoxRLe = 0
 R("MmetsoxSLb_Transp"): 0 M10fthfc = 0
 R("R_EX_metsox_S_L_e_"): - MmetsoxSLe = 0
 R("R_EX_metsox_S_L_e__Rev"): + MmetsoxSLe = 0
 R("R_METSOX1tex"): + MmetsoxSLp - MmetsoxSLe = 0
 R("R_METSOX1tex_Rev"): - MmetsoxSLp + MmetsoxSLe = 0
 R("Mmg2b_Transp"): 0 M10fthfc = 0
 R("R_EX_mg2_e_"): - Mmg2e = 0
 R("R_EX_mg2_e__Rev"): + Mmg2e = 0
 R("R_MG2tpp"): - Mmg2p + Mmg2c = 0
 R("R_MG2tex"): + Mmg2p - Mmg2e = 0
 R("R_MG2tex_Rev"): - Mmg2p + Mmg2e = 0
 R("R_MCITL2"): + Msuccc + Mpyrc - Mmicitc = 0
 R("R_MCITL2_Rev"): - Msuccc - Mpyrc + Mmicitc = 0
 R("Mminohpb_Transp"): 0 M10fthfc = 0
 R("R_EX_minohp_e_"): - Mminohpe = 0
 R("R_EX_minohp_e__Rev"): + Mminohpe = 0
 R("R_MINOHPtexi"): + Mminohpp - Mminohpe = 0
 R("R_MME"): + MmmcoaSc - MmmcoaRc = 0
 R("R_MME_Rev"): - MmmcoaSc + MmmcoaRc = 0
 R("R_MMM2"): - Msuccoac + MmmcoaRc = 0
 R("Mmmetb_Transp"): 0 M10fthfc = 0
 R("R_EX_mmet_e_"): - Mmmete = 0
 R("R_EX_mmet_e__Rev"): + Mmmete = 0
 R("R_MMETtex"): + Mmmetp - Mmmete = 0
 R("R_MMETtex_Rev"): - Mmmetp + Mmmete = 0
 R("Mmn2b_Transp"): 0 M10fthfc = 0
 R("R_EX_mn2_e_"): - Mmn2e = 0
 R("R_EX_mn2_e__Rev"): + Mmn2e = 0
 R("R_MN2tpp"): - Mmn2p + Mmn2c = 0
 R("R_MNtex"): + Mmn2p - Mmn2e = 0
 R("R_MNtex_Rev"): - Mmn2p + Mmn2e = 0
 R("R_MNLptspp"): + Mpyrc - Mpepc - Mmnlp + Mmnl1pc = 0
 R("Mmnlb_Transp"): 0 M10fthfc = 0
 R("R_EX_mnl_e_"): - Mmnle = 0
 R("R_EX_mnl_e__Rev"): + Mmnle = 0
 R("R_MNLtex"): + Mmnlp - Mmnle = 0
 R("R_MNLtex_Rev"): - Mmnlp + Mmnle = 0
 R("Mmobdb_Transp"): 0 M10fthfc = 0
 R("R_EX_mobd_e_"): - Mmobde = 0
 R("R_EX_mobd_e__Rev"): + Mmobde = 0
 R("R_MOBDtex"): + Mmobdp - Mmobde = 0
 R("R_MOBDtex_Rev"): - Mmobdp + Mmobde = 0
 R("Mmso3b_Transp"): 0 M10fthfc = 0
 R("R_EX_mso3_e_"): - Mmso3e = 0
 R("R_EX_mso3_e__Rev"): + Mmso3e = 0
 R("R_MSO3tex"): + Mmso3p - Mmso3e = 0
 R("R_MSO3tex_Rev"): - Mmso3p + Mmso3e = 0
 R("Mn2ob_Transp"): 0 M10fthfc = 0
 R("R_EX_n2o_e_"): - Mn2oe = 0
 R("R_EX_n2o_e__Rev"): + Mn2oe = 0
 R("R_N2Otpp"): - Mn2op + Mn2oc = 0
 R("R_N2Otpp_Rev"): + Mn2op - Mn2oc = 0
 R("R_N2Otex"): + Mn2op - Mn2oe = 0
 R("R_N2Otex_Rev"): - Mn2op + Mn2oe = 0
 R("Mna1b_Transp"): 0 M10fthfc = 0
 R("R_EX_na1_e_"): - Mna1e = 0
 R("R_EX_na1_e__Rev"): + Mna1e = 0
 R("R_PNTOt4pp"): - MpntoRp + MpntoRc - Mna1p + Mna1c = 0
 R("R_PPAt4pp"): - Mppap + Mppac - Mna1p + Mna1c = 0
 R("R_PROt4pp"): - MproLp + MproLc - Mna1p + Mna1c = 0
 R("R_SERt4pp"): - MserLp + MserLc - Mna1p + Mna1c = 0
 R("R_THRt4pp"): - MthrLp + MthrLc - Mna1p + Mna1c = 0
 R("R_NAtex"): + Mna1p - Mna1e = 0
 R("R_NAtex_Rev"): - Mna1p + Mna1e = 0
 R("Mnacb_Transp"): 0 M10fthfc = 0
 R("R_EX_nac_e_"): - Mnace = 0
 R("R_EX_nac_e__Rev"): + Mnace = 0
 R("R_NACtpp"): - Mnacp + Mnacc = 0
 R("R_NACtex"): + Mnacp - Mnace = 0
 R("R_NACtex_Rev"): - Mnacp + Mnace = 0
 R("R_NADTRHD"): - Mnadphc + Mnadpc + Mnadhc - Mnadc = 0
 R("Mnh4b_Transp"): 0 M10fthfc = 0
 R("R_EX_nh4_e_"): - Mnh4e = 0
 R("R_EX_nh4_e__Rev"): + Mnh4e = 0
 R("R_NH4tpp"): - Mnh4p + Mnh4c = 0
 R("R_NH4tpp_Rev"): + Mnh4p - Mnh4c = 0
 R("R_SERD_D"): - MserDc + Mpyrc + Mnh4c = 0
 R("R_SERD_L"): - MserLc + Mpyrc + Mnh4c = 0
 R("R_NH4tex"): + Mnh4p - Mnh4e = 0
 R("R_NH4tex_Rev"): - Mnh4p + Mnh4e = 0
 R("Mni2b_Transp"): 0 M10fthfc = 0
 R("R_EX_ni2_e_"): - Mni2e = 0
 R("R_EX_ni2_e__Rev"): + Mni2e = 0
 R("R_NI2tpp"): - Mni2p + Mni2c = 0
 R("R_NI2tex"): + Mni2p - Mni2e = 0
 R("R_NI2tex_Rev"): - Mni2p + Mni2e = 0
 R("Mnmnb_Transp"): 0 M10fthfc = 0
 R("R_EX_nmn_e_"): - Mnmne = 0
 R("R_EX_nmn_e__Rev"): + Mnmne = 0
 R("R_NMNPtpp"): - Mnmnp + Mnmnc = 0
 R("R_NMNtex"): + Mnmnp - Mnmne = 0
 R("R_NMNtex_Rev"): - Mnmnp + Mnmne = 0
 R("Mno2b_Transp"): 0 M10fthfc = 0
 R("R_EX_no2_e_"): - Mno2e = 0
 R("R_EX_no2_e__Rev"): + Mno2e = 0
 R("R_NO3t7pp"): - Mno3p + Mno3c + Mno2p - Mno2c = 0
 R("R_NO2tex"): + Mno2p - Mno2e = 0
 R("R_NO2tex_Rev"): - Mno2p + Mno2e = 0
 R("Mno3b_Transp"): 0 M10fthfc = 0
 R("R_EX_no3_e_"): - Mno3e = 0
 R("R_EX_no3_e__Rev"): + Mno3e = 0
 R("R_NO3tex"): + Mno3p - Mno3e = 0
 R("R_NO3tex_Rev"): - Mno3p + Mno3e = 0
 R("Mnob_Transp"): 0 M10fthfc = 0
 R("R_EX_no_e_"): - Mnoe = 0
 R("R_EX_no_e__Rev"): + Mnoe = 0
 R("R_NOtpp"): - Mnop + Mnoc = 0
 R("R_NOtpp_Rev"): + Mnop - Mnoc = 0
 R("R_NOtex"): + Mnop - Mnoe = 0
 R("R_NOtex_Rev"): - Mnop + Mnoe = 0
 R("Mo16a4colipab_Transp"): 0 M10fthfc = 0
 R("R_EX_o16a4colipa_e_"): - Mo16a4colipae = 0
 R("R_EX_o16a4colipa_e__Rev"): + Mo16a4colipae = 0
 R("R_O16A4COLIPAtex"): - Mo16a4colipap + Mo16a4colipae = 0
 R("R_O16AUNDtpp"): + Mo16aundp - Mo16aundc = 0
 R("Mo2b_Transp"): 0 M10fthfc = 0
 R("R_EX_o2_e_"): - Mo2e = 0
 R("R_EX_o2_e__Rev"): + Mo2e = 0
 R("R_O2tpp"): - Mo2p + Mo2c = 0
 R("R_O2tpp_Rev"): + Mo2p - Mo2c = 0
 R("R_O2tex"): + Mo2p - Mo2e = 0
 R("R_O2tex_Rev"): - Mo2p + Mo2e = 0
 R("Mo2sb_Transp"): 0 M10fthfc = 0
 R("R_EX_o2s_e_"): - Mo2se = 0
 R("R_EX_o2s_e__Rev"): + Mo2se = 0
 R("R_O2Stex"): + Mo2sp - Mo2se = 0
 R("R_O2Stex_Rev"): - Mo2sp + Mo2se = 0
 R("Mocdcab_Transp"): 0 M10fthfc = 0
 R("R_EX_ocdca_e_"): - Mocdcae = 0
 R("R_EX_ocdca_e__Rev"): + Mocdcae = 0
 R("R_OCDCAtexi"): + Mocdcap - Mocdcae = 0
 R("Mocdceab_Transp"): 0 M10fthfc = 0
 R("R_EX_ocdcea_e_"): - Mocdceae = 0
 R("R_EX_ocdcea_e__Rev"): + Mocdceae = 0
 R("R_OCDCEAtexi"): + Mocdceap - Mocdceae = 0
 R("Moctab_Transp"): 0 M10fthfc = 0
 R("R_EX_octa_e_"): - Moctae = 0
 R("R_EX_octa_e__Rev"): + Moctae = 0
 R("R_OCTAtex"): + Moctap - Moctae = 0
 R("R_OCTAtex_Rev"): - Moctap + Moctae = 0
 R("R_ODECOAI"): - Modecoac + Mod2coac = 0
 R("Mornb_Transp"): 0 M10fthfc = 0
 R("R_EX_orn_e_"): - Morne = 0
 R("R_EX_orn_e__Rev"): + Morne = 0
 R("R_PTRCORNt7pp"): - Mptrcp + Mptrcc + Mornp - Mornc = 0
 R("R_PTRCORNt7pp_Rev"): + Mptrcp - Mptrcc - Mornp + Mornc = 0
 R("R_ORNtex"): + Mornp - Morne = 0
 R("R_ORNtex_Rev"): - Mornp + Morne = 0
 R("R_ORPT"): + Mprppc - Mppic + Morotc - Morot5pc = 0
 R("R_ORPT_Rev"): - Mprppc + Mppic - Morotc + Morot5pc = 0
 R("Morotb_Transp"): 0 M10fthfc = 0
 R("R_EX_orot_e_"): - Morote = 0
 R("R_EX_orot_e__Rev"): + Morote = 0
 R("R_OROTtex"): + Morotp - Morote = 0
 R("R_OROTtex_Rev"): - Morotp + Morote = 0
 R("Moxamb_Transp"): 0 M10fthfc = 0
 R("R_DM_OXAM"): - Moxamc = 0
 R("Mpacaldb_Transp"): 0 M10fthfc = 0
 R("R_EX_pacald_e_"): - Mpacalde = 0
 R("R_EX_pacald_e__Rev"): + Mpacalde = 0
 R("R_PACALDtex"): + Mpacaldp - Mpacalde = 0
 R("R_PACALDtex_Rev"): - Mpacaldp + Mpacalde = 0
 R("Mpeamnb_Transp"): 0 M10fthfc = 0
 R("R_EX_peamn_e_"): - Mpeamne = 0
 R("R_EX_peamn_e__Rev"): + Mpeamne = 0
 R("R_PEAMNtex"): + Mpeamnp - Mpeamne = 0
 R("R_PEAMNtex_Rev"): - Mpeamnp + Mpeamne = 0
 R("R_SBTptspp"): - MsbtDp + Msbt6pc + Mpyrc - Mpepc = 0
 R("R_SUCptspp"): - Msucrp + Msuc6pc + Mpyrc - Mpepc = 0
 R("R_TREptspp"): - Mtrep + Mtre6pc + Mpyrc - Mpepc = 0
 R("R_UAGCVT"): - Muacgamc + Muaccgc + Mpic - Mpepc = 0
 R("MpheLb_Transp"): 0 M10fthfc = 0
 R("R_EX_phe_L_e_"): - MpheLe = 0
 R("R_EX_phe_L_e__Rev"): + MpheLe = 0
 R("R_PHEtex"): + MpheLp - MpheLe = 0
 R("R_PHEtex_Rev"): - MpheLp + MpheLe = 0
 R("Mphemeb_Transp"): 0 M10fthfc = 0
 R("R_EX_pheme_e_"): - Mphemee = 0
 R("R_EX_pheme_e__Rev"): + Mphemee = 0
 R("R_PHEMEtiex"): - Mphemep + Mphemee = 0
 R("Mpib_Transp"): 0 M10fthfc = 0
 R("R_EX_pi_e_"): - Mpie = 0
 R("R_EX_pi_e__Rev"): + Mpie = 0
 R("R_PUNP7"): - Mxtsnc + Mxanc + Mr1pc - Mpic = 0
 R("R_PUNP7_Rev"): + Mxtsnc - Mxanc - Mr1pc + Mpic = 0
 R("R_PYNP2r"): - Muric + Murac + Mr1pc - Mpic = 0
 R("R_PYNP2r_Rev"): + Muric - Murac - Mr1pc + Mpic = 0
 R("R_PItex"): + Mpip - Mpie = 0
 R("R_PItex_Rev"): - Mpip + Mpie = 0
 R("MpntoRb_Transp"): 0 M10fthfc = 0
 R("R_EX_pnto_R_e_"): - MpntoRe = 0
 R("R_EX_pnto_R_e__Rev"): + MpntoRe = 0
 R("R_PNTOtex"): + MpntoRp - MpntoRe = 0
 R("R_PNTOtex_Rev"): - MpntoRp + MpntoRe = 0
 R("Mppab_Transp"): 0 M10fthfc = 0
 R("R_EX_ppa_e_"): - Mppae = 0
 R("R_EX_ppa_e__Rev"): + Mppae = 0
 R("R_PPCSCT"): + Msuccoac - Msuccc - Mppcoac + Mppac = 0
 R("R_PPAtex"): + Mppap - Mppae = 0
 R("R_PPAtex_Rev"): - Mppap + Mppae = 0
 R("Mppalb_Transp"): 0 M10fthfc = 0
 R("R_EX_ppal_e_"): - Mppale = 0
 R("R_EX_ppal_e__Rev"): + Mppale = 0
 R("R_PPALtpp"): + Mppalp - Mppalc = 0
 R("R_PPALtpp_Rev"): - Mppalp + Mppalc = 0
 R("R_PPALtex"): + Mppalp - Mppale = 0
 R("R_PPALtex_Rev"): - Mppalp + Mppale = 0
 R("R_UPPRT"): - Murac + Mumpc - Mprppc + Mppic = 0
 R("R_XPPT"): + Mxmpc - Mxanc - Mprppc + Mppic = 0
 R("Mpppnb_Transp"): 0 M10fthfc = 0
 R("R_EX_pppn_e_"): - Mpppne = 0
 R("R_EX_pppn_e__Rev"): + Mpppne = 0
 R("R_PPPNtex"): + Mpppnp - Mpppne = 0
 R("R_PPPNtex_Rev"): - Mpppnp + Mpppne = 0
 R("Mpptb_Transp"): 0 M10fthfc = 0
 R("R_EX_ppt_e_"): - Mppte = 0
 R("R_EX_ppt_e__Rev"): + Mppte = 0
 R("R_PPTtex"): + Mpptp - Mppte = 0
 R("R_PPTtex_Rev"): - Mpptp + Mppte = 0
 R("R_PRMICI"): + Mprlpc - Mprfpc = 0
 R("R_PRMICI_Rev"): - Mprlpc + Mprfpc = 0
 R("Mproglyb_Transp"): 0 M10fthfc = 0
 R("R_EX_progly_e_"): - Mproglye = 0
 R("R_EX_progly_e__Rev"): + Mproglye = 0
 R("R_PROGLYtex"): + Mproglyp - Mproglye = 0
 R("R_PROGLYtex_Rev"): - Mproglyp + Mproglye = 0
 R("MproLb_Transp"): 0 M10fthfc = 0
 R("R_EX_pro_L_e_"): - MproLe = 0
 R("R_EX_pro_L_e__Rev"): + MproLe = 0
 R("R_PROtex"): + MproLp - MproLe = 0
 R("R_PROtex_Rev"): - MproLp + MproLe = 0
 R("Mpsclysb_Transp"): 0 M10fthfc = 0
 R("R_EX_psclys_e_"): - Mpsclyse = 0
 R("R_EX_psclys_e__Rev"): + Mpsclyse = 0
 R("R_PSCLYStex"): + Mpsclysp - Mpsclyse = 0
 R("R_PSCLYStex_Rev"): - Mpsclysp + Mpsclyse = 0
 R("MpserLb_Transp"): 0 M10fthfc = 0
 R("R_EX_pser_L_e_"): - MpserLe = 0
 R("R_EX_pser_L_e__Rev"): + MpserLe = 0
 R("R_PSERtex"): + MpserLp - MpserLe = 0
 R("R_PSERtex_Rev"): - MpserLp + MpserLe = 0
 R("Mptrcb_Transp"): 0 M10fthfc = 0
 R("R_EX_ptrc_e_"): - Mptrce = 0
 R("R_EX_ptrc_e__Rev"): + Mptrce = 0
 R("R_PTRCtex"): + Mptrcp - Mptrce = 0
 R("R_PTRCtex_Rev"): - Mptrcp + Mptrce = 0
 R("Mpyrb_Transp"): 0 M10fthfc = 0
 R("R_EX_pyr_e_"): - Mpyre = 0
 R("R_EX_pyr_e__Rev"): + Mpyre = 0
 R("R_PYRtex"): + Mpyrp - Mpyre = 0
 R("R_PYRtex_Rev"): - Mpyrp + Mpyre = 0
 R("R_PPM"): + Mr5pc - Mr1pc = 0
 R("R_PPM_Rev"): - Mr5pc + Mr1pc = 0
 R("Mr5pb_Transp"): 0 M10fthfc = 0
 R("R_EX_r5p_e_"): - Mr5pe = 0
 R("R_EX_r5p_e__Rev"): + Mr5pe = 0
 R("R_RPI"): + Mru5pDc - Mr5pc = 0
 R("R_RPI_Rev"): - Mru5pDc + Mr5pc = 0
 R("R_R5Ptex"): + Mr5pp - Mr5pe = 0
 R("R_R5Ptex_Rev"): - Mr5pp + Mr5pe = 0
 R("MribDb_Transp"): 0 M10fthfc = 0
 R("R_EX_rib_D_e_"): - MribDe = 0
 R("R_EX_rib_D_e__Rev"): + MribDe = 0
 R("R_RIBtex"): + MribDp - MribDe = 0
 R("R_RIBtex_Rev"): - MribDp + MribDe = 0
 R("R_RMI"): - Mrmnc + Mrmlc = 0
 R("R_RMI_Rev"): + Mrmnc - Mrmlc = 0
 R("Mrmnb_Transp"): 0 M10fthfc = 0
 R("R_EX_rmn_e_"): - Mrmne = 0
 R("R_EX_rmn_e__Rev"): + Mrmne = 0
 R("R_RMNtex"): + Mrmnp - Mrmne = 0
 R("R_RMNtex_Rev"): - Mrmnp + Mrmne = 0
 R("R_RPE"): + Mxu5pDc - Mru5pDc = 0
 R("R_RPE_Rev"): - Mxu5pDc + Mru5pDc = 0
 R("R_RBP4E"): + Mxu5pDc - Mru5pLc = 0
 R("R_RBP4E_Rev"): - Mxu5pDc + Mru5pLc = 0
 R("R_X5PL3E"): - Mxu5pLc + Mru5pLc = 0
 R("MsbtDb_Transp"): 0 M10fthfc = 0
 R("R_EX_sbt_D_e_"): - MsbtDe = 0
 R("R_EX_sbt_D_e__Rev"): + MsbtDe = 0
 R("R_SBTtex"): + MsbtDp - MsbtDe = 0
 R("R_SBTtex_Rev"): - MsbtDp + MsbtDe = 0
 R("MserDb_Transp"): 0 M10fthfc = 0
 R("R_EX_ser_D_e_"): - MserDe = 0
 R("R_EX_ser_D_e__Rev"): + MserDe = 0
 R("R_DSERtex"): + MserDp - MserDe = 0
 R("R_DSERtex_Rev"): - MserDp + MserDe = 0
 R("MserLb_Transp"): 0 M10fthfc = 0
 R("R_EX_ser_L_e_"): - MserLe = 0
 R("R_EX_ser_L_e__Rev"): + MserLe = 0
 R("R_SERtex"): + MserLp - MserLe = 0
 R("R_SERtex_Rev"): - MserLp + MserLe = 0
 R("Mskmb_Transp"): 0 M10fthfc = 0
 R("R_EX_skm_e_"): - Mskme = 0
 R("R_EX_skm_e__Rev"): + Mskme = 0
 R("R_SKMtex"): + Mskmp - Mskme = 0
 R("R_SKMtex_Rev"): - Mskmp + Mskme = 0
 R("Mso2b_Transp"): 0 M10fthfc = 0
 R("R_EX_so2_e_"): - Mso2e = 0
 R("R_EX_so2_e__Rev"): + Mso2e = 0
 R("R_SO2tpp"): - Mso2p + Mso2c = 0
 R("R_SO2tpp_Rev"): + Mso2p - Mso2c = 0
 R("R_SO2tex"): + Mso2p - Mso2e = 0
 R("R_SO2tex_Rev"): - Mso2p + Mso2e = 0
 R("Mso3b_Transp"): 0 M10fthfc = 0
 R("R_EX_so3_e_"): - Mso3e = 0
 R("R_EX_so3_e__Rev"): + Mso3e = 0
 R("R_SO3tex"): + Mso3p - Mso3e = 0
 R("R_SO3tex_Rev"): - Mso3p + Mso3e = 0
 R("Mso4b_Transp"): 0 M10fthfc = 0
 R("R_EX_so4_e_"): - Mso4e = 0
 R("R_EX_so4_e__Rev"): + Mso4e = 0
 R("R_SO4tex"): + Mso4p - Mso4e = 0
 R("R_SO4tex_Rev"): - Mso4p + Mso4e = 0
 R("Mspmdb_Transp"): 0 M10fthfc = 0
 R("R_EX_spmd_e_"): - Mspmde = 0
 R("R_EX_spmd_e__Rev"): + Mspmde = 0
 R("R_SPMDtex"): + Mspmdp - Mspmde = 0
 R("R_SPMDtex_Rev"): - Mspmdp + Mspmde = 0
 R("Msuccb_Transp"): 0 M10fthfc = 0
 R("R_EX_succ_e_"): - Msucce = 0
 R("R_EX_succ_e__Rev"): + Msucce = 0
 R("R_TARTRt7pp"): - MtartrLp + MtartrLc + Msuccp - Msuccc = 0
 R("R_TARTRt7pp_Rev"): + MtartrLp - MtartrLc - Msuccp + Msuccc = 0
 R("R_SUCCtex"): + Msuccp - Msucce = 0
 R("R_SUCCtex_Rev"): - Msuccp + Msucce = 0
 R("Msucrb_Transp"): 0 M10fthfc = 0
 R("R_EX_sucr_e_"): - Msucre = 0
 R("R_EX_sucr_e__Rev"): + Msucre = 0
 R("R_SUCRtex"): + Msucrp - Msucre = 0
 R("R_SUCRtex_Rev"): - Msucrp + Msucre = 0
 R("Msulfacb_Transp"): 0 M10fthfc = 0
 R("R_EX_sulfac_e_"): - Msulface = 0
 R("R_EX_sulfac_e__Rev"): + Msulface = 0
 R("R_SULFACtex"): + Msulfacp - Msulface = 0
 R("R_SULFACtex_Rev"): - Msulfacp + Msulface = 0
 R("MtartrLb_Transp"): 0 M10fthfc = 0
 R("R_EX_tartr_L_e_"): - MtartrLe = 0
 R("R_EX_tartr_L_e__Rev"): + MtartrLe = 0
 R("R_TARTRtex"): + MtartrLp - MtartrLe = 0
 R("R_TARTRtex_Rev"): - MtartrLp + MtartrLe = 0
 R("Mtaurb_Transp"): 0 M10fthfc = 0
 R("R_EX_taur_e_"): - Mtaure = 0
 R("R_EX_taur_e__Rev"): + Mtaure = 0
 R("R_TAURtex"): + Mtaurp - Mtaure = 0
 R("R_TAURtex_Rev"): - Mtaurp + Mtaure = 0
 R("Mtcyntb_Transp"): 0 M10fthfc = 0
 R("R_EX_tcynt_e_"): - Mtcynte = 0
 R("R_EX_tcynt_e__Rev"): + Mtcynte = 0
 R("R_TCYNTtex"): + Mtcyntp - Mtcynte = 0
 R("R_TCYNTtex_Rev"): - Mtcyntp + Mtcynte = 0
 R("R_TDECOAI"): - Mtdecoac + Mtd2coac = 0
 R("Mthmb_Transp"): 0 M10fthfc = 0
 R("R_EX_thm_e_"): - Mthme = 0
 R("R_EX_thm_e__Rev"): + Mthme = 0
 R("R_THMtex"): + Mthmp - Mthme = 0
 R("R_THMtex_Rev"): - Mthmp + Mthme = 0
 R("MthrLb_Transp"): 0 M10fthfc = 0
 R("R_EX_thr_L_e_"): - MthrLe = 0
 R("R_EX_thr_L_e__Rev"): + MthrLe = 0
 R("R_THRtex"): + MthrLp - MthrLe = 0
 R("R_THRtex_Rev"): - MthrLp + MthrLe = 0
 R("Mthrpb_Transp"): 0 M10fthfc = 0
 R("R_EX_thrp_e_"): - Mthrpe = 0
 R("R_EX_thrp_e__Rev"): + Mthrpe = 0
 R("R_THRPtex"): + Mthrpp - Mthrpe = 0
 R("R_THRPtex_Rev"): - Mthrpp + Mthrpe = 0
 R("Mthymb_Transp"): 0 M10fthfc = 0
 R("R_EX_thym_e_"): - Mthyme = 0
 R("R_EX_thym_e__Rev"): + Mthyme = 0
 R("Mthymdb_Transp"): 0 M10fthfc = 0
 R("R_EX_thymd_e_"): - Mthymde = 0
 R("R_EX_thymd_e__Rev"): + Mthymde = 0
 R("R_THMDtex"): + Mthymdp - Mthymde = 0
 R("R_THMDtex_Rev"): - Mthymdp + Mthymde = 0
 R("R_THYMtex"): + Mthymp - Mthyme = 0
 R("R_THYMtex_Rev"): - Mthymp + Mthyme = 0
 R("Mtmab_Transp"): 0 M10fthfc = 0
 R("R_EX_tma_e_"): - Mtmae = 0
 R("R_EX_tma_e__Rev"): + Mtmae = 0
 R("R_TMAtex"): + Mtmap - Mtmae = 0
 R("R_TMAtex_Rev"): - Mtmap + Mtmae = 0
 R("Mtmaob_Transp"): 0 M10fthfc = 0
 R("R_EX_tmao_e_"): - Mtmaoe = 0
 R("R_EX_tmao_e__Rev"): + Mtmaoe = 0
 R("R_TMAOtex"): + Mtmaop - Mtmaoe = 0
 R("R_TMAOtex_Rev"): - Mtmaop + Mtmaoe = 0
 R("Mtreb_Transp"): 0 M10fthfc = 0
 R("R_EX_tre_e_"): - Mtree = 0
 R("R_EX_tre_e__Rev"): + Mtree = 0
 R("R_TREtex"): + Mtrep - Mtree = 0
 R("R_TREtex_Rev"): - Mtrep + Mtree = 0
 R("MtrpLb_Transp"): 0 M10fthfc = 0
 R("R_EX_trp_L_e_"): - MtrpLe = 0
 R("R_EX_trp_L_e__Rev"): + MtrpLe = 0
 R("R_TRPtex"): + MtrpLp - MtrpLe = 0
 R("R_TRPtex_Rev"): - MtrpLp + MtrpLe = 0
 R("Mtsulb_Transp"): 0 M10fthfc = 0
 R("R_EX_tsul_e_"): - Mtsule = 0
 R("R_EX_tsul_e__Rev"): + Mtsule = 0
 R("R_TSULtex"): + Mtsulp - Mtsule = 0
 R("R_TSULtex_Rev"): - Mtsulp + Mtsule = 0
 R("Mttdcab_Transp"): 0 M10fthfc = 0
 R("R_EX_ttdca_e_"): - Mttdcae = 0
 R("R_EX_ttdca_e__Rev"): + Mttdcae = 0
 R("R_TTDCAtexi"): + Mttdcap - Mttdcae = 0
 R("Mttdceab_Transp"): 0 M10fthfc = 0
 R("R_EX_ttdcea_e_"): - Mttdceae = 0
 R("R_EX_ttdcea_e__Rev"): + Mttdceae = 0
 R("R_TTDCEAtexi"): + Mttdceap - Mttdceae = 0
 R("Mtungsb_Transp"): 0 M10fthfc = 0
 R("R_EX_tungs_e_"): - Mtungse = 0
 R("R_EX_tungs_e__Rev"): + Mtungse = 0
 R("R_TUNGStex"): + Mtungsp - Mtungse = 0
 R("R_TUNGStex_Rev"): - Mtungsp + Mtungse = 0
 R("Mtymb_Transp"): 0 M10fthfc = 0
 R("R_EX_tym_e_"): - Mtyme = 0
 R("R_EX_tym_e__Rev"): + Mtyme = 0
 R("R_TYMtex"): + Mtymp - Mtyme = 0
 R("R_TYMtex_Rev"): - Mtymp + Mtyme = 0
 R("MtyrLb_Transp"): 0 M10fthfc = 0
 R("R_EX_tyr_L_e_"): - MtyrLe = 0
 R("R_EX_tyr_L_e__Rev"): + MtyrLe = 0
 R("R_TYRtex"): + MtyrLp - MtyrLe = 0
 R("R_TYRtex_Rev"): - MtyrLp + MtyrLe = 0
 R("Mtyrpb_Transp"): 0 M10fthfc = 0
 R("R_EX_tyrp_e_"): - Mtyrpe = 0
 R("R_EX_tyrp_e__Rev"): + Mtyrpe = 0
 R("R_TYRPtex"): + Mtyrpp - Mtyrpe = 0
 R("R_TYRPtex_Rev"): - Mtyrpp + Mtyrpe = 0
 R("Muacgamb_Transp"): 0 M10fthfc = 0
 R("R_EX_uacgam_e_"): - Muacgame = 0
 R("R_EX_uacgam_e__Rev"): + Muacgame = 0
 R("R_ACGAMT"): + Munagac + Mumpc - Mudcppc - Muacgamc = 0
 R("R_UAG2E"): + Muacmamc - Muacgamc = 0
 R("R_UAG2E_Rev"): - Muacmamc + Muacgamc = 0
 R("R_UACGAMtex"): + Muacgamp - Muacgame = 0
 R("R_UACGAMtex_Rev"): - Muacgamp + Muacgame = 0
 R("R_PAPPT3"): + Mumpc - Mugmdac - Mudcppc + Muagmdac = 0
 R("R_UDCPPtppi"): - Mudcppp + Mudcppc = 0
 R("R_UPLA4FNT"): + MuLa4fnc - MudpLa4fnc + Mudpc - Mudcppc = 0
 R("Mudpacgalb_Transp"): 0 M10fthfc = 0
 R("R_EX_udpacgal_e_"): - Mudpacgale = 0
 R("R_EX_udpacgal_e__Rev"): + Mudpacgale = 0
 R("R_UDPACGALtex"): + Mudpacgalp - Mudpacgale = 0
 R("R_UDPACGALtex_Rev"): - Mudpacgalp + Mudpacgale = 0
 R("Mudpgalb_Transp"): 0 M10fthfc = 0
 R("R_EX_udpgal_e_"): - Mudpgale = 0
 R("R_EX_udpgal_e__Rev"): + Mudpgale = 0
 R("R_UDPG4E"): - Mudpgc + Mudpgalc = 0
 R("R_UDPG4E_Rev"): + Mudpgc - Mudpgalc = 0
 R("R_UDPGALM"): + Mudpgalfurc - Mudpgalc = 0
 R("R_UDPGALtex"): + Mudpgalp - Mudpgale = 0
 R("R_UDPGALtex_Rev"): - Mudpgalp + Mudpgale = 0
 R("Mudpgb_Transp"): 0 M10fthfc = 0
 R("R_EX_udpg_e_"): - Mudpge = 0
 R("R_EX_udpg_e__Rev"): + Mudpge = 0
 R("R_UDPGtex"): + Mudpgp - Mudpge = 0
 R("R_UDPGtex_Rev"): - Mudpgp + Mudpge = 0
 R("Mudpglcurb_Transp"): 0 M10fthfc = 0
 R("R_EX_udpglcur_e_"): - Mudpglcure = 0
 R("R_EX_udpglcur_e__Rev"): + Mudpglcure = 0
 R("R_UDPGLCURtex"): + Mudpglcurp - Mudpglcure = 0
 R("R_UDPGLCURtex_Rev"): - Mudpglcurp + Mudpglcure = 0
 R("R_ULA4Ntppi"): + MuLa4np - MuLa4nc = 0
 R("Mumpb_Transp"): 0 M10fthfc = 0
 R("R_EX_ump_e_"): - Mumpe = 0
 R("R_EX_ump_e__Rev"): + Mumpe = 0
 R("R_UMPtex"): + Mumpp - Mumpe = 0
 R("R_UMPtex_Rev"): - Mumpp + Mumpe = 0
 R("R_ECAtpp"): + Munagamufp - Munagamufc = 0
 R("Murab_Transp"): 0 M10fthfc = 0
 R("R_EX_ura_e_"): - Murae = 0
 R("R_EX_ura_e__Rev"): + Murae = 0
 R("R_URAtex"): + Murap - Murae = 0
 R("R_URAtex_Rev"): - Murap + Murae = 0
 R("Mureab_Transp"): 0 M10fthfc = 0
 R("R_EX_urea_e_"): - Mureae = 0
 R("R_EX_urea_e__Rev"): + Mureae = 0
 R("R_UREAtpp"): - Mureap + Mureac = 0
 R("R_UREAtpp_Rev"): + Mureap - Mureac = 0
 R("R_UREAtex"): + Mureap - Mureae = 0
 R("R_UREAtex_Rev"): - Mureap + Mureae = 0
 R("Murib_Transp"): 0 M10fthfc = 0
 R("R_EX_uri_e_"): - Murie = 0
 R("R_EX_uri_e__Rev"): + Murie = 0
 R("R_URItex"): + Murip - Murie = 0
 R("R_URItex_Rev"): - Murip + Murie = 0
 R("MvalLb_Transp"): 0 M10fthfc = 0
 R("R_EX_val_L_e_"): - MvalLe = 0
 R("R_EX_val_L_e__Rev"): + MvalLe = 0
 R("R_VALtex"): + MvalLp - MvalLe = 0
 R("R_VALtex_Rev"): - MvalLp + MvalLe = 0
 R("Mxanb_Transp"): 0 M10fthfc = 0
 R("R_EX_xan_e_"): - Mxane = 0
 R("R_EX_xan_e__Rev"): + Mxane = 0
 R("R_XANtpp"): - Mxanp + Mxanc = 0
 R("R_XANtpp_Rev"): + Mxanp - Mxanc = 0
 R("R_XANtex"): + Mxanp - Mxane = 0
 R("R_XANtex_Rev"): - Mxanp + Mxane = 0
 R("Mxmpb_Transp"): 0 M10fthfc = 0
 R("R_EX_xmp_e_"): - Mxmpe = 0
 R("R_EX_xmp_e__Rev"): + Mxmpe = 0
 R("R_XMPtex"): + Mxmpp - Mxmpe = 0
 R("R_XMPtex_Rev"): - Mxmpp + Mxmpe = 0
 R("Mxtsnb_Transp"): 0 M10fthfc = 0
 R("R_EX_xtsn_e_"): - Mxtsne = 0
 R("R_EX_xtsn_e__Rev"): + Mxtsne = 0
 R("R_XTSNtex"): + Mxtsnp - Mxtsne = 0
 R("R_XTSNtex_Rev"): - Mxtsnp + Mxtsne = 0
 R("MxylDb_Transp"): 0 M10fthfc = 0
 R("R_EX_xyl_D_e_"): - MxylDe = 0
 R("R_EX_xyl_D_e__Rev"): + MxylDe = 0
 R("R_XYLI1"): + MxyluDc - MxylDc = 0
 R("R_XYLI1_Rev"): - MxyluDc + MxylDc = 0
 R("R_XYLtex"): + MxylDp - MxylDe = 0
 R("R_XYLtex_Rev"): - MxylDp + MxylDe = 0
 R("MxyluLb_Transp"): 0 M10fthfc = 0
 R("R_EX_xylu_L_e_"): - MxyluLe = 0
 R("R_EX_xylu_L_e__Rev"): + MxyluLe = 0
 R("R_XYLUtex"): + MxyluLp - MxyluLe = 0
 R("R_XYLUtex_Rev"): - MxyluLp + MxyluLe = 0
 R("Mzn2b_Transp"): 0 M10fthfc = 0
 R("R_EX_zn2_e_"): - Mzn2e = 0
 R("R_EX_zn2_e__Rev"): + Mzn2e = 0
 R("R_ZN2tpp"): - Mzn2p + Mzn2c = 0
 R("R_Zn2tex"): + Mzn2p - Mzn2e = 0
 R("R_Zn2tex_Rev"): - Mzn2p + Mzn2e = 0

Bounds
 0 <= M10fthfc <= 10000
 0 <= M12dgr120c <= 10000
 0 <= M12dgr120p <= 10000
 0 <= M12dgr140c <= 10000
 0 <= M12dgr140p <= 10000
 0 <= M12dgr141c <= 10000
 0 <= M12dgr141p <= 10000
 0 <= M12dgr160c <= 10000
 0 <= M12dgr160p <= 10000
 0 <= M12dgr161c <= 10000
 0 <= M12dgr161p <= 10000
 0 <= M12dgr180c <= 10000
 0 <= M12dgr180p <= 10000
 0 <= M12dgr181c <= 10000
 0 <= M12dgr181p <= 10000
 0 <= M12ppdRc <= 10000
 0 <= M12ppdRe <= 10000
 0 <= M12ppdRp <= 10000
 0 <= M12ppdSc <= 10000
 0 <= M12ppdSe <= 10000
 0 <= M12ppdSp <= 10000
 0 <= M13dpgc <= 10000
 0 <= M14glucanc <= 10000
 0 <= M14glucane <= 10000
 0 <= M14glucanp <= 10000
 0 <= M15dapc <= 10000
 0 <= M15dape <= 10000
 0 <= M15dapp <= 10000
 0 <= M1agpe120p <= 10000
 0 <= M1agpe140p <= 10000
 0 <= M1agpe141p <= 10000
 0 <= M1agpe160p <= 10000
 0 <= M1agpe161p <= 10000
 0 <= M1agpe180p <= 10000
 0 <= M1agpe181p <= 10000
 0 <= M1agpg120p <= 10000
 0 <= M1agpg140p <= 10000
 0 <= M1agpg141p <= 10000
 0 <= M1agpg160p <= 10000
 0 <= M1agpg161p <= 10000
 0 <= M1agpg180p <= 10000
 0 <= M1agpg181p <= 10000
 0 <= M1ddecg3pc <= 10000
 0 <= M1ddecg3pp <= 10000
 0 <= M1hdec9eg3pc <= 10000
 0 <= M1hdec9eg3pp <= 10000
 0 <= M1hdecg3pc <= 10000
 0 <= M1hdecg3pp <= 10000
 0 <= M1odec11eg3pc <= 10000
 0 <= M1odec11eg3pp <= 10000
 0 <= M1odecg3pc <= 10000
 0 <= M1odecg3pp <= 10000
 0 <= M1pyr5cc <= 10000
 0 <= M1tdec7eg3pc <= 10000
 0 <= M1tdec7eg3pp <= 10000
 0 <= M1tdecg3pc <= 10000
 0 <= M1tdecg3pp <= 10000
 0 <= M23campe <= 10000
 0 <= M23campp <= 10000
 0 <= M23ccmpe <= 10000
 0 <= M23ccmpp <= 10000
 0 <= M23cgmpe <= 10000
 0 <= M23cgmpp <= 10000
 0 <= M23cumpe <= 10000
 0 <= M23cumpp <= 10000
 0 <= M23dappac <= 10000
 0 <= M23dappae <= 10000
 0 <= M23dappap <= 10000
 0 <= M23ddhbc <= 10000
 0 <= M23dhbac <= 10000
 0 <= M23dhbc <= 10000
 0 <= M23dhbzsc <= 10000
 0 <= M23dhdpc <= 10000
 0 <= M23dhmbc <= 10000
 0 <= M23dhmpc <= 10000
 0 <= M23dogulnc <= 10000
 0 <= M25aicsc <= 10000
 0 <= M25dkglcnc <= 10000
 0 <= M25drappc <= 10000
 0 <= M26dapLLc <= 10000
 0 <= M26dapMc <= 10000
 0 <= M26dapMe <= 10000
 0 <= M26dapMp <= 10000
 0 <= M2agpe120c <= 10000
 0 <= M2agpe120p <= 10000
 0 <= M2agpe140c <= 10000
 0 <= M2agpe140p <= 10000
 0 <= M2agpe141c <= 10000
 0 <= M2agpe141p <= 10000
 0 <= M2agpe160c <= 10000
 0 <= M2agpe160p <= 10000
 0 <= M2agpe161c <= 10000
 0 <= M2agpe161p <= 10000
 0 <= M2agpe180c <= 10000
 0 <= M2agpe180p <= 10000
 0 <= M2agpe181c <= 10000
 0 <= M2agpe181p <= 10000
 0 <= M2agpg120c <= 10000
 0 <= M2agpg120p <= 10000
 0 <= M2agpg140c <= 10000
 0 <= M2agpg140p <= 10000
 0 <= M2agpg141c <= 10000
 0 <= M2agpg141p <= 10000
 0 <= M2agpg160c <= 10000
 0 <= M2agpg160p <= 10000
 0 <= M2agpg161c <= 10000
 0 <= M2agpg161p <= 10000
 0 <= M2agpg180c <= 10000
 0 <= M2agpg180p <= 10000
 0 <= M2agpg181c <= 10000
 0 <= M2agpg181p <= 10000
 0 <= M2ahbutc <= 10000
 0 <= M2amsac <= 10000
 0 <= M2aobutc <= 10000
 0 <= M2cpr5pc <= 10000
 0 <= M2dda7pc <= 10000
 0 <= M2ddecg3pc <= 10000
 0 <= M2ddecg3pp <= 10000
 0 <= M2ddg6pc <= 10000
 0 <= M2ddglcnc <= 10000
 0 <= M2ddglcne <= 10000
 0 <= M2ddglcnp <= 10000
 0 <= M2dh3dgal6pc <= 10000
 0 <= M2dh3dgalc <= 10000
 0 <= M2dhglcnc <= 10000
 0 <= M2dhgulnc <= 10000
 0 <= M2dhpc <= 10000
 0 <= M2dmmq8c <= 10000
 0 <= M2dmmql8c <= 10000
 0 <= M2dr1pc <= 10000
 0 <= M2dr5pc <= 10000
 0 <= M2h3oppanc <= 10000
 0 <= M2hdec9eg3pc <= 10000
 0 <= M2hdec9eg3pp <= 10000
 0 <= M2hdecg3pc <= 10000
 0 <= M2hdecg3pp <= 10000
 0 <= M2ippmc <= 10000
 0 <= M2kmbc <= 10000
 0 <= M2mahmpc <= 10000
 0 <= M2mcacnc <= 10000
 0 <= M2mcitc <= 10000
 0 <= M2me4pc <= 10000
 0 <= M2mecdpc <= 10000
 0 <= M2obutc <= 10000
 0 <= M2odec11eg3pc <= 10000
 0 <= M2odec11eg3pp <= 10000
 0 <= M2odecg3pc <= 10000
 0 <= M2odecg3pp <= 10000
 0 <= M2ohphc <= 10000
 0 <= M2ombzlc <= 10000
 0 <= M2omhmblc <= 10000
 0 <= M2ommblc <= 10000
 0 <= M2omphc <= 10000
 0 <= M2ophc <= 10000
 0 <= M2p4c2mec <= 10000
 0 <= M2pgc <= 10000
 0 <= M2pglycc <= 10000
 0 <= M2shchcc <= 10000
 0 <= M2tdec7eg3pc <= 10000
 0 <= M2tdec7eg3pp <= 10000
 0 <= M2tdecg3pc <= 10000
 0 <= M2tdecg3pp <= 10000
 0 <= M2tpr3dpcoac <= 10000
 0 <= M34dhpace <= 10000
 0 <= M34dhpacp <= 10000
 0 <= M34hppc <= 10000
 0 <= M3ampe <= 10000
 0 <= M3ampp <= 10000
 0 <= M3c2hmpc <= 10000
 0 <= M3c3hmpc <= 10000
 0 <= M3c4mopc <= 10000
 0 <= M3cmpe <= 10000
 0 <= M3cmpp <= 10000
 0 <= M3dhgulnc <= 10000
 0 <= M3dhgulnpc <= 10000
 0 <= M3dhqc <= 10000
 0 <= M3dhskc <= 10000
 0 <= M3gmpe <= 10000
 0 <= M3gmpp <= 10000
 0 <= M3haACPc <= 10000
 0 <= M3hbcoac <= 10000
 0 <= M3hcddec5eACPc <= 10000
 0 <= M3hcinnmc <= 10000
 0 <= M3hcinnme <= 10000
 0 <= M3hcinnmp <= 10000
 0 <= M3hcmrs7eACPc <= 10000
 0 <= M3hcpalm9eACPc <= 10000
 0 <= M3hcvac11eACPc <= 10000
 0 <= M3hdcoac <= 10000
 0 <= M3hddcoac <= 10000
 0 <= M3hddecACPc <= 10000
 0 <= M3hdecACPc <= 10000
 0 <= M3hhcoac <= 10000
 0 <= M3hhdcoac <= 10000
 0 <= M3hhexACPc <= 10000
 0 <= M3hmrsACPc <= 10000
 0 <= M3hocoac <= 10000
 0 <= M3hoctaACPc <= 10000
 0 <= M3hoctACPc <= 10000
 0 <= M3hodcoac <= 10000
 0 <= M3hpalmACPc <= 10000
 0 <= M3hpppnc <= 10000
 0 <= M3hpppne <= 10000
 0 <= M3hpppnp <= 10000
 0 <= M3htdcoac <= 10000
 0 <= M3ig3pc <= 10000
 0 <= M3mobc <= 10000
 0 <= M3mopc <= 10000
 0 <= M3ocddec5eACPc <= 10000
 0 <= M3ocmrs7eACPc <= 10000
 0 <= M3ocpalm9eACPc <= 10000
 0 <= M3ocvac11eACPc <= 10000
 0 <= M3odcoac <= 10000
 0 <= M3oddcoac <= 10000
 0 <= M3oddecACPc <= 10000
 0 <= M3odecACPc <= 10000
 0 <= M3ohcoac <= 10000
 0 <= M3ohdcoac <= 10000
 0 <= M3ohexACPc <= 10000
 0 <= M3omrsACPc <= 10000
 0 <= M3oocoac <= 10000
 0 <= M3ooctACPc <= 10000
 0 <= M3ooctdACPc <= 10000
 0 <= M3oodcoac <= 10000
 0 <= M3opalmACPc <= 10000
 0 <= M3ophbc <= 10000
 0 <= M3otdcoac <= 10000
 0 <= M3pgc <= 10000
 0 <= M3phpc <= 10000
 0 <= M3psmec <= 10000
 0 <= M3salac <= 10000
 0 <= M3umpe <= 10000
 0 <= M3umpp <= 10000
 0 <= M4abutc <= 10000
 0 <= M4abute <= 10000
 0 <= M4abutnc <= 10000
 0 <= M4abutp <= 10000
 0 <= M4abzc <= 10000
 0 <= M4adchoc <= 10000
 0 <= M4ahmmpc <= 10000
 0 <= M4ampmc <= 10000
 0 <= M4c2mec <= 10000
 0 <= M4h2opntnc <= 10000
 0 <= M4hbac <= 10000
 0 <= M4hbzc <= 10000
 0 <= M4hoxpacde <= 10000
 0 <= M4hoxpacdp <= 10000
 0 <= M4hthrc <= 10000
 0 <= M4mhetzc <= 10000
 0 <= M4mopc <= 10000
 0 <= M4mpetzc <= 10000
 0 <= M4paspc <= 10000
 0 <= M4perc <= 10000
 0 <= M4ppanc <= 10000
 0 <= M4ppcysc <= 10000
 0 <= M4r5auc <= 10000
 0 <= M5aizcc <= 10000
 0 <= M5aopc <= 10000
 0 <= M5aprbuc <= 10000
 0 <= M5apruc <= 10000
 0 <= M5caizc <= 10000
 0 <= M5dglcnc <= 10000
 0 <= M5dglcne <= 10000
 0 <= M5dglcnp <= 10000
 0 <= M5dh4dglcc <= 10000
 0 <= M5dribc <= 10000
 0 <= M5mdr1pc <= 10000
 0 <= M5mdru1pc <= 10000
 0 <= M5mtac <= 10000
 0 <= M5mthfc <= 10000
 0 <= M5mtrc <= 10000
 0 <= M5prdmbzc <= 10000
 0 <= M6hmhptc <= 10000
 0 <= M6hmhptppc <= 10000
 0 <= M6pgcc <= 10000
 0 <= M6pglc <= 10000
 0 <= M8aonnc <= 10000
 0 <= Maacaldc <= 10000
 0 <= Maacoac <= 10000
 0 <= Maactc <= 10000
 0 <= Macacc <= 10000
 0 <= Macace <= 10000
 0 <= Macacp <= 10000
 0 <= MacACPc <= 10000
 0 <= Macaldc <= 10000
 0 <= Macalde <= 10000
 0 <= Macaldp <= 10000
 0 <= Macanthc <= 10000
 0 <= Macc <= 10000
 0 <= Maccoac <= 10000
 0 <= Mace <= 10000
 0 <= Macetolc <= 10000
 0 <= Macg5pc <= 10000
 0 <= Macg5sac <= 10000
 0 <= Macgal1pe <= 10000
 0 <= Macgal1pp <= 10000
 0 <= Macgale <= 10000
 0 <= Macgalp <= 10000
 0 <= Macgam1pc <= 10000
 0 <= Macgam1pe <= 10000
 0 <= Macgam1pp <= 10000
 0 <= Macgam6pc <= 10000
 0 <= Macgamc <= 10000
 0 <= Macgame <= 10000
 0 <= Macgamp <= 10000
 0 <= MacglcDc <= 10000
 0 <= Macgluc <= 10000
 0 <= Macmaltc <= 10000
 0 <= Macmanac <= 10000
 0 <= Macmanae <= 10000
 0 <= Macmanap <= 10000
 0 <= Macmanapc <= 10000
 0 <= Macmum6pc <= 10000
 0 <= Macmume <= 10000
 0 <= Macmump <= 10000
 0 <= Macnamc <= 10000
 0 <= Macname <= 10000
 0 <= Macnamp <= 10000
 0 <= Macolipae <= 10000
 0 <= Macolipap <= 10000
 0 <= MaconCc <= 10000
 0 <= Maconmc <= 10000
 0 <= MaconTc <= 10000
 0 <= Macornc <= 10000
 0 <= Macp <= 10000
 0 <= MACPc <= 10000
 0 <= Macserc <= 10000
 0 <= Macsere <= 10000
 0 <= Macserp <= 10000
 0 <= MactACPc <= 10000
 0 <= Mactpc <= 10000
 0 <= Madec <= 10000
 0 <= Madee <= 10000
 0 <= Madep <= 10000
 0 <= Madnc <= 10000
 0 <= Madne <= 10000
 0 <= Madnp <= 10000
 0 <= Madocbic <= 10000
 0 <= Madocbipc <= 10000
 0 <= Madocblc <= 10000
 0 <= Madocble <= 10000
 0 <= Madocblp <= 10000
 0 <= Madpc <= 10000
 0 <= Madpglcc <= 10000
 0 <= MadphepDDc <= 10000
 0 <= MadphepLDc <= 10000
 0 <= Madpribc <= 10000
 0 <= Magc <= 10000
 0 <= Magdpcbic <= 10000
 0 <= Mage <= 10000
 0 <= Magmc <= 10000
 0 <= Magme <= 10000
 0 <= Magmp <= 10000
 0 <= Mahcysc <= 10000
 0 <= Mahdtc <= 10000
 0 <= Maicarc <= 10000
 0 <= Mairc <= 10000
 0 <= Makgc <= 10000
 0 <= Makge <= 10000
 0 <= Makgp <= 10000
 0 <= Malaalac <= 10000
 0 <= Malaalae <= 10000
 0 <= Malaalap <= 10000
 0 <= MalaBc <= 10000
 0 <= MalaBe <= 10000
 0 <= MalaBp <= 10000
 0 <= MalacSc <= 10000
 0 <= MalaDc <= 10000
 0 <= MalaDe <= 10000
 0 <= MalaDp <= 10000
 0 <= MalaLc <= 10000
 0 <= MalaLe <= 10000
 0 <= MalaLp <= 10000
 0 <= Malatrnac <= 10000
 0 <= Mall6pc <= 10000
 0 <= MallDc <= 10000
 0 <= MallDe <= 10000
 0 <= MallDp <= 10000
 0 <= Malltnc <= 10000
 0 <= Malltne <= 10000
 0 <= Malltnp <= 10000
 0 <= Mallttc <= 10000
 0 <= Mallul6pc <= 10000
 0 <= Malppp <= 10000
 0 <= Maltrnc <= 10000
 0 <= Mametamc <= 10000
 0 <= Mametc <= 10000
 0 <= Mamobc <= 10000
 0 <= Mampc <= 10000
 0 <= Mampe <= 10000
 0 <= Mampp <= 10000
 0 <= Manhgm3pc <= 10000
 0 <= Manhgm3pp <= 10000
 0 <= Manhgm4pc <= 10000
 0 <= Manhgm4pp <= 10000
 0 <= Manhgmc <= 10000
 0 <= Manhgme <= 10000
 0 <= Manhgmp <= 10000
 0 <= Manhm3pc <= 10000
 0 <= Manhm4pc <= 10000
 0 <= Manhmc <= 10000
 0 <= Manthc <= 10000
 0 <= Map4ac <= 10000
 0 <= Map5ac <= 10000
 0 <= Mapg120c <= 10000
 0 <= Mapg140c <= 10000
 0 <= Mapg141c <= 10000
 0 <= Mapg160c <= 10000
 0 <= Mapg161c <= 10000
 0 <= Mapg180c <= 10000
 0 <= Mapg181c <= 10000
 0 <= MapoACPc <= 10000
 0 <= Mapsc <= 10000
 0 <= Mara5pc <= 10000
 0 <= MarabLc <= 10000
 0 <= MarabLe <= 10000
 0 <= MarabLp <= 10000
 0 <= Maragundc <= 10000
 0 <= Marbt6pc <= 10000
 0 <= Marbtnc <= 10000
 0 <= Marbtne <= 10000
 0 <= Marbtnfe3c <= 10000
 0 <= Marbtnfe3e <= 10000
 0 <= Marbtnfe3p <= 10000
 0 <= Marbtnp <= 10000
 0 <= MargLc <= 10000
 0 <= MargLe <= 10000
 0 <= MargLp <= 10000
 0 <= Margsucc <= 10000
 0 <= Margtrnac <= 10000
 0 <= Mascb6pc <= 10000
 0 <= MascbLe <= 10000
 0 <= MascbLp <= 10000
 0 <= MasnLc <= 10000
 0 <= MasnLe <= 10000
 0 <= MasnLp <= 10000
 0 <= Masntrnac <= 10000
 0 <= Maso3c <= 10000
 0 <= Maso3e <= 10000
 0 <= Maso3p <= 10000
 0 <= Maso4c <= 10000
 0 <= MaspLc <= 10000
 0 <= MaspLe <= 10000
 0 <= MaspLp <= 10000
 0 <= Maspsac <= 10000
 0 <= Masptrnac <= 10000
 0 <= MathrLc <= 10000
 0 <= Matpc <= 10000
 0 <= Mb2coac <= 10000
 0 <= Mbbtcoac <= 10000
 0 <= Mbetaldc <= 10000
 0 <= Mbglycogenc <= 10000
 0 <= Mbtalc <= 10000
 0 <= Mbtcoac <= 10000
 0 <= Mbtnc <= 10000
 0 <= Mbtnsoc <= 10000
 0 <= Mbut2eACPc <= 10000
 0 <= MbutACPc <= 10000
 0 <= Mbutc <= 10000
 0 <= Mbute <= 10000
 0 <= Mbutp <= 10000
 0 <= Mbutso3c <= 10000
 0 <= Mbutso3e <= 10000
 0 <= Mbutso3p <= 10000
 0 <= Mca2c <= 10000
 0 <= Mca2e <= 10000
 0 <= Mca2p <= 10000
 0 <= Mcampc <= 10000
 0 <= Mcbaspc <= 10000
 0 <= Mcbic <= 10000
 0 <= Mcbie <= 10000
 0 <= Mcbip <= 10000
 0 <= Mcbl1c <= 10000
 0 <= Mcbl1e <= 10000
 0 <= Mcbl1p <= 10000
 0 <= Mcbpc <= 10000
 0 <= Mcd2c <= 10000
 0 <= Mcd2e <= 10000
 0 <= Mcd2p <= 10000
 0 <= Mcddec5eACPc <= 10000
 0 <= Mcdec3eACPc <= 10000
 0 <= Mcdpc <= 10000
 0 <= Mcdpdddecgc <= 10000
 0 <= Mcdpdhdec9egc <= 10000
 0 <= Mcdpdhdecgc <= 10000
 0 <= Mcdpdodec11egc <= 10000
 0 <= Mcdpdodecgc <= 10000
 0 <= Mcdpdtdec7egc <= 10000
 0 <= Mcdpdtdecgc <= 10000
 0 <= Mcechdddc <= 10000
 0 <= Mcenchdddc <= 10000
 0 <= Mcglyc <= 10000
 0 <= Mcglye <= 10000
 0 <= Mcglyp <= 10000
 0 <= Mcholc <= 10000
 0 <= Mchole <= 10000
 0 <= Mcholp <= 10000
 0 <= Mchorc <= 10000
 0 <= Mcinnmc <= 10000
 0 <= Mcitc <= 10000
 0 <= Mcite <= 10000
 0 <= Mcitp <= 10000
 0 <= McitrLc <= 10000
 0 <= Mckdoc <= 10000
 0 <= Mclc <= 10000
 0 <= Mcle <= 10000
 0 <= Mclp <= 10000
 0 <= Mclpn120p <= 10000
 0 <= Mclpn140p <= 10000
 0 <= Mclpn141p <= 10000
 0 <= Mclpn160p <= 10000
 0 <= Mclpn161p <= 10000
 0 <= Mclpn180p <= 10000
 0 <= Mclpn181p <= 10000
 0 <= Mcmpc <= 10000
 0 <= Mcmpe <= 10000
 0 <= Mcmpp <= 10000
 0 <= Mco2c <= 10000
 0 <= Mco2e <= 10000
 0 <= Mco2p <= 10000
 0 <= Mcoac <= 10000
 0 <= Mcobalt2c <= 10000
 0 <= Mcobalt2e <= 10000
 0 <= Mcobalt2p <= 10000
 0 <= Mcolipac <= 10000
 0 <= Mcolipae <= 10000
 0 <= Mcolipap <= 10000
 0 <= Mcpe160c <= 10000
 0 <= Mcpe180c <= 10000
 0 <= Mcpg160c <= 10000
 0 <= Mcpg180c <= 10000
 0 <= Mcpgnc <= 10000
 0 <= Mcpgne <= 10000
 0 <= Mcpgnp <= 10000
 0 <= Mcpgnunc <= 10000
 0 <= Mcpgnune <= 10000
 0 <= Mcpgnunp <= 10000
 0 <= Mcpppg3c <= 10000
 0 <= Mcrnc <= 10000
 0 <= Mcrncoac <= 10000
 0 <= McrnDc <= 10000
 0 <= McrnDcoac <= 10000
 0 <= McrnDp <= 10000
 0 <= Mcrne <= 10000
 0 <= Mcrnp <= 10000
 0 <= Mcsnc <= 10000
 0 <= Mcsne <= 10000
 0 <= Mcsnp <= 10000
 0 <= Mctbtc <= 10000
 0 <= Mctbtcoac <= 10000
 0 <= Mctbtp <= 10000
 0 <= Mctpc <= 10000
 0 <= Mcu2c <= 10000
 0 <= Mcu2e <= 10000
 0 <= Mcu2p <= 10000
 0 <= Mcuc <= 10000
 0 <= Mcue <= 10000
 0 <= Mcup <= 10000
 0 <= Mcyanc <= 10000
 0 <= Mcyane <= 10000
 0 <= Mcyanp <= 10000
 0 <= Mcyntc <= 10000
 0 <= Mcynte <= 10000
 0 <= Mcyntp <= 10000
 0 <= McysDc <= 10000
 0 <= McysDe <= 10000
 0 <= McysDp <= 10000
 0 <= McysLc <= 10000
 0 <= McysLe <= 10000
 0 <= McysLp <= 10000
 0 <= McystLc <= 10000
 0 <= Mcystrnac <= 10000
 0 <= Mcytdc <= 10000
 0 <= Mcytde <= 10000
 0 <= Mcytdp <= 10000
 0 <= Mdad2c <= 10000
 0 <= Mdad2e <= 10000
 0 <= Mdad2p <= 10000
 0 <= Mdad5c <= 10000
 0 <= Mdadpc <= 10000
 0 <= Mdampc <= 10000
 0 <= Mdampe <= 10000
 0 <= Mdampp <= 10000
 0 <= Mdannc <= 10000
 0 <= Mdatpc <= 10000
 0 <= Mdb4pc <= 10000
 0 <= Mdc2coac <= 10000
 0 <= MdcaACPc <= 10000
 0 <= Mdcac <= 10000
 0 <= Mdcacoac <= 10000
 0 <= Mdcae <= 10000
 0 <= Mdcampc <= 10000
 0 <= Mdcap <= 10000
 0 <= Mdcdpc <= 10000
 0 <= Mdcmpc <= 10000
 0 <= Mdcmpe <= 10000
 0 <= Mdcmpp <= 10000
 0 <= Mdctpc <= 10000
 0 <= Mdcytc <= 10000
 0 <= Mdcyte <= 10000
 0 <= Mdcytp <= 10000
 0 <= Mdd2coac <= 10000
 0 <= MddcaACPc <= 10000
 0 <= Mddcac <= 10000
 0 <= Mddcacoac <= 10000
 0 <= Mddcae <= 10000
 0 <= Mddcap <= 10000
 0 <= Mdgdpc <= 10000
 0 <= Mdgmpc <= 10000
 0 <= Mdgmpe <= 10000
 0 <= Mdgmpp <= 10000
 0 <= Mdgsnc <= 10000
 0 <= Mdgsne <= 10000
 0 <= Mdgsnp <= 10000
 0 <= Mdgtpc <= 10000
 0 <= Mdhac <= 10000
 0 <= Mdhae <= 10000
 0 <= Mdhap <= 10000
 0 <= Mdhapc <= 10000
 0 <= Mdhcinnmc <= 10000
 0 <= Mdhfc <= 10000
 0 <= Mdhmptpc <= 10000
 0 <= Mdhnac <= 10000
 0 <= Mdhnptc <= 10000
 0 <= MdhorSc <= 10000
 0 <= Mdhpmpc <= 10000
 0 <= Mdhpppnc <= 10000
 0 <= Mdhptc <= 10000
 0 <= Mdhptdc <= 10000
 0 <= Mdidpc <= 10000
 0 <= Mdimpc <= 10000
 0 <= Mdimpe <= 10000
 0 <= Mdimpp <= 10000
 0 <= Mdinc <= 10000
 0 <= Mdine <= 10000
 0 <= Mdinp <= 10000
 0 <= Mditpc <= 10000
 0 <= Mdkmppc <= 10000
 0 <= Mdmbzidc <= 10000
 0 <= Mdmlzc <= 10000
 0 <= Mdmppc <= 10000
 0 <= Mdmsc <= 10000
 0 <= Mdmse <= 10000
 0 <= Mdmsoc <= 10000
 0 <= Mdmsoe <= 10000
 0 <= Mdmsop <= 10000
 0 <= Mdmsp <= 10000
 0 <= Mdnadc <= 10000
 0 <= Mdopae <= 10000
 0 <= Mdopap <= 10000
 0 <= Mdpcoac <= 10000
 0 <= Mdsbaoxp <= 10000
 0 <= Mdsbardp <= 10000
 0 <= Mdsbcoxp <= 10000
 0 <= Mdsbcrdp <= 10000
 0 <= Mdsbdoxc <= 10000
 0 <= Mdsbdrdc <= 10000
 0 <= Mdsbgoxp <= 10000
 0 <= Mdsbgrdp <= 10000
 0 <= Mdsclc <= 10000
 0 <= Mdtbtc <= 10000
 0 <= Mdtdp4aaddgc <= 10000
 0 <= Mdtdp4addgc <= 10000
 0 <= Mdtdp4d6dgc <= 10000
 0 <= Mdtdp4d6dmc <= 10000
 0 <= Mdtdpc <= 10000
 0 <= Mdtdpgluc <= 10000
 0 <= Mdtdprmnc <= 10000
 0 <= Mdtmpc <= 10000
 0 <= Mdtmpe <= 10000
 0 <= Mdtmpp <= 10000
 0 <= Mdttpc <= 10000
 0 <= Mdudpc <= 10000
 0 <= Mdumpc <= 10000
 0 <= Mdumpe <= 10000
 0 <= Mdumpp <= 10000
 0 <= Mduric <= 10000
 0 <= Mdurie <= 10000
 0 <= Mdurip <= 10000
 0 <= Mdutpc <= 10000
 0 <= Mdxyl5pc <= 10000
 0 <= Mdxylc <= 10000
 0 <= Me4pc <= 10000
 0 <= Meca2undp <= 10000
 0 <= Meca3undp <= 10000
 0 <= Meca4colipae <= 10000
 0 <= Meca4colipap <= 10000
 0 <= Meca4undp <= 10000
 0 <= Meig3pc <= 10000
 0 <= Menlipae <= 10000
 0 <= Menlipap <= 10000
 0 <= Menterc <= 10000
 0 <= Mentere <= 10000
 0 <= Menterp <= 10000
 0 <= Methac <= 10000
 0 <= Methae <= 10000
 0 <= Methap <= 10000
 0 <= Methso3c <= 10000
 0 <= Methso3e <= 10000
 0 <= Methso3p <= 10000
 0 <= Metohc <= 10000
 0 <= Metohe <= 10000
 0 <= Metohp <= 10000
 0 <= Mf1pc <= 10000
 0 <= Mf6pc <= 10000
 0 <= Mf6pe <= 10000
 0 <= Mf6pp <= 10000
 0 <= Mfadc <= 10000
 0 <= Mfadh2c <= 10000
 0 <= Mfaldc <= 10000
 0 <= Mfalde <= 10000
 0 <= Mfaldp <= 10000
 0 <= Mfc1pc <= 10000
 0 <= MfclLc <= 10000
 0 <= Mfdpc <= 10000
 0 <= Mfe2c <= 10000
 0 <= Mfe2e <= 10000
 0 <= Mfe2p <= 10000
 0 <= Mfe3c <= 10000
 0 <= Mfe3dcite <= 10000
 0 <= Mfe3dcitp <= 10000
 0 <= Mfe3dhbzsc <= 10000
 0 <= Mfe3dhbzse <= 10000
 0 <= Mfe3dhbzsp <= 10000
 0 <= Mfe3e <= 10000
 0 <= Mfe3hoxc <= 10000
 0 <= Mfe3hoxe <= 10000
 0 <= Mfe3hoxp <= 10000
 0 <= Mfe3hoxunc <= 10000
 0 <= Mfe3hoxune <= 10000
 0 <= Mfe3hoxunp <= 10000
 0 <= Mfe3p <= 10000
 0 <= Mfecrmc <= 10000
 0 <= Mfecrme <= 10000
 0 <= Mfecrmp <= 10000
 0 <= Mfecrmunc <= 10000
 0 <= Mfecrmune <= 10000
 0 <= Mfecrmunp <= 10000
 0 <= Mfeenterc <= 10000
 0 <= Mfeentere <= 10000
 0 <= Mfeenterp <= 10000
 0 <= Mfeoxamc <= 10000
 0 <= Mfeoxame <= 10000
 0 <= Mfeoxamp <= 10000
 0 <= Mfeoxamunc <= 10000
 0 <= Mfeoxamune <= 10000
 0 <= Mfeoxamunp <= 10000
 0 <= Mfgamc <= 10000
 0 <= Mfldoxc <= 10000
 0 <= Mfldrdc <= 10000
 0 <= Mfmettrnac <= 10000
 0 <= Mfmnc <= 10000
 0 <= Mfmnh2c <= 10000
 0 <= Mforc <= 10000
 0 <= Mfore <= 10000
 0 <= Mforp <= 10000
 0 <= Mfpramc <= 10000
 0 <= Mfpricac <= 10000
 0 <= Mfrdpc <= 10000
 0 <= Mfruc <= 10000
 0 <= Mfrue <= 10000
 0 <= Mfrulysc <= 10000
 0 <= Mfrulyse <= 10000
 0 <= Mfrulysp <= 10000
 0 <= Mfrulyspc <= 10000
 0 <= Mfrup <= 10000
 0 <= Mfruurc <= 10000
 0 <= Mfruure <= 10000
 0 <= Mfruurp <= 10000
 0 <= MfucLc <= 10000
 0 <= MfucLe <= 10000
 0 <= MfucLp <= 10000
 0 <= Mfumc <= 10000
 0 <= Mfume <= 10000
 0 <= Mfump <= 10000
 0 <= Mg1pc <= 10000
 0 <= Mg1pe <= 10000
 0 <= Mg1pp <= 10000
 0 <= Mg3pc <= 10000
 0 <= Mg3pcc <= 10000
 0 <= Mg3pce <= 10000
 0 <= Mg3pcp <= 10000
 0 <= Mg3pec <= 10000
 0 <= Mg3pee <= 10000
 0 <= Mg3pep <= 10000
 0 <= Mg3pgc <= 10000
 0 <= Mg3pge <= 10000
 0 <= Mg3pgp <= 10000
 0 <= Mg3pic <= 10000
 0 <= Mg3pie <= 10000
 0 <= Mg3pip <= 10000
 0 <= Mg3psc <= 10000
 0 <= Mg3pse <= 10000
 0 <= Mg3psp <= 10000
 0 <= Mg6pc <= 10000
 0 <= Mg6pe <= 10000
 0 <= Mg6pp <= 10000
 0 <= Mgagicolipac <= 10000
 0 <= Mgal1pc <= 10000
 0 <= Mgal1pe <= 10000
 0 <= Mgal1pp <= 10000
 0 <= MgalbDe <= 10000
 0 <= MgalbDp <= 10000
 0 <= Mgalc <= 10000
 0 <= MgalctDc <= 10000
 0 <= MgalctDe <= 10000
 0 <= MgalctDp <= 10000
 0 <= MgalctnDc <= 10000
 0 <= MgalctnDe <= 10000
 0 <= MgalctnDp <= 10000
 0 <= MgalctnLc <= 10000
 0 <= MgalctnLe <= 10000
 0 <= MgalctnLp <= 10000
 0 <= Mgale <= 10000
 0 <= Mgalp <= 10000
 0 <= Mgalt1pc <= 10000
 0 <= Mgalte <= 10000
 0 <= Mgaltp <= 10000
 0 <= Mgalurc <= 10000
 0 <= Mgalure <= 10000
 0 <= Mgalurp <= 10000
 0 <= Mgam1pc <= 10000
 0 <= Mgam6pc <= 10000
 0 <= Mgam6pe <= 10000
 0 <= Mgam6pp <= 10000
 0 <= Mgame <= 10000
 0 <= Mgamp <= 10000
 0 <= Mgaragundc <= 10000
 0 <= Mgarc <= 10000
 0 <= Mgbbtnc <= 10000
 0 <= Mgbbtne <= 10000
 0 <= Mgbbtnp <= 10000
 0 <= Mgcaldc <= 10000
 0 <= Mgdpc <= 10000
 0 <= Mgdpddmanc <= 10000
 0 <= Mgdpe <= 10000
 0 <= Mgdpfucc <= 10000
 0 <= Mgdpmannc <= 10000
 0 <= Mgdpofucc <= 10000
 0 <= Mgdpp <= 10000
 0 <= Mgdptpc <= 10000
 0 <= Mgfgaragundc <= 10000
 0 <= Mgg4abutc <= 10000
 0 <= Mggagicolipac <= 10000
 0 <= Mggbutalc <= 10000
 0 <= Mgggagicolipac <= 10000
 0 <= Mggptrcc <= 10000
 0 <= Mgicolipac <= 10000
 0 <= MglcDc <= 10000
 0 <= MglcDe <= 10000
 0 <= MglcDp <= 10000
 0 <= Mglcnc <= 10000
 0 <= Mglcne <= 10000
 0 <= Mglcnp <= 10000
 0 <= Mglcrc <= 10000
 0 <= Mglcre <= 10000
 0 <= Mglcrp <= 10000
 0 <= Mglcur1pe <= 10000
 0 <= Mglcur1pp <= 10000
 0 <= Mglcurc <= 10000
 0 <= Mglcure <= 10000
 0 <= Mglcurp <= 10000
 0 <= MglnLc <= 10000
 0 <= MglnLe <= 10000
 0 <= MglnLp <= 10000
 0 <= Mglntrnac <= 10000
 0 <= Mglu1sac <= 10000
 0 <= Mglu5pc <= 10000
 0 <= Mglu5sac <= 10000
 0 <= Mglucysc <= 10000
 0 <= MgluDc <= 10000
 0 <= MgluLc <= 10000
 0 <= MgluLe <= 10000
 0 <= MgluLp <= 10000
 0 <= Mglutrnac <= 10000
 0 <= Mglxc <= 10000
 0 <= Mglyaldc <= 10000
 0 <= Mglyalde <= 10000
 0 <= Mglyaldp <= 10000
 0 <= Mglybc <= 10000
 0 <= Mglybe <= 10000
 0 <= Mglybp <= 10000
 0 <= Mglyc <= 10000
 0 <= Mglyc2pe <= 10000
 0 <= Mglyc2pp <= 10000
 0 <= Mglyc3pc <= 10000
 0 <= Mglyc3pe <= 10000
 0 <= Mglyc3pp <= 10000
 0 <= Mglycc <= 10000
 0 <= Mglyce <= 10000
 0 <= Mglycltc <= 10000
 0 <= Mglyclte <= 10000
 0 <= Mglycltp <= 10000
 0 <= Mglycogenc <= 10000
 0 <= Mglycp <= 10000
 0 <= MglycRc <= 10000
 0 <= MglycRe <= 10000
 0 <= MglycRp <= 10000
 0 <= Mglye <= 10000
 0 <= Mglyp <= 10000
 0 <= Mglytrnac <= 10000
 0 <= Mgmhep17bpc <= 10000
 0 <= Mgmhep1pc <= 10000
 0 <= Mgmhep7pc <= 10000
 0 <= Mgmpc <= 10000
 0 <= Mgmpe <= 10000
 0 <= Mgmpp <= 10000
 0 <= Mgp4gc <= 10000
 0 <= Mgrdpc <= 10000
 0 <= Mgrxoxc <= 10000
 0 <= Mgrxrdc <= 10000
 0 <= Mgsnc <= 10000
 0 <= Mgsne <= 10000
 0 <= Mgsnp <= 10000
 0 <= Mgthoxc <= 10000
 0 <= Mgthoxe <= 10000
 0 <= Mgthoxp <= 10000
 0 <= Mgthrdc <= 10000
 0 <= Mgthrde <= 10000
 0 <= Mgthrdp <= 10000
 0 <= Mgtpc <= 10000
 0 <= Mgtpe <= 10000
 0 <= Mgtpp <= 10000
 0 <= Mgtspmdc <= 10000
 0 <= Mguac <= 10000
 0 <= Mguae <= 10000
 0 <= Mguap <= 10000
 0 <= Mh2c <= 10000
 0 <= Mh2e <= 10000
 0 <= Mh2mb4pc <= 10000
 0 <= Mh2o2c <= 10000
 0 <= Mh2o2e <= 10000
 0 <= Mh2o2p <= 10000
 0 <= Mh2oc <= 10000
 0 <= Mh2oe <= 10000
 0 <= Mh2op <= 10000
 0 <= Mh2p <= 10000
 0 <= Mh2sc <= 10000
 0 <= Mh2se <= 10000
 0 <= Mh2sp <= 10000
 0 <= Mhacolipae <= 10000
 0 <= Mhalipae <= 10000
 0 <= Mhc <= 10000
 0 <= Mhco3c <= 10000
 0 <= MhcysLc <= 10000
 0 <= Mhdcac <= 10000
 0 <= Mhdcae <= 10000
 0 <= Mhdcap <= 10000
 0 <= Mhdceac <= 10000
 0 <= Mhdceae <= 10000
 0 <= Mhdceap <= 10000
 0 <= Mhdcoac <= 10000
 0 <= Mhdd2coac <= 10000
 0 <= MhdeACPc <= 10000
 0 <= Mhe <= 10000
 0 <= MhemeOc <= 10000
 0 <= MhexACPc <= 10000
 0 <= Mhg2c <= 10000
 0 <= Mhg2e <= 10000
 0 <= Mhg2p <= 10000
 0 <= Mhhlipac <= 10000
 0 <= MhisLc <= 10000
 0 <= MhisLe <= 10000
 0 <= MhisLp <= 10000
 0 <= Mhispc <= 10000
 0 <= Mhistdc <= 10000
 0 <= Mhistrnac <= 10000
 0 <= Mhknddc <= 10000
 0 <= Mhkntdc <= 10000
 0 <= Mhlipac <= 10000
 0 <= Mhmbilc <= 10000
 0 <= Mhmfurnc <= 10000
 0 <= Mhmgthc <= 10000
 0 <= MhomLc <= 10000
 0 <= MhomLe <= 10000
 0 <= MhomLp <= 10000
 0 <= Mhp <= 10000
 0 <= Mhphhlipac <= 10000
 0 <= Mhpyrc <= 10000
 0 <= Mhqnc <= 10000
 0 <= Mhx2coac <= 10000
 0 <= Mhxac <= 10000
 0 <= Mhxae <= 10000
 0 <= Mhxanc <= 10000
 0 <= Mhxane <= 10000
 0 <= Mhxanp <= 10000
 0 <= Mhxap <= 10000
 0 <= Mhxcoac <= 10000
 0 <= Miaspc <= 10000
 0 <= Michorc <= 10000
 0 <= Micitc <= 10000
 0 <= Micolipac <= 10000
 0 <= MidonLc <= 10000
 0 <= MidonLe <= 10000
 0 <= MidonLp <= 10000
 0 <= Midpc <= 10000
 0 <= MileLc <= 10000
 0 <= MileLe <= 10000
 0 <= MileLp <= 10000
 0 <= Miletrnac <= 10000
 0 <= Mimacpc <= 10000
 0 <= Mimpc <= 10000
 0 <= Mimpe <= 10000
 0 <= Mimpp <= 10000
 0 <= Mindolec <= 10000
 0 <= Mindolee <= 10000
 0 <= Mindolep <= 10000
 0 <= Minostc <= 10000
 0 <= Minoste <= 10000
 0 <= Minostp <= 10000
 0 <= Minsc <= 10000
 0 <= Minse <= 10000
 0 <= Minsp <= 10000
 0 <= Mipdpc <= 10000
 0 <= Misetacc <= 10000
 0 <= Misetace <= 10000
 0 <= Misetacp <= 10000
 0 <= Mitpc <= 10000
 0 <= Mkc <= 10000
 0 <= Mkdo2lipid4c <= 10000
 0 <= Mkdo2lipid4e <= 10000
 0 <= Mkdo2lipid4Lc <= 10000
 0 <= Mkdo2lipid4p <= 10000
 0 <= Mkdo2lipid4pc <= 10000
 0 <= Mkdo8pc <= 10000
 0 <= Mkdoc <= 10000
 0 <= Mkdolipid4c <= 10000
 0 <= Mke <= 10000
 0 <= Mkp <= 10000
 0 <= Mkphphhlipac <= 10000
 0 <= MlacDc <= 10000
 0 <= MlacDe <= 10000
 0 <= MlacDp <= 10000
 0 <= MlacLc <= 10000
 0 <= MlacLe <= 10000
 0 <= MlacLp <= 10000
 0 <= MLalaDgluc <= 10000
 0 <= MLalaDgluMdapc <= 10000
 0 <= MLalaDgluMdapDalac <= 10000
 0 <= MLalaDgluMdapDalae <= 10000
 0 <= MLalaDgluMdapDalap <= 10000
 0 <= MLalaDgluMdape <= 10000
 0 <= MLalaDgluMdapp <= 10000
 0 <= MLalaLgluc <= 10000
 0 <= MlaldDc <= 10000
 0 <= MlaldLc <= 10000
 0 <= Mlctsc <= 10000
 0 <= Mlctse <= 10000
 0 <= Mlctsp <= 10000
 0 <= MleuLc <= 10000
 0 <= MleuLe <= 10000
 0 <= MleuLp <= 10000
 0 <= Mleutrnac <= 10000
 0 <= MlgtSc <= 10000
 0 <= Mlipac <= 10000
 0 <= Mlipacoldc <= 10000
 0 <= Mlipacolde <= 10000
 0 <= Mlipacoldp <= 10000
 0 <= Mlipae <= 10000
 0 <= Mlipap <= 10000
 0 <= MlipidAc <= 10000
 0 <= MlipidAdsc <= 10000
 0 <= MlipidXc <= 10000
 0 <= Mlppp <= 10000
 0 <= MlysLc <= 10000
 0 <= MlysLe <= 10000
 0 <= MlysLp <= 10000
 0 <= Mlystrnac <= 10000
 0 <= MlyxLc <= 10000
 0 <= MlyxLe <= 10000
 0 <= MlyxLp <= 10000
 0 <= MmalACPc <= 10000
 0 <= Mmalcoac <= 10000
 0 <= MmalDc <= 10000
 0 <= MmalDe <= 10000
 0 <= MmalDp <= 10000
 0 <= MmalLc <= 10000
 0 <= MmalLe <= 10000
 0 <= MmalLp <= 10000
 0 <= Mmalt6pc <= 10000
 0 <= Mmaltc <= 10000
 0 <= Mmalte <= 10000
 0 <= Mmalthpc <= 10000
 0 <= Mmalthxc <= 10000
 0 <= Mmalthxe <= 10000
 0 <= Mmalthxp <= 10000
 0 <= Mmaltp <= 10000
 0 <= Mmaltptc <= 10000
 0 <= Mmaltpte <= 10000
 0 <= Mmaltptp <= 10000
 0 <= Mmalttrc <= 10000
 0 <= Mmalttre <= 10000
 0 <= Mmalttrp <= 10000
 0 <= Mmaltttrc <= 10000
 0 <= Mmaltttre <= 10000
 0 <= Mmaltttrp <= 10000
 0 <= Mman1pc <= 10000
 0 <= Mman6pc <= 10000
 0 <= Mman6pe <= 10000
 0 <= Mman6pglycc <= 10000
 0 <= Mman6pp <= 10000
 0 <= Mmanac <= 10000
 0 <= Mmanc <= 10000
 0 <= Mmane <= 10000
 0 <= Mmanglyce <= 10000
 0 <= Mmanglycp <= 10000
 0 <= Mmanp <= 10000
 0 <= Mmelibc <= 10000
 0 <= Mmelibe <= 10000
 0 <= Mmelibp <= 10000
 0 <= Mmercppyrc <= 10000
 0 <= MmetDc <= 10000
 0 <= MmetDe <= 10000
 0 <= MmetDp <= 10000
 0 <= Mmethfc <= 10000
 0 <= MmetLc <= 10000
 0 <= MmetLe <= 10000
 0 <= MmetLp <= 10000
 0 <= MmetsoxRLc <= 10000
 0 <= MmetsoxRLe <= 10000
 0 <= MmetsoxRLp <= 10000
 0 <= MmetsoxSLc <= 10000
 0 <= MmetsoxSLe <= 10000
 0 <= MmetsoxSLp <= 10000
 0 <= Mmettrnac <= 10000
 0 <= Mmg2c <= 10000
 0 <= Mmg2e <= 10000
 0 <= Mmg2p <= 10000
 0 <= Mmi1pDc <= 10000
 0 <= Mmicitc <= 10000
 0 <= Mminohpe <= 10000
 0 <= Mminohpp <= 10000
 0 <= Mmlthfc <= 10000
 0 <= MmmcoaRc <= 10000
 0 <= MmmcoaSc <= 10000
 0 <= Mmmetc <= 10000
 0 <= Mmmete <= 10000
 0 <= Mmmetp <= 10000
 0 <= Mmn2c <= 10000
 0 <= Mmn2e <= 10000
 0 <= Mmn2p <= 10000
 0 <= Mmnl1pc <= 10000
 0 <= Mmnle <= 10000
 0 <= Mmnlp <= 10000
 0 <= Mmobdc <= 10000
 0 <= Mmobde <= 10000
 0 <= Mmobdp <= 10000
 0 <= Mmql8c <= 10000
 0 <= Mmqn8c <= 10000
 0 <= Mmso3c <= 10000
 0 <= Mmso3e <= 10000
 0 <= Mmso3p <= 10000
 0 <= Mmthgxlc <= 10000
 0 <= Mmurein3p3pp <= 10000
 0 <= Mmurein3px3pp <= 10000
 0 <= Mmurein3px4pp <= 10000
 0 <= Mmurein4p3pp <= 10000
 0 <= Mmurein4p4pp <= 10000
 0 <= Mmurein4px4p4pp <= 10000
 0 <= Mmurein4px4pp <= 10000
 0 <= Mmurein4px4px4pp <= 10000
 0 <= Mmurein5p3pp <= 10000
 0 <= Mmurein5p4pp <= 10000
 0 <= Mmurein5p5p5pp <= 10000
 0 <= Mmurein5p5pp <= 10000
 0 <= Mmurein5px3pp <= 10000
 0 <= Mmurein5px4pp <= 10000
 0 <= Mmurein5px4px4pp <= 10000
 0 <= MmyrsACPc <= 10000
 0 <= MN1aspmdc <= 10000
 0 <= Mn2oc <= 10000
 0 <= Mn2oe <= 10000
 0 <= Mn2op <= 10000
 0 <= Mn8aspmdc <= 10000
 0 <= Mna1c <= 10000
 0 <= Mna1e <= 10000
 0 <= Mna1p <= 10000
 0 <= Mnacc <= 10000
 0 <= Mnace <= 10000
 0 <= Mnacp <= 10000
 0 <= Mnadc <= 10000
 0 <= Mnadhc <= 10000
 0 <= Mnadpc <= 10000
 0 <= Mnadphc <= 10000
 0 <= Mncamc <= 10000
 0 <= Mnh4c <= 10000
 0 <= Mnh4e <= 10000
 0 <= Mnh4p <= 10000
 0 <= Mni2c <= 10000
 0 <= Mni2e <= 10000
 0 <= Mni2p <= 10000
 0 <= Mnicrntc <= 10000
 0 <= Mnmnc <= 10000
 0 <= Mnmne <= 10000
 0 <= Mnmnp <= 10000
 0 <= MNmtrpc <= 10000
 0 <= Mno2c <= 10000
 0 <= Mno2e <= 10000
 0 <= Mno2p <= 10000
 0 <= Mno3c <= 10000
 0 <= Mno3e <= 10000
 0 <= Mno3p <= 10000
 0 <= Mnoc <= 10000
 0 <= Mnoe <= 10000
 0 <= Mnop <= 10000
 0 <= Mo16a2undp <= 10000
 0 <= Mo16a3undp <= 10000
 0 <= Mo16a4colipae <= 10000
 0 <= Mo16a4colipap <= 10000
 0 <= Mo16a4undp <= 10000
 0 <= Mo16aundc <= 10000
 0 <= Mo16aundp <= 10000
 0 <= Mo2c <= 10000
 0 <= Mo2e <= 10000
 0 <= Mo2p <= 10000
 0 <= Mo2sc <= 10000
 0 <= Mo2se <= 10000
 0 <= Mo2sp <= 10000
 0 <= Moaac <= 10000
 0 <= Moc2coac <= 10000
 0 <= MocACPc <= 10000
 0 <= Moccoac <= 10000
 0 <= MocdcaACPc <= 10000
 0 <= Mocdcac <= 10000
 0 <= Mocdcae <= 10000
 0 <= Mocdcap <= 10000
 0 <= Mocdceac <= 10000
 0 <= Mocdceae <= 10000
 0 <= Mocdceap <= 10000
 0 <= Moctac <= 10000
 0 <= Moctae <= 10000
 0 <= Moctap <= 10000
 0 <= Moctdpc <= 10000
 0 <= MocteACPc <= 10000
 0 <= Mod2coac <= 10000
 0 <= Modecoac <= 10000
 0 <= Mohpbc <= 10000
 0 <= Mop4enc <= 10000
 0 <= Mornc <= 10000
 0 <= Morne <= 10000
 0 <= Mornp <= 10000
 0 <= Morot5pc <= 10000
 0 <= Morotc <= 10000
 0 <= Morote <= 10000
 0 <= Morotp <= 10000
 0 <= Moxamc <= 10000
 0 <= Moxurc <= 10000
 0 <= Mpa120c <= 10000
 0 <= Mpa120p <= 10000
 0 <= Mpa140c <= 10000
 0 <= Mpa140p <= 10000
 0 <= Mpa141c <= 10000
 0 <= Mpa141p <= 10000
 0 <= Mpa160c <= 10000
 0 <= Mpa160p <= 10000
 0 <= Mpa161c <= 10000
 0 <= Mpa161p <= 10000
 0 <= Mpa180c <= 10000
 0 <= Mpa180p <= 10000
 0 <= Mpa181c <= 10000
 0 <= Mpa181p <= 10000
 0 <= Mpacaldc <= 10000
 0 <= Mpacalde <= 10000
 0 <= Mpacaldp <= 10000
 0 <= Mpacc <= 10000
 0 <= MpalmACPc <= 10000
 0 <= Mpan4pc <= 10000
 0 <= MpantRc <= 10000
 0 <= Mpapc <= 10000
 0 <= Mpapsc <= 10000
 0 <= Mpdx5pc <= 10000
 0 <= Mpe120c <= 10000
 0 <= Mpe120p <= 10000
 0 <= Mpe140c <= 10000
 0 <= Mpe140p <= 10000
 0 <= Mpe141c <= 10000
 0 <= Mpe141p <= 10000
 0 <= Mpe160c <= 10000
 0 <= Mpe160p <= 10000
 0 <= Mpe161c <= 10000
 0 <= Mpe161p <= 10000
 0 <= Mpe180c <= 10000
 0 <= Mpe180p <= 10000
 0 <= Mpe181c <= 10000
 0 <= Mpe181p <= 10000
 0 <= Mpeamne <= 10000
 0 <= Mpeamnp <= 10000
 0 <= Mpepc <= 10000
 0 <= Mpg120c <= 10000
 0 <= Mpg120p <= 10000
 0 <= Mpg140c <= 10000
 0 <= Mpg140p <= 10000
 0 <= Mpg141c <= 10000
 0 <= Mpg141p <= 10000
 0 <= Mpg160c <= 10000
 0 <= Mpg160p <= 10000
 0 <= Mpg161c <= 10000
 0 <= Mpg161p <= 10000
 0 <= Mpg180c <= 10000
 0 <= Mpg180p <= 10000
 0 <= Mpg181c <= 10000
 0 <= Mpg181p <= 10000
 0 <= Mpgp120c <= 10000
 0 <= Mpgp120p <= 10000
 0 <= Mpgp140c <= 10000
 0 <= Mpgp140p <= 10000
 0 <= Mpgp141c <= 10000
 0 <= Mpgp141p <= 10000
 0 <= Mpgp160c <= 10000
 0 <= Mpgp160p <= 10000
 0 <= Mpgp161c <= 10000
 0 <= Mpgp161p <= 10000
 0 <= Mpgp180c <= 10000
 0 <= Mpgp180p <= 10000
 0 <= Mpgp181c <= 10000
 0 <= Mpgp181p <= 10000
 0 <= Mphaccoac <= 10000
 0 <= MpheLc <= 10000
 0 <= MpheLe <= 10000
 0 <= MpheLp <= 10000
 0 <= Mphemec <= 10000
 0 <= Mphemee <= 10000
 0 <= Mphemep <= 10000
 0 <= Mphetrnac <= 10000
 0 <= Mphhlipac <= 10000
 0 <= Mphomc <= 10000
 0 <= Mphphhlipac <= 10000
 0 <= Mphpyrc <= 10000
 0 <= Mphthrc <= 10000
 0 <= Mpic <= 10000
 0 <= Mpie <= 10000
 0 <= Mpip <= 10000
 0 <= Mpmcoac <= 10000
 0 <= Mpmtcoac <= 10000
 0 <= MpntoRc <= 10000
 0 <= MpntoRe <= 10000
 0 <= MpntoRp <= 10000
 0 <= Mppac <= 10000
 0 <= Mppae <= 10000
 0 <= Mppalc <= 10000
 0 <= Mppale <= 10000
 0 <= Mppalp <= 10000
 0 <= Mppap <= 10000
 0 <= Mppapc <= 10000
 0 <= Mppbngc <= 10000
 0 <= Mppcoac <= 10000
 0 <= Mppgppc <= 10000
 0 <= Mpphnc <= 10000
 0 <= Mppic <= 10000
 0 <= Mppp9c <= 10000
 0 <= Mpppg9c <= 10000
 0 <= Mpppic <= 10000
 0 <= Mpppnc <= 10000
 0 <= Mpppne <= 10000
 0 <= Mpppnp <= 10000
 0 <= Mppte <= 10000
 0 <= Mpptp <= 10000
 0 <= Mpramc <= 10000
 0 <= Mpranc <= 10000
 0 <= Mprbampc <= 10000
 0 <= Mprbatpc <= 10000
 0 <= Mprfpc <= 10000
 0 <= Mprlpc <= 10000
 0 <= Mproglyc <= 10000
 0 <= Mproglye <= 10000
 0 <= Mproglyp <= 10000
 0 <= MproLc <= 10000
 0 <= MproLe <= 10000
 0 <= MproLp <= 10000
 0 <= Mprotrnac <= 10000
 0 <= Mprppc <= 10000
 0 <= Mps120c <= 10000
 0 <= Mps140c <= 10000
 0 <= Mps141c <= 10000
 0 <= Mps160c <= 10000
 0 <= Mps161c <= 10000
 0 <= Mps180c <= 10000
 0 <= Mps181c <= 10000
 0 <= Mpsclysc <= 10000
 0 <= Mpsclyse <= 10000
 0 <= Mpsclysp <= 10000
 0 <= MpserLc <= 10000
 0 <= MpserLe <= 10000
 0 <= MpserLp <= 10000
 0 <= Mptrcc <= 10000
 0 <= Mptrce <= 10000
 0 <= Mptrcp <= 10000
 0 <= Mpyam5pc <= 10000
 0 <= Mpydamc <= 10000
 0 <= Mpydx5pc <= 10000
 0 <= Mpydxc <= 10000
 0 <= Mpydxnc <= 10000
 0 <= Mpyrc <= 10000
 0 <= Mpyre <= 10000
 0 <= Mpyrp <= 10000
 0 <= Mq8c <= 10000
 0 <= Mq8h2c <= 10000
 0 <= Mqulnc <= 10000
 0 <= Mr15bpc <= 10000
 0 <= Mr1pc <= 10000
 0 <= Mr5pc <= 10000
 0 <= Mr5pe <= 10000
 0 <= Mr5pp <= 10000
 0 <= Mragundc <= 10000
 0 <= Mrbflvrdc <= 10000
 0 <= MrblLc <= 10000
 0 <= Mrdmbzic <= 10000
 0 <= Mrhcysc <= 10000
 0 <= MribDc <= 10000
 0 <= MribDe <= 10000
 0 <= MribDp <= 10000
 0 <= Mribflvc <= 10000
 0 <= Mrml1pc <= 10000
 0 <= Mrmlc <= 10000
 0 <= Mrmnc <= 10000
 0 <= Mrmne <= 10000
 0 <= Mrmnp <= 10000
 0 <= Mru5pDc <= 10000
 0 <= Mru5pLc <= 10000
 0 <= Ms7pc <= 10000
 0 <= Msarcsc <= 10000
 0 <= Msbt6pc <= 10000
 0 <= MsbtDe <= 10000
 0 <= MsbtDp <= 10000
 0 <= Msbzcoac <= 10000
 0 <= Msc <= 10000
 0 <= Msclc <= 10000
 0 <= Msectrnac <= 10000
 0 <= Mselnc <= 10000
 0 <= Mselnpc <= 10000
 0 <= Mserampc <= 10000
 0 <= MserDc <= 10000
 0 <= MserDe <= 10000
 0 <= MserDp <= 10000
 0 <= MserLc <= 10000
 0 <= MserLe <= 10000
 0 <= MserLp <= 10000
 0 <= Msertrnac <= 10000
 0 <= Msertrnasecoooc <= 10000
 0 <= MSfglutthc <= 10000
 0 <= Mshemec <= 10000
 0 <= Mskm5pc <= 10000
 0 <= Mskmc <= 10000
 0 <= Mskme <= 10000
 0 <= Mskmp <= 10000
 0 <= Msl26dac <= 10000
 0 <= Msl2a6oc <= 10000
 0 <= Mso2c <= 10000
 0 <= Mso2e <= 10000
 0 <= Mso2p <= 10000
 0 <= Mso3c <= 10000
 0 <= Mso3e <= 10000
 0 <= Mso3p <= 10000
 0 <= Mso4c <= 10000
 0 <= Mso4e <= 10000
 0 <= Mso4p <= 10000
 0 <= Mspmdc <= 10000
 0 <= Mspmde <= 10000
 0 <= Mspmdp <= 10000
 0 <= Mssaltppc <= 10000
 0 <= Mstcoac <= 10000
 0 <= Msuc6pc <= 10000
 0 <= Msucargc <= 10000
 0 <= Msucbzc <= 10000
 0 <= Msuccc <= 10000
 0 <= Msucce <= 10000
 0 <= Msuccoac <= 10000
 0 <= Msuccp <= 10000
 0 <= Msucgluc <= 10000
 0 <= Msucgsac <= 10000
 0 <= Msuchmsc <= 10000
 0 <= Msucornc <= 10000
 0 <= Msucre <= 10000
 0 <= Msucrp <= 10000
 0 <= Msucsalc <= 10000
 0 <= Msulfacc <= 10000
 0 <= Msulface <= 10000
 0 <= Msulfacp <= 10000
 0 <= Mt3c11vaceACPc <= 10000
 0 <= Mt3c5ddeceACPc <= 10000
 0 <= Mt3c7mrseACPc <= 10000
 0 <= Mt3c9palmeACPc <= 10000
 0 <= Mtag6pDc <= 10000
 0 <= MtagdpDc <= 10000
 0 <= Mtagurc <= 10000
 0 <= MtartrLc <= 10000
 0 <= MtartrLe <= 10000
 0 <= MtartrLp <= 10000
 0 <= Mtaurc <= 10000
 0 <= Mtaure <= 10000
 0 <= Mtaurp <= 10000
 0 <= Mtcyntc <= 10000
 0 <= Mtcynte <= 10000
 0 <= Mtcyntp <= 10000
 0 <= Mtd2coac <= 10000
 0 <= Mtdcoac <= 10000
 0 <= Mtddec2eACPc <= 10000
 0 <= MtdeACPc <= 10000
 0 <= Mtdec2eACPc <= 10000
 0 <= Mtdecoac <= 10000
 0 <= Mthdpc <= 10000
 0 <= Mthex2eACPc <= 10000
 0 <= Mthfc <= 10000
 0 <= Mthmc <= 10000
 0 <= Mthme <= 10000
 0 <= Mthmmpc <= 10000
 0 <= Mthmp <= 10000
 0 <= Mthmppc <= 10000
 0 <= MthrLc <= 10000
 0 <= MthrLe <= 10000
 0 <= MthrLp <= 10000
 0 <= Mthrpe <= 10000
 0 <= Mthrpp <= 10000
 0 <= Mthrtrnac <= 10000
 0 <= Mthymc <= 10000
 0 <= Mthymdc <= 10000
 0 <= Mthymde <= 10000
 0 <= Mthymdp <= 10000
 0 <= Mthyme <= 10000
 0 <= Mthymp <= 10000
 0 <= Mtmac <= 10000
 0 <= Mtmae <= 10000
 0 <= Mtmaoc <= 10000
 0 <= Mtmaoe <= 10000
 0 <= Mtmaop <= 10000
 0 <= Mtmap <= 10000
 0 <= Mtmrs2eACPc <= 10000
 0 <= Mtoct2eACPc <= 10000
 0 <= Mtoctd2eACPc <= 10000
 0 <= Mtpalm2eACPc <= 10000
 0 <= Mtrdoxc <= 10000
 0 <= Mtrdrdc <= 10000
 0 <= Mtre6pc <= 10000
 0 <= Mtrec <= 10000
 0 <= Mtree <= 10000
 0 <= Mtrep <= 10000
 0 <= Mtrnaalac <= 10000
 0 <= Mtrnaargc <= 10000
 0 <= Mtrnaasnc <= 10000
 0 <= Mtrnaaspc <= 10000
 0 <= Mtrnacysc <= 10000
 0 <= Mtrnaglnc <= 10000
 0 <= Mtrnagluc <= 10000
 0 <= Mtrnaglyc <= 10000
 0 <= Mtrnahisc <= 10000
 0 <= Mtrnailec <= 10000
 0 <= Mtrnaleuc <= 10000
 0 <= Mtrnalysc <= 10000
 0 <= Mtrnametc <= 10000
 0 <= Mtrnaphec <= 10000
 0 <= Mtrnaproc <= 10000
 0 <= Mtrnasecysc <= 10000
 0 <= Mtrnaserc <= 10000
 0 <= Mtrnathrc <= 10000
 0 <= Mtrnatrpc <= 10000
 0 <= Mtrnatyrc <= 10000
 0 <= Mtrnavalc <= 10000
 0 <= MtrpLc <= 10000
 0 <= MtrpLe <= 10000
 0 <= MtrpLp <= 10000
 0 <= Mtrptrnac <= 10000
 0 <= Mtsulc <= 10000
 0 <= Mtsule <= 10000
 0 <= Mtsulp <= 10000
 0 <= Mttdcac <= 10000
 0 <= Mttdcae <= 10000
 0 <= Mttdcap <= 10000
 0 <= Mttdceac <= 10000
 0 <= Mttdceae <= 10000
 0 <= Mttdceap <= 10000
 0 <= Mtungsc <= 10000
 0 <= Mtungse <= 10000
 0 <= Mtungsp <= 10000
 0 <= Mtyme <= 10000
 0 <= Mtymp <= 10000
 0 <= MtyrLc <= 10000
 0 <= MtyrLe <= 10000
 0 <= MtyrLp <= 10000
 0 <= Mtyrpe <= 10000
 0 <= Mtyrpp <= 10000
 0 <= Mtyrtrnac <= 10000
 0 <= Mu23gac <= 10000
 0 <= Mu3agac <= 10000
 0 <= Mu3hgac <= 10000
 0 <= Muaagmdac <= 10000
 0 <= Muaccgc <= 10000
 0 <= Muacgamc <= 10000
 0 <= Muacgame <= 10000
 0 <= Muacgamp <= 10000
 0 <= Muacmamc <= 10000
 0 <= Muacmamuc <= 10000
 0 <= Muagmdac <= 10000
 0 <= Muamac <= 10000
 0 <= Muamagc <= 10000
 0 <= Muamrc <= 10000
 0 <= Mudcpdpc <= 10000
 0 <= Mudcpdpp <= 10000
 0 <= Mudcppc <= 10000
 0 <= Mudcppp <= 10000
 0 <= Mudpacgale <= 10000
 0 <= Mudpacgalp <= 10000
 0 <= Mudpc <= 10000
 0 <= Mudpgalc <= 10000
 0 <= Mudpgale <= 10000
 0 <= Mudpgalfurc <= 10000
 0 <= Mudpgalp <= 10000
 0 <= Mudpgc <= 10000
 0 <= Mudpge <= 10000
 0 <= Mudpglcurc <= 10000
 0 <= Mudpglcure <= 10000
 0 <= Mudpglcurp <= 10000
 0 <= Mudpgp <= 10000
 0 <= MudpLa4fnc <= 10000
 0 <= MudpLa4nc <= 10000
 0 <= MudpLa4oc <= 10000
 0 <= Mugmdac <= 10000
 0 <= Mugmdc <= 10000
 0 <= MuLa4fnc <= 10000
 0 <= MuLa4nc <= 10000
 0 <= MuLa4np <= 10000
 0 <= Mum4pc <= 10000
 0 <= Mumpc <= 10000
 0 <= Mumpe <= 10000
 0 <= Mumpp <= 10000
 0 <= Munagac <= 10000
 0 <= Munagamuc <= 10000
 0 <= Munagamufc <= 10000
 0 <= Munagamufp <= 10000
 0 <= Muppg3c <= 10000
 0 <= Murac <= 10000
 0 <= Murae <= 10000
 0 <= Murap <= 10000
 0 <= Muratec <= 10000
 0 <= Murdglycc <= 10000
 0 <= Mureac <= 10000
 0 <= Mureae <= 10000
 0 <= Mureap <= 10000
 0 <= Muric <= 10000
 0 <= Murie <= 10000
 0 <= Murip <= 10000
 0 <= Mutpc <= 10000
 0 <= MvalLc <= 10000
 0 <= MvalLe <= 10000
 0 <= MvalLp <= 10000
 0 <= Mvaltrnac <= 10000
 0 <= Mxanc <= 10000
 0 <= Mxane <= 10000
 0 <= Mxanp <= 10000
 0 <= Mxdpc <= 10000
 0 <= Mxmpc <= 10000
 0 <= Mxmpe <= 10000
 0 <= Mxmpp <= 10000
 0 <= Mxtpc <= 10000
 0 <= Mxtsnc <= 10000
 0 <= Mxtsne <= 10000
 0 <= Mxtsnp <= 10000
 0 <= Mxu5pDc <= 10000
 0 <= Mxu5pLc <= 10000
 0 <= MxylDc <= 10000
 0 <= MxylDe <= 10000
 0 <= MxylDp <= 10000
 0 <= MxyluDc <= 10000
 0 <= MxyluLc <= 10000
 0 <= MxyluLe <= 10000
 0 <= MxyluLp <= 10000
 0 <= Mzn2c <= 10000
 0 <= Mzn2e <= 10000
 0 <= Mzn2p <= 10000

End
