\* Problem: Unknown *\

Maximize
 obj: + Mr5pc

Subject To
 R("R_Biomass_Ecoli_core_N__w_GAM_"): - 0.8977 Mr5pc - 2.8328 Mpyrc
 + 55.703 Mpic - 0.5191 Mpepc - 1.7867 Moaac - 18.225 Mnadphc
 + 18.225 Mnadpc + 3.547 Mnadhc - 3.547 Mnadc + 41.025 Mhc
 - 55.703 Mh2oc - 0.205 Mg6pc - 0.129 Mg3pc - 0.0709 Mf6pc - 0.361 Me4pc
 + 3.7478 Mcoac - 55.703 Matpc - 1.0789 Makgc + 55.703 Madpc
 - 3.7478 Maccoac - 1.496 M3pgc = 0
 R("R_GAPD"): - Mpic + Mnadhc - Mnadc + Mhc - Mg3pc + M13dpgc = 0
 R("R_PGK"): - Matpc + Madpc - M3pgc + M13dpgc = 0
 R("R_ENO"): + Mpepc + Mh2oc - M2pgc = 0
 R("R_PGM"): + M3pgc - M2pgc = 0
 R("R_GND"): + Mru5pDc + Mnadphc - Mnadpc + Mco2c - M6pgcc = 0
 R("R_PGL"): + Mhc - Mh2oc - M6pglc + M6pgcc = 0
 R("R_G6PDH2r"): + Mnadphc - Mnadpc + Mhc - Mg6pc + M6pglc = 0
 R("Macb_Transp"): + Macb = 0
 R("R_EX_ac_e_"): - Mace + Macb = 0
 R("R_ACKr"): - Matpc + Madpc + Mactpc - Macc = 0
 R("R_ACt2r"): - Mhe + Mhc - Mace + Macc = 0
 R("R_ADHEr"): - 2 Mnadhc + 2 Mnadc - 2 Mhc + Metohc + Mcoac - Maccoac
 = 0
 R("R_CS"): - Moaac + Mhc - Mh2oc + Mcoac + Mcitc - Maccoac = 0
 R("R_MALS"): + MmalLc + Mhc - Mh2oc - Mglxc + Mcoac - Maccoac = 0
 R("R_PDH"): - Mpyrc + Mnadhc - Mnadc - Mcoac + Mco2c + Maccoac = 0
 R("R_PFL"): - Mpyrc + Mforc - Mcoac + Maccoac = 0
 R("R_PTAr"): - Mpic + Mcoac + Mactpc - Maccoac = 0
 R("Mace_Transp"): + 2 Mace = 0
 R("R_ADK1"): - Matpc - Mampc + 2 Madpc = 0
 R("R_ATPM"): + Mpic + Mhc - Mh2oc - Matpc + Madpc = 0
 R("R_ATPS4r"): - Mpic - 4 Mhe + 3 Mhc + Mh2oc + Matpc - Madpc = 0
 R("R_PFK"): + Mhc + Mfdpc - Mf6pc - Matpc + Madpc = 0
 R("R_PPCK"): + Mpepc - Moaac + Mco2c - Matpc + Madpc = 0
 R("R_PYK"): + Mpyrc - Mpepc - Mhc + Matpc - Madpc = 0
 R("R_SUCOAS"): + Msuccoac - Msuccc + Mpic - Mcoac - Matpc + Madpc = 0
 R("Makgb_Transp"): + Makgb = 0
 R("R_EX_akg_e_"): - Makge + Makgb = 0
 R("R_AKGDH"): + Msuccoac + Mnadhc - Mnadc - Mcoac + Mco2c - Makgc = 0
 R("R_AKGt2r"): - Mhe + Mhc - Makge + Makgc = 0
 R("R_ICDHyr"): + Mnadphc - Mnadpc - Micitc + Mco2c + Makgc = 0
 R("Makge_Transp"): + 2 Makge = 0
 R("R_PPS"): - Mpyrc + Mpic + Mpepc + 2 Mhc - Mh2oc - Matpc + Mampc = 0
 R("R_ACONT"): + Micitc - Mcitc = 0
 R("Mco2b_Transp"): + Mco2b = 0
 R("R_EX_co2_e_"): - Mco2e + Mco2b = 0
 R("R_CO2t"): - Mco2e + Mco2c = 0
 R("R_ME1"): + Mpyrc + Mnadhc - Mnadc - MmalLc + Mco2c = 0
 R("R_ME2"): + Mpyrc + Mnadphc - Mnadpc - MmalLc + Mco2c = 0
 R("R_PPC"): + Mpic - Mpepc + Moaac + Mhc - Mh2oc - Mco2c = 0
 R("Mco2e_Transp"): + 2 Mco2e = 0
 R("R_FBA"): + Mg3pc - Mfdpc + Mdhapc = 0
 R("R_TPI"): + Mg3pc - Mdhapc = 0
 R("R_TALA"): - Ms7pc - Mg3pc + Mf6pc + Me4pc = 0
 R("R_TKT2"): - Mxu5pDc + Mg3pc + Mf6pc - Me4pc = 0
 R("Metohb_Transp"): + Metohb = 0
 R("R_EX_etoh_e_"): - Metohe + Metohb = 0
 R("R_ETOHt2r"): - Mhe + Mhc - Metohe + Metohc = 0
 R("Metohe_Transp"): + 2 Metohe = 0
 R("R_FBP"): + Mpic - Mh2oc - Mfdpc + Mf6pc = 0
 R("R_PGI"): - Mg6pc + Mf6pc = 0
 R("R_FRD"): + Msuccc - Mfumc - Mfadh2c + Mfadc = 0
 R("R_SUCD1i"): - Msuccc + Mfumc + Mfadh2c - Mfadc = 0
 R("R_SUCD4"): + Mq8h2c - Mq8c - Mfadh2c + Mfadc = 0
 R("Mforb_Transp"): + Mforb = 0
 R("R_EX_for_e_"): - Mfore + Mforb = 0
 R("R_FORt"): - Mfore + Mforc = 0
 R("Mfore_Transp"): + 2 Mfore = 0
 R("Mfumb_Transp"): + Mfumb = 0
 R("R_EX_fum_e_"): - Mfume + Mfumb = 0
 R("R_FUM"): + MmalLc - Mh2oc - Mfumc = 0
 R("R_FUMt2_2"): - 2 Mhe + 2 Mhc - Mfume + Mfumc = 0
 R("Mfume_Transp"): + 2 Mfume = 0
 R("R_TKT1"): - Mxu5pDc + Ms7pc - Mr5pc + Mg3pc = 0
 R("R_GLCpts"): + Mpyrc - Mpepc - MglcDe + Mg6pc = 0
 R("MglcDb_Transp"): + MglcDb = 0
 R("R_EX_glc_e_"): - MglcDe + MglcDb = 0
 R("MglcDe_Transp"): + 2 MglcDe = 0
 R("R_ICL"): + Msuccc - Micitc + Mglxc = 0
 R("Mh2ob_Transp"): + Mh2ob = 0
 R("R_EX_h2o_e_"): - Mh2oe + Mh2ob = 0
 R("R_CYTBD"): - Mq8h2c + Mq8c - 0.5 Mo2c + 2 Mhe - 2 Mhc + Mh2oc = 0
 R("R_H2Ot"): - Mh2oe + Mh2oc = 0
 R("Mh2oe_Transp"): + 2 Mh2oe = 0
 R("Mhb_Transp"): + Mhb = 0
 R("R_EX_h_e_"): - Mhe + Mhb = 0
 R("R_D_LACt2"): - MlacDe + MlacDc - Mhe + Mhc = 0
 R("R_LDH_D"): + Mpyrc + Mnadhc - Mnadc - MlacDc + Mhc = 0
 R("R_MDH"): + Moaac + Mnadhc - Mnadc - MmalLc + Mhc = 0
 R("R_NADH11"): + Mq8h2c - Mq8c - Mnadhc + Mnadc + 2 Mhe - 3 Mhc = 0
 R("R_PYRt2r"): - Mpyre + Mpyrc - Mhe + Mhc = 0
 R("R_SUCCt2_2"): - Msucce + Msuccc - 2 Mhe + 2 Mhc = 0
 R("R_SUCCt2b"): + Msucce - Msuccc + Mhe - Mhc = 0
 R("R_THD2"): + Mnadphc - Mnadpc - Mnadhc + Mnadc - 2 Mhe + 2 Mhc = 0
 R("Mhe_Transp"): + 13 Mhe = 0
 R("MlacDb_Transp"): + MlacDb = 0
 R("R_EX_lac_D_e_"): - MlacDe + MlacDb = 0
 R("MlacDe_Transp"): + 2 MlacDe = 0
 R("R_NADTRHD"): - Mnadphc + Mnadpc + Mnadhc - Mnadc = 0
 R("Mo2b_Transp"): + Mo2b = 0
 R("R_EX_o2_e_"): - Mo2e + Mo2b = 0
 R("R_O2t"): - Mo2e + Mo2c = 0
 R("Mo2e_Transp"): + 2 Mo2e = 0
 R("Mpib_Transp"): + Mpib = 0
 R("R_EX_pi_e_"): - Mpie + Mpib = 0
 R("R_PIt"): + Mpie - Mpic = 0
 R("Mpie_Transp"): + 2 Mpie = 0
 R("Mpyrb_Transp"): + Mpyrb = 0
 R("R_EX_pyr_e_"): - Mpyre + Mpyrb = 0
 R("Mpyre_Transp"): + 2 Mpyre = 0
 R("R_RPI"): + Mru5pDc - Mr5pc = 0
 R("R_RPE"): + Mxu5pDc - Mru5pDc = 0
 R("Msuccb_Transp"): + Msuccb = 0
 R("R_EX_succ_e_"): - Msucce + Msuccb = 0
 R("Msucce_Transp"): + 3 Msucce = 0

Bounds
 0 <= M13dpgc <= 10000
 0 <= M2pgc <= 10000
 0 <= M3pgc <= 10000
 0 <= M6pgcc <= 10000
 0 <= M6pglc <= 10000
 0 <= Macb <= 10000
 0 <= Macc <= 10000
 0 <= Maccoac <= 10000
 0 <= Mace <= 10000
 0 <= Mactpc <= 10000
 0 <= Madpc <= 10000
 0 <= Makgb <= 10000
 0 <= Makgc <= 10000
 0 <= Makge <= 10000
 0 <= Mampc <= 10000
 0 <= Matpc <= 10000
 0 <= Mcitc <= 10000
 0 <= Mco2b <= 10000
 0 <= Mco2c <= 10000
 0 <= Mco2e <= 10000
 0 <= Mcoac <= 10000
 0 <= Mdhapc <= 10000
 0 <= Me4pc <= 10000
 0 <= Metohb <= 10000
 0 <= Metohc <= 10000
 0 <= Metohe <= 10000
 0 <= Mf6pc <= 10000
 0 <= Mfadc <= 10000
 0 <= Mfadh2c <= 10000
 0 <= Mfdpc <= 10000
 0 <= Mforb <= 10000
 0 <= Mforc <= 10000
 0 <= Mfore <= 10000
 0 <= Mfumb <= 10000
 0 <= Mfumc <= 10000
 0 <= Mfume <= 10000
 0 <= Mg3pc <= 10000
 0 <= Mg6pc <= 10000
 0 <= MglcDb <= 10000
 0 <= MglcDe <= 10000
 0 <= Mglxc <= 10000
 0 <= Mh2ob <= 10000
 0 <= Mh2oc <= 10000
 0 <= Mh2oe <= 10000
 0 <= Mhb <= 10000
 0 <= Mhc <= 10000
 0 <= Mhe <= 10000
 0 <= Micitc <= 10000
 0 <= MlacDb <= 10000
 0 <= MlacDc <= 10000
 0 <= MlacDe <= 10000
 0 <= MmalLc <= 10000
 0 <= Mnadc <= 10000
 0 <= Mnadhc <= 10000
 0 <= Mnadpc <= 10000
 0 <= Mnadphc <= 10000
 0 <= Mo2b <= 10000
 0 <= Mo2c <= 10000
 0 <= Mo2e <= 10000
 0 <= Moaac <= 10000
 0 <= Mpepc <= 10000
 0 <= Mpib <= 10000
 0 <= Mpic <= 10000
 0 <= Mpie <= 10000
 0 <= Mpyrb <= 10000
 0 <= Mpyrc <= 10000
 0 <= Mpyre <= 10000
 0 <= Mq8c <= 10000
 0 <= Mq8h2c <= 10000
 0 <= Mr5pc <= 10000
 0 <= Mru5pDc <= 10000
 0 <= Ms7pc <= 10000
 0 <= Msuccb <= 10000
 0 <= Msuccc <= 10000
 0 <= Msucce <= 10000
 0 <= Msuccoac <= 10000
 0 <= Mxu5pDc <= 10000

End
