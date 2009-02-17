Maximize
Z: R("R_Biomass_Ecoli_core_N__w_GAM_")

Subject To
M13dpgc: R("R_GAPD") - R("R_GAPD_Rev") + R("R_PGK") - R("R_PGK_Rev")
  = 0.

M2pgc: -R("R_ENO") + R("R_ENO_Rev") - R("R_PGM") + R("R_PGM_Rev")  =
 0.

M3pgc: -1.496 R("R_Biomass_Ecoli_core_N__w_GAM_") - R("R_PGK") +
 R("R_PGK_Rev") + R("R_PGM") - R("R_PGM_Rev")  = 0.

M6pgcc: -R("R_GND") + R("R_PGL")  = 0.

M6pglc: R("R_G6PDH2r") - R("R_G6PDH2r_Rev") - R("R_PGL")  = 0.

Macb: R("Macb_Transp") + R("R_EX_ac_e_") - R("R_EX_ac_e__Rev")  = 0.

Macc: -R("R_ACKr") + R("R_ACKr_Rev") + R("R_ACt2r") -
 R("R_ACt2r_Rev")  = 0.

Maccoac: -R("R_ADHEr") + R("R_ADHEr_Rev") - 3.7478
 R("R_Biomass_Ecoli_core_N__w_GAM_") - R("R_CS") - R("R_MALS") + R("R_PDH") + R("R_PFL")
 - R("R_PTAr") + R("R_PTAr_Rev")  = 0.

Mace: 2 R("Mace_Transp") - R("R_ACt2r") + R("R_ACt2r_Rev") -
 R("R_EX_ac_e_") + R("R_EX_ac_e__Rev")  = 0.

Mactpc: R("R_ACKr") - R("R_ACKr_Rev") + R("R_PTAr") -
 R("R_PTAr_Rev")  = 0.

Madpc: R("R_ACKr") - R("R_ACKr_Rev") + 2 R("R_ADK1") - 2
 R("R_ADK1_Rev") + R("R_ATPM") - R("R_ATPS4r") + R("R_ATPS4r_Rev") + 55.703
 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_PFK") + R("R_PGK") -
 R("R_PGK_Rev") + R("R_PPCK") - R("R_PYK") + R("R_SUCOAS") - R("R_SUCOAS_Rev")  =
 0.

Makgb: R("Makgb_Transp") + R("R_EX_akg_e_") - R("R_EX_akg_e__Rev") 
 = 0.

Makgc: -R("R_AKGDH") + R("R_AKGt2r") - R("R_AKGt2r_Rev") - 1.0789
 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_ICDHyr") - R("R_ICDHyr_Rev") 
 = 0.

Makge: 2 R("Makge_Transp") - R("R_AKGt2r") + R("R_AKGt2r_Rev") -
 R("R_EX_akg_e_") + R("R_EX_akg_e__Rev")  = 0.

Mampc: -R("R_ADK1") + R("R_ADK1_Rev") + R("R_PPS")  = 0.

Matpc: -R("R_ACKr") + R("R_ACKr_Rev") - R("R_ADK1") +
 R("R_ADK1_Rev") - R("R_ATPM") + R("R_ATPS4r") - R("R_ATPS4r_Rev") - 55.703
 R("R_Biomass_Ecoli_core_N__w_GAM_") - R("R_PFK") - R("R_PGK") + R("R_PGK_Rev")
 - R("R_PPCK") - R("R_PPS") + R("R_PYK") - R("R_SUCOAS") +
 R("R_SUCOAS_Rev")  = 0.

Mcitc: -R("R_ACONT") + R("R_ACONT_Rev") + R("R_CS")  = 0.

Mco2b: R("Mco2b_Transp") + R("R_EX_co2_e_") - R("R_EX_co2_e__Rev") 
 = 0.

Mco2c: R("R_AKGDH") + R("R_CO2t") - R("R_CO2t_Rev") + R("R_GND") +
 R("R_ICDHyr") - R("R_ICDHyr_Rev") + R("R_ME1") + R("R_ME2") +
 R("R_PDH") - R("R_PPC") + R("R_PPCK")  = 0.

Mco2e: 2 R("Mco2e_Transp") - R("R_CO2t") + R("R_CO2t_Rev") -
 R("R_EX_co2_e_") + R("R_EX_co2_e__Rev")  = 0.

Mcoac: R("R_ADHEr") - R("R_ADHEr_Rev") - R("R_AKGDH") + 3.7478
 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_CS") + R("R_MALS") - R("R_PDH")
 - R("R_PFL") + R("R_PTAr") - R("R_PTAr_Rev") - R("R_SUCOAS") +
 R("R_SUCOAS_Rev")  = 0.

Mdhapc: R("R_FBA") - R("R_FBA_Rev") - R("R_TPI") + R("R_TPI_Rev")  =
 0.

Me4pc: -0.361 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_TALA") -
 R("R_TALA_Rev") - R("R_TKT2") + R("R_TKT2_Rev")  = 0.

Metohb: R("Metohb_Transp") + R("R_EX_etoh_e_") -
 R("R_EX_etoh_e__Rev")  = 0.

Metohc: R("R_ADHEr") - R("R_ADHEr_Rev") + R("R_ETOHt2r") -
 R("R_ETOHt2r_Rev")  = 0.

Metohe: 2 R("Metohe_Transp") - R("R_ETOHt2r") + R("R_ETOHt2r_Rev") -
 R("R_EX_etoh_e_") + R("R_EX_etoh_e__Rev")  = 0.

Mf6pc: -0.0709 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_FBP") -
 R("R_PFK") + R("R_PGI") - R("R_PGI_Rev") + R("R_TALA") - R("R_TALA_Rev")
 + R("R_TKT2") - R("R_TKT2_Rev")  = 0.

Mfadc: R("R_FRD") - R("R_SUCD1i") + R("R_SUCD4") - R("R_SUCD4_Rev") 
 = 0.

Mfadh2c: -R("R_FRD") + R("R_SUCD1i") - R("R_SUCD4") +
 R("R_SUCD4_Rev")  = 0.

Mfdpc: -R("R_FBA") + R("R_FBA_Rev") - R("R_FBP") + R("R_PFK")  = 0.

Mforb: R("Mforb_Transp") + R("R_EX_for_e_") - R("R_EX_for_e__Rev") 
 = 0.

Mforc: R("R_FORt") - R("R_FORt_Rev") + R("R_PFL")  = 0.

Mfore: 2 R("Mfore_Transp") - R("R_EX_for_e_") + R("R_EX_for_e__Rev")
 - R("R_FORt") + R("R_FORt_Rev")  = 0.

Mfumb: R("Mfumb_Transp") + R("R_EX_fum_e_") - R("R_EX_fum_e__Rev") 
 = 0.

Mfumc: -R("R_FRD") - R("R_FUM") + R("R_FUM_Rev") + R("R_FUMt2_2") +
 R("R_SUCD1i")  = 0.

Mfume: 2 R("Mfume_Transp") - R("R_EX_fum_e_") + R("R_EX_fum_e__Rev")
 - R("R_FUMt2_2")  = 0.

Mg3pc: -0.129 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_FBA") -
 R("R_FBA_Rev") - R("R_GAPD") + R("R_GAPD_Rev") - R("R_TALA") +
 R("R_TALA_Rev") + R("R_TKT1") - R("R_TKT1_Rev") + R("R_TKT2") - R("R_TKT2_Rev")
 + R("R_TPI") - R("R_TPI_Rev")  = 0.

Mg6pc: -0.205 R("R_Biomass_Ecoli_core_N__w_GAM_") - R("R_G6PDH2r") +
 R("R_G6PDH2r_Rev") + R("R_GLCpts") - R("R_PGI") + R("R_PGI_Rev")  =
 0.

MglcDb: R("MglcDb_Transp") + R("R_EX_glc_e_") - R("R_EX_glc_e__Rev")
  = 0.

MglcDe: 2 R("MglcDe_Transp") - R("R_EX_glc_e_") +
 R("R_EX_glc_e__Rev") - R("R_GLCpts")  = 0.

Mglxc: R("R_ICL") - R("R_MALS")  = 0.

Mh2ob: R("Mh2ob_Transp") + R("R_EX_h2o_e_") - R("R_EX_h2o_e__Rev") 
 = 0.

Mh2oc: -R("R_ATPM") + R("R_ATPS4r") - R("R_ATPS4r_Rev") - 55.703
 R("R_Biomass_Ecoli_core_N__w_GAM_") - R("R_CS") + R("R_CYTBD") +
 R("R_ENO") - R("R_ENO_Rev") - R("R_FBP") - R("R_FUM") + R("R_FUM_Rev") +
 R("R_H2Ot") - R("R_H2Ot_Rev") - R("R_MALS") - R("R_PGL") - R("R_PPC") -
 R("R_PPS")  = 0.

Mh2oe: 2 R("Mh2oe_Transp") - R("R_EX_h2o_e_") + R("R_EX_h2o_e__Rev")
 - R("R_H2Ot") + R("R_H2Ot_Rev")  = 0.

Mhb: R("Mhb_Transp") + R("R_EX_h_e_") - R("R_EX_h_e__Rev")  = 0.

Mhc: R("R_ACt2r") - R("R_ACt2r_Rev") - 2 R("R_ADHEr") + 2
 R("R_ADHEr_Rev") + R("R_AKGt2r") - R("R_AKGt2r_Rev") + R("R_ATPM") + 3
 R("R_ATPS4r") - 3 R("R_ATPS4r_Rev") + 41.025
 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_CS") - 2 R("R_CYTBD") + R("R_D_LACt2") - R("R_D_LACt2_Rev") +
 R("R_ETOHt2r") - R("R_ETOHt2r_Rev") + 2 R("R_FUMt2_2") +
 R("R_G6PDH2r") - R("R_G6PDH2r_Rev") + R("R_GAPD") - R("R_GAPD_Rev") + R("R_LDH_D")
 - R("R_LDH_D_Rev") + R("R_MALS") + R("R_MDH") - R("R_MDH_Rev") - 3
 R("R_NADH11") + R("R_PFK") + R("R_PGL") + R("R_PPC") + 2 R("R_PPS") -
 R("R_PYK") + R("R_PYRt2r") - R("R_PYRt2r_Rev") + 2 R("R_SUCCt2_2") -
 R("R_SUCCt2b") + 2 R("R_THD2")  = 0.

Mhe: 13 R("Mhe_Transp") - R("R_ACt2r") + R("R_ACt2r_Rev") -
 R("R_AKGt2r") + R("R_AKGt2r_Rev") - 4 R("R_ATPS4r") + 4 R("R_ATPS4r_Rev") + 2
 R("R_CYTBD") - R("R_D_LACt2") + R("R_D_LACt2_Rev") - R("R_ETOHt2r") +
 R("R_ETOHt2r_Rev") - R("R_EX_h_e_") + R("R_EX_h_e__Rev") - 2
 R("R_FUMt2_2") + 2 R("R_NADH11") - R("R_PYRt2r") + R("R_PYRt2r_Rev") - 2
 R("R_SUCCt2_2") + R("R_SUCCt2b") - 2 R("R_THD2")  = 0.

Micitc: R("R_ACONT") - R("R_ACONT_Rev") - R("R_ICDHyr") +
 R("R_ICDHyr_Rev") - R("R_ICL")  = 0.

MlacDb: R("MlacDb_Transp") + R("R_EX_lac_D_e_") -
 R("R_EX_lac_D_e__Rev")  = 0.

MlacDc: R("R_D_LACt2") - R("R_D_LACt2_Rev") - R("R_LDH_D") +
 R("R_LDH_D_Rev")  = 0.

MlacDe: 2 R("MlacDe_Transp") - R("R_D_LACt2") + R("R_D_LACt2_Rev") -
 R("R_EX_lac_D_e_") + R("R_EX_lac_D_e__Rev")  = 0.

MmalLc: R("R_FUM") - R("R_FUM_Rev") + R("R_MALS") - R("R_MDH") +
 R("R_MDH_Rev") - R("R_ME1") - R("R_ME2")  = 0.

Mnadc: 2 R("R_ADHEr") - 2 R("R_ADHEr_Rev") - R("R_AKGDH") - 3.547
 R("R_Biomass_Ecoli_core_N__w_GAM_") - R("R_GAPD") + R("R_GAPD_Rev") -
 R("R_LDH_D") + R("R_LDH_D_Rev") - R("R_MDH") + R("R_MDH_Rev") -
 R("R_ME1") + R("R_NADH11") - R("R_NADTRHD") - R("R_PDH") + R("R_THD2")  = 0.

Mnadhc: -2 R("R_ADHEr") + 2 R("R_ADHEr_Rev") + R("R_AKGDH") + 3.547
 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_GAPD") - R("R_GAPD_Rev") +
 R("R_LDH_D") - R("R_LDH_D_Rev") + R("R_MDH") - R("R_MDH_Rev") +
 R("R_ME1") - R("R_NADH11") + R("R_NADTRHD") + R("R_PDH") - R("R_THD2")  =
 0.

Mnadpc: 18.225 R("R_Biomass_Ecoli_core_N__w_GAM_") - R("R_G6PDH2r")
 + R("R_G6PDH2r_Rev") - R("R_GND") - R("R_ICDHyr") + R("R_ICDHyr_Rev")
 - R("R_ME2") + R("R_NADTRHD") - R("R_THD2")  = 0.

Mnadphc: -18.225 R("R_Biomass_Ecoli_core_N__w_GAM_") +
 R("R_G6PDH2r") - R("R_G6PDH2r_Rev") + R("R_GND") + R("R_ICDHyr") -
 R("R_ICDHyr_Rev") + R("R_ME2") - R("R_NADTRHD") + R("R_THD2")  = 0.

Mo2b: R("Mo2b_Transp") + R("R_EX_o2_e_") - R("R_EX_o2_e__Rev")  = 0.

Mo2c: -0.5 R("R_CYTBD") + R("R_O2t") - R("R_O2t_Rev")  = 0.

Mo2e: 2 R("Mo2e_Transp") - R("R_EX_o2_e_") + R("R_EX_o2_e__Rev") -
 R("R_O2t") + R("R_O2t_Rev")  = 0.

Moaac: -1.7867 R("R_Biomass_Ecoli_core_N__w_GAM_") - R("R_CS") +
 R("R_MDH") - R("R_MDH_Rev") + R("R_PPC") - R("R_PPCK")  = 0.

Mpepc: -0.5191 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_ENO") -
 R("R_ENO_Rev") - R("R_GLCpts") - R("R_PPC") + R("R_PPCK") + R("R_PPS") -
 R("R_PYK")  = 0.

Mpib: R("Mpib_Transp") + R("R_EX_pi_e_") - R("R_EX_pi_e__Rev")  = 0.

Mpic: R("R_ATPM") - R("R_ATPS4r") + R("R_ATPS4r_Rev") + 55.703
 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_FBP") - R("R_GAPD") +
 R("R_GAPD_Rev") - R("R_PIt") + R("R_PIt_Rev") + R("R_PPC") + R("R_PPS") -
 R("R_PTAr") + R("R_PTAr_Rev") + R("R_SUCOAS") - R("R_SUCOAS_Rev")  = 0.

Mpie: 2 R("Mpie_Transp") - R("R_EX_pi_e_") + R("R_EX_pi_e__Rev") +
 R("R_PIt") - R("R_PIt_Rev")  = 0.

Mpyrb: R("Mpyrb_Transp") + R("R_EX_pyr_e_") - R("R_EX_pyr_e__Rev") 
 = 0.

Mpyrc: -2.8328 R("R_Biomass_Ecoli_core_N__w_GAM_") + R("R_GLCpts") +
 R("R_LDH_D") - R("R_LDH_D_Rev") + R("R_ME1") + R("R_ME2") -
 R("R_PDH") - R("R_PFL") - R("R_PPS") + R("R_PYK") + R("R_PYRt2r") -
 R("R_PYRt2r_Rev")  = 0.

Mpyre: 2 R("Mpyre_Transp") - R("R_EX_pyr_e_") + R("R_EX_pyr_e__Rev")
 - R("R_PYRt2r") + R("R_PYRt2r_Rev")  = 0.

Mq8c: R("R_CYTBD") - R("R_NADH11") - R("R_SUCD4") + R("R_SUCD4_Rev")
  = 0.

Mq8h2c: -R("R_CYTBD") + R("R_NADH11") + R("R_SUCD4") -
 R("R_SUCD4_Rev")  = 0.

Mr5pc: -0.8977 R("R_Biomass_Ecoli_core_N__w_GAM_") - R("R_RPI") +
 R("R_RPI_Rev") - R("R_TKT1") + R("R_TKT1_Rev")  = 0.

Mru5pDc: R("R_GND") - R("R_RPE") + R("R_RPE_Rev") + R("R_RPI") -
 R("R_RPI_Rev")  = 0.

Ms7pc: -R("R_TALA") + R("R_TALA_Rev") + R("R_TKT1") -
 R("R_TKT1_Rev")  = 0.

Msuccb: R("Msuccb_Transp") + R("R_EX_succ_e_") -
 R("R_EX_succ_e__Rev")  = 0.

Msuccc: R("R_FRD") + R("R_ICL") + R("R_SUCCt2_2") - R("R_SUCCt2b") -
 R("R_SUCD1i") - R("R_SUCOAS") + R("R_SUCOAS_Rev")  = 0.

Msucce: 3 R("Msucce_Transp") - R("R_EX_succ_e_") +
 R("R_EX_succ_e__Rev") - R("R_SUCCt2_2") + R("R_SUCCt2b")  = 0.

Msuccoac: R("R_AKGDH") + R("R_SUCOAS") - R("R_SUCOAS_Rev")  = 0.

Mxu5pDc: R("R_RPE") - R("R_RPE_Rev") - R("R_TKT1") + R("R_TKT1_Rev")
 - R("R_TKT2") + R("R_TKT2_Rev")  = 0.

Bounds

0 <= R("R_PGM") <= 999999
0 <= R("R_PGM_Rev") <= 999999
0 <= R("R_ENO") <= 999999
0 <= R("R_ENO_Rev") <= 999999
0 <= R("R_EX_ac_e_") <= 999999
0 <= R("R_EX_ac_e__Rev") <= 999999
0 <= R("R_EX_akg_e_") <= 999999
0 <= R("R_EX_akg_e__Rev") <= 999999
0 <= R("R_PGK") <= 999999
0 <= R("R_PGK_Rev") <= 999999
0 <= R("R_ACKr") <= 999999
0 <= R("R_ACKr_Rev") <= 999999
0 <= R("R_ADK1") <= 999999
0 <= R("R_ADK1_Rev") <= 999999
0 <= R("R_ACONT") <= 999999
0 <= R("R_ACONT_Rev") <= 999999
0 <= R("R_EX_co2_e_") <= 999999
0 <= R("R_EX_co2_e__Rev") <= 999999
0 <= R("R_CO2t") <= 999999
0 <= R("R_CO2t_Rev") <= 999999
0 <= R("R_TPI") <= 999999
0 <= R("R_TPI_Rev") <= 999999
0 <= R("R_EX_etoh_e_") <= 999999
0 <= R("R_EX_etoh_e__Rev") <= 999999
0 <= R("R_FBA") <= 999999
0 <= R("R_FBA_Rev") <= 999999
0 <= R("R_EX_for_e_") <= 999999
0 <= R("R_EX_for_e__Rev") <= 999999
0 <= R("R_FORt") <= 999999
0 <= R("R_FORt_Rev") <= 999999
0 <= R("R_EX_fum_e_") <= 999999
0 <= R("R_EX_fum_e__Rev") <= 999999
0 <= R("R_PGI") <= 999999
0 <= R("R_PGI_Rev") <= 999999
0 <= R("R_EX_glc_e_") <= 999999
0 <= R("R_EX_glc_e__Rev") <= 999999
0 <= R("R_FUM") <= 999999
0 <= R("R_FUM_Rev") <= 999999
0 <= R("R_EX_h2o_e_") <= 999999
0 <= R("R_EX_h2o_e__Rev") <= 999999
0 <= R("R_H2Ot") <= 999999
0 <= R("R_H2Ot_Rev") <= 999999
0 <= R("R_EX_h_e_") <= 999999
0 <= R("R_EX_h_e__Rev") <= 999999
0 <= R("R_ACt2r") <= 999999
0 <= R("R_ACt2r_Rev") <= 999999
0 <= R("R_AKGt2r") <= 999999
0 <= R("R_AKGt2r_Rev") <= 999999
0 <= R("R_ETOHt2r") <= 999999
0 <= R("R_ETOHt2r_Rev") <= 999999
0 <= R("R_EX_lac_D_e_") <= 999999
0 <= R("R_EX_lac_D_e__Rev") <= 999999
0 <= R("R_D_LACt2") <= 999999
0 <= R("R_D_LACt2_Rev") <= 999999
0 <= R("R_LDH_D") <= 999999
0 <= R("R_LDH_D_Rev") <= 999999
0 <= R("R_MDH") <= 999999
0 <= R("R_MDH_Rev") <= 999999
0 <= R("R_ADHEr") <= 999999
0 <= R("R_ADHEr_Rev") <= 999999
0 <= R("R_G6PDH2r") <= 999999
0 <= R("R_G6PDH2r_Rev") <= 999999
0 <= R("R_ICDHyr") <= 999999
0 <= R("R_ICDHyr_Rev") <= 999999
0 <= R("R_EX_o2_e_") <= 999999
0 <= R("R_EX_o2_e__Rev") <= 999999
0 <= R("R_O2t") <= 999999
0 <= R("R_O2t_Rev") <= 999999
0 <= R("R_PIt") <= 999999
0 <= R("R_PIt_Rev") <= 999999
0 <= R("R_PTAr") <= 999999
0 <= R("R_PTAr_Rev") <= 999999
0 <= R("R_ATPS4r") <= 999999
0 <= R("R_ATPS4r_Rev") <= 999999
0 <= R("R_GAPD") <= 999999
0 <= R("R_GAPD_Rev") <= 999999
0 <= R("R_EX_pi_e_") <= 999999
0 <= R("R_EX_pi_e__Rev") <= 999999
0 <= R("R_EX_pyr_e_") <= 999999
0 <= R("R_EX_pyr_e__Rev") <= 999999
0 <= R("R_PYRt2r") <= 999999
0 <= R("R_PYRt2r_Rev") <= 999999
0 <= R("R_SUCD4") <= 999999
0 <= R("R_SUCD4_Rev") <= 999999
0 <= R("R_RPI") <= 999999
0 <= R("R_RPI_Rev") <= 999999
0 <= R("R_RPE") <= 999999
0 <= R("R_RPE_Rev") <= 999999
0 <= R("R_TALA") <= 999999
0 <= R("R_TALA_Rev") <= 999999
0 <= R("R_SUCOAS") <= 999999
0 <= R("R_SUCOAS_Rev") <= 999999
0 <= R("R_EX_succ_e_") <= 999999
0 <= R("R_EX_succ_e__Rev") <= 999999
0 <= R("R_TKT2") <= 999999
0 <= R("R_TKT2_Rev") <= 999999
0 <= R("R_TKT1") <= 999999
0 <= R("R_TKT1_Rev") <= 999999
0 <= R("R_PFK") <= 999999
0 <= R("R_FRD") <= 999999
0 <= R("R_PGL") <= 999999
0 <= R("R_ATPM") <= 999999
0 <= R("R_FBP") <= 999999
0 <= R("R_MALS") <= 999999
0 <= R("R_FUMt2_2") <= 999999
0 <= R("R_ICL") <= 999999
0 <= R("R_AKGDH") <= 999999
0 <= R("R_ME1") <= 999999
0 <= R("R_GND") <= 999999
0 <= R("R_ME2") <= 999999
0 <= R("R_THD2") <= 999999
0 <= R("R_NADTRHD") <= 999999
0 <= R("R_PPCK") <= 999999
0 <= R("R_CS") <= 999999
0 <= R("R_GLCpts") <= 999999
0 <= R("R_PPC") <= 999999
0 <= R("R_PYK") <= 999999
0 <= R("R_PFL") <= 999999
0 <= R("R_PPS") <= 999999
0 <= R("R_PDH") <= 999999
0 <= R("R_NADH11") <= 999999
0 <= R("R_CYTBD") <= 999999
0 <= R("R_Biomass_Ecoli_core_N__w_GAM_") <= 999999
0 <= R("R_SUCD1i") <= 999999
0 <= R("R_SUCCt2b") <= 999999
0 <= R("R_SUCCt2_2") <= 999999


\ Some model modifications follow

-999999 <= R("Mace_Transp") <= 0
-999999 <= R("Macb_Transp") <= 0
-999999 <= R("Makge_Transp") <= 0
-999999 <= R("Makgb_Transp") <= 0
-999999 <= R("Mco2e_Transp") <= 0
-999999 <= R("Mco2b_Transp") <= 0
-999999 <= R("Mco2e_Transp") <= 0
-999999 <= R("Metohe_Transp") <= 0
-999999 <= R("Metohb_Transp") <= 0
-999999 <= R("Mfore_Transp") <= 0
-999999 <= R("Mforb_Transp") <= 0
-999999 <= R("Mfore_Transp") <= 0
-999999 <= R("Mfume_Transp") <= 0
-999999 <= R("Mfumb_Transp") <= 0
-999999 <= R("MglcDe_Transp") <= 0
-999999 <= R("MglcDb_Transp") <= 0
-999999 <= R("Mh2oe_Transp") <= 0
-999999 <= R("Mh2ob_Transp") <= 0
-999999 <= R("Mh2oe_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("Mhb_Transp") <= 0
-999999 <= R("Mace_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("Makge_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("Metohe_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("MlacDe_Transp") <= 0
-999999 <= R("MlacDb_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("MlacDe_Transp") <= 0
-999999 <= R("Mo2e_Transp") <= 0
-999999 <= R("Mo2b_Transp") <= 0
-999999 <= R("Mo2e_Transp") <= 0
-999999 <= R("Mpie_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("Mpie_Transp") <= 0
-999999 <= R("Mpib_Transp") <= 0
-999999 <= R("Mpyre_Transp") <= 0
-999999 <= R("Mpyrb_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("Mpyre_Transp") <= 0
-999999 <= R("Msucce_Transp") <= 0
-999999 <= R("Msuccb_Transp") <= 0
-999999 <= R("Mfume_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("MglcDe_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("Msucce_Transp") <= 0
-999999 <= R("Mhe_Transp") <= 0
-999999 <= R("Msucce_Transp") <= 0
