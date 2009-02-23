#!/usr/bin/env python
# encoding: utf-8
"""
RndFluxesForMCanalysis.py

Created by Nikolaus Sonnenschein on 2009-01-23.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

from ifba.GlpkWrap import metabolism, randomMedia, util
import copy


def main(path2template, resultsPath, runs):
    import gzip
    tmpLP = metabolism.Metabolism(util.ImportCplex(path2template))
    tmpLP.setReactionObjectiveMinimizeRest('R("R_Ec_biomass_iAF1260_core_59p81M")')
    include = ('R("Mhb_Transp")',
             'R("Mna1b_Transp")',
             'R("Mkb_Transp")',
             'R("Mca2b_Transp")',
             'R("Mcu2b_Transp")',
             'R("Mmg2b_Transp")',
             'R("Mzn2b_Transp")',
             'R("Mmobdb_Transp")',
             'R("Mfe2b_Transp")',
             'R("Mfe3b_Transp")',
             'R("Mcobalt2b_Transp")',
             'R("Mmn2b_Transp")',
             'R("Mclb_Transp")')
    lp = randomMedia.Almaas(copy.copy(tmpLP), alwaysInc=include)
    print lp.lp
    f = lp.generateFluxdist()
    for item in range(runs):
        f = lp.generateFluxdist()
        stringDump = randomMedia.dict2tsv(lp.currDict) + "\n" + f.tsv()
        print stringDump
        path = resultsPath + "iAf1260_fluxDist_" + str(item) + ".tsv.gz"
        lp.lp.initialize()
        gzip.open(path, 'w').write(stringDump)

if __name__ == '__main__':
	main("../models/iAF1260template2.lp", 
	"/Users/niko/arbeit/Publishing/MetabolicControlPaper2/data/RandomFluxDistributions/", 500)

