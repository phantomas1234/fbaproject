#!/usr/bin/env python
# encoding: utf-8
"""
metaboliteKOs.py

Created by Nikolaus Sonnenschein on 2009-06-17.
Copyright (c) 2009 . All rights reserved.
"""

import sys
import os
from ifba.glpki import glpki
from ifba.GlpkWrap import util, metabolism, randomMedia
from ifba import general
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
    lp = Almaas(copy.copy(tmpLP), alwaysInc=include)
    med.generateFluxdist(minGrowth=0.5)
    conditions = med.currDict()
    print conditions
    # for item in range(runs):
    #     f = lp.generateFluxdist()
    #     print f
    #     # stringDump = general.util.dict2tsv(lp.currDict) + "\n" + f.tsv()
    #     # general.util.dict2tsv(metKOresult)
    #     # print stringDump
    #     # path = resultsPath + "iAf1260_metaboliteEssentialities_" + str(item) + ".tsv.gz"
    #     # lp.lp.initialize()
    #     # gzip.open(path, 'w').write(stringDump)



if __name__ == '__main__':
    path = '../../ifba/models/iAF1260template2.lp' # no carb
    resultsPath = '/Users/arbeit/Programs/BottleNecksChemicalConstraints/MetaboliteEssentialities'
    main(path, resultsPath, runs)

    lp = metabolism.Metabolism(util.ImportCplex(path))
    dir(randomMedia)
    # print lp.modifyColumnBounds({'R("MglcDb_Transp")':(0,20), 'R("Mo2b_Transp")':(0,20)})
    mets = lp.fba().getActiveMetabolites()
    print 'active metabolties: ', mets
    # metKOresult = dict()
    # print lp
    # for m in mets:
    #     lp.deleteMetabolites([m])
    #     try:
    #         objVal = lp.fba()['R("R_Ec_biomass_iAF1260_core_59p81M")']
    #     except:
    #         print 'hey'
    #         objVal = 0.
    #     print m, objVal
    #     metKOresult[m] = objVal
    #     lp.undo()
    # open('metaboliteKOresult.tsv', 'w').write(general.util.dict2tsv(metKOresult))
