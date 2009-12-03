#!/usr/bin/env python
# encoding: utf-8
"""
AlternativeTRN.py

Created by Nikolaus Sonnenschein on 2009-11-12.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import copy
import re
from ifba.GlpkWrap import util, metabolism, randomMedia
from ifba.GlpkWrap.fluxdist import FluxDist

def filterTransporters(dictionary, transpRegEx=".*_Transp"):
    tmpDict = dict()
    for i in dictionary.keys():
        if re.search(transpRegEx, i):
            tmpDict[i] = dictionary[i]
    return tmpDict
        

if __name__ == '__main__':
    oxygenUptake = 8
    carbonUptake = 4
    path = '../../ifba/models/iAF1260templateMinMax.lp'
    transp = ['R("Mo2b_Transp")', 'R("MglcDb_Transp")', 'R("MlacDb_Transp")','R("Mmaltb_Transp")']
    lp = metabolism.Metabolism(util.ImportCplex(path))
    # Allow secretion for all transporters
    # allTransp = lp.getTransporters()
    # colBounds = lp.getColumnBounds()
    # tmpDict = dict()
    # for t in allTransp:
    #     if not colBounds[t][1] > 0:
    #         tmpDict[t] = (-999999, 0)
    # lp.modifyColumnBounds(tmpDict)
    # lp.eraseHistory()
    # Remove Glucose and Oxygen from the medium
    lp.modifyColumnBounds({'R("Mo2b_Transp")':(0,0), 'R("MglcDb_Transp")':(0,0)})
    lp.eraseHistory()
    # Test glucose conditions
    lp.modifyColumnBounds({'R("Mo2b_Transp")':(0,oxygenUptake), 'R("MglcDb_Transp")':(0,carbonUptake)})
    fluxdist = lp.fba()
    print fluxdist['R("R_Ec_biomass_iAF1260_core_59p81M")']
    print [(i, fluxdist[i]) for i in transp]
    print filterTransporters(dict(fluxdist.getActiveFluxDist()))
    lp.initialize()
    # Test lactose conditions
    # print [(elem, lp.getColumnBounds()[elem]) for elem in transp]
    lp.modifyColumnBounds({'R("Mo2b_Transp")':(0,oxygenUptake), 'R("MlacDb_Transp")':(0,carbonUptake)})
    fluxdist = lp.fba()
    print fluxdist['R("R_Ec_biomass_iAF1260_core_59p81M")']
    print [(i, fluxdist[i]) for i in transp]
    print filterTransporters(dict(fluxdist.getActiveFluxDist()))
    lp.initialize()
    # Test lactose and glucose conditions
    lp.modifyColumnBounds({'R("Mo2b_Transp")':(0,oxygenUptake), 'R("MlacDb_Transp")':(0,carbonUptake), 'R("MglcDb_Transp")':(0,carbonUptake)})
    fluxdist = lp.fba()
    print fluxdist['R("R_Ec_biomass_iAF1260_core_59p81M")']
    print [(i, fluxdist[i]) for i in transp]
    print filterTransporters(dict(fluxdist.getActiveFluxDist()))
    lp.initialize()
    # Test maltose conditions
    lp.modifyColumnBounds({'R("Mo2b_Transp")':(0,oxygenUptake), 'R("Mmaltb_Transp")':(0,carbonUptake)})
    fluxdist = lp.fba()
    print fluxdist['R("R_Ec_biomass_iAF1260_core_59p81M")']
    print [(i, fluxdist[i]) for i in transp]
    print filterTransporters(dict(fluxdist.getActiveFluxDist()))
    lp.initialize()
    