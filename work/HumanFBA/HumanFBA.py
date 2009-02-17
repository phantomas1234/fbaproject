#!/usr/bin/env python
# encoding: utf-8
"""
HumanFBA.py

Created by Nikolaus Sonnenschein on 2008-10-06.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from ifba.GlpkWrap import glpk, util, fluxdist, metabolism
import Bio.Geo
import numpy


        











def humanFBA():
    """docstring for main"""
    lp = metabolism.Metabolism(util.ImportCplex('validHumanFBAmodel.lp'))
    print dir(lp)
    print lp
    lp.setReactionObjectiveMinimizeRest('R("R_Obj")')
    # lp.deleteReactions(['R("R_BILGLCURte")', 'R("R_GARFT")'])
    f = lp.fba()
    fluxdist = f.getActiveFluxDist()
    print fluxdist
    print len(fluxdist)
    

def joinSamplesFromExperiment(expDict):
    """docstring for joinSamplesFromExperiment"""
    pass
    
    
    
if __name__ == '__main__':
    humanFBA()


    threshold1 = .8
    threshold2 = .8



    surrogatData = dict([(r, random.uniform(0., 1.)) for r in lp.getReactions()])
    print surrogatData
    for key, value in surrogatData.items():
        if value < .8:
            surrogatData[key] = 1. - value
        else:
            surrogatData.pop(key)
            print surrogatData

    # lp.setReactionObjectiveMinimizeRest('R("R_Ec_biomass_iAF1260_core_59p81M")')
    # lp.modifyColumnBounds({'R("MglcDb_Transp")':(0, 20), 'R("R_ATPM")':(8.39, 8.39), 'R("Mo2b_Transp")':(0, 18.5)})
    # lp2 = copy.copy(lp)
    # fd = lp2.fba()
    # growth = fd['R("R_Ec_biomass_iAF1260_core_59p81M")']
    # print growth
    # lp3 = copy.copy(lp2)
    # lp3.modifyColumnBounds({'R("R_Ec_biomass_iAF1260_core_59p81M")':(growth*threshold1, growth)})
    # 
    # lp3.setObjectiveFunction(surrogatData)
    # 
    # fd = lp3.fba()
    # growth = fd['R("R_Ec_biomass_iAF1260_core_59p81M")']
    # print growth
    # print lp3.getObjectiveFunction()
