#!/usr/bin/env python
# encoding: utf-8
"""
moritzTest.py

Created by Nikolaus Sonnenschein on 2008-09-24.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from ifba.GlpkWrap import util, metabolism, randomMedia, knockouts, fluxdist
import copy

def readMinMed(path):
    f = open(path, 'r')
    med = eval(f.read())
    f.close()
    return dict(med)

def sampleFBA():
    """docstring for main"""
    lp = metabolism.Metabolism(util.ImportCplex('../models/iAF1260template2.lp'))
    lp.setReactionObjectiveMinimizeRest('R("R_Ec_biomass_iAF1260_core_59p81M")')
    # print lp
    # lp.simplex()
    # print "The current objective value: ", lp.getObjVal()
    lp.modifyColumnBounds({'R("MglcDb_Transp")':(0, 20), 'R("R_ATPM")':(8.39, 8.39), 'R("Mo2b_Transp")':(0, 18.5)})
    lp.simplex()
    lp.simplex()
    print "The current objective value after some media modifications: ", lp.getObjVal()
    data = fluxdist.FluxDist(lp)
    print data.getActiveReactions()
    print data.getActiveFluxDist()

    lp.modifyColumnBounds({'R("R_PGK")':(0,0), 'R("R_PGK_Rev")':(0,0), 'R("R_PGK_Rev")':(0,0)})
    lp.simplex()
    print "The current objective value after a Reaction Knockout: ", lp.getObjVal()
    lp.undo()
    # lp.initialize()
    lp.simplex()
    print "The current objective value after undoing the Knockout: ", lp.getObjVal()
    
    print dir(lp)

    # print data
    # print lp.getObjVal()
    # print lp.getReactions()
    # # print lp.getTransporters()
    # transp = lp.getTransporters()
    # print lp.translateColumnNames(transp)
    # lp.modifyColumnBounds({'R("R_PGK")':(0, 20), 'R("R_ENO")':(0, 20)})
    # print lp.history

if __name__ == '__main__':
	sampleFBA()

