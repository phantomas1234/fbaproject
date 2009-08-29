#!/usr/bin/env python
# encoding: utf-8
"""
CheckEssential.py

Created by Nikolaus Sonnenschein on 2009-08-20.
Copyright (c) 2009 . All rights reserved.
"""

import sys
import os
import copy
from ifba.GlpkWrap import util, metabolism

def essentialQ(model, reac, threshold=0.005):
    model.deleteReactions((reac,))
    model.simplex()
    objVal = model.getObjVal()
    print objVal
    if objVal <= threshold:
        returnValue = True
    else:
        returnValue = False
    model.initialize()
    return returnValue

def checkEssential(model_path, reacs=None):
    lp = metabolism.Metabolism(util.ImportCplex(model_path))
    if reacs:
        result = analyseBlockedReactions(lp, reacs)
    else:
        result = analyseBlockedReactions(lp)
    return result

if __name__ == '__main__':
    model_path = '/Users/niko/arbeit/Data/SBMLmodels/HomoSapiens/humanFBA_ATPObjective_MartinsMedium.lp'
    lp = metabolism.Metabolism(util.ImportCplex(model_path))
    activeReacs1 = set(lp.fba().getActiveReactions())
    print lp.getObjVal()
    lp.deleteReactions(('R("R_GHMT2r")','R("R_GHMT2r_Rev")', 'R("R_GHMT2rm")', 'R("R_GHMT2rm_Rev")'))
    activeReacs2 = set(lp.fba().getActiveReactions())
    print lp.getObjVal()
    print activeReacs1.difference(activeReacs2)
    
    