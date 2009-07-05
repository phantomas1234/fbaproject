#!/usr/bin/env python
# encoding: utf-8
"""
newApproach.py

Data:
Set1: 5 affected environmental conditions
Set2: 35 no effect environmental conditions

Algorithm:
1. Remove irrelevant reactions:
Keep only
Complement[(Not Blocked Reactions in Set1), (Not blocked Reactions in Set2)]
2.


Created by Nikolaus Sonnenschein on 2009-06-12.
Copyright (c) 2009 . All rights reserved.
"""

import sys
import os
import copy
import pickle
from ifba.glpki import glpki
from ifba.GlpkWrap import util, metabolism
from ifba.blockedReactions import blockedReactions as bc


biom = 'R("R_Ec_biomass_iAF1260_core_59p81M")'

def loadCarbonSources(path):
    """docstring for loadCarbonSources"""
    try:
        file = open(path, 'r')
        content = file.read().split('\n')
        file.close()
        return content
    except:
        pass

def prepareModel():
    """docstring for prepareModel"""
    path = '../../ifba/models/iAF1260template_minimalMed_noCarb.lp' # no carb    
    lp = metabolism.Metabolism(util.ImportCplex(path))
    lp.smcp.presolve = glpki.GLP_ON
    return lp

def determineFreeReactions(lp, setCarbonSources, freeReacsSet=None):
    if not freeReacsSet:
        allReactions = lp.getReactions()
        freeReacsSet = set(allReactions)
    for carbSource in setCarbonSources:
        print 'length of freeReactionSet', len(freeReacsSet)
        print carbSource
        lp.modifyColumnBounds({carbSource:(0,20)})
        print 'normal FBA: ', lp.fba()[biom]
        blockedReactions = set(bc.analyseBlockedReactions(lp, freeReacsSet))
        lp.initialize()
        freeReacsSet.difference_update(blockedReactions)
        print 5*'\n'
    return freeReacsSet

if __name__ == '__main__':
    set1 = loadCarbonSources('targetCarbonSources.txt')
    set2 = loadCarbonSources('viableCarbonSources.txt')
    lp = prepareModel()
    # freeReacsSet1 = determineFreeReactions(lp, set1)
    # pickle.dump(freeReacsSet1, open('set1freeReactions.pickle', 'w'))
    freeReacsSet1 = pickle.load(open('set1freeReactions.pickle'))
    print len(freeReacsSet1)
    freeReacsSet2 = determineFreeReactions(lp, set2, freeReacsSet1)
    pickle.dump(freeReacsSet2, open('set2freeReactions.pickle', 'w'))    


