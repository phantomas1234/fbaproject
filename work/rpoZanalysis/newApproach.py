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
from ifba.glpki import glpki
from ifba.GlpkWrap import util, metabolism
from ifba.blockedReactions import blockedReactions as bc
import copy

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

def main():
    pass


if __name__ == '__main__':
    path = '../../ifba/models/iAF1260template_minimalMed_noCarb.lp' # no carb    
    lp = metabolism.Metabolism(util.ImportCplex(path))
    lp.smcp.presolve = glpki.GLP_ON
    print lp
    set1 = loadCarbonSources('targetCarbonSources.txt')
    set2 = loadCarbonSources('viableCarbonSources.txt')
    print set1
    carbSource = set1[1]
    print carbSource
    lp.modifyColumnBounds({carbSource:(0,10)})
    lp.simplex()
    print 'asdf'
    print lp.translateColumnNames(['R("R_PGK")'])
    print 'asfd'
    print lp.translateColumnNames(['R("R_Ec_biomass_iAF1260_core_59p81M")'])
    print lp.getObjVal()
    blockedResults = dict()
    for reac in set1:
        print reac
        blockedReactions = bc.analyseBlockedReactions(lp)
        blockedResults[reac] = blockedReactions
    print blockedResults
    open('debug.txt', 'w').write(str(blockedResults))
    
    # util.WriteCplex(openLP)
    # blockedReactions = analyseBlockedReactions(openLP)
    # print blockedReactions
    # open(sys.argv[2], 'w').write("\n".join(blockedReactions))

    
    main()

