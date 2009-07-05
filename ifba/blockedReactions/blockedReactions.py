#!/usr/bin/env python
# encoding: utf-8
"""
blockedReactions.py

Created by Nikolaus Sonnenschein on 2008-03-13.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from ifba.GlpkWrap import util, metabolism, randomMedia, knockouts, glpk
import copy

def openAllDoors(lp, defaulBound=1000.):
    transp = lp.getTransporters()
    bounds = [(-defaulBound, defaulBound) for i in transp]
    boundsDict = dict(zip(transp, bounds))
    lp.modifyColumnBounds(boundsDict)
    return copy.copy(lp)

def blockedQ(lp, reaction):
    """A predicate function that returns true if a reactions is blocked."""
    returnValue = None
    print reaction
    lp.setReactionObjective(reaction)
    lp.getObjectiveFunction()
    try:
        lp.simplex()
        objVal = lp.getObjVal()
        print objVal
        if -1e-6 < objVal < 1e-6:
            lp.setOptFlag('MiN')
            lp.simplex()
            objVal = lp.getObjVal()
            print objVal
            if -1e-6 < objVal < 1e-6:
                returnValue = True
            else:
                returnValue = False
        else:
            returnValue = False
    except:
        returnValue = True
    lp.initialize()
    return returnValue

def analyseBlockedReactions(lp, reactions=None):
    """Returns a list of blocked reactions in the lp model."""
    blockedReactions = list()
    if not reactions:
        reactions = lp.getColumnIDs()
    for reaction in reactions:
        if blockedQ(lp, reaction):
            blockedReactions.append(reaction)
            print "Hit ya! A blocked one dude!", reaction
        else:
            print "A normal reaction!", reaction
    return blockedReactions

if __name__ == '__main__':

    # path = sys.argv[1]
    path = '../models/iAF1260template.lp'
    lp = metabolism.Metabolism(util.ImportCplex(path))
    openLP = openAllDoors(lp)
    openLP.smcp.presolve = glpk.GLP_ON
    print openLP
    util.WriteCplex(openLP)
    blockedReactions = analyseBlockedReactions(openLP)
    print blockedReactions
    open(sys.argv[2], 'w').write("\n".join(blockedReactions))

