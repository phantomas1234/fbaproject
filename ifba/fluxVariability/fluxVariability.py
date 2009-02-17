#!/usr/bin/env python
# encoding: utf-8
"""
fluxVariability.py

Created by Nikolaus Sonnenschein on 2008-02-12.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import re
import copy
import time
from ifba.GlpkWrap import metabolism, util, randomMedia, fluxdist, glpk

def varibilityQ(lp, reaction):
    """Returns a tuple of the lower and upper limit of the reactions's flux
    variability."""
    pass

def checkVariability(lp, reactions):
    """Returns a dictionary"""
    pass

class Variability(object):
    """docstring for Variability"""
    def __init__(self, lp, noise=1.):
        super(Variability, self).__init__()
        self.lp = lp
        self.noise = noise
        # self.lp.smcp.presolve = glpk.GLP_ON
        self._fixObjVal()
        self.activeR = self.lp.fba().getActiveReactions()
    
    def _fixObjVal(self):
        """Fixes the bounds of the initial objective function to the objective
        value determined by FBA."""
        objVal = self.lp.getObjVal()
        print objVal
        objFunc = 'R("R_Ec_biomass_iAF1260_core_59p81M")'
        print objFunc
        self.lp.modifyColumnBounds({objFunc : (objVal*self.noise, objVal)})
    
    def _lowerLimit(self):
        """Returns the lower flux limit of reaction."""
        self.lp.setOptFlag('min')
        self.lp.simplex()
        self.lp.undo()
        return self.lp.getObjVal()
    
    def _upperLimit(self):
        """Returns the upper flux limit of reaction."""
        self.lp.setOptFlag('max')
        self.lp.simplex()
        self.lp.undo()
        return self.lp.getObjVal()
    
    def variabilityGenerator(self):
        """Returns a generator that iterates over all active reactions in the
        initial flux distribution. Every iteration produces a tuple like this:
        (analysedReactionID, lowerLimit, upperLimit)."""
        for reaction in self.activeR:
            self.lp.setReactionObjective(reaction)
            yield (reaction, (self._lowerLimit(), self._upperLimit()))

class FluxCoupling(object):
    """A Class implementing the FCF (Flux Coupling Finder) Framework described
    in Burgard et al. Flux coupling analysis of genome-scale metabolic network
    reconstructions. Genome Res (2004) vol. 14 (2) pp. 301-312"""# TODO Its not finished yet
    def __init__(self, lp):
        self.lp = lp
        self.lp = self.lp
    
    def computeFluxRatio(self, v_i, v_j):
        """Computes the the ration of v1/v2."""
        fishyFlag = 0
        self.lp.setReactionObjective(v_i)
        self.lp.modifyColumnBounds({v_j : (1., 1.)})
        self.lp.setOptFlag('min')
        try:
            self.lp.simplex()
            minRatio = self.lp.getObjVal()
        except Exception, msg:
            print msg
            print "solving exact ..."
            try:
                self.lp.exact()
                minRatio = self.lp.getObjVal()
            except Exception, msg:
                print msg
                minRatio = self.lp.getObjVal()
                fishyFlag = 1
        self.lp.setOptFlag('max')
        try:
            self.lp.simplex()
            maxRatio = self.lp.getObjVal()
        except Exception, msg:
            print msg
            print "solving exact ..."
            try:
                self.lp.exact()
                maxRatio = self.lp.getObjVal()
            except Exception, msg:
                print msg
                maxRatio = self.lp.getObjVal()
                fishyFlag = 1
        self.lp.initialize()
        return (minRatio, maxRatio, fishyFlag)
    
    def optimalRatioGenerator(self, reacList=None):
        """Returns a generator object the computes in every iteration the
        optimal minimal and maximal ratios of v_i and v_j."""
        if reacList:
            pass
        else:
            reacList = self.lp.getReactions()
        for v_i in reacList:
            for v_j in self.reacList:
                if v_i == v_j:
                    break
                (minRatio, maxRatio, fishyFlag) = self.computeFluxRatio(v_i, v_j)
                yield (v_i, v_j, minRatio, maxRatio, fishyFlag)

def openAllDoors(lp, defaulBound=1000.):
    transp = lp.getTransporters()
    bounds = [(-defaulBound, defaulBound) for i in transp]
    boundsDict = dict(zip(transp, bounds))
    lp.modifyColumnBounds(boundsDict)
    return copy.copy(lp)

def fluxVariabilityMain():
    path = 'iAF1260templateFluxCoupling.lp'
    lp = metabolism.Metabolism(util.ImportCplex(path))
    lp.modifyColumnBounds({'R("MglcDb_Transp")':(0, 20), 'R("Mo2b_Transp")':(0, 18.5)})
    lp.simplex()
    print lp.getObjVal()
    varGenerator = Variability(lp).variabilityGenerator()
    vResult = []
    for var in varGenerator:
        print var
        vResult.append(var)
    print vResult
    open('VariabilityAnalysis.txt', 'w').write(repr(vResult))


def fluxCouplingMain():
    # path = 'iAF1260templateFluxCoupling.lp'
    # path = '../models/iAF1260template2.lp'
    path = '../models/Ecoli_Core_template.lp'
    lp = metabolism.Metabolism(util.ImportCplex(path, terminal="ON"))
    # reacs = [r for r in lp.getReactions() if re.match('.*EX.*', r)]
    # print reacs
    # lp.modifyColumnBounds({'R("R_ATPM")':(0,9999999)})
    print lp.getObjectiveFunction()
    openAllDoors(lp, 999999)
    print lp.getObjectiveFunction()
    # fluxDist = lp.fba()
    # print lp.getObjVal()
    # print fluxDist
    # reactions = set(lp.getReactions()).intersection(set(fluxDist.getActiveReactions()))
    lp2 = copy.copy(lp)
    logFile = open('fluxCoupling.log', 'w')
    logFile.write('')
    logFile.close()
    # lp2.toggleVerbosity()
    allreadyCoupled = list()
    # reactions = ('R("R_ENO")','R("R_PGK")', 'R("R_PGK_Rev")', 'R("R_PGK_Rev")', 'R("R_TPI")','R("R_TPI_Rev")')
    reactions = lp2.getReactions()
    for i in reactions:
    # for i in ('R("R_Ec_biomass_iAF1260_core_59p81M")',):
        for j in reactions:
            if i == j:
                continue
            optRatio = FluxCoupling(lp2).computeFluxRatio(i, j)
            optRatioStr = tuple([repr(elem) for elem in list(optRatio)])
            print str((i, j, optRatio))
            logFile = open('fluxCoupling.log', 'a')
            logFile.write('\t'.join((i,j) + optRatioStr) + "\n")
            logFile.close()

if __name__ == '__main__':
    
    # fluxVariabilityMain()
    fluxCouplingMain()




# old stuff

# def fluxCouplingSingleReaction():
#     # path = 'iAF1260templateFluxCoupling.lp'
#     path = '../models/iAF1260template2.lp'
#     lp = metabolism.Metabolism(util.ImportCplex(path, terminal="ON"))
#     lp.modifyColumnBounds({'R("R_ATPM")':(0,99999999)})
#     openAllDoors(lp,99999999)
#     print lp
#     # fac = 100000
#     # lp.modifyColumnBounds({'R("MglcDb_Transp")':(0, 20*fac), 'R("Mo2b_Transp")':(0, 18*fac)})
#     lp.toggleVerbosity()
#     print lp
#     lp.simplex()
#     fluxDist =  lp.fba()
#     print fluxDist.tsv()
#     reactions = fluxDist.getActiveReactions()
#     lp2 = copy.copy(lp)
#     lp2.toggleVerbosity()
#     logFile = open('fluxCoupling.log', 'w')
#     logFile.write('')
#     logFile.close()
#     
#     # optRatio = FluxCoupling(lp2).computeFluxRatio('R("R_LPADSS")', 'R("R_EDTXS1")')
#     # print str(('R("R_LPADSS")', 'R("R_EDTXS1")', optRatio))
#     optRatio = FluxCoupling(lp2).computeFluxRatio('R("R_MTHFC")', 'R("R_GART")')
#     print optRatio
#     optRatio = FluxCoupling(lp2).computeFluxRatio('R("R_GART")', 'R("R_MTHFC")')
#     print optRatio
# 
#     optRatio = FluxCoupling(lp2).computeFluxRatio('R("R_MTHFC")', 'R("R_ACKr")')
#     print optRatio
#     optRatio = FluxCoupling(lp2).computeFluxRatio('R("R_ACKr")', 'R("R_MTHFC")')
#     print optRatio
# 
#     optRatio = FluxCoupling(lp2).computeFluxRatio('R("R_MTHFC")', 'R("R_ACKr_Rev")')
#     print optRatio
#     optRatio = FluxCoupling(lp2).computeFluxRatio('R("R_ACKr_Rev")', 'R("R_MTHFC")')
#     print optRatio
