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
from ifba.GlpkWrap import metabolism, util, randomMedia, fluxdist, glpk, exceptions
from termcolor import colored

class Variability(object):
    """docstring for Variability"""
    def __init__(self, lp, threshold=10**-16):
        super(Variability, self).__init__()
        self.lp = lp
        self.theshold = threshold
    
    def _lowerLimit(self):
        """Returns the lower flux limit of reaction."""
        self.lp.setOptFlag('min')
        try:
            self.lp.simplex()
        except:
            print "Running with presolver"
            self.lp.smcp.presolve = glpk.GLP_ON
            self.lp.simplex()
            self.lp.smcp.presolve = glpk.GLP_OFF
        self.lp.undo()
        return self.lp.getObjVal()
    
    def _upperLimit(self):
        """Returns the upper flux limit of reaction."""
        self.lp.setOptFlag('max')
        try:
            self.lp.simplex()
        except:
            print "Running with presolver"
            self.lp.smcp.presolve = glpk.GLP_ON
            self.lp.simplex()
            self.lp.smcp.presolve = glpk.GLP_OFF
        self.lp.undo()
        return self.lp.getObjVal()
    
    def variabilityQ(self, reaction):
        """Returns a tuple of the lower and upper limit of the reactions's flux
        variability."""
        self.lp.setReactionObjective(reaction)
        var = (self._lowerLimit(), self._upperLimit())
        self.lp.initialize()
        return var
    
    def variabilityGenerator(self):
        """Returns a generator that iterates over all active reactions in the
        initial flux distribution. Every iteration produces a tuple like this:
        (analysedReactionID, lowerLimit, upperLimit)."""
        for reaction in self.lp.getReactions():
            print reaction
            revFlag = 0
            if re.search(".*_Rev", reaction):
                revFlag = 1
                self.lp.modifyColumnBounds({reaction.replace('_Rev',''):(0., 0.)})
            if not revFlag:
                try:
                    self.lp.modifyColumnBounds({reaction.split('")')[0] + "_Rev" + '")':(0., 0.)})
                except:
                    pass
            self.lp.setReactionObjective(reaction)
            lowLimit = self._lowerLimit()
            upperLimit = self._upperLimit()
            self.lp.initialize()
            print reaction, lowLimit, upperLimit
            yield (reaction, (lowLimit, upperLimit))
    
    def getBlocked(self):
        """Returns a list of blocked reactions"""
        gen = self.variabilityGenerator()
        blocked = list()
        for var in gen:
            if var[1][1] < self.theshold:
                print 'Blocked!!!', var[0]
                blocked.append(var[0])
        return blocked
    


class FluxCoupling(object):
    """A Class implementing the FCF (Flux Coupling Finder) Framework described
    in Burgard et al. Flux coupling analysis of genome-scale metabolic network
    reconstructions. Genome Res (2004) vol. 14 (2) pp. 301-312"""# TODO Its not finished yet
    def __init__(self, lp):
        self.lp = lp
    
    def computeFluxRatio(self, v_i, v_j):
        """Computes the the ration of v1/v2."""
        fishyFlag = 0
        self.lp.setReactionObjective(v_i)
        # print lp.fba().getActiveFluxDist()
        self.lp.modifyColumnBounds({v_j : (1., 1.)})
        self.lp.setOptFlag('min')
        try:
            self.lp.simplex()
            minRatio = self.lp.getObjVal()
            dualDict = lp.getShadowPriceDict()
            (rxn1minShadow, rxn2minShadow) = (dualDict[v_i], dualDict[v_j])
        except Exception, msg:
            minRatio = self.lp.getObjVal()
            fishyFlag = 1
        self.lp.setOptFlag('max')
        try:
            self.lp.simplex()
            maxRatio = self.lp.getObjVal()
            dualDict = lp.getShadowPriceDict()
            (rxn1maxShadow, rxn2maxShadow) = (dualDict[v_i], dualDict[v_j])
        except Exception, msg:
            maxRatio = self.lp.getObjVal()
            fishyFlag = 1
        self.lp.initialize()
        return (minRatio, maxRatio, rxn1maxShadow, rxn1minShadow, rxn2maxShadow, rxn2minShadow, fishyFlag)

    def fluxCouplingFinder(self, reacs=None):
        if not reacs:
            reacs=self.lp.getReactions()
        reacs = [r for r in reacs if not re.search('R_EX.*', r)]
        print len(reacs)
        bounds = self.lp.getColumnBounds()
        uncoupled = list()
        blocked = list()
        directionallyCoupled = list()
        fullCoupled= list()
        partiallyCoupled = list()
        allreadyCoupled = list()
        for i, rxn1 in enumerate(reacs):
            if i in allreadyCoupled:
                continue
            if i in blocked:
                continue
            print rxn1
            self.lp.setReactionObjective(rxn1)
            self.lp.setOptFlag('max')
            self.lp.simplex()
            maxObjVal = self.lp.getObjVal()
            if maxObjVal <= 0:
                print colored(rxn1, 'red'), " is blocked"
                blocked.append(i)
                self.lp.initialize
                continue
            self.lp.initialize()
            print maxObjVal
            print "still ", len(reacs) - i - len(allreadyCoupled) - len(blocked), " to check"
            print "directionallyCoupled: ", directionallyCoupled
            print "fullCoupled: ", fullCoupled
            print "partiallyCoupled: ", partiallyCoupled
            for j, rxn2 in enumerate(reacs):
                if j <= i:
                    continue
                if j in blocked:
                    continue
                self.lp.setReactionObjective(rxn2)
                self.lp.simplex()
                if self.lp.getObjVal() < 1.:
                    print colored(rxn2, 'red'), " is blocked"
                    blocked.append(j)
                    self.lp.initialize
                    continue
                self.lp.initialize
                (r_min, r_max, rxn1maxShadow, rxn1minShadow, rxn2maxShadow, rxn2minShadow, fishy) = self.computeFluxRatio(rxn1, rxn2)
                info = ' '.join(("------> minRatio: ", str(chop(r_min)), "maxRatio: ", str(chop(r_max)), "rxn1maxShadow: ", str(chop(rxn1maxShadow)), "rxn1minShadow: ", str(chop(rxn1minShadow)), "rxn2maxShadow: ", str(chop(rxn2maxShadow)), "rxn2minShadow: ", str(chop(rxn2minShadow))))
                if fishy:
                    print "fishy: ", rxn1, " ", rxn2, " ", fishy
                    self.lp.smcp.presolve = glpk.GLP_ON
                    (r_min, r_max, rxn1maxShadow, rxn1minShadow, rxn2maxShadow, rxn2minShadow, fishy) = self.computeFluxRatio(rxn1, rxn2)
                    info = ' '.join(("------> minRatio: ", str(chop(r_min)), "maxRatio: ", str(chop(r_max)), "rxn1maxShadow: ", str(chop(rxn1maxShadow)), "rxn1minShadow: ", str(chop(rxn1minShadow)), "rxn2maxShadow: ", str(chop(rxn2maxShadow)), "rxn2minShadow: ", str(chop(rxn2minShadow))))
                    self.lp.smcp.presolve = glpk.GLP_OFF
                    print "fishy: ", rxn1, " ", rxn2, " ", fishy
                    if fishy:
                        continue
                    # sys.exit(-1)
                    continue
                if chop(rxn2maxShadow) <= 0. and chop(rxn2minShadow) == 0.:
                    print "uncoupled: ", rxn1, rxn2
                    print info
                    uncoupled.append((rxn1, rxn2))
                    continue
                elif chop(rxn2maxShadow) > 0. and chop(rxn2minShadow) == 0.:
                    directionallyCoupled.append((rxn1, rxn2))
                    print colored("directionallyCoupled: ", 'red'), rxn1, rxn2
                    print info
                    continue
                elif chop(rxn2maxShadow) == 0. and chop(rxn2minShadow) > 0.:
                    directionallyCoupled.append((rxn1,rxn2))
                    print colored("reverse directionallyCoupled: ", 'red'), rxn1, rxn2
                    print info
                    continue
                elif chop(rxn2maxShadow) == chop(rxn2minShadow):
                    allreadyCoupled.append(j)
                    fullCoupled.append((rxn1, rxn2))
                    print colored("fullycoupled: ", 'red'), rxn1, rxn2
                    print info
                    continue
                elif chop(rxn2maxShadow) != chop(rxn2minShadow):
                    allreadyCoupled.append(j)
                    partiallyCoupled.append((rxn1, rxn2))
                    print colored("partiallycoupled: ", 'red'), rxn1, rxn2
                    print info
                else:
                    raise Exception, ' '.join('Coupling for', rxn1, 'and', rxn2, 'could not be determined!')

def chop(val, cutoff=1e-10):
    if (val <= cutoff) and (val >= -cutoff):
        return 0.
    else:
        return val

def openAllDoors(lp, defaulBound=1000.):
    transp = lp.getTransporters()
    bounds = [(-defaulBound, defaulBound) for i in transp]
    boundsDict = dict(zip(transp, bounds))
    lp.modifyColumnBounds(boundsDict)
    return copy.copy(lp)

def fluxVariabilityMain(lp):
    varGenerator = Variability(lp).variabilityGenerator()
    vResult = []
    for var in varGenerator:
        print var
        vResult.append(var)
    # open('VariabilityAnalysis.txt', 'w').write(repr(vResult))
    return vResult

def freeBounds(lp):
    lp.modifyColumnBounds(dict([(r, ('-inf', 1000000.))for r in lp.getTransporters()]))
    lp.modifyColumnBounds(dict([(r, (0., 1000000.))for r in lp.getReactions()]))
    # util.WriteCplex(lp, 'iAF1260_template_for_FluxCouplingAnalysis.lp')
    return copy.copy(lp)
    
def scaleColumnBounds(lp, scale=1.):
    rxns = lp.getColumnIDs()
    columnBounds = lp.getColumnBounds()
    newBounds = dict()
    for r in rxns:
        tmpBounds = columnBounds[r]
        newBounds[r] = (tmpBounds[0] * scale, tmpBounds[1] * scale)
    lp.modifyColumnBounds(newBounds)
    return copy.copy(lp)

def fluxCouplingMain2(lp):
    blocked = Variability(lp).getBlocked()
    lp.initialize()
    print len(lp.getMetabolites())
    lp.deleteReactions(['R("R_Ec_biomass_iAF1260_core_59p81M")'])
    (sub, prod) = lp.getSubstratesAndProducts('R("R_Ec_biomass_iAF1260_core_59p81M")')
    list(sub).extend(list(prod))
    bioMets = sub + prod
    print bioMets
    print len(lp.getMetabolites())
    lp.addMetaboliteDrains(bioMets)
    print len(lp.getMetabolites())
    lp.eraseHistory()
    blocked = open('blocked.tsv', 'r').read().split('\n')
    reactions = list(set(lp.getReactions()) - set(blocked))
    varObject = Variability(lp)
    for r in reactions[0:20]:
        var = varObject.variabilityQ(r)
        if var[1] <= 1.:
            print r, var
    lp.initialize()
    FluxCoupling(lp).fluxCouplingFinder(reacs=reactions[0:1000])

if __name__ == '__main__':
    path = 'iAF1260templateFluxCoupling.lp'
    lp = freeBounds(metabolism.Metabolism(util.ImportCplex(path, terminal="OFF")))
    lp = copy.copy(lp)
    fluxCouplingMain2(lp)

