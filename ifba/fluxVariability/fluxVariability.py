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
    def __init__(self, lp, threshold=1e-10, presolve=True):
        super(Variability, self).__init__()
        self.lp = lp
        self.threshold = threshold
        if presolve:
            self.lp.smcp.presolve = glpk.GLP_ON
    
    def _lowerLimit(self):
        """Returns the lower flux limit of reaction."""
        self.lp.setOptFlag('min')
        try:
            self.lp.simplex()
        except Exception, msg:
            print msg
            # print "Running with presolver"
            # self.lp.smcp.presolve = glpk.GLP_ON
            # self.lp.simplex()
            # self.lp.smcp.presolve = glpk.GLP_OFF
        self.lp.undo()
        return self.lp.getObjVal()
    
    def _upperLimit(self):
        """Returns the upper flux limit of reaction."""
        self.lp.setOptFlag('max')
        try:
            self.lp.simplex()
        except Exception, msg:
            print msg
            # print "Running with presolver"
            # self.lp.smcp.presolve = glpk.GLP_ON
            # self.lp.simplex()
            # self.lp.smcp.presolve = glpk.GLP_OFF
        self.lp.undo()
        return self.lp.getObjVal()
    
    def variabilityQ(self, reaction):
        """Returns a tuple of the lower and upper limit of the reactions's flux
        variability."""
        self.lp.setReactionObjective(reaction)
        var = (self._lowerLimit(), self._upperLimit())
        self.lp.initialize()
        return var
    
    def variabilityGenerator(self, reactions2check=None):
        """Returns a generator that iterates over all active reactions in the
        initial flux distribution. Every iteration produces a tuple like this:
        (analysedReactionID, lowerLimit, upperLimit)."""
        if not reactions2check:
            reactions2check = self.lp.getReactions()
        for reaction in reactions2check:
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
    
    def getBlocked(self, reactions2check=None):
        """Returns a list of blocked reactions"""
        if not reactions2check:
            reactions2check = self.lp.getReactions()
        gen = self.variabilityGenerator(reactions2check=reactions2check)
        blocked = list()
        for var in gen:
            if var[1][0] > -self.threshold and var[1][1] < self.threshold:
                print 'Blocked!!!', var[0]
                blocked.append(var[0])
        return blocked


class FluxCoupling(object):
    """A Class implementing the FCF (Flux Coupling Finder) Framework described
    in Burgard et al. Flux coupling analysis of genome-scale metabolic network
    reconstructions. Genome Res (2004) vol. 14 (2) pp. 301-312"""# TODO Its not finished yet
    def __init__(self, lp):
        self.lp = lp
        print self.lp
    
    def computeFluxRatio(self, v_i, v_j):
        """Computes the the ration of v1/v2."""
        fishyFlag = 0
        self.lp.setReactionObjective(v_i)
        self.lp.modifyColumnBounds({v_j : (1., 1.)})
        # print self.lp.fba().getActiveFluxDist()
        self.lp.setOptFlag('min')
        revFlag = 0
        if re.search(".*_Rev", v_i):
            revFlag = 1
            try:
                self.lp.modifyColumnBounds({v_i.replace('_Rev',''):(0., 0.)})
            except:
                pass
        if not revFlag:
            try:
                self.lp.modifyColumnBounds({v_i.split('")')[0] + "_Rev" + '")':(0., 0.)}) # TODO _Rev does not imply not _Rev
            except:
                pass
        revFlag = 0
        if re.search(".*_Rev", v_j):
            revFlag = 1
            try:
                self.lp.modifyColumnBounds({v_j.replace('_Rev',''):(0., 0.)})
            except:
                pass
        if not revFlag:
            try:
                self.lp.modifyColumnBounds({v_j.split('")')[0] + "_Rev" + '")':(0., 0.)})
            except:
                pass
        try:
            self.lp.simplex()
            minRatio = self.lp.getObjVal()
            dualDict = self.lp.getShadowPriceDict()
            (rxn1minShadow, rxn2minShadow) = (dualDict[v_i], dualDict[v_j])
        except Exception, msg:
            minRatio = self.lp.getObjVal()            
            dualDict = self.lp.getShadowPriceDict()
            (rxn1minShadow, rxn2minShadow) = (dualDict[v_i], dualDict[v_j])
            fishyFlag = 1
        # fdTmp = fluxdist.FluxDist(self.lp).getActiveFluxDist()
        # print fdTmp
        # print '\n'
        # for elem in fdTmp:
        #     (sub, prod) = self.lp.getSubstratesAndProducts(elem[0])
        #     if "Mppic" in sub:
        #         print "sub", elem[0], elem[1],sub
        #     if "Mppic" in prod:
        #         print "prod", elem[0], elem[1],prod
        self.lp.setOptFlag('max')
        try:
            self.lp.simplex()
        except exceptions.SolutionUnbounded:
            maxRatio = 'inf'
        except Exception, msg:
            print msg
            maxRatio = self.lp.getObjVal()
            fishyFlag = 1
        else:
            maxRatio = self.lp.getObjVal()
        finally:
            dualDict = self.lp.getShadowPriceDict()
            (rxn1maxShadow, rxn2maxShadow) = (dualDict[v_i], dualDict[v_j])
        # fdTmp = fluxdist.FluxDist(self.lp).getActiveFluxDist()
        # print fdTmp
        # print '\n'
        # for elem in fdTmp:
        #     (sub, prod) = self.lp.getSubstratesAndProducts(elem[0])
        #     if "Mppic" in sub:
        #         print "sub", elem[0], elem[1], sub
        #     if "Mppic" in prod:
        #         print "prod", elem[0], elem[1], prod
        self.lp.initialize()
        return (minRatio, maxRatio, rxn1maxShadow, rxn1minShadow, rxn2maxShadow, rxn2minShadow, fishyFlag)

    def fluxCouplingFinder(self, rxns=None, printUncoupled=False):
        if not rxns:
            print "No reactions specified. Checking all reactions!"
            rxns=self.lp.getReactions()
        reacs = [r for r in rxns if not re.search('R_EX.*', r)]
        print "No. of reactions to check: ", len(reacs)
        bounds = self.lp.getColumnBounds()
        uncoupled = list()
        blocked = list()
        directionallyCoupled = list()
        fullyCoupled= list()
        partiallyCoupled = list()
        allreadyCoupled = list()
        for i, rxn1 in enumerate(reacs):
            if i in allreadyCoupled:
                continue
            if i in blocked:
                continue
            print "Reaction ", rxn1, " is being checked!"
            self.lp.setReactionObjective(rxn1)
            self.lp.setOptFlag('max')
            revFlag = 0
            if re.search(".*_Rev", rxn1):
                revFlag = 1
                try:
                    self.lp.modifyColumnBounds({rxn1.replace('_Rev',''):(0., 0.)})
                except:
                    pass
            if not revFlag:
                try:
                    self.lp.modifyColumnBounds({rxn1.split('")')[0] + "_Rev" + '")':(0., 0.)}) # TODO _Rev does not imply not _Rev
                except:
                    pass
            revFlag = 0
            self.lp.simplex()
            maxObjVal = self.lp.getObjVal()
            if maxObjVal <= 0:
                print colored(rxn1, 'red'), " is blocked"
                blocked.append(i)
                self.lp.initialize()
                continue
            self.lp.initialize()
            print maxObjVal
            print "still ", len(reacs) - i - len(blocked), " to check"
            print "directionallyCoupled: ", len(directionallyCoupled)
            print "fullyCoupled: ", len(fullyCoupled)
            print "partiallyCoupled: ", len(partiallyCoupled)
            for j, rxn2 in enumerate(reacs):
                if j <= i:
                    continue
                if j in blocked:
                    continue
                self.lp.setReactionObjective(rxn2)
                self.lp.setOptFlag('max')
                revFlag = 0
                if re.search(".*_Rev", rxn2):
                    revFlag = 1
                    try:
                        self.lp.modifyColumnBounds({rxn2.replace('_Rev',''):(0., 0.)})
                        if rxn2 == 'R("R_ADSL1r_Rev")':
                            print self.lp.getColumnBounds()['R("R_ADSL1r")']
                    except:
                        pass
                if not revFlag:
                    try:
                        self.lp.modifyColumnBounds({rxn2.split('")')[0] + "_Rev" + '")':(0., 0.)}) # TODO _Rev does not imply not _Rev
                    except:
                        pass
                revFlag = 0
                try:
                    self.lp.simplex()
                except:
                    continue
                if rxn2 == 'R("R_ADSL1r_Rev")':
                    print self.lp.getColumnBounds()['R("R_ADSL1r_Rev")']
                    print self.lp.getColumnBounds()['R("R_ADSL1r")']
                    print self.lp.getObjVal()

                if self.lp.getObjVal() < 1.:
                    print colored(rxn2, 'red'), " is blocked"
                    blocked.append(j)
                    self.lp.initialize()
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
                    if printUncoupled:
                        print "uncoupled: ", rxn1, rxn2
                        print info
                    uncoupled.append((rxn1, rxn2))
                    continue
                elif chop(rxn2maxShadow) > 0. and chop(rxn2minShadow) == 0.:
                    directionallyCoupled.append((rxn1, rxn2))
                    print colored("directionallyCoupled: ", 'red'), rxn1, rxn2
                    print info
                    continue
                elif chop(rxn2maxShadow) <= 0. and chop(rxn2minShadow) > 0.:
                    directionallyCoupled.append((rxn2,rxn1))
                    print colored("reverse directionallyCoupled: ", 'red'), rxn1, rxn2
                    print info
                    continue
                elif chop(rxn2maxShadow - rxn2minShadow) == 0. and chop(r_min - r_max) == 0.:
                    allreadyCoupled.append(j)
                    fullyCoupled.append((rxn1, rxn2))
                    print colored("fullycoupled: ", 'red'), rxn1, rxn2
                    print info
                    continue
                elif chop(r_min - r_max) != 0.:
                    print colored("partiallycoupled: ", 'red'), rxn1, rxn2
                    print info
                    (r_min, r_max, rxn1maxShadow, rxn1minShadow, rxn2maxShadow, rxn2minShadow, fishy) = self.computeFluxRatio(rxn2, rxn1)
                    info = ' '.join(("------> minRatio: ", str(chop(r_min)), "maxRatio: ", str(chop(r_max)), "rxn1maxShadow: ", str(chop(rxn1maxShadow)), "rxn1minShadow: ", str(chop(rxn1minShadow)), "rxn2maxShadow: ", str(chop(rxn2maxShadow)), "rxn2minShadow: ", str(chop(rxn2minShadow))))
                    print info
                    if chop(rxn2maxShadow) > 0. and chop(rxn2minShadow) == 0.:
                        print colored("!!!directionallyCoupled: ", 'red'), rxn2, rxn1
                        directionallyCoupled.append((rxn2, rxn1))
                    else:
                        allreadyCoupled.append(j)
                        partiallyCoupled.append((rxn1, rxn2))
                else:
                    raise Exception, ' '.join('Coupling for', rxn1, 'and', rxn2, 'could not be determined!')
        return {"fullyCoupled":tuple(fullyCoupled), "partiallyCoupled":tuple(partiallyCoupled), "directionallyCoupled":tuple(directionallyCoupled)}

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

def scaleColumnBounds(lp, scale=1.):
    rxns = lp.getColumnIDs()
    columnBounds = lp.getColumnBounds()
    newBounds = dict()
    for r in rxns:
        tmpBounds = columnBounds[r]
        newBounds[r] = (tmpBounds[0] * scale, tmpBounds[1] * scale)
    lp.modifyColumnBounds(newBounds)
    return copy.copy(lp)

def fluxVariabilityMain(lp):
    varGenerator = Variability(lp).variabilityGenerator()
    vResult = []
    for var in varGenerator:
        print var
        vResult.append(var)
    # open('VariabilityAnalysis.txt', 'w').write(repr(vResult))
    return vResult

def prepareLP(lp, biomRxn=None):
    if biomRxn:
        lp.deleteReactions([biomRxn])
        (sub, prod) = lp.getSubstratesAndProducts(biomRxn)
        bioMets = sub
        lp.addMetaboliteDrains(bioMets)
    lp.modifyColumnBounds(dict([(r, ('-inf', 10000000.))for r in lp.getTransporters()]))
    lp.modifyColumnBounds(dict([(r, (0., 10000000.))for r in lp.getReactions()]))
    lp.eraseHistory()
    return lp

def fluxCouplingMain(lp, rxns=None):
    if not rxns:
        rxns = lp.getReactions()
    print lp
    # blocked = Variability(lp).getBlocked(rxns)
    blocked = Variability(lp).getBlocked()
    print blocked
    reactions = list(set(rxns) - set(blocked))
    print reactions
    lp.deleteReactionsFromStoich(blocked)
    lp.eraseHistory()
    varObject = Variability(lp)
    print "\n"
    print "Checking if every reaction can reach a flux of 1!"
    for r in reactions:
        var = varObject.variabilityQ(r)
        if var[1] <= 1.:
            print r, var
            raise Exception, "Reaction " + r + " cannot exceed a flux >=1"
    lp.initialize()
    print '\n'
    print 'Starting flux coupling analysis!'
    print '\n'
    lp.smcp.presolve = glpk.GLP_OFF
    print lp
    
    return FluxCoupling(lp).fluxCouplingFinder(rxns=reactions)

if __name__ == '__main__':
    path = 'iAF1260templateFluxCoupling.lp'
    lp = prepareLP(metabolism.Metabolism(util.ImportCplex(path, terminal="OFF")))
    # print fluxCouplingMain(lp, rxns=lp.getReactions())
    print fluxCouplingMain(lp, rxns=('R("R_AACPS7")', 'R("R_GMHEPK")'))

