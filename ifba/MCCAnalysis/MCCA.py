#!/usr/bin/env python
# encoding: utf-8
"""
MCCA.py

Created by Nikolaus Sonnenschein on 2010-02-01.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import re
import copy
import time
from ifba.GlpkWrap import util, exceptions
from ifba.GlpkWrap.metabolism import Metabolism
from ifba.glpki.glpki import *
from ifba.general.util import chop
from termcolor import colored


class AnalysisOfConservedPools(Metabolism):
    """Base class for MCCA and MCPI
    
    Nikolaev et al. Elucidation and structural analysis of conserved pools for 
    genome-scale metabolic reconstructions. Biophys J (2005) vol. 88 (1) pp. 37-49
    """
    def __init__(self, lp):
        super(AnalysisOfConservedPools, self).__init__(self.__prepareProblem(lp))

    def __prepareProblem(self, lp):
        """Prepares the lp construct for analysis
        
        Transposes the stoichiometric matrix S_ij -> S^T
        Columns(Reactions) become rows and rows(Metabolites) become columns
        Previous column, now row bounds are set to 0
        Previous row, now column bounds are set to > 0
        """
        oldProb = Metabolism(lp)
        biomR = [r for r in oldProb.getColumnIDs() if re.search(".*.*", r, re.IGNORECASE)][0]
        print biomR
        oldProb.deleteReactionsFromStoich([biomR])
        oldProb.eraseHistory()
        newProb = Metabolism(glp_create_prob())
        for elem in range(1, oldProb.getNumCols() + 1):
            columnsDict = dict()
            identifier = glp_get_col_name(oldProb.lp, elem)
            colCoefList = oldProb.getColumnCoef(elem)
            columnsDict[identifier] = (0, 0, colCoefList)
            newProb.addRows(columnsDict)
        for elem in range(1, oldProb.getNumRows() + 1):
            rowsDict = dict()
            identifier = glp_get_row_name(oldProb.lp, elem)
            rowCoefList = oldProb.getRowCoef(elem)
            rowsDict[identifier] = (0, 10000, rowCoefList)
            newProb.addColumns(rowsDict)
        newProb.eraseHistory()
        return newProb.lp


class MCCA(AnalysisOfConservedPools):
    """Implements metabolite concentration coupling analysis (MCCA)
    
    Nikolaev et al. Elucidation and structural analysis of conserved pools for 
    genome-scale metabolic reconstructions. Biophys J (2005) vol. 88 (1) pp. 37-49
    """
    def __init__(self, lp):
        super(MCCA, self).__init__(lp)
        
    def computeMinMaxRatio(self, m_i, m_j):
        """Computes the the ration of v1/v2."""
        fishyFlag = 0
        self.setReactionObjective(m_i)
        self.modifyColumnBounds({m_j : (1., 1.)})
        self.setOptFlag('min')
        try:
            self.simplex()
            minRatio = self.getObjVal()
            dualDict = self.getShadowPriceDict()
            (rxn1minShadow, rxn2minShadow) = (dualDict[m_i], dualDict[m_j])
        except Exception, msg:
            minRatio = self.getObjVal()            
            dualDict = self.getShadowPriceDict()
            (rxn1minShadow, rxn2minShadow) = (dualDict[m_i], dualDict[m_j])
            fishyFlag = 1
        self.setOptFlag('max')
        try:
            self.simplex()
        except exceptions.SolutionUnbounded:
            maxRatio = 'inf'
        except Exception, msg:
            print msg
            maxRatio = self.getObjVal()
            fishyFlag = 1
        else:
            maxRatio = self.getObjVal()
        finally:
            dualDict = self.getShadowPriceDict()
            (rxn1maxShadow, rxn2maxShadow) = (dualDict[m_i], dualDict[m_j])
        self.initialize()
        return (minRatio, maxRatio, rxn1maxShadow, rxn1minShadow, rxn2maxShadow, rxn2minShadow, fishyFlag)

    def mcca(self, metabolites=None, printUncoupled=True):
        if not metabolites:
            print "No metabolites specified. Checking all metabolites!"
            metabolites=self.getColumnIDs()
        bounds = self.getColumnBounds()
        uncoupled = list()
        blocked = list()
        directionallyCoupled = list()
        fullyCoupled= list()
        partiallyCoupled = list()
        allreadyCoupled = list()
        for i, met1 in enumerate(metabolites):
            if i in allreadyCoupled:
                continue
            if i in blocked:
                continue
            print "Metabolite ", met1, " is being checked!"
            self.setReactionObjective(met1)
            self.setOptFlag('max')
            try:
                self.simplex()
                maxObjVal = self.getObjVal()
            except exceptions.SolutionUnbounded:
                maxObjVal = 'inf'
            if maxObjVal <= 0:
                print colored(met1, 'red'), " is blocked"
                blocked.append(i)
                self.initialize()
                continue
            self.initialize()
            print "still ", len(metabolites) - i , " to check"
            print "directionallyCoupled: ", len(directionallyCoupled)
            print "fullyCoupled: ", len(fullyCoupled)
            print "partiallyCoupled: ", len(partiallyCoupled)
            for j, met2 in enumerate(metabolites):
                if j <= i:
                    continue
                if j in blocked:
                    continue
                self.setReactionObjective(met2)
                self.setOptFlag('max')
                util.WriteCplex(self, 'debug.lp')
                try:
                    self.simplex()
                    maxObjVal = self.getObjVal()
                except:
                    continue
                if self.getObjVal() < 1.:
                    print colored(met2, 'red'), " is blocked"
                    blocked.append(j)
                    self.initialize()
                    continue
                self.initialize()
                (r_min, r_max, met1maxShadow, met1minShadow, met2maxShadow, met2minShadow, fishy) = self.computeMinMaxRatio(met1, met2)
                info = ' '.join(("------> minRatio: ", str(chop(r_min)), "maxRatio: ", str(chop(r_max)), "met1maxShadow: ", str(chop(met1maxShadow)), "met1minShadow: ", str(chop(met1minShadow)), "met2maxShadow: ", str(chop(met2maxShadow)), "met2minShadow: ", str(chop(met2minShadow))))
                if fishy:
                    print "fishy: ", met1, " ", met2, " ", fishy
                    self.smcp.presolve = glpk.GLP_ON
                    (r_min, r_max, met1maxShadow, met1minShadow, met2maxShadow, met2minShadow, fishy) = self.computeMinMaxRatio(met1, met2)
                    info = ' '.join(("------> minRatio: ", str(chop(r_min)), "maxRatio: ", str(chop(r_max)), "met1maxShadow: ", str(chop(met1maxShadow)), "met1minShadow: ", str(chop(met1minShadow)), "met2maxShadow: ", str(chop(met2maxShadow)), "met2minShadow: ", str(chop(met2minShadow))))
                    self.smcp.presolve = glpk.GLP_OFF
                    print "fishy: ", met1, " ", met2, " ", fishy
                    if fishy:
                        continue
                    # sys.exit(-1)
                    continue
                if chop(met2maxShadow) <= 0. and chop(met2minShadow) == 0.:
                    if printUncoupled:
                        print "uncoupled: ", met1, met2
                        print info
                    uncoupled.append((met1, met2))
                    continue
                elif chop(met2maxShadow) > 0. and chop(met2minShadow) == 0.:
                    directionallyCoupled.append((met1, met2))
                    print colored("directionallyCoupled: ", 'red'), met1, met2
                    print info
                    continue
                elif chop(met2maxShadow) <= 0. and chop(met2minShadow) > 0.:
                    directionallyCoupled.append((met2,met1))
                    print colored("reverse directionallyCoupled: ", 'red'), met1, met2
                    print info
                    continue
                elif chop(met2maxShadow - met2minShadow) == 0. and chop(r_min - r_max) == 0.:
                    allreadyCoupled.append(j)
                    fullyCoupled.append((met1, met2))
                    print colored("fullycoupled: ", 'red'), met1, met2
                    print info
                    continue
                elif chop(r_min - r_max) != 0.:
                    print colored("partiallycoupled: ", 'red'), met1, met2
                    print info
                    (r_min, r_max, met1maxShadow, met1minShadow, met2maxShadow, met2minShadow, fishy) = self.computeMinMaxRatio(met2, met1)
                    info = ' '.join(("------> minRatio: ", str(chop(r_min)), "maxRatio: ", str(chop(r_max)), "met1maxShadow: ", str(chop(met1maxShadow)), "met1minShadow: ", str(chop(met1minShadow)), "met2maxShadow: ", str(chop(met2maxShadow)), "met2minShadow: ", str(chop(met2minShadow))))
                    print info
                    if chop(met2maxShadow) > 0. and chop(met2minShadow) == 0.:
                        print colored("!!!directionallyCoupled: ", 'red'), met2, met1
                        directionallyCoupled.append((met1, met2))
                    else:
                        allreadyCoupled.append(j)
                        partiallyCoupled.append((met1, met2))
                else:
                    raise Exception, ' '.join('Coupling for', met1, 'and', met2, 'could not be determined!')
        return {"fullyCoupled":tuple(fullyCoupled), "partiallyCoupled":tuple(partiallyCoupled), "directionallyCoupled":tuple(directionallyCoupled)}
    



class MCPI(AnalysisOfConservedPools):
    """Implements minimal conserved pools identification (MCPI)
    
    Nikolaev et al. Elucidation and structural analysis of conserved pools for 
    genome-scale metabolic reconstructions. Biophys J (2005) vol. 88 (1) pp. 37-49    
    """
    def __init__(self, lp):
        super(MCPI, self).__init__(lp)
        

if __name__ == '__main__':
    import pickle
    # mcpiObj = MCCA(util.ImportCplex("../models/Ecoli_Core_template.lp"))
    mcpiObj = MCCA(util.ImportCplex("../models/Ecoli_Core_template_rev.lp"))
    # mcpiObj = MCCA(util.ImportCplex("../models/iAF1260template.lp"))
    rxns2remove = [r for r in mcpiObj.getRowIDs() if re.search('.*_Transp.*', r)]
    mcpiObj.deleteMetabolitesFromStoich(rxns2remove)
    mcpiObj.eraseHistory()
    util.WriteCplex(mcpiObj, 'debug.lp')
    result = mcpiObj.mcca(printUncoupled=True)
    pickle.dump(result, open('EcoliCoreMetaboliteCouplings.pcl', 'w'))
    
    
    
    