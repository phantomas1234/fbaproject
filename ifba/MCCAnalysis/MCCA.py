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
        # oldProb = Metabolism(lp)
        # newProb = Metabolism(glp_create_prob())
        # numCols = oldProb.getNumCols()
        # numRows = oldProb.getNumRows()
        # print numCols, numRows
        # rowsDict = dict()
        # for elem in range(1, oldProb.getNumRows() + 1):
        #     identifier = glp_get_row_name(oldProb.lp, elem)
        #     rowCoefList = oldProb.getRowCoef(elem)
        #     rowsDict[identifier] = (0, 'inf', rowCoefList)
        # # print rowsDict
        # 
        # glp_add_cols(newProb.lp, numRows)
        # # glp_add_rows(newProb.lp, numCols)
        # print newProb
        # for i, item in enumerate(rowsDict):
        #     index = numCols+1+i
        #     lb = rowsDict[item][0]
        #     ub = rowsDict[item][1]
        #     ia = intArray(numCols+1)
        #     da = doubleArray(numCols+1)
        #     sparseCoef = rowsDict[item][2]
        #     print sparseCoef
        #     for j in range(1, numCols + 1):
        #         if j in sparseCoef:
        #             print 'sadf'
        #             ia[j] = j
        #             da[j] = sparseCoef[j]
        #         else:
        #             ia[j] = j
        #             da[j] = 0.
        #     print ia
        #     glp_set_col_name(newProb.lp, index, item)
        #     print "asdf"
        #     glp_set_obj_coef(newProb.lp, index, 0.)
        #     newProb.lp._setColumnBound(index, lb, ub)
        #     glp_set_mat_col(self.reciever.lp, index, numRows, ia, da)

        oldProb = Metabolism(lp)
        newProb = Metabolism(glp_create_prob())
        rowsDict = dict()
        for elem in range(1, oldProb.getNumRows() + 1):
            identifier = glp_get_row_name(oldProb.lp, elem)
            rowCoefList = oldProb.getRowCoef(elem)
            rowsDict[identifier] = (0, 'inf', rowCoefList)
        print rowsDict
        newProb.addColumns(rowsDict)
        
        columnsDict = dict()
        for elem in range(1, oldProb.getNumCols() + 1):
            identifier = glp_get_col_name(oldProb.lp, elem)
            colCoefList = oldProb.getColumnCoef(elem)
            print identifier, colCoefList
            columnsDict[identifier] = (0, 0, colCoefList)
        print columnsDict
        newProb.addRows(columnsDict)
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
        print m_i
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

    def mcca(metabolites=None):
        if not metabolites:
            print "No metabolites specified. Checking all reactions!"
            metabolites=self.lp.getColumnIDs()
        for i, met1 in enumerate(metabolites):
            if i in allreadyCoupled:
                continue
            if i in blocked:
                continue
            print "Metabolite ", met1, " is being checked!"
            self.lp.setReactionObjective(met1)
            self.lp.setOptFlag('max')
            self.lp.simplex()
            maxObjVal = self.lp.getObjVal()
            if maxObjVal <= 0:
                print colored(met1, 'red'), " is blocked"
                blocked.append(i)
                self.lp.initialize
                continue
            self.lp.initialize()
            print maxObjVal
            print "still ", len(reacs) - i - len(allreadyCoupled) - len(blocked), " to check"
            print "directionallyCoupled: ", len(directionallyCoupled)
            print "fullyCoupled: ", len(fullyCoupled)
            print "partiallyCoupled: ", len(partiallyCoupled)
            for j, met2 in enumerate(metabolites):
                if j <= i:
                    continue
                if j in blocked:
                    continue
                self.lp.setReactionObjective(met2)
                try:
                    self.lp.simplex()
                except:
                    continue
                if self.lp.getObjVal() < 1.:
                    print colored(met2, 'red'), " is blocked"
                    blocked.append(j)
                    self.lp.initialize
                    continue
                self.lp.initialize
                (r_min, r_max, met1maxShadow, met1minShadow, met2maxShadow, met2minShadow, fishy) = self.minMaxRatio(met1, met2)
                info = ' '.join(("------> minRatio: ", str(chop(r_min)), "maxRatio: ", str(chop(r_max)), "met1maxShadow: ", str(chop(met1maxShadow)), "met1minShadow: ", str(chop(met1minShadow)), "met2maxShadow: ", str(chop(met2maxShadow)), "met2minShadow: ", str(chop(met2minShadow))))
                if fishy:
                    print "fishy: ", met1, " ", met2, " ", fishy
                    self.lp.smcp.presolve = glpk.GLP_ON
                    (r_min, r_max, met1maxShadow, met1minShadow, met2maxShadow, met2minShadow, fishy) = self.minMaxRatio(met1, met2)
                    info = ' '.join(("------> minRatio: ", str(chop(r_min)), "maxRatio: ", str(chop(r_max)), "met1maxShadow: ", str(chop(met1maxShadow)), "met1minShadow: ", str(chop(met1minShadow)), "met2maxShadow: ", str(chop(met2maxShadow)), "met2minShadow: ", str(chop(met2minShadow))))
                    self.lp.smcp.presolve = glpk.GLP_OFF
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
                    (r_min, r_max, met1maxShadow, met1minShadow, met2maxShadow, met2minShadow, fishy) = self.computeFluxRatio(met2, met1)
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
    print Metabolism(util.ImportCplex("../models/Ecoli_Core_template.lp"))
    mcpiObj = MCCA(util.ImportCplex("../models/Ecoli_Core_template.lp"))
    print mcpiObj
    util.WriteCplex(mcpiObj, 'debug.lp')
    print mcpiObj.getColumnIDs()
    print mcpiObj.getRowIDs()
    print mcpiObj.getColumnBounds()
    print mcpiObj.getRowBounds()
    print mcpiObj.computeMinMaxRatio('Matpc', 'Madpc')

