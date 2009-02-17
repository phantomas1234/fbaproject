#!/usr/bin/env python
# encoding: utf-8
"""
ContextFBA.py

Created by Nikolaus Sonnenschein on 2009-01-13.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from ifba.GlpkWrap import util, metabolism, randomMedia, fluxdist, glpk
from ifba.general.util import sumDicts, filterDict, dict2mathematica, dict2tsv
from ifba.combinatorics.combinatorics import SetCombine
import copy
import random


class ContextFBA(metabolism.Metabolism):
    """docstring for ClassName"""
    def __init__(self, lp, rmfID=None, level=.9):
        metabolism.Metabolism.__init__(self, lp)
        self.rmfID = rmfID
        self.level = level
        self._partA(level=self.level)
        self.history = []
    
    def _partA(self, level=None):
        """docstring for _partA"""
        self.setObjectiveFunction({self.rmfID : 1.})
        print self.getObjectiveFunction()
        rmfFlux = self.fba()[self.rmfID]
        leveledRmf = rmfFlux * level
        self.modifyColumnBounds({self.rmfID:(leveledRmf, rmfFlux)})
    
    def _generateContextObjective(self, expData, cutoff):
        """docstring for _generateContextObjective"""
        for iD, val in expData.items():
            if val > cutoff:
                expData[iD] = 0.
            else:
                expData[iD] = cutoff - val
        return expData
    
    def computeInconsistency(self, contxtObjective, contextFluxDist):
        """docstring for _computeInconsistency"""
        inconsistency = 0.
        for key, val in contxtObjective.items():
            inconsistency += contxtObjective[key] * contextFluxDist[key]
        return inconsistency
    
    def contextFBA(self, expData, cutoff=None):
        contxtObj = self._generateContextObjective(expData, cutoff)
        self.setObjectiveFunction(contxtObj)
        self.setOptFlag("Min")
        contextFluxDist = self.fba()
        return (contextFluxDist, self.computeInconsistency(contxtObj, contextFluxDist))


def loadReactionData(path):
    """docstring for loadExpData"""
    file = open(path, 'r')
    data = []
    for line in file:
        data.append(tuple(line.rstrip().split('\t')))
    return dict([(reac.replace('[', '(').replace(']', ')'), float(val)) for (reac, val) in data])

def main():
    dataPaths = tuple([file.rstrip() for file in os.popen("ls *.tsv")])
    
    model_path = '../HumanFBA/validHumanFBAmodel.lp'
    stdFBAfluxDist = metabolism.Metabolism(util.ImportCplex(model_path)).fba()
    print "Standard FBA fluxDist: ", stdFBAfluxDist.getActiveFluxDist()
    
    print 2*"\n"
    reactionSets = list()
    for path in dataPaths:
        cntxtFBA = ContextFBA(util.ImportCplex(model_path), rmfID='R("R_Obj")', level=.5)
        # print cntxtFBA.getColumnBounds()['R("R_Obj")']
        (cntxtFluxDist, incon) = cntxtFBA.contextFBA(loadReactionData(path), cutoff=.9)
        print "results for: ", path
        print "ContextFBA fluxDist: ", cntxtFluxDist.getActiveFluxDist()
        print "Inconsistency Score: ", incon
        print 2*'\n'
        reactionSets.append(set(cntxtFluxDist.getActiveReactions()))
    for i in range(0, 3):
        for j in range(0,3):
            print dataPaths[i], dataPaths[j]
            print float(len(reactionSets[i].intersection(reactionSets[j])))/float(max(len(reactionSets[i]),len(reactionSets[j])))

def main2():
    model_path = '../HumanFBA/validHumanFBAmodel.lp'
    path = "pUC-18_EBER2-L2.tsv"
    cntxtFBA = ContextFBA(util.ImportCplex(model_path, terminal="ON"), rmfID='R("R_MCDm")', level=.8)
    reacs = cntxtFBA.getReactions()
    for r in reacs:
        print "checking ", r
        cntxtFBA = ContextFBA(util.ImportCplex(model_path, terminal="OFF"), rmfID=r, level=.8)
        (cntxtFluxDist, incon) = cntxtFBA.contextFBA(loadReactionData(path), cutoff=.9)
        print "results for: ", r, " for ", path
        # print "ContextFBA fluxDist: ", cntxtFluxDist.getActiveFluxDist()
        print "Inconsistency Score: ", incon


if __name__ == '__main__':
    main2()
    



# Old Stuff

# if __name__ == '__main__':
#
#     def posNegSplit(real):
#         if real > 0.:
#             return 1
#         else:
#             return -1
#
#     threshold1 = .8
#
#     dat = open('/Users/niko/arbeit/Publishing/MetabolicControlPaper2/code/test2.tsv', 'r')
#     datDict = dict()
#     for line in dat:
#         tmp = line.rstrip().split('\t')
#         datDict[tmp[0]] = float(tmp[1])
#     dat.close()
#
#     datDict2 = dict([(key, -1 * value) for key, value in datDict.items()])
#     print datDict2
#
#
#     path = '../models/iAF1260template2.lp'
#     lp = metabolism.Metabolism(util.ImportCplex(path))
#
#
#     lp.setReactionObjectiveMinimizeRest('R("R_Ec_biomass_iAF1260_core_59p81M")')
#
#     transp = lp.getTransporters()
#     superMedium = dict()
#     for i in transp:
#         superMedium[i] = (-20., 20.)
#     # lp.modifyColumnBounds({'R("MglcDb_Transp")':(0, 20), 'R("R_ATPM")':(8.39, 8.39), 'R("Mo2b_Transp")':(0, 18.5)})
#     lp.modifyColumnBounds(superMedium)
#     lp2 = copy.copy(lp)
#     fd = lp2.fba()
#     growth = fd['R("R_Ec_biomass_iAF1260_core_59p81M")']
#     print growth
#
#     lp3 = copy.copy(lp2)
#     lp3.modifyColumnBounds({'R("R_Ec_biomass_iAF1260_core_59p81M")':(growth*threshold1, growth)})
#
#     lp3.setObjectiveFunction(datDict)
#
#     fd = lp3.fba()
#     open("/Users/niko/arbeit/Data/fluxOne.tsv", 'w').write(fd.tsv())
#     growth = fd['R("R_Ec_biomass_iAF1260_core_59p81M")']
#     print growth
#     tmp1= fd.getActiveReactions()
#     tmp1.sort()
#
#
#
#     lp3.setObjectiveFunction(datDict2)
#
#     fd = lp3.fba()
#     open("/Users/niko/arbeit/Data/fluxTwo.tsv", 'w').write(fd.tsv())
#     growth = fd['R("R_Ec_biomass_iAF1260_core_59p81M")']
#     print growth
#     tmp2 = fd.getActiveReactions()
#     tmp2.sort()
#
#
#     print
#     print set(tmp1).difference(set(tmp2)), "\n"
#     print set(tmp2).difference(set(tmp1)),

    
    
    # context-specific fba
    #
    # print help(filterDict)
    # print filterDict({1:.1,2:.99},.3)
    #
    # threshold1 = .8
    # threshold2 = .8
    # path = '../models/iAF1260template2.lp'
    # lp = metabolism.Metabolism(util.ImportCplex(path))
    # surrogatData = dict([(r, random.uniform(0., 1.)) for r in lp.getReactions()])
    # print surrogatData
    # for key, value in surrogatData.items():
    #     if value < .8:
    #         surrogatData[key] = 1. - value
    #     else:
    #         surrogatData.pop(key)
    # print surrogatData
    #
    #
    # lp.setReactionObjectiveMinimizeRest('R("R_Ec_biomass_iAF1260_core_59p81M")')
    # lp.modifyColumnBounds({'R("MglcDb_Transp")':(0, 20), 'R("R_ATPM")':(8.39, 8.39), 'R("Mo2b_Transp")':(0, 18.5)})
    # lp2 = copy.copy(lp)
    # fd = lp2.fba()
    # growth = fd['R("R_Ec_biomass_iAF1260_core_59p81M")']
    # print growth
    # lp3 = copy.copy(lp2)
    # lp3.modifyColumnBounds({'R("R_Ec_biomass_iAF1260_core_59p81M")':(growth*threshold1, growth)})
    #
    # lp3.setObjectiveFunction(surrogatData)
    #
    # fd = lp3.fba()
    # growth = fd['R("R_Ec_biomass_iAF1260_core_59p81M")']
    # print growth
    # print lp3.getObjectiveFunction()
