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
from ifba.general.combinatorics import SetCombine
import copy
import random


class ContextFBA(metabolism.Metabolism):
    """Implementation of Context-specific FBA.
    
    lp -> glpk problem struct
    rmfID -> Specification of the cell objective
    level -> The percentage of the normal FBA objective flux that should be
        fixed in the following contextFBA simulations
    
    Becker und Palsson. Context-specific metabolic networks are consistent
    with experiments. PLoS Comput Biol (2008) vol. 4 (5) pp. e1000082
    """
    def __init__(self, lp, rmfID=None, level=.9):
        metabolism.Metabolism.__init__(self, lp)
        self.rmfID = rmfID
        self.level = level
        self._partA(level=self.level)
        self.history = []
    
    def _partA(self, level=None):
        """docstring for _partA"""
        self.setObjectiveFunction({self.rmfID : 1.})
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
        # inconsistency = 0.
        # for key, val in contxtObjective.items():
        #     inconsistency += contxtObjective[key] * contextFluxDist[key]
        # return inconsistency
        iterator = (contxtObjective[key] * contextFluxDist[key] \
                    for key, val in contxtObjective.items())
        return sum(iterator)
    
    def contextFBA(self, expData, cutoff=None):
        """Uses an weighted reaction set in form of a dictionary """
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
    file.close()
    return dict([(reac.replace('[', '(').replace(']', ')'), float(val)) for (reac, val) in data])

def main():
    dataPaths = tuple([file.rstrip() for file in os.popen("ls *.tsv")])
    
    model_path = '../../work/HumanFBA/validHumanFBAmodel.lp'
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
    model_path = '../../work/HumanFBA/validHumanFBAmodel.lp'
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
        
def cntxtAnalysis(path, outpath):
    dataPaths = tuple([file.rstrip() for file in os.popen("ls "+path+"*.tsv")])
    print dataPaths
    outputFile = open(outpath, 'w').close()
    outputFile = open(outpath, 'a')
    model_path = '../../work/HumanFBA/validHumanFBAmodel.lp'
    stdFBAfluxDist = metabolism.Metabolism(util.ImportCplex(model_path)).fba()
    print "Standard FBA fluxDist: ", stdFBAfluxDist.getActiveFluxDist()
    print 2*"\n"
    reactionSets = list()
    lp = util.ImportCplex(model_path)
    cntxtFBA = ContextFBA(lp, rmfID='R("R_Obj")', level=.8)
    for path in dataPaths:
        (cntxtFluxDist, incon) = cntxtFBA.contextFBA(loadReactionData(path), \
        cutoff=.8)
        print "results for: ", path
        print "ContextFBA fluxDist: ", cntxtFluxDist.getActiveFluxDist()
        print "Inconsistency Score: ", incon
        print 2*'\n'
        reactionSets.append(set(cntxtFluxDist.getActiveReactions()))
        outputFile.write(path + '\t'+ str(incon) + "\n")
    outputFile.close()


if __name__ == '__main__':
    print cntxtAnalysis("/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/", \
    "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/eberResults.tsv")
    print cntxtAnalysis("/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/WT/", \
    "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/wtResults.tsv")
    

