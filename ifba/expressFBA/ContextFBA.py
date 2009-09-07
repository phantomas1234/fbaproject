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
import numpy


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
        # print self.getObjective()
        # rmfFlux = self.fba()[self.rmfID]
        self.simplex()
        rmfFlux = self.getObjVal()
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
        
    def _generate_milp_context_objective(self):
        """docstring for _generate_milp_context_objective"""
        pass
    
    def _prep_for_switches(self):
        """docstring for _prep_for_switches"""
        pass

    def _add_switches(self, reacs):
        """Adds new columns to the lp.
        columns is a dictionary of the form {columnID : (lb, ub, coeffList)}
        """
        cols = dict()
        switches = ["switch_"+r for r in reacs]
        for r in switches:
            cols[r] = (0., 1., [0. for m in xrange(self.getNumCols())])
        self.addColumns(cols)
        colKinds = dict()
        # for r in switches
        # self.setColumnKinds()
    
    # def computeInconsistency(self, contxtObjective, contextFluxDist):
    #     """docstring for _computeInconsistency"""
    #     stuff = contxtObjective.items()
    #     iterator = (val * contextFluxDist[key] for key, val in stuff)
    #     # print contxtObjective
    #     return sum(iterator)

    def computeInconsistency(self, contxtObjective, contextFluxDist):
        """docstring for _computeInconsistency"""
        return sum(numpy.array(contxtObjective) * numpy.array(contextFluxDist.getFluxArray()))
    
    def contextFBA(self, expData, cutoff=None):
        """Uses a weighted reaction set in form of a dictionary """
        contxtObj = self._generateContextObjective(copy.copy(expData), cutoff)
        self.setObjectiveFunction(contxtObj)
        self.setOptFlag("Min")
        contextFluxDist = self.fba()
        return (contextFluxDist, self.computeInconsistency(self.getObjective(), contextFluxDist))        
        # return (contextFluxDist, self.computeInconsistency(self.getObjectiveFunction(), contextFluxDist))
        
    def milp_context_fba(self):
        """docstring for milp_context_fba"""
        contxtObj = self._generateContextObjective(copy.copy(expData), cutoff)
        self.setObjectiveFunction(contxtObj)



def loadReactionData(path):
    """docstring for loadExpData"""
    file = open(path, 'r')
    data = []
    for line in file:
        data.append(tuple(line.rstrip().split('\t')))
    file.close()
    return dict([(reac.replace('[', '(').replace(']', ')'), float(val)) for (reac, val) in data])


if __name__ == '__main__':
    model_path = '../../work/HumanFBA/validHumanFBAmodel.lp'
    lp = util.ImportCplex(model_path)
    cntxtFBA = ContextFBA(lp, rmfID='R("R_Obj")', level=.8)
    import random
    reacs = random.sample(cntxtFBA.getReactions(), 7)
    print reacs
    # cntxtFBA._add_switches()
    # print loadReactionData("/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/HB004PO_EBER0_6.tsv")
    # (cntxtFluxDist, incon) = cntxtFBA.contextFBA(loadReactionData("/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/HB004PO_EBER0_6.tsv"), .5)
    # print incon
