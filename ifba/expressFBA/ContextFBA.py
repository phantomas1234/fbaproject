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


if __name__ == '__main__':
    pass