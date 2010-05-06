#!/usr/bin/env python
# encoding: utf-8
"""
ContextFBA.py

Created by Nikolaus Sonnenschein on 2009-01-13.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

import re
import copy
import numpy
import pickle
from ifba.GlpkWrap import util, metabolism, randomMedia, fluxdist, glpk
from ifba.general.util import sumDicts, filterDict, dict2mathematica, dict2tsv
from ifba.general.combinatorics import SetCombine
from ifba.glpki import glpki


class ContextFBA(metabolism.Metabolism):
    """Implementation of Context-specific FBA.
    
    lp -> glpk problem struct
    rmfID -> Specification of the cell objective
    level -> The percentage of the normal FBA objective flux that should be
        fixed in the following contextFBA simulations
    
    Becker und Palsson. Context-specific metabolic networks are consistent
    with experiments. PLoS Comput Biol (2008) vol. 4 (5) pp. e1000082
    
    Furthermore it implements a new approach where the number of (and not the flux
    through) unexpressed reactions is minimized.
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
        colIDs = self.getColumnIDs()
        for iD, val in expData.items():
            if iD in colIDs:
                if val > cutoff:
                    expData[iD] = 0.
                else:
                    expData[iD] = cutoff - val
                revId = iD.split('")')[0] + "_Rev" + '")'
                if revId in colIDs:
                    if val > cutoff:
                        expData[revId] = 0.
                    else:
                        expData[revId] = cutoff - val
        return expData
    
    def computeInconsistency(self, contxtObjective, contextFluxDist):
        """docstring for _computeInconsistency"""
        return sum(numpy.array(contxtObjective) * numpy.array(contextFluxDist.getFluxArray()))
        
    def getContributors(self, fluxDist):
        obj = self.getObjectiveFunction()
        contribs = dict()
        for reac in obj:
            flux = fluxDist[reac]
            if flux != 0:
                contribs[reac] = fluxDist[reac]
        return contribs
    
    def contextFBA(self, expData, cutoff=None):
        """Uses a weighted reaction set in form of a dictionary """
        contxtObj = self._generateContextObjective(copy.copy(expData), cutoff)
        self.setObjectiveFunction(contxtObj)
        self.setOptFlag("Min")
        contextFluxDist = self.fba()
        # return (contextFluxDist, self.computeInconsistency(self.getObjective(), contextFluxDist), self.getContributors(contextFluxDist))
        return (contextFluxDist, self.computeInconsistency(self.getObjective(), contextFluxDist))


class MilpContextFBA(ContextFBA):
    """docstring for MilpContextFBA"""
    def __init__(self, lp, rmfID=None, level=.9):
        super(MilpContextFBA, self).__init__(lp, rmfID=rmfID, level=level)
        self._add_switches()
    
    def _generate_milp_context_objective(self, reacs):
        """docstring for _generate_milp_context_objective"""
        milpObjective = dict()
        for r in reacs:
            milpObjective["switch_"+r] = 1.
        return milpObjective
    
    def _add_switches(self):
        """Adds new columns to the lp.
        columns is a dictionary of the form {columnID : (lb, ub, coeffList)}
        """
        reacs = self.getColumnIDs()
        cols = dict()
        for r in reacs:
            cols["switch_"+r] = (0., 1., dict())
        print 'adding columns'
        self.addColumns(cols)
        print 'added columns'
        colKinds = dict()
        for r in reacs:
            colKinds[self.translateColumnNames(("switch_"+r,))[0]] = glpki.GLP_BV
        self.setColumnKinds(colKinds)
        colBounds = self.getColumnBounds()
        newRows = dict()
        for r in reacs:
            # print "processing ", r
            viPos = self.translateColumnNames((r,))[0]
            yPos = self.translateColumnNames(("switch_"+r,))[0]
            coefUb = {viPos:1., yPos:-colBounds[r][1]}
            coefLb = {viPos:1., yPos:-colBounds[r][0]}
            newRows["Ub_switch_"+r] = ('inf', 0., coefUb)
            newRows["Lb_switch_"+r] = (0., 'inf', coefLb)
        print 'adding rows'
        self.addRows(newRows)
        print 'added rows'

    def milpContextFBA(self, expData, cutoff=None):
        """docstring for milp_context_fba"""
        print expData['R("R_GLNS")']
        bannedReacs = self._generateContextObjective(copy.copy(expData), cutoff).keys()
        print bannedReacs
        milpContxtObj = self._generate_milp_context_objective(bannedReacs)
        self.setObjectiveFunction(milpContxtObj)
        self.setOptFlag("Min")
        contextFluxDist = self.fba(method='intopt')
        return (contextFluxDist, dict([(id, val) for id, val in self.mipValues().items() if val == 1. and id in milpContxtObj]))

def loadReactionData(path):
    """docstring for loadExpData"""
    file = open(path, 'r')
    data = []
    for line in file:
        data.append(tuple(line.rstrip().split('\t')))
    file.close()
    return dict([(reac.replace('[', '(').replace(']', ')'), float(val)) for (reac, val) in data])

def loadMultipleReactionData(path, separator="///"):
    """docstring for loadExpData"""
    file = open(path, 'r')
    data = []
    tmpData = []
    for line in file:
        if re.search(separator, line):
            data.append(dict([(reac.replace('[', '(').replace(']', ')'), float(val)) for (reac, val) in tmpData]))
            tmpData = []
            continue
        tmpData.append(tuple(line.rstrip().split('\t')))
    data.append(dict([(reac.replace('[', '(').replace(']', ')'), float(val)) for (reac, val) in tmpData]))
    file.close()
    return data


if __name__ == '__main__':
    testReacs = ('R("R_GHMT2r")','R("R_GHMT2r_Rev")', 'R("R_GHMT2rm")', 'R("R_GHMT2rm_Rev")')
    model_path = '../../work/HumanFBA/validHumanFBAmodel.lp'
    # lp = util.ImportCplex(model_path)
    # cntxtFBA = ContextFBA(lp, rmfID='R("R_Obj")', level=.8)
    # import random
    # reacs = random.sample(cntxtFBA.getReactions(), 1)
    # print reacs
    # # cntxtFBA._add_switches(reacs)
    # cntxtFBA.milpContextFBAdebug(testReacs)
    # util.WriteCplex(cntxtFBA, 'withSwitches.lp')

    # lp = util.ImportCplex(model_path)
    # cntxtFBA = ContextFBA(lp, rmfID='R("R_Obj")', level=.8)
    # print cntxtFBA
    # import random
    # reacs = random.sample(cntxtFBA.getReactions(), 1)
    # print reacs
    # # cntxtFBA._add_switches(reacs)
    # print cntxtFBA.milpContextFBAdebug(('R("R_Obj")',))[0]
    # util.WriteCplex(cntxtFBA, 'withSwitches.lp')
    # print cntxtFBA.mipValues()
    import glob
    data = glob.glob("/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/*.tsv")
    print data
    lp = util.ImportCplex(model_path)
    milpCntxtFBA = MilpContextFBA(lp, rmfID='R("R_Obj")', level=.8)
    util.WriteCplex(milpCntxtFBA, file_path='milpCntxtFBA.lp')
    # f = open('milpCntxtFBA.pickle', 'w')
    # pickle.dump(milpCntxtFBA, f)
    # f.close()
    # f = open('milpCntxtFBA.pickle', 'r')
    # milpCntxtFBA = pickle.load(f)
    # f.close()
    lp = util.ImportCplex(model_path)
    cntxtFBA = ContextFBA(lp, rmfID='R("R_Obj")', level=.8)
    for dat in data:
        print dat
        tmpDat = loadReactionData(dat)
        print milpCntxtFBA.milpContextFBA(tmpDat, cutoff=0.)[-1]
        print cntxtFBA.contextFBA(tmpDat, cutoff=0.)[-1]

    # print result[0]
    # print result[1]
    # util.WriteCplex(cntxtFBA, 'withSwitches.lp')


    # print loadReactionData("/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/HB004PO_EBER0_6.tsv")
    # (cntxtFluxDist, incon) = cntxtFBA.contextFBA(loadReactionData("/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/HB004PO_EBER0_6.tsv"), .5)
    # print incon
