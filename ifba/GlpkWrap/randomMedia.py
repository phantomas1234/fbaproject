#!/usr/bin/env python
# encoding: utf-8
"""
randomMedia.py

Created by Nikolaus Sonnenschein on 2008-02-18.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
import random
import string
import time
import copy
import re

from ifba.GlpkWrap import fluxdist, metabolism, util
from ifba.GlpkWrap.fluxdist import FBAsimulationResult
from ifba.glpki.glpki import *
from ifba.general.util import sumDicts, randomString

# random.seed(100)
# random.seed(101)
# TODO: Remove if not needed

INCLUDE = ('R("Mhb_Transp")',
         'R("Mna1b_Transp")',
         'R("Mkb_Transp")',
         'R("Mca2b_Transp")',
         'R("Mcu2b_Transp")',
         'R("Mmg2b_Transp")',
         'R("Mzn2b_Transp")',
         'R("Mmobdb_Transp")',
         'R("Mfe2b_Transp")',
         'R("Mfe3b_Transp")',
         'R("Mcobalt2b_Transp")',
         'R("Mmn2b_Transp")',
         'R("Mclb_Transp")')

class RandomMediaSimulations(object):
    """A class simplifying random media simulations."""
    def __init__(self, path2template, objective, includes, description, minimizeRest=True, optimizationRoutine='fba', koQ=True):
        self.path2template = path2template
        self.objective = objective
        self.includes = includes
        self.koQ = koQ
        self.descr = description
        self.tmpLP = metabolism.Metabolism(util.ImportCplex(self.path2template))
        # print tmpLP.pFBA().getActiveFluxes()
        if minimizeRest:
            self.tmpLP.setReactionObjectiveMinimizeRest(self.objective)
        else:
            self.tmpLP.setReactionObjective(self.objective)
        self.optimizationRoutine = optimizationRoutine
        self.tmpLP.eraseHistory()
        self.almaas = Almaas(self.tmpLP, alwaysInc=self.includes, optimizationRoutine=self.optimizationRoutine)
    
    def run(self, *args, **kwargs):
        """Run one random medium simulation and return a FBAsimulationResult object."""
        f = self.almaas.generateFluxdist()
        knockoutEffects = dict()
        self.koQ = False
        if self.koQ:
            self.almaas.lp.modifyColumnBounds(self.almaas.lastBounds)
            knockoutEffects = self.almaas.lp.singleKoAnalysis(f.getActiveReactions())
            wt = f[self.objective]
            for k in knockoutEffects:
                knockoutEffects[k] = knockoutEffects[k] / wt
        self.almaas.lp.initialize()
        # print dict([(k, v) for k, v in self.almaas.lastBounds.items() if re.search('.*_Transp.*', k)])
        return FBAsimulationResult(f, knockoutEffects, self.almaas.lastBounds, self.almaas.lp.getObjectiveFunction(), time.time(), self.path2template, self.descr)

    def __del__(self):
        """docstring for __del__"""
        glp_delete_prob(self.almaas.lp.lp) # FIXME this is a dirty hack
        del self

def dict2tsv(condDict):
    """Convert a dict into TSV format."""
    string = str()
    for i in condDict:
        string += i + "\t" + "{%f, %f}" % condDict[i] + "\n"
    return string

class Almaas(object):
    """Randomized media analysis.
    
    Adds the necessary functionality to the metabolism class.
    """
    def __init__(self, lp, default_bound=20, percRange=(5, 100), alwaysInc=set(), optimizationRoutine='fba'):
        super(Almaas, self).__init__()
        self.lp = lp
        self._preconditioning()
        self.def_bnd = default_bound
        self.currDict = {}
        self.percRange = percRange
        self.debugDict = {}
        self.alwaysInc = alwaysInc
        self.optimizationRoutine = optimizationRoutine
        self.lastBounds = None

    def _preconditioning(self):
        """docstring for _preconditioning"""
        transporters = self.lp.getTransporters()
        # columnBounds = self.lp.getColumnBounds()
        # for t in transporters:
        #     print t, columnBounds[t]
        d = dict()
        for t in transporters:
            d[t] = (-99999,0)
        self.lp.modifyColumnBounds(d)
        self.lp.eraseHistory()
        # transporters = self.lp.getTransporters()
        # columnBounds = self.lp.getColumnBounds()
        # print 10*'\n'
        # for t in transporters:
        #     print t, columnBounds[t]

    def _rndBndDict(self, list):
        boundDict = {}
        rnd = random.uniform
        for i in list:
            # boundDict[i] = (-rnd(0, self.def_bnd), rnd(0, self.def_bnd))
            boundDict[i] = (-99999, rnd(0, self.def_bnd))
        return boundDict
    
    def _randomPercantage(self):
        # inc = 10
        num = len(self.lp.getTransporters())
        # val = range(self.percRange[0], self.percRange[1] + inc , inc)
        # return (num/100) * random.choice(val)
        return (num/100) * random.randrange(self.percRange[0], self.percRange[1])
    
    def generateRandCond(self, num):
        transp = self.lp.getTransporters()
        sample = set(random.sample(transp, num)).union(self.alwaysInc)
        self.currDict = self._rndBndDict(sample)
        self.lp.modifyColumnBounds(self.currDict)
        self.lastBounds = self.lp.getColumnBounds()
        # tmpDict = {}
        # for i in self.currDict:
        #     tmpDict[i] = 1
        # self.debugDict = sumDicts(self.debugDict, tmpDict)
    
    def generateFluxdist(self, minGrowth=0.2):
        try:
            getattr(self.lp, self.optimizationRoutine)
        except AttributeError, msg:
            print "\nThe optimization routine seems to wrongly specified!"
            raise AttributeError, msg
        growth = 0.
        num = self._randomPercantage()
        while growth <= minGrowth:
            self.generateRandCond(num)
            try:
                fluxdist = getattr(self.lp, self.optimizationRoutine)()
                growth = fluxdist[self.lp.getReactionObjective()]
                self.lp.initialize()
            except Exception, msg:
                # print msg
                growth = 0.
                self.lp.initialize()
        # print self.currDict
        return fluxdist

def main(path2template, resultsPath, runs):
    import gzip
    tmpLP = metabolism.Metabolism(util.ImportCplex(path2template))
    tmpLP.setReactionObjectiveMinimizeRest('R("R_Ec_biomass_iAF1260_core_59p81M")')
    include = ('R("Mhb_Transp")',
             'R("Mna1b_Transp")',
             'R("Mkb_Transp")',
             'R("Mca2b_Transp")',
             'R("Mcu2b_Transp")',
             'R("Mmg2b_Transp")',
             'R("Mzn2b_Transp")',
             'R("Mmobdb_Transp")',
             'R("Mfe2b_Transp")',
             'R("Mfe3b_Transp")',
             'R("Mcobalt2b_Transp")',
             'R("Mmn2b_Transp")',
             'R("Mclb_Transp")')
    lp = Almaas(copy.copy(tmpLP), alwaysInc=include)
    print lp.lp
    f = lp.generateFluxdist()
    for item in range(runs):
        f = lp.generateFluxdist()
        stringDump = dict2tsv(lp.currDict) + "\n" + f.tsv()
        print stringDump
        path = resultsPath + "iAf1260_fluxDist_" + str(item) + ".tsv.gz"
        lp.lp.initialize()
        gzip.open(path, 'w').write(stringDump)

def main2(path2template, resultsPath, runs):
    import gzip
    tmpLP = metabolism.Metabolism(util.ImportCplex(path2template))
    tmpLP.setReactionObjectiveMinimizeRest('R("R_BiomassEcoli")')
    lp = Almaas(copy.copy(tmpLP))
    print lp.lp
    f = lp.generateFluxdist()
    print f
    # for item in range(runs):
    #     f = lp.generateFluxdist()
    #     stringDump = dict2tsv(lp.currDict) + "\n" + f.tsv()
    #     print stringDump
    #     path = resultsPath + "iAf1260_fluxDist_" + str(item) + ".tsv.gz"
    #     lp.lp.initialize()
    #     gzip.open(path, 'w').write(stringDump)


if __name__ == '__main__':
    # main('../models/iAF1260template2.lp', "/Volumes/Rudi1/Niko/20080926_NewRandomMediaSimulations_iAF1260/", 100000)
    # main2('../models/iJR904_rev_template.lp', '', 1)
    # lp = RandomMediaSimulations('../models/iAF1260template2.lp', 'R("R_Ec_biomass_iAF1260_core_59p81M")', [], minimizeRest=True)
    # print lp.almaas.lp.getObjective()
    obj = RandomMediaSimulations('../models/iAF1260template2.lp', 'R("R_Ec_biomass_iAF1260_core_59p81M")', INCLUDE, minimizeRest=False, optimizationRoutine='fba', koQ=True)
    result = obj.run()
    print result.fluxactivity
    print result.knockoutEffects

    obj = RandomMediaSimulations('../models/iAF1260template2.lp', 'R("R_Ec_biomass_iAF1260_core_59p81M")', INCLUDE, minimizeRest=False, optimizationRoutine='pFBA', koQ=True)
    result = obj.run()
    print result.fluxactivity
    print result.knockoutEffects

    
