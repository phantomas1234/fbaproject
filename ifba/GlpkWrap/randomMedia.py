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

from ifba.GlpkWrap import fluxdist, metabolism, util
from ifba.GlpkWrap.fluxdist import FBAsimulationResult
from ifba.glpki.glpki import *
from ifba.general.util import sumDicts, randomString

# random.seed(100)
# random.seed(101)
# TODO: Remove if not needed

class RandomMediaSimulations(object):
    """A class simplifying random media simulations."""
    def __init__(self, path2template, objective, includes, minimizeRest=True):
        self.path2template = path2template
        self.objective = objective
        self.includes = includes
        tmpLP = metabolism.Metabolism(util.ImportCplex(self.path2template))
        if minimizeRest:
            tmpLP.setReactionObjectiveMinimizeRest(self.objective)
        else:
            tmpLP.setReactionObjective(self.objective)
        self.almaas = Almaas(copy.copy(tmpLP), alwaysInc=self.includes)
    
    def run(self, *args, **kwargs):
        """Run one random medium simulation and return a FBAsimulationResult object."""
        f = self.almaas.generateFluxdist()
        return FBAsimulationResult(f, self.almaas.lp.getColumnBounds(), self.almaas.lp.getObjectiveFunction(), time.time(), self.path2template)


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
    def __init__(self, lp, default_bound=20, percRange=(10, 100), alwaysInc=set()):
        super(Almaas, self).__init__(lp)
        self.lp = lp
        self.def_bnd = default_bound
        # self.lp.smcp.presolve = GLP_ON #TODO: I don't know if this is correct
        self.currDict = {}
        self.percRange = percRange
        self.debugDict = {}
        self.alwaysInc = alwaysInc
    
    def _rndBndDict(self, list):
        boundDict = {}
        rnd = random.uniform
        for i in list:
            boundDict[i] = (-rnd(0, self.def_bnd), rnd(0, self.def_bnd)) # TODO Fix it
            # boundDict[i] = (-rnd(0, self.def_bnd), rnd(0, self.def_bnd))
        return boundDict
    
    def _randomPercantage(self):
        # inc = 10
        num = len(self.lp.transporters)
        # val = range(self.percRange[0], self.percRange[1] + inc , inc)
        # return (num/100) * random.choice(val)
        return (num/100) * random.randrange(self.percRange[0], self.percRange[1])
    
    def generateRandCond(self, num):
        transp = self.lp.transporters
        sample = set(random.sample(transp, num)).union(self.alwaysInc)
        dict = self._rndBndDict(sample)
        self.lp.modifyColumnBounds(dict)
        self.currDict = dict
        # TODO: Caution! Debug stuff follows
        tmpDict = {}
        for i in self.currDict:
            tmpDict[i] = 1
        self.debugDict = sumDicts(self.debugDict, tmpDict)
    
    def generateFluxdist(self, minGrowth=0.05):
        growth = 0.
        num = self._randomPercantage()
        while growth <= minGrowth:
            self.generateRandCond(num)
            try:
                self.lp.simplex()
                growth = self.lp.getObjVal()
                self.lp.initialize()
            except Exception, msg:
                # print msg
                growth = 0.
                self.lp.initialize()
        return fluxdist.FluxDist(self.lp)


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
    main2('../models/iJR904_rev_template.lp', '', 1)
    
