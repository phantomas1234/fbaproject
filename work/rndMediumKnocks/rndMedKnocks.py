#!/usr/bin/env python
# encoding: utf-8
"""
rndMedKnocks.py

Created by Nikolaus Sonnenschein on 2008-05-02.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from ifba.GlpkWrap import util, metabolism, randomMedia, fluxdist
import copy
from ifba.general.util import sumDicts, filterDict, dict2mathematica, dict2tsv
import copy

class rndMedKnocks(object):
    """A class that defines an object that can be passed to distFBA client.
    
    It implements a run method that returns a string represenation of a result
    that can be send to a distFBA server by a distFBA client."""
    def __init__(self, lp, threshold=0.05):
        super(rndMedKnocks, self).__init__()
        self.lp = lp
        self.activityDict = dict()
        self.wtKnocks = dict()
        self.threshold = threshold
        
    def run(self, *args, **kwargs):
        """docstring for run"""
        med = randomMedia.Almaas(self.lp)
        fluxDist = med.generateFluxdist(minGrowth=0.5)
        self.lp.modifyColumnBounds(med.currDict)
        cond = med.currDict
        print cond
        self.lp.simplex()
        wt = self.lp.getObjVal()
        targets = fluxDist.getActiveReactions()
        for reac in targets:
            self.activityDict[reac] = 1
        print len(targets)
        # print targets
        print "Growth: ", wt
        for reac in targets:
            print "reaction", reac, "gets knocked out"
            self.lp.modifyColumnBounds({reac : (0,0)})
            self.lp.simplex()
            objVal = self.lp.getObjVal()
            print "Growth: ", objVal
            if (objVal / wt) < self.threshold:
                print 'mutant: ', objVal, ' wildtyp: ', wt, ' threshold: ', (wt * self.threshold)
                self.wtKnocks[reac] = 1
            self.lp.undo()
        self.lp.initialize()
        print self.lp.getObjectiveFunction()
        return (self.wtKnocks, self.activityDict, cond)


def test_without_networking(lp):
    d = dict()
    activityDict = dict()
    counter = 0
    while True:
        counter += 1
        test = rndMedKnocks(lp, threshold=0.05)
        data = test.run('Sensless Crap Argument Stub')
        d = sumDicts(d, data[0])
        activityDict = sumDicts(activityDict, data[1])
        print d, activityDict
        # prefix = "/Users/niko/tmp2newRndMedKnocks/"
        # path = prefix + "dict" + str(counter) + ".txt"
        # print "writing to " + path
        # open(path, 'w').write(dict2tsv(activityDict) + "\n" + dict2tsv(d) + "\n" + dict2tsv(data[2]))


def extractCondition(path):
    """docstring for extractCondition"""
    res2 = dict()
    f = open(path, 'r')
    res = f.readlines()
    res.reverse()
    f.close()
    for line in res:
        if line != '\n':
            tmp = line.rstrip().replace('[', '(').replace(']', ')').split('\t')
            res2[tmp[0]] = eval(tmp[1])
        else:
            break
    return res2
    

if __name__ == '__main__':
    path = '../models/iAF1260template.lp'
    path = '../models/iJR904template.lp'
    lp = metabolism.Metabolism(util.ImportCplex(path))
    lp.setReactionObjectiveMinimizeRest('R("R_BiomassEcoli")')
    #     # test_without_networking(copy.copy(lp))
    # paths = ['/Volumes/Rudi1/Niko/20081016_relativeEssentiality_0.05/dict'+str(i)+".txt" for i in (38, 39)]
    # print paths
    # conditions = [extractCondition(path) for path in paths]
    test_without_networking(copy.copy(lp))
    