#!/usr/bin/env python
# encoding: utf-8
"""
rpoZanalysis.py

Created by Nikolaus Sonnenschein on 2008-01-02.
Copyright (c) 2008 Jacobs University of Bremen. All rights reserved.
"""

import sys
import copy
import pprint
import yaml

from ifba import distFBA
from ifba.glpki.glpki import *
from combinatorics import combinatorics as comb
from ifba.distFBA import server

from ifba.GlpkWrap import metabolism, util, fluxdist, knockouts


def loadSource(path):
    """Loads a list of transportReactions. Format:
    R("Macgamb_Transp")
    R("Madnb_Transp")
    R("MalaDb_Transp")..."""
    file = open(path, 'r')
    sources = [line.strip() for line in file]
    file.close()
    return sources
    
def loadMinimalMed(path):
    """Loads a list of transportReactions. Format:
    R["R_ATPM"] -> {8.39, 8.39}
    R["Mo2b_Transp"] -> {0, 18.5}..."""
    file = open(path, 'r')
    content = eval(file.read())
    file.close()
    return dict(content)
    # return eval("{" + content + "}")

def init(path):
    """The iAF1260 template gets initialized."""
    struct = util.ImportCplex(path)
    blub = metabolism.Metabolism(struct)
    blub.toggleVerbosity()
    return blub

def main():
	g = comb.SetCombine(range(1,200),range(1,300)).generate()
	s = server.Server(g)
	s.run()
	
def insertTranspAndEvaluate(lp, transp):
    lp.modifyColumnBounds({transp : (0, 20)})
    lp.simplex()
    return lp.getObjVal()

def insertTranspAndKOall(lp, transp):
    lp.modifyColumnBounds({transp : (0, 20)})
    lp.simplex()
    lp.undo()
    return lp.getObjVal()


if __name__ == '__main__':

    # loading the targets
    targets = loadSource('./targetCarbonSources.txt')
    # loading the viable sources
    viable = loadSource('./viableCarbonSources.txt')
    # loading minimal medium conditions
    minMed = loadMinimalMed('./minimalMedium.txt')
    
    lp = init('../GlpkWrap/test_data/iAF1260template.lp')
    lp.modifyColumnBounds(minMed)
    lpready = copy.copy(lp)
    lpready.toggleVerbosity()
    
    bigCollector = dict()
    # e.g. the targets
    print "show of the targets obj values"
    print targets[0:2]
    for transp in targets[0::1]:
        lpready = copy.copy(lp)
        lpready.toggleVerbosity()
        oVal = insertTranspAndEvaluate(lpready, transp)
        fluxd = fluxdist.FluxDist(lpready)
        koTargets = fluxd.getActiveReactions()
        print len(koTargets)
        tmp = knockouts.KnockOut(lpready)
        smallCollector = list()
        for ko in koTargets:
            tmp.knockOut(ko)
            lpready.simplex()
            if lpready.getObjVal() > 0.02:
                smallCollector.append(ko)
            lpready.undo()
        print smallCollector, len(smallCollector)
        bigCollector[transp] = smallCollector
    print bigCollector
    open('result.txt', 'w').write(str(bigCollector))
        

    

