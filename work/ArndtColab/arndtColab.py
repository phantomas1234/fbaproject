#!/usr/bin/env python
# encoding: utf-8
"""
arndtColab.py

Created by Nikolaus Sonnenschein on 2009-03-10.
Copyright (c) 2009 Jacobs University Bremen. All rights reserved.
"""

import sys
import os
from ifba.GlpkWrap import util, metabolism
from ifba.expressFBA.ContextFBA import ContextFBA, loadReactionData

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
    cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/", \
    "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/eberResults2.tsv")

