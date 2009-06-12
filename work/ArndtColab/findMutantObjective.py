#!/usr/bin/env python
# encoding: utf-8
"""
findMutantObjective.py

Created by Nikolaus Sonnenschein on 2009-05-14.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

import sys
import os
from ifba.GlpkWrap import util, metabolism
from ifba.expressFBA.ContextFBA import ContextFBA, loadReactionData

def checkInconsistencyForAllReactions(data, level=.8, cutoff=.8):
    model_path = '../../work/HumanFBA/validHumanFBAmodel.lp'
    lp = util.ImportCplex(model_path)
    reacs = metabolism.Metabolism(lp).getReactions()
    print reacs
    print 2*"\n"
    testReacs =reacs[0:5]
    testReacs.append('R("R_Obj")')
    for r in testReacs:
        print r
        lp = util.ImportCplex(model_path)
        cntxtFBA = ContextFBA(lp, rmfID=r, level=.8)
        if cntxtFBA.getObjVal() > 0:
            for path in data:
                print path
                (cntxtFluxDist, incon) = cntxtFBA.contextFBA(loadReactionData(path), \
                cutoff=cutoff)
                print "results for: ", path
                print "ContextFBA fluxDist: ", cntxtFluxDist.getActiveFluxDist()
                print "Inconsistency Score: ", incon
                print 2*'\n'
                stuff = "\t".join(cntxtFluxDist.getActiveReactions())
        else:
            print "no rmfID flux!"



if __name__ == '__main__':
    data = ('/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/HB00ECR_EBER0_5.tsv', 
    "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/HB00ECR_EBER0_6.tsv")
    checkInconsistencyForAllReactions(data)
