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
import copy

def openAllDoors(lp, defaulBound=20.):
    transp = lp.getTransporters()
    bounds = [(-999999, defaulBound) for i in transp]
    boundsDict = dict(zip(transp, bounds))
    lp.modifyColumnBounds(boundsDict)
    return copy.copy(lp)

def cntxtAnalysis(path, outpath, level=.8, cutoff=.8):
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
    cntxtFBA = ContextFBA(lp, rmfID='R("R_Obj")', level=level)
    for path in dataPaths:
        (cntxtFluxDist, incon) = cntxtFBA.contextFBA(loadReactionData(path), \
        cutoff=cutoff)
        print "results for: ", path
        print "ContextFBA fluxDist: ", cntxtFluxDist.getActiveFluxDist()
        print "Inconsistency Score: ", incon
        print 2*'\n'
        stuff = "\t".join(cntxtFluxDist.getActiveReactions())
        outputFile.write(path + '\t'+ str(incon) + "\t" + stuff + "\n")
    outputFile.close()
    
def cntxtAnalysisFullMedium(path, outpath, level=.8, cutoff=.8):
    dataPaths = tuple([file.rstrip() for file in os.popen("ls "+path+"*.tsv")])
    print dataPaths
    outputFile = open(outpath, 'w').close()
    outputFile = open(outpath, 'a')
    model_path = '../../work/HumanFBA/validHumanFBAmodel.lp'
    openModel = openAllDoors(metabolism.Metabolism(util.ImportCplex(model_path)))
    stdFBAfluxDist = openModel.fba()
    print "Standard FBA fluxDist: ", stdFBAfluxDist.getActiveFluxDist()
    print 2*"\n"
    reactionSets = list()
    openModel2 = copy.copy(openModel)
    cntxtFBA = ContextFBA(openModel2.lp, rmfID='R("R_Obj")', level=level)
    for path in dataPaths:
        (cntxtFluxDist, incon) = cntxtFBA.contextFBA(loadReactionData(path), \
        cutoff=cutoff)
        print "results for: ", path
        print "ContextFBA fluxDist: ", cntxtFluxDist.getActiveFluxDist()
        print "Inconsistency Score: ", incon
        print 2*'\n'
        stuff = "\t".join(cntxtFluxDist.getActiveReactions())
        outputFile.write(path + '\t'+ str(incon) + "\t" + stuff + "\n")
    outputFile.close()    

def analyzeData(paths, outputPath, method=cntxtAnalysis, cutoff=.5, level=.8):
    """docstring for analyzeData"""
    if outputPath[-1] is not "/":
        outputPath = outputPath + "/"
    for path in paths:
        outPutFilePath = "".join((outputPath, path.split("/")[-2], "_level", str(level), "_cutoff", str(cutoff), "_method.", method.__name__, ".tsv"))
        method(path, outPutFilePath, cutoff=cutoff, level=level)



if __name__ == '__main__':
    paths = ["/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/96cell-human/",
            "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/aids-human/",
            "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/apoptosis-human/",
            "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/lgcs-mouse/",
            "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/mcz-human/",
            "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/monkey-human/",
            "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/scd-human/",
            "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/pancreas-mouse/"
    ]
    print paths
    analyzeData(paths, "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/", method=cntxtAnalysisFullMedium)

    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/96cell-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/96cellResults_level.8_cutoff.5.tsv", cutoff=.5, level=.8)
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/96cell-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/96cellResults_level.8_cutoff.8.tsv", cutoff=.8, level=.8)
    # 
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/aids-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/aidsResults_level.8_cutoff.5.tsv", cutoff=.5, level=.8)
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/aids-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/aidsResults_level.8_cutoff.8.tsv", cutoff=.8, level=.8)
    # 
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/apoptosis-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/apoptosis-human_Results_level.8_cutoff.5.tsv", cutoff=.5, level=.8)
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/apoptosis-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/apoptosis-human_Results_level.8_cutoff.8.tsv", cutoff=.8, level=.8)
    # 
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/lgcs-mouse/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/lgcs-mouse_MouseResults_level.8_cutoff.5.tsv", cutoff=.5, level=.8)
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/lgcs-mouse/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/lgcs-mouse_MouseResults_level.8_cutoff.8.tsv", cutoff=.8, level=.8)
    # 
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/mcz-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/mcz-human_Results_level.8_cutoff.5.tsv", cutoff=.5, level=.8)
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/mcz-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/mcz-human_Results_level.8_cutoff.8.tsv", cutoff=.8, level=.8)
    # 
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/monkey-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/monkey-human_Results_level.8_cutoff.5.tsv", cutoff=.5, level=.8)
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/monkey-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/monkey-human_Results_level.8_cutoff.8.tsv", cutoff=.8, level=.8)
    # 
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/pancreas-mouse/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/pancreas-mouse_Results_level.8_cutoff.5.tsv", cutoff=.5, level=.8)
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/pancreas-mouse/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/pancreas-mouse_Results_level.8_cutoff.8.tsv", cutoff=.8, level=.8)
    # 
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/scd-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/scd-human_Results_level.8_cutoff.5.tsv", cutoff=.5, level=.8)
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/scd-human/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/scd-human_Results_level.8_cutoff.8.tsv", cutoff=.8, level=.8)





    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/96cell/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/96cellResults_level.8_cutoff.8.tsv", cutoff=.8, level=.8)

    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/96cell/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/96cellResults_level.8_cutoff.5.tsv", cutoff=.5, level=.8)


    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/EBER/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/eberResults_level.8_cutoff.5.tsv", cutoff=.5)
    # cntxtAnalysis("~/arbeit/Publishing/Arndt_Colab/data/contextData/WT/", \
    # "/Users/niko/arbeit/Publishing/Arndt_Colab/data/contextData/wtResults_level.8_cutoff.5.tsv", cutoff=.5)
