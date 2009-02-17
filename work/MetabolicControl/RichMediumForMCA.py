#!/usr/bin/env python
# encoding: utf-8
"""
RichMediumForMCA.py

Created by Nikolaus Sonnenschein on 2009-01-27.
Copyright (c) 2009 Jacobs University of Bremen. All rights reserved.
"""

from ifba.GlpkWrap import metabolism, randomMedia, util
import copy


def generateRichMedium(transporters):
    richMed = dict()
    for t in transporters:
        richMed[t] = (-20., 20)
    return richMed

def main(path2template, resultsPath, runs):
    import gzip
    tmpLP = metabolism.Metabolism(util.ImportCplex(path2template))
    tmpLP.setReactionObjectiveMinimizeRest('R("R_Ec_biomass_iAF1260_core_59p81M")')
    lp = copy.copy(tmpLP)
    richMed = generateRichMedium(lp.getTransporters())
    print richMed
    lp.modifyColumnBounds(richMed)    
    richFluxDist = lp.fba()
    stringDump = randomMedia.dict2tsv(richMed) + "\n" + richFluxDist.tsv()
    print stringDump
    path = resultsPath + "iAf1260_fluxDist_RichMedium.tsv.gz"
    gzip.open(path, 'w').write(stringDump)

if __name__ == '__main__':
	main("../models/iAF1260template2.lp", 
	"/Users/niko/arbeit/Publishing/MetabolicControlPaper2/data/RandomFluxDistributions/", 500)

